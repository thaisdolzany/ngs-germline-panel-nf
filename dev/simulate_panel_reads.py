#!/usr/bin/env python3
"""Simulate a diploid targeted-sequencing sample with a known truth set.

Truth variants are drawn from real population catalogues (1000 Genomes SNPs and
Mills/1000G indels) that fall inside the panel targets. Each variant receives a
genotype (0/1 or 1/1), two haplotypes are built per target region, and
paired-end reads are sampled from those haplotypes with a simple error model
(higher substitution rate and lower base quality at the 3' end of each read).
Capture efficiency varies between regions (log-normal), and a small fraction of
regions is simulated at very low depth to mimic capture dropouts.

Only the Python standard library is used, so the output is reproducible from
the seed alone.

Outputs:
    <out>/<sample>_R1.fastq.gz
    <out>/<sample>_R2.fastq.gz
    <out>/<sample>.truth.vcf   (sorted, 1-based, with GT)
"""
import argparse
import gzip
import os
import random


COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def read_bed(path):
    with open(path) as fh:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.split("\t") for l in fh if l.strip())]


def in_targets(chrom, pos0, end0, targets):
    return any(c == chrom and s <= pos0 and end0 <= e for c, s, e in targets)


def load_candidates(vcf_gz, targets, ref, kind, max_indel):
    out = []
    with gzip.open(vcf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 6)
            chrom, pos, ref_a, alt = f[0], int(f[1]), f[3].upper(), f[4].upper()
            if "," in alt or "N" in ref_a + alt or "<" in alt or "*" in alt:
                continue
            is_snv = len(ref_a) == 1 and len(alt) == 1
            if kind == "snv" and not is_snv:
                continue
            if kind == "indel":
                if is_snv or len(ref_a) == len(alt) or abs(len(ref_a) - len(alt)) > max_indel:
                    continue
                if ref_a[0] != alt[0]:
                    continue  # keep left-anchored indels only
            p0 = pos - 1
            if chrom not in ref or ref[chrom][p0:p0 + len(ref_a)] != ref_a:
                continue
            if in_targets(chrom, p0, p0 + len(ref_a), targets):
                out.append((chrom, pos, ref_a, alt))
    return out


def pick_variants(cands, n, min_gap, rng, taken):
    rng.shuffle(cands)
    chosen = []
    for v in cands:
        if len(chosen) >= n:
            break
        chrom, pos, ref_a, _ = v
        if any(c == chrom and abs(p - pos) < min_gap + len(r) for c, p, r, _ in taken + chosen):
            continue
        chosen.append(v)
    return chosen


def build_haplotype(seq, start0, end0, variants):
    """Apply variants (list of (pos1, ref, alt)) to seq[start0:end0]."""
    out, cur = [], start0
    for pos1, ref_a, alt in sorted(variants):
        p0 = pos1 - 1
        if p0 < cur or p0 + len(ref_a) > end0:
            continue
        out.append(seq[cur:p0])
        out.append(alt)
        cur = p0 + len(ref_a)
    out.append(seq[cur:end0])
    return "".join(out)


def mutate(read, rng, err_body, err_tail, tail_len):
    bases, quals = list(read), []
    n = len(bases)
    for i in range(n):
        tail = i >= n - tail_len
        err = err_tail if tail else err_body
        q = 20 if tail else 36
        if rng.random() < err:
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
            q = rng.choice([8, 12, 15]) if tail else 25
        quals.append(chr(q + 33))
    return "".join(bases), "".join(quals)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--snps", required=True, help="bgzipped VCF of candidate SNVs")
    ap.add_argument("--indels", required=True, help="bgzipped VCF of candidate indels")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-snv", type=int, default=60)
    ap.add_argument("--n-indel", type=int, default=15)
    ap.add_argument("--max-indel", type=int, default=12)
    ap.add_argument("--depth", type=float, default=120.0, help="mean diploid depth")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--frag-mean", type=int, default=300)
    ap.add_argument("--frag-sd", type=int, default=30)
    ap.add_argument("--flank", type=int, default=300, help="bases simulated beyond each target")
    ap.add_argument("--capture-sd", type=float, default=0.35, help="log-normal sigma of per-region depth")
    ap.add_argument("--dropout", type=float, default=0.03, help="fraction of regions simulated at 8%% depth")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    ref = read_fasta(args.fasta)
    targets = read_bed(args.bed)

    snv_c = load_candidates(args.snps, targets, ref, "snv", args.max_indel)
    ind_c = load_candidates(args.indels, targets, ref, "indel", args.max_indel)
    indels = pick_variants(ind_c, args.n_indel, 20, rng, [])
    snvs = pick_variants(snv_c, args.n_snv, 20, rng, indels)
    truth = []
    for v in indels + snvs:
        gt = rng.choice(["0/1", "0/1", "1/1"])
        truth.append(v + (gt,))
    truth.sort(key=lambda t: (t[0], t[1]))

    # simulation regions: targets plus flanks, merged
    regs = sorted((c, max(0, s - args.flank), min(len(ref[c]), e + args.flank)) for c, s, e in targets)
    merged = [list(regs[0])]
    for c, s, e in regs[1:]:
        if c == merged[-1][0] and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([c, s, e])

    os.makedirs(args.outdir, exist_ok=True)
    r1 = gzip.open(os.path.join(args.outdir, f"{args.sample}_R1.fastq.gz"), "wt", compresslevel=9)
    r2 = gzip.open(os.path.join(args.outdir, f"{args.sample}_R2.fastq.gz"), "wt", compresslevel=9)
    rl, n_pairs = args.read_len, 0
    for chrom, s, e in merged:
        eff = rng.lognormvariate(0, args.capture_sd)
        if rng.random() < args.dropout:
            eff = 0.08
        hap_vars = [[], []]
        for c, pos, ra, alt, gt in truth:
            if c != chrom or not (s <= pos - 1 < e):
                continue
            if gt == "1/1":
                hap_vars[0].append((pos, ra, alt))
                hap_vars[1].append((pos, ra, alt))
            else:
                hap_vars[rng.randrange(2)].append((pos, ra, alt))
        for h in (0, 1):
            hap = build_haplotype(ref[chrom], s, e, hap_vars[h])
            if len(hap) < args.frag_mean + 3 * args.frag_sd:
                continue
            pairs = int(round(eff * (args.depth / 2.0) * len(hap) / (2 * rl)))
            for _ in range(pairs):
                flen = max(rl + 10, min(len(hap), int(rng.gauss(args.frag_mean, args.frag_sd))))
                st = rng.randrange(0, len(hap) - flen + 1)
                frag = hap[st:st + flen]
                if rng.random() < 0.5:
                    frag = revcomp(frag)
                a, qa = mutate(frag[:rl], rng, 0.001, 0.01, 20)
                b, qb = mutate(revcomp(frag)[:rl], rng, 0.001, 0.01, 20)
                n_pairs += 1
                name = f"{args.sample}:{n_pairs}"
                r1.write(f"@{name}/1\n{a}\n+\n{qa}\n")
                r2.write(f"@{name}/2\n{b}\n+\n{qb}\n")
    r1.close()
    r2.close()

    with open(os.path.join(args.outdir, f"{args.sample}.truth.vcf"), "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c in sorted({t[0] for t in truth}):
            fh.write(f"##contig=<ID={c},length={len(ref[c])}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write(f"##source=simulate_panel_reads.py seed={args.seed}\n")
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")
        for c, pos, ra, alt, gt in truth:
            fh.write(f"{c}\t{pos}\t.\t{ra}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")

    n_ind = sum(1 for t in truth if len(t[2]) != len(t[3]))
    print(f"{args.sample}: {n_pairs} read pairs, {len(truth)} truth variants "
          f"({len(truth) - n_ind} SNV, {n_ind} indel)")


if __name__ == "__main__":
    main()
