#!/usr/bin/env python3
"""Build a target BED from a GENCODE GTF.

For each requested gene, takes the CDS of the MANE Select transcript, pads each
interval by a fixed number of bases (to cover splice regions) and merges
overlapping intervals. Output is a 4-column BED (0-based, half-open) with the
gene symbol in the name column.

Usage:
    make_panel_bed.py --gtf chr21_gencode.gtf --genes APP,CXADR --pad 20 > panel.bed
"""
import argparse
import re
import sys
from collections import defaultdict


def parse_attrs(field):
    return dict(re.findall(r'(\S+) "([^"]*)"', field))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genes", required=True, help="comma-separated gene symbols")
    ap.add_argument("--pad", type=int, default=20)
    args = ap.parse_args()

    wanted = [g.strip() for g in args.genes.split(",") if g.strip()]
    cds = defaultdict(list)
    chrom_of = {}
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            if "MANE_Select" not in f[8]:
                continue
            a = parse_attrs(f[8])
            gene = a.get("gene_name")
            if gene not in wanted:
                continue
            cds[gene].append((int(f[3]) - 1, int(f[4])))  # GTF 1-based closed -> BED
            chrom_of[gene] = f[0]

    missing = [g for g in wanted if g not in cds]
    if missing:
        sys.exit(f"No MANE Select CDS found for: {', '.join(missing)}")

    rows = []
    for gene in wanted:
        ivs = sorted((max(0, s - args.pad), e + args.pad) for s, e in cds[gene])
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        rows += [(chrom_of[gene], s, e, gene) for s, e in merged]

    for chrom, s, e, gene in sorted(rows, key=lambda r: (r[0], r[1])):
        print(f"{chrom}\t{s}\t{e}\t{gene}")


if __name__ == "__main__":
    main()
