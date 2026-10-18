# Building the test data

These scripts are not part of the pipeline. They document how the files in
`assets/` were produced so the test data can be rebuilt or changed.

Inputs, all from [nf-core/test-datasets](https://github.com/nf-core/test-datasets)
(branch `modules`, folder `data/genomics/homo_sapiens/genome/chr21`):

- `sequence/genome.fasta` (GRCh38 chr21)
- `sequence/chr21_gencode.gtf`
- `germlineresources/1000G_phase1.snps.hg38.vcf.gz`
- `germlineresources/mills_and_1000G.indels.hg38.vcf.gz`

```bash
python3 dev/make_panel_bed.py \
    --gtf chr21_gencode.gtf \
    --genes APP,CXADR,NRIP1,BACH1,NCAM2,JAM2,ADAMTS1,CCT8,LTN1 \
    --pad 20 > assets/panel_chr21_demo.bed

python3 dev/simulate_panel_reads.py --fasta genome.fasta --bed assets/panel_chr21_demo.bed \
    --snps 1000G_phase1.snps.hg38.vcf.gz --indels mills_and_1000G.indels.hg38.vcf.gz \
    --sample SAMPLE01 --seed 11 --depth 150 --outdir assets/test

python3 dev/simulate_panel_reads.py --fasta genome.fasta --bed assets/panel_chr21_demo.bed \
    --snps 1000G_phase1.snps.hg38.vcf.gz --indels mills_and_1000G.indels.hg38.vcf.gz \
    --sample SAMPLE02 --seed 22 --depth 60 --outdir assets/test
```

`--depth` is the nominal mean before per-region capture variation is
applied, so the observed depth on target ends up higher (about 198x and 75x).
Both scripts use only the Python standard library.
