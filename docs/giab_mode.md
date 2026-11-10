# GIAB benchmark mode (realistic data)

The synthetic dataset proves the machinery works; GIAB shows how the pipeline does on a real human genome. Run it on your own machine; it is not part of CI because of download size.

## Steps

```bash
micromamba activate ngs-vv-demo
bash scripts/fetch_giab.sh                                   # ~5–15 min, mostly network
snakemake --cores 4 --configfile config/giab_hg002_chr20.yaml
cat results_giab/HG002/eval/metrics.tsv
```

`fetch_giab.sh` downloads GRCh38 chr20 (UCSC), the GIAB v4.2.1 HG002 benchmark VCF and confident BED (sliced to the region), and streams only the reads in the region from a GIAB BAM. Default region: `chr20:10,000,000–12,000,000` (2 Mb).

## What to expect

bcftools on real data scores below the synthetic 1.0. Typical results for bcftools on high-coverage Illumina HG002 are high-0.9s for SNVs and noticeably lower for indels. Record what you get, then decide thresholds for this dataset in a separate locked file; do not reuse the synthetic thresholds without evidence.

## Haplotype-aware comparison with hap.py (recommended for GIAB)

The built-in comparator matches exact normalised alleles. On real data, complex variants can be written differently and get counted as FP + FN. hap.py resolves this:

```bash
docker run --rm -v "$PWD:/data" jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  /data/data/giab/truth.vcf.gz /data/results_giab/HG002/HG002.vcf.gz \
  -f /data/data/giab/confident.bed -r /data/data/giab/chr20.fa \
  -T /data/data/giab/region.bed -o /data/results_giab/HG002/happy
```

Comparing the built-in metrics with hap.py's is itself a useful finding to discuss: it validates (or exposes the limits of) the test tool.

## Stratification (closes risk H-09)

GIAB publishes stratification BEDs (homopolymers, low mappability, segmental duplications, GC extremes). Passing them to hap.py with `--stratification` reports accuracy per region type, which is how a real validation shows where the pipeline is weak.

The BAM URL in the script could not be verified in the build environment; if it fails, set `BAM_URL` to a current GRCh38 Illumina HG002 BAM from the GIAB FTP.
