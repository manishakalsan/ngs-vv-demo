#!/usr/bin/env bash
# Fetch a small GIAB HG002 (GRCh38) chr20 slice for the realistic benchmark.
#
# Run on YOUR machine (needs internet; ~25 MB reference + a few hundred MB of
# streamed BAM slice). Output goes to data/giab/ (git-ignored).
#
#   bash scripts/fetch_giab.sh                 # default 2 Mb region
#   REGION=chr20:10000000-11000000 bash scripts/fetch_giab.sh
#
# Then:  snakemake --cores 4 --configfile config/giab_hg002_chr20.yaml
#
# URLs: truth VCF/BED paths confirmed against public documentation on
# 2026-11-10. The BAM path is the commonly used NIST 2x250 novoalign BAM; if it
# 404s, browse the GIAB FTP (data/AshkenazimTrio/HG002_NA24385_son/) for a
# current GRCh38 Illumina BAM and set BAM_URL.
set -euo pipefail
cd "$(dirname "$0")/.."

REGION=${REGION:-chr20:10000000-12000000}
OUT=data/giab
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp
BAM_URL=${BAM_URL:-$GIAB/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam}
TRUTH_DIRS=("$GIAB/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
            "$GIAB/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38")
TRUTH_VCF=HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
TRUTH_BED=HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
REF_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz

for t in samtools bcftools curl; do command -v "$t" >/dev/null || { echo "missing $t (activate the conda env)"; exit 1; }; done
mkdir -p "$OUT"
CHR=${REGION%%:*}; RANGE=${REGION#*:}; START=${RANGE%-*}; END=${RANGE#*-}
echo "$CHR	$START	$END" > "$OUT/region.bed"

echo "[1/4] reference $CHR (UCSC hg38 = GRCh38 primary sequence)"
if [ ! -s "$OUT/chr20.fa" ]; then
  curl -fsSL "$REF_URL" | gunzip -c | tr 'acgtn' 'ACGTN' > "$OUT/chr20.fa"   # unmask soft-masked bases
fi
samtools faidx "$OUT/chr20.fa"
cut -f1,2 "$OUT/chr20.fa.fai"

echo "[2/4] truth VCF + confident BED (v4.2.1)"
got=""
for d in "${TRUTH_DIRS[@]}"; do
  if curl -fsSL -o "$OUT/$TRUTH_BED" "$d/$TRUTH_BED"; then
    bcftools view -R "$OUT/region.bed" -Oz -o "$OUT/truth.vcf.gz" "$d/$TRUTH_VCF"
    got=$d; break
  fi
done
[ -n "$got" ] || { echo "could not download truth set from GIAB FTP"; exit 1; }
bcftools index -t "$OUT/truth.vcf.gz"
# confident regions clipped to the analysed region
awk -v c="$CHR" -v s="$START" -v e="$END" 'BEGIN{OFS="\t"} $1==c && $3>s && $2<e {
  a=($2<s)?s:$2; b=($3>e)?e:$3; print $1,a,b }' "$OUT/$TRUTH_BED" > "$OUT/confident.bed"
echo "  truth variants in region: $(bcftools view -H "$OUT/truth.vcf.gz" | wc -l)"
echo "  confident bp in region:   $(awk '{s+=$3-$2} END{print s}' "$OUT/confident.bed")"

echo "[3/4] stream reads overlapping $REGION from the GIAB BAM (no full download)"
samtools view -b -o "$OUT/slice.bam" "$BAM_URL" "$REGION"
echo "[4/4] BAM slice -> paired FASTQ (pairs with a mate outside the region are dropped)"
samtools collate -u -O "$OUT/slice.bam" "$OUT/tmp_collate" \
  | samtools fastq -F 0x900 -1 "$OUT/HG002_R1.fastq.gz" -2 "$OUT/HG002_R2.fastq.gz" -0 /dev/null -s /dev/null -n -
rm -f "$OUT/slice.bam" ./*.bai
echo "  read pairs: $(( $(zcat "$OUT/HG002_R1.fastq.gz" | wc -l) / 4 ))"
echo
echo "Done. Run: snakemake --cores 4 --configfile config/giab_hg002_chr20.yaml"
echo "Optional haplotype-aware benchmark: see docs/giab_mode.md (hap.py)"
