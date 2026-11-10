#!/usr/bin/env python3
"""Generate the tiny synthetic test dataset with a known truth set.

Why synthetic: CI must run in seconds and without downloading GIAB data.
Truth is known exactly because we inserted the variants ourselves.
GIAB is used for the realistic benchmark (scripts/fetch_giab.sh).

Outputs (deterministic for a given --seed):
  ref.fa                 30 kb single-contig reference 'synth1'
  sample_R1.fastq.gz     paired 2x100 bp reads, ~40x, 0.2% substitution errors
  sample_R2.fastq.gz
  truth.vcf              inserted SNVs and indels with genotypes (normalised)
  confident.bed          confident region (excludes contig edges)
"""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

CONTIG = "synth1"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def make_reference(rng: random.Random, length: int) -> str:
    # ~45% GC, no long homopolymers (keeps indel placement unambiguous)
    seq = []
    for _ in range(length):
        while True:
            b = rng.choices("ACGT", weights=[0.275, 0.225, 0.225, 0.275])[0]
            if len(seq) < 3 or not (seq[-1] == seq[-2] == seq[-3] == b):
                break
        seq.append(b)
    return "".join(seq)


def make_variants(rng: random.Random, ref: str, n_snv: int, n_indel: int,
                  margin: int, spacing: int):
    """Return list of (pos1, ref_allele, alt_allele, gt). Left-aligned indels."""
    L = len(ref)
    used: list[int] = []
    variants = []

    def free(p):
        return all(abs(p - u) >= spacing for u in used)

    # one variant deliberately OUTSIDE the confident region (tests BED handling)
    p = 200
    alt = rng.choice([b for b in "ACGT" if b != ref[p - 1]])
    variants.append((p, ref[p - 1], alt, "0/1"))
    used.append(p)

    kinds = ["SNV"] * n_snv + ["INS"] * (n_indel // 2) + ["DEL"] * (n_indel - n_indel // 2)
    rng.shuffle(kinds)
    for kind in kinds:
        for _attempt in range(10_000):
            p = rng.randint(margin, L - margin)  # 1-based anchor position
            if not free(p):
                continue
            gt = rng.choice(["0/1", "0/1", "1/1"])
            anchor = ref[p - 1]
            if kind == "SNV":
                alt = rng.choice([b for b in "ACGT" if b != anchor])
                variants.append((p, anchor, alt, gt))
            elif kind == "DEL":
                k = rng.randint(1, 6)
                deleted = ref[p:p + k]
                # left-aligned iff anchor base != last deleted base
                if anchor == deleted[-1]:
                    continue
                variants.append((p, anchor + deleted, anchor, gt))
            else:  # INS
                k = rng.randint(1, 6)
                ins = "".join(rng.choice("ACGT") for _ in range(k))
                # left-aligned iff last inserted base != anchor; also avoid
                # the insertion merely repeating following reference bases
                if ins[-1] == anchor or ins[0] == ref[p]:
                    continue
                variants.append((p, anchor, anchor + ins, gt))
            used.append(p)
            break
        else:
            raise RuntimeError("could not place variant; reduce counts")
    return sorted(variants)


def apply_haplotype(ref: str, variants, hap_index: int) -> str:
    out, cursor = [], 0
    for pos, r, a, gt in variants:
        alleles = gt.split("/")
        if alleles[hap_index] != "1":
            continue
        start = pos - 1
        out.append(ref[cursor:start])
        out.append(a)
        cursor = start + len(r)
    out.append(ref[cursor:])
    return "".join(out)


def add_errors(rng, seq, rate):
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def simulate_reads(rng, haps, n_pairs, read_len, frag_mean, frag_sd, err):
    r1, r2 = [], []
    for i in range(n_pairs):
        hap = haps[rng.randrange(2)]
        frag = max(read_len + 20, int(rng.gauss(frag_mean, frag_sd)))
        start = rng.randint(0, len(hap) - frag)
        f = hap[start:start + frag]
        if rng.random() < 0.5:
            f = revcomp(f)
        a = add_errors(rng, f[:read_len], err)
        b = add_errors(rng, revcomp(f)[:read_len], err)
        q = "F" * read_len  # Phred 37
        r1.append(f"@synth_{i}/1\n{a}\n+\n{q}\n")
        r2.append(f"@synth_{i}/2\n{b}\n+\n{q}\n")
    return r1, r2


def write_gz(path: Path, lines):
    # mtime=0 makes the .gz byte-identical across runs (reproducible data)
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write("".join(lines).encode())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="data/synthetic")
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--length", type=int, default=30_000)
    ap.add_argument("--snvs", type=int, default=40)
    ap.add_argument("--indels", type=int, default=12)
    ap.add_argument("--coverage", type=float, default=40)
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--error-rate", type=float, default=0.002)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    ref = make_reference(rng, args.length)
    variants = make_variants(rng, ref, args.snvs, args.indels, margin=800, spacing=150)
    haps = [apply_haplotype(ref, variants, 0), apply_haplotype(ref, variants, 1)]
    n_pairs = int(args.length * args.coverage / (2 * args.read_len))
    r1, r2 = simulate_reads(rng, haps, n_pairs, args.read_len, 300, 30, args.error_rate)

    with open(out / "ref.fa", "w") as fh:
        fh.write(f">{CONTIG}\n")
        for i in range(0, len(ref), 60):
            fh.write(ref[i:i + 60] + "\n")
    write_gz(out / "sample_R1.fastq.gz", r1)
    write_gz(out / "sample_R2.fastq.gz", r2)

    with open(out / "truth.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={CONTIG},length={len(ref)}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write(f"##source=simulate_test_data.py seed={args.seed}\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n")
        for pos, r, a, gt in variants:
            fh.write(f"{CONTIG}\t{pos}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{gt}\n")

    with open(out / "confident.bed", "w") as fh:
        fh.write(f"{CONTIG}\t500\t{len(ref) - 500}\n")

    n_snv = sum(1 for v in variants if len(v[1]) == len(v[2]) == 1)
    print(f"wrote {out}: {len(variants)} variants ({n_snv} SNV, {len(variants) - n_snv} indel), "
          f"{n_pairs} read pairs")


if __name__ == "__main__":
    main()
