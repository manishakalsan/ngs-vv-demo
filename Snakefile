# ngs-vv-demo: a small germline variant-calling pipeline built to be verified.
#
#   FASTQ -> validate -> fastp -> BWA-MEM -> samtools sort/index
#         -> bcftools mpileup/call -> norm -> filter -> VCF -> evaluate vs truth
#
# Run:  snakemake --cores 2
# Another config:   snakemake --cores 2 --configfile config/changes/v1.1_skip_indels.yaml

from pathlib import Path

configfile: str(Path(workflow.basedir) / "config" / "config.yaml")

REPO = Path(workflow.basedir)
OUT = config["outdir"]
# --configfile is merged onto config.yaml, so a run can name the samples it processes
SAMPLES = config.get("run_samples") or list(config["samples"])
C = config["caller"]
F = config["fastp"]
T = int(config.get("threads", 2))


def p(path):
    """Resolve config paths relative to the repo, so tests can run from tmp dirs."""
    path = Path(path)
    return str(path if path.is_absolute() else REPO / path)


REF = f"{OUT}/reference/ref.fa"


rule all:
    input:
        expand(f"{OUT}/{{s}}/{{s}}.vcf.gz.tbi", s=SAMPLES),
        expand(f"{OUT}/{{s}}/eval/metrics.json", s=SAMPLES),
        f"{OUT}/software_versions.txt",


rule validate_inputs:
    """REQ-001/002/003: fail fast, before any compute, on bad inputs."""
    input:
        r1=lambda w: p(config["samples"][w.s]["r1"]),
        r2=lambda w: p(config["samples"][w.s]["r2"]),
        ref=p(config["reference"]),
    output:
        f"{OUT}/{{s}}/qc/input_validation.json",
    params:
        expected=" ".join(f"--expected-contig {c}" for c in config["expected_contigs"]),
        beds=" ".join(f"--bed {p(b)}" for b in [config["confident_bed"], config.get("target_bed")] if b),
        script=REPO / "scripts/validate_inputs.py",
    shell:
        "python {params.script} --r1 {input.r1} --r2 {input.r2} --reference {input.ref} "
        "{params.expected} {params.beds} --out {output}"


rule prepare_reference:
    input:
        p(config["reference"]),
    output:
        fa=REF,
        fai=REF + ".fai",
        bwt=REF + ".bwt",
    log:
        f"{OUT}/logs/prepare_reference.log",
    shell:
        "cp {input} {output.fa} && samtools faidx {output.fa} && bwa index {output.fa} 2> {log}"


rule fastp:
    input:
        ok=f"{OUT}/{{s}}/qc/input_validation.json",
        r1=lambda w: p(config["samples"][w.s]["r1"]),
        r2=lambda w: p(config["samples"][w.s]["r2"]),
    output:
        r1=temp(f"{OUT}/{{s}}/trimmed/{{s}}_R1.fastq.gz"),
        r2=temp(f"{OUT}/{{s}}/trimmed/{{s}}_R2.fastq.gz"),
        json=f"{OUT}/{{s}}/qc/{{s}}.fastp.json",
        html=f"{OUT}/{{s}}/qc/{{s}}.fastp.html",
    log:
        f"{OUT}/logs/{{s}}.fastp.log",
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-q {F[qualified_quality]} -l {F[min_length]} -w 1 "
        "-j {output.json} -h {output.html} 2> {log}"


rule align:
    """-K fixes the batch size so output does not depend on thread count (REQ-006)."""
    input:
        r1=f"{OUT}/{{s}}/trimmed/{{s}}_R1.fastq.gz",
        r2=f"{OUT}/{{s}}/trimmed/{{s}}_R2.fastq.gz",
        ref=REF,
        bwt=REF + ".bwt",
    output:
        bam=f"{OUT}/{{s}}/{{s}}.bam",
        bai=f"{OUT}/{{s}}/{{s}}.bam.bai",
    threads: T
    log:
        f"{OUT}/logs/{{s}}.align.log",
    shell:
        "(bwa mem -t {threads} -K 10000000 -R '@RG\\tID:{wildcards.s}\\tSM:{wildcards.s}\\tPL:ILLUMINA' "
        "{input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output.bam} - "
        "&& samtools index {output.bam}) 2> {log}"


rule call_variants:
    input:
        bam=f"{OUT}/{{s}}/{{s}}.bam",
        bai=f"{OUT}/{{s}}/{{s}}.bam.bai",
        ref=REF,
        fai=REF + ".fai",
    output:
        vcf=f"{OUT}/{{s}}/{{s}}.vcf.gz",
        tbi=f"{OUT}/{{s}}/{{s}}.vcf.gz.tbi",
    params:
        indels="" if C["call_indels"] else "-I",
        regions=f"-R {p(config['target_bed'])}" if config.get("target_bed") else "",
    log:
        f"{OUT}/logs/{{s}}.call.log",
    shell:
        "(bcftools mpileup -f {input.ref} -q {C[min_mapq]} -Q {C[min_baseq]} -a AD,DP "
        "{params.indels} {params.regions} {input.bam} "
        "| bcftools call -mv "
        "| bcftools norm -m- -f {input.ref} "
        "| bcftools filter -s LowQual -e 'QUAL<{C[min_qual]}' -Oz -o {output.vcf} "
        "&& bcftools index -t {output.vcf}) 2> {log}"


rule evaluate:
    """Compare to truth in confident regions -> PPA/PPV per variant class."""
    input:
        vcf=f"{OUT}/{{s}}/{{s}}.vcf.gz",
        truth=p(config["truth_vcf"]),
        bed=p(config["confident_bed"]),
    output:
        json=f"{OUT}/{{s}}/eval/metrics.json",
        tsv=f"{OUT}/{{s}}/eval/metrics.tsv",
    params:
        script=REPO / "scripts/evaluate.py",
    shell:
        "python {params.script} --calls {input.vcf} --truth {input.truth} --bed {input.bed} "
        "--json {output.json} --tsv {output.tsv}"


rule software_versions:
    """Configuration management: record exactly which tool versions produced the results."""
    output:
        f"{OUT}/software_versions.txt",
    params:
        v=config["pipeline_version"],
    shell:
        "( echo 'ngs-vv-demo {params.v}'; "
        "echo \"snakemake $(snakemake --version)\"; "
        "fastp --version 2>&1; "
        "echo \"bwa $(bwa 2>&1 | grep Version | cut -d' ' -f2)\"; "
        "samtools --version | head -1; "
        "bcftools --version | head -1 ) > {output}"
