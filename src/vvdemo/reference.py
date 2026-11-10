"""Reference / genome-build consistency checks (REQ-003).

A wrong-build reference is a classic silent failure: alignment and calling
run fine, but coordinates do not match the truth set or downstream
annotation. The check compares contig names (and optionally lengths) in the
reference .fai index against what the run expects.
"""

from __future__ import annotations

from pathlib import Path

from .fastq import InputValidationError


def read_fai(fai: str | Path) -> dict[str, int]:
    contigs = {}
    with open(fai) as fh:
        for line in fh:
            if not line.strip():
                continue
            name, length = line.split("\t")[:2]
            contigs[name] = int(length)
    return contigs


def read_fasta_contigs(fasta: str | Path) -> dict[str, int]:
    """Contig lengths straight from a FASTA (used when no .fai exists yet)."""
    contigs: dict[str, int] = {}
    name = None
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                contigs[name] = 0
            elif name is not None:
                contigs[name] += len(line.strip())
    if not contigs:
        raise InputValidationError(f"{fasta}: no FASTA records found")
    return contigs


def bed_contigs(bed: str | Path) -> set[str]:
    out = set()
    with open(bed) as fh:
        for line in fh:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                out.add(line.split("\t")[0])
    return out


def check_reference(
    fasta: str | Path,
    expected_contigs: dict[str, int] | None = None,
    bed: str | Path | None = None,
) -> dict[str, int]:
    """Fail if the reference does not contain the expected contigs/lengths,
    or if the region BED refers to contigs the reference lacks."""
    fasta = Path(fasta)
    if not fasta.exists():
        raise InputValidationError(f"{fasta}: reference not found")
    contigs = read_fasta_contigs(fasta)

    if expected_contigs:
        for name, length in expected_contigs.items():
            if name not in contigs:
                raise InputValidationError(
                    f"{fasta}: expected contig '{name}' is missing "
                    f"(found {sorted(contigs)[:5]}); wrong genome build or naming (chr20 vs 20)?"
                )
            if length and contigs[name] != length:
                raise InputValidationError(
                    f"{fasta}: contig '{name}' length {contigs[name]} != expected {length}; wrong genome build?"
                )
    if bed:
        missing = bed_contigs(bed) - set(contigs)
        if missing:
            raise InputValidationError(
                f"{bed}: contigs {sorted(missing)} not in reference {fasta}; build or naming mismatch"
            )
    return contigs
