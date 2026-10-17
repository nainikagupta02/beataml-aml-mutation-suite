#!/usr/bin/env python3
"""Extract a four-gene GENCODE target transcriptome and matching GTF."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


GENES = {"RUNX1", "RUNX1T1", "PML", "RARA"}
GENE_RE = re.compile(r'(?:^|;\s*)gene_name "([^"]+)";')
TRANSCRIPT_RE = re.compile(r'(?:^|;\s*)transcript_id "([^"]+)";')


def opener(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--out-gtf", type=Path, required=True)
    parser.add_argument("--out-fasta", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    transcripts: set[str] = set()
    genes_seen: set[str] = set()
    gtf_records = 0
    with opener(args.gtf) as src, args.out_gtf.open("w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            gene_match = GENE_RE.search(line)
            if not gene_match or gene_match.group(1) not in GENES:
                continue
            dst.write(line)
            genes_seen.add(gene_match.group(1))
            transcript_match = TRANSCRIPT_RE.search(line)
            if transcript_match:
                transcripts.add(transcript_match.group(1))
            gtf_records += 1

    fasta_transcripts: set[str] = set()
    write = False
    with opener(args.fasta) as src, args.out_fasta.open("w") as dst:
        for line in src:
            if line.startswith(">"):
                transcript = line[1:].split("|", 1)[0]
                write = transcript in transcripts
                if write:
                    fasta_transcripts.add(transcript)
            if write:
                dst.write(line)

    missing_genes = GENES - genes_seen
    missing_transcripts = transcripts - fasta_transcripts
    if missing_genes or missing_transcripts:
        raise RuntimeError(
            f"reference mismatch: missing genes={sorted(missing_genes)}, "
            f"missing transcripts={len(missing_transcripts)}"
        )
    args.report.write_text(
        f"genes\t{','.join(sorted(genes_seen))}\n"
        f"transcripts\t{len(transcripts)}\n"
        f"gtf_records\t{gtf_records}\n"
        f"fasta_transcripts\t{len(fasta_transcripts)}\n"
    )


if __name__ == "__main__":
    main()
