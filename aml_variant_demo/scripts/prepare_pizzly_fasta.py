#!/usr/bin/env python3
"""Make a GENCODE transcript FASTA compatible with pizzly.

GENCODE places pipe-delimited metadata in each FASTA record name.  Pizzly
expects the first whitespace-delimited token to be the transcript_id used in
the GTF, so pipes in header lines must be replaced before building the
kallisto index.  Sequence lines are copied byte-for-byte.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    records = 0
    changed = 0
    with args.input.open("rt") as source, args.output.open("wt") as target:
        for line in source:
            if line.startswith(">"):
                records += 1
                fixed = line.replace("|", " ")
                changed += fixed != line
                target.write(fixed)
            else:
                target.write(line)

    if not records or changed != records:
        raise RuntimeError(
            f"expected every FASTA header to require conversion; "
            f"records={records}, changed={changed}"
        )

    with args.output.open() as handle:
        first_header = next(line.rstrip() for line in handle if line.startswith(">"))
    first_token = first_header[1:].split()[0]
    if "|" in first_token or not first_token.startswith("ENST"):
        raise RuntimeError(f"invalid corrected FASTA identifier: {first_token}")

    args.report.write_text(
        "metric\tvalue\n"
        f"records\t{records}\n"
        f"headers_changed\t{changed}\n"
        f"first_identifier\t{first_token}\n"
        f"input_sha256\t{sha256(args.input)}\n"
        f"output_sha256\t{sha256(args.output)}\n"
    )


if __name__ == "__main__":
    main()
