#!/usr/bin/env python3
"""Calculate exact-allele precision/recall/F1 from bcftools isec output."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def records(path: Path) -> int:
    opener = __import__("gzip").open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return sum(1 for line in handle if not line.startswith("#"))


def metrics(kind: str, directory: Path) -> dict[str, int | float | str]:
    fn = records(directory / "0000.vcf")
    fp = records(directory / "0001.vcf")
    tp = records(directory / "0002.vcf")
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "type": kind,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "F1": f1,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snv-dir", type=Path, required=True)
    parser.add_argument("--indel-dir", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-txt", type=Path, required=True)
    args = parser.parse_args()
    rows = [metrics("SNV", args.snv_dir), metrics("INDEL", args.indel_dir)]
    fields = ["type", "TP", "FP", "FN", "precision", "recall", "F1"]
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    lines = ["type\tTP\tFP\tFN\tprecision\trecall\tF1"]
    for row in rows:
        lines.append(
            f"{row['type']}\t{row['TP']}\t{row['FP']}\t{row['FN']}\t"
            f"{row['precision']:.4f}\t{row['recall']:.4f}\t{row['F1']:.4f}"
        )
    args.output_txt.write_text("\n".join(lines) + "\n")
    print(args.output_txt.read_text(), end="")


if __name__ == "__main__":
    main()

