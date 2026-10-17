#!/usr/bin/env python3
"""Create exact-symbol exon BEDs and apply the preregistered truth-count gate."""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
from pathlib import Path


GENE_RE = re.compile(r'(?:^|;\s*)gene_name "([^"]+)";')
GENE_TYPE_RE = re.compile(r'(?:^|;\s*)gene_type "([^"]+)";')


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def write_exons(
    gtf: Path,
    output: Path,
    selected: set[str] | None,
    protein_coding_only: bool = False,
) -> tuple[set[str], int]:
    found: set[str] = set()
    count = 0
    with open_text(gtf) as src, output.open("w") as dst:
        for line in src:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            match = GENE_RE.search(fields[8])
            if not match:
                continue
            gene = match.group(1)
            if selected is not None and gene not in selected:
                continue
            if protein_coding_only:
                type_match = GENE_TYPE_RE.search(fields[8])
                if not type_match or type_match.group(1) != "protein_coding":
                    continue
            chrom = fields[0]
            if not re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)", chrom):
                continue
            dst.write(f"{chrom}\t{int(fields[3]) - 1}\t{fields[4]}\t{gene}\n")
            found.add(gene)
            count += 1
    return found, count


def run(argv: list[str], stdout: Path | None = None) -> None:
    if stdout is None:
        subprocess.run(argv, check=True)
    else:
        with stdout.open("w") as handle:
            subprocess.run(argv, check=True, stdout=handle)


def sort_merge(source: Path, output: Path) -> None:
    sorted_path = output.with_suffix(".sorted.tmp")
    run(["bedtools", "sort", "-i", str(source)], sorted_path)
    run(["bedtools", "merge", "-i", str(sorted_path)], output)
    sorted_path.unlink()


def intersect(a: Path, b: Path, output: Path) -> None:
    run(["bedtools", "intersect", "-a", str(a), "-b", str(b)], output)


def count_truth(vcf: Path, bed: Path) -> int:
    result = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", str(vcf), "-b", str(bed)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return sum(1 for line in result.stdout.splitlines() if line and not line.startswith("#"))


def bed_stats(path: Path) -> tuple[int, int]:
    intervals = bases = 0
    with path.open() as handle:
        for line in handle:
            chrom, start, end, *_ = line.rstrip().split("\t")
            intervals += 1
            bases += int(end) - int(start)
    return intervals, bases


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--truth-bed", type=Path, required=True)
    parser.add_argument("--snv", type=Path, required=True)
    parser.add_argument("--indel", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--min-truth", type=int, default=200)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    genes = {line.strip() for line in args.genes.read_text().splitlines() if line.strip()}
    raw_panel = args.outdir / "impact_exons.raw.bed"
    found, raw_count = write_exons(args.gtf, raw_panel, genes)
    missing = sorted(genes - found)
    (args.outdir / "impact_genes_missing_from_gencode.txt").write_text(
        "\n".join(missing) + ("\n" if missing else "")
    )
    panel_exons = args.outdir / "impact_exons.hg38.bed"
    sort_merge(raw_panel, panel_exons)
    raw_panel.unlink()
    panel_validation = args.outdir / "impact_exons.in_truth.bed"
    intersect(panel_exons, args.truth_bed, panel_validation)

    panel_snv = count_truth(args.snv, panel_validation)
    panel_indel = count_truth(args.indel, panel_validation)
    selected = panel_validation
    region_name = "IMPACT468 high-confidence exons"

    if panel_snv + panel_indel < args.min_truth:
        raw_all = args.outdir / "all_basic_exons.raw.bed"
        _, all_raw_count = write_exons(args.gtf, raw_all, None, protein_coding_only=True)
        all_exons = args.outdir / "all_basic_exons.hg38.bed"
        sort_merge(raw_all, all_exons)
        raw_all.unlink()
        selected = args.outdir / "all_basic_exons.in_truth.bed"
        intersect(all_exons, args.truth_bed, selected)
        region_name = "all GENCODE v50 basic protein-coding high-confidence exons"
    else:
        all_raw_count = 0

    final_bed = args.outdir / "track1_validation.bed"
    final_bed.write_bytes(selected.read_bytes())
    final_snv = count_truth(args.snv, final_bed)
    final_indel = count_truth(args.indel, final_bed)
    intervals, bases = bed_stats(final_bed)

    genome = args.outdir / "reference.genome"
    with args.fai.open() as src, genome.open("w") as dst:
        for line in src:
            chrom, length, *_ = line.split("\t")
            dst.write(f"{chrom}\t{length}\n")
    padded = args.outdir / "track1_slice.padded200.bed"
    padded_raw = args.outdir / "padded.tmp.bed"
    run(["bedtools", "slop", "-i", str(final_bed), "-g", str(genome), "-b", "200"], padded_raw)
    sort_merge(padded_raw, padded)
    padded_raw.unlink()

    (args.outdir / "region_stats.tsv").write_text(
        "region_set\tintervals\tbases\tsnv_truth\tindel_truth\tcombined_truth\n"
        f"impact_panel\t{bed_stats(panel_validation)[0]}\t{bed_stats(panel_validation)[1]}\t"
        f"{panel_snv}\t{panel_indel}\t{panel_snv + panel_indel}\n"
        f"final\t{intervals}\t{bases}\t{final_snv}\t{final_indel}\t{final_snv + final_indel}\n"
    )
    (args.outdir / "gate_report.txt").write_text(
        f"Configured minimum combined truth variants: {args.min_truth}\n"
        f"IMPACT source symbols: {len(genes)}\n"
        f"IMPACT symbols found in GENCODE v50: {len(found)}\n"
        f"IMPACT symbols missing from GENCODE v50: {len(missing)}\n"
        f"Raw IMPACT exon records: {raw_count}\n"
        f"IMPACT gate: SNV={panel_snv}, INDEL={panel_indel}, combined={panel_snv + panel_indel}\n"
        f"Gate decision: {'PASS' if panel_snv + panel_indel >= args.min_truth else 'WIDEN'}\n"
        f"Final region: {region_name}\n"
        f"Final BED: intervals={intervals}, bases={bases}\n"
        f"Final truth: SNV={final_snv}, INDEL={final_indel}, combined={final_snv + final_indel}\n"
        f"All-exon raw records parsed when widening: {all_raw_count}\n"
    )
    print((args.outdir / "gate_report.txt").read_text(), end="")
    if final_snv + final_indel < args.min_truth:
        raise RuntimeError(
            f"Final region still has only {final_snv + final_indel} truth variants; refusing to call"
        )


if __name__ == "__main__":
    main()
