#!/usr/bin/env python3
"""Poll the local process table to capture peak RSS for Track 1 on macOS."""

from __future__ import annotations

import csv
import subprocess
import time
from pathlib import Path


OUTPUT = Path("results/track1/process_monitor.tsv")
PROCESS_NAMES = ("SLICE_REMOTE_BAMS", "MUTECT2_CALL", "FILTER_CALLS", "BENCHMARK_CALLS")


def belongs(name: str, command: str) -> bool:
    if name == "SLICE_REMOTE_BAMS":
        return any(token in command for token in ("samtools view", "samtools index", "samtools quickcheck"))
    if name == "MUTECT2_CALL":
        return "gatk-package" in command and "Mutect2" in command and "FilterMutectCalls" not in command
    if name == "FILTER_CALLS":
        return "gatk-package" in command and "FilterMutectCalls" in command
    return any(token in command for token in ("bcftools view", "bcftools norm", "bcftools isec", "calculate_metrics.py"))


def snapshot():
    result = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,rss=,command="],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    rows = []
    for line in result.stdout.splitlines():
        parts = line.strip().split(None, 3)
        if len(parts) == 4:
            rows.append((int(parts[0]), int(parts[1]), int(parts[2]), parts[3]))
    return rows


def write(stats):
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["process", "observed_seconds", "peak_rss_mib"])
        for name, values in stats.items():
            elapsed = max(0.0, values["last"] - values["first"]) if values["first"] else 0.0
            writer.writerow([name, f"{elapsed:.1f}", f"{values['peak_kib'] / 1024:.1f}"])


def main():
    stats = {name: {"first": 0.0, "last": 0.0, "peak_kib": 0} for name in PROCESS_NAMES}
    saw_pipeline = False
    absent = 0
    while True:
        now = time.monotonic()
        rows = snapshot()
        pipeline_alive = any("nextflow" in cmd and "run main.nf" in cmd for _, _, _, cmd in rows)
        saw_pipeline = saw_pipeline or pipeline_alive
        absent = 0 if pipeline_alive else absent + 1
        for name in PROCESS_NAMES:
            matching = [rss for _, _, rss, cmd in rows if belongs(name, cmd)]
            if matching:
                values = stats[name]
                values["first"] = values["first"] or now
                values["last"] = now
                values["peak_kib"] = max(values["peak_kib"], sum(matching))
        write(stats)
        if saw_pipeline and absent >= 10:
            break
        time.sleep(0.5)


if __name__ == "__main__":
    main()
