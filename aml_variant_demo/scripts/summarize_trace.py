#!/usr/bin/env python3
"""Combine Nextflow resource trace with retained work-directory disk use."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def main() -> None:
    trace = Path(sys.argv[1])
    output = Path(sys.argv[2])
    with trace.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = ["process", "status", "realtime", "%cpu", "peak_rss", "peak_vmem", "rchar", "wchar", "workdir_disk"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            workdir = row.get("workdir", "")
            disk = "NA"
            if workdir and Path(workdir).exists():
                disk = subprocess.run(
                    ["du", "-sh", workdir], check=True, text=True, stdout=subprocess.PIPE
                ).stdout.split()[0]
            writer.writerow({field: (disk if field == "workdir_disk" else row.get(field, "")) for field in fields})


if __name__ == "__main__":
    main()

