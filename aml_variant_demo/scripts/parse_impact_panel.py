#!/usr/bin/env python3
"""Parse the official cBioPortal IMPACT468 response into HGNC symbols."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> None:
    source, output, report = map(Path, sys.argv[1:4])
    payload = json.loads(source.read_text())
    if payload.get("genePanelId") != "IMPACT468":
        raise ValueError(f"Unexpected panel: {payload.get('genePanelId')!r}")
    symbols = sorted({g["hugoGeneSymbol"].strip() for g in payload["genes"]})
    if not 450 <= len(symbols) <= 480:
        raise ValueError(f"Implausible IMPACT468 gene count: {len(symbols)}")
    output.write_text("\n".join(symbols) + "\n")
    report.write_text(
        f"panel_id\t{payload['genePanelId']}\n"
        f"description\t{payload.get('description', '')}\n"
        f"unique_symbols\t{len(symbols)}\n"
        f"contains_GTF2I\t{'yes' if 'GTF2I' in symbols else 'no'}\n"
    )


if __name__ == "__main__":
    main()

