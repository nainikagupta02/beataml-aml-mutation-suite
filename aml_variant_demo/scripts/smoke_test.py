#!/usr/bin/env python3
"""Run installation smoke tests without downloading biological data."""

from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path


PREFIX = Path(__file__).resolve().parents[1] / ".conda" / "envs" / "aml-variant-demo"
BIN = PREFIX / "bin"
ENV = os.environ.copy()
ENV.update(
    {
        "PATH": f"{BIN}:{ENV.get('PATH', '')}",
        "CONDA_PREFIX": str(PREFIX),
        "JAVA_HOME": str(PREFIX / "lib" / "jvm"),
        "NXF_HOME": str(PREFIX.parents[2] / ".nextflow"),
        "CLOUDSDK_CONFIG": str(PREFIX.parents[2] / ".gcloud"),
    }
)


def run(name: str, argv: list[str], allowed: set[int] = {0}) -> bool:
    try:
        result = subprocess.run(
            argv,
            env=ENV,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=60,
            check=False,
        )
    except Exception as exc:
        print(f"FAIL\t{name}\t{type(exc).__name__}: {exc}")
        return False
    output = " | ".join(line.strip() for line in result.stdout.splitlines() if line.strip())
    output = output[:500]
    status = "PASS" if result.returncode in allowed else "FAIL"
    print(f"{status}\t{name}\texit={result.returncode}\t{output}")
    return status == "PASS"


def main() -> int:
    print(f"host\t{platform.platform()}\t{platform.machine()}")
    checks = [
        ("python", [str(BIN / "python"), "--version"], {0}),
        ("java", [str(PREFIX / "lib" / "jvm" / "bin" / "java"), "-version"], {0}),
        ("nextflow", [str(BIN / "nextflow"), "-version"], {0}),
        ("samtools", [str(BIN / "samtools"), "--version"], {0}),
        ("bcftools", [str(BIN / "bcftools"), "--version"], {0}),
        ("htslib/bgzip", [str(BIN / "bgzip"), "--version"], {0}),
        ("htslib/tabix", [str(BIN / "tabix"), "--version"], {0}),
        ("bedtools", [str(BIN / "bedtools"), "--version"], {0}),
        ("gatk", [str(BIN / "gatk"), "--version"], {0}),
        ("kallisto", [str(BIN / "kallisto"), "version"], {0}),
        # These tools use a nonzero status to display usage when no input is supplied.
        ("pizzly", [str(BIN / "pizzly"), "--help"], {0, 1}),
        ("getITD", [str(BIN / "getitd"), "--help"], {0}),
        ("SRA prefetch", [str(BIN / "prefetch"), "--version"], {0}),
        ("SRA fasterq-dump", [str(BIN / "fasterq-dump"), "--version"], {0}),
        ("gsutil", [str(BIN / "gsutil"), "version", "-l"], {0}),
        ("pandas", [str(BIN / "python"), "-c", "import pandas; print(pandas.__version__)"], {0}),
    ]
    ok = [run(*check) for check in checks]
    return 0 if all(ok) else 1


if __name__ == "__main__":
    raise SystemExit(main())
