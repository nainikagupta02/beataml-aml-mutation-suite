# Installation notes

Platform checked before installation: native macOS on Apple silicon (`arm64`).
Windows would require WSL2; that does not apply to this machine. Kallisto 0.48
and Pizzly 0.37.3 are absent from the `osx-arm64` Bioconda channel, while both
are available for `osx-64`. Rosetta 2 was tested successfully, so this project
uses an isolated `osx-64` environment on the Apple-silicon host.

`google-cloud-sdk` is the Conda package that supplies the requested `gsutil`
executable. PyPI currently has no distribution named `getitd`, so the requested
`pip install getitd` cannot succeed. Bioconda now provides getITD 1.5.17 as a
noarch package; that verified package is included in `environment.yml` and was
tested. ITDseek is therefore not needed. Pizzly is pinned against Kallisto 0.48
because newer Kallisto index formats are incompatible with Pizzly.

The environment is kept project-local to avoid modifying the existing Miniforge
installation, whose Conda Rattler plugin currently crashes on this macOS version.
Mamba itself works. Create it with:

```bash
export CONDA_ENVS_PATH="$PWD/.conda/envs"
export CONDA_PKGS_DIRS="$PWD/.conda/pkgs_x86"
export CONDA_SUBDIR=osx-64
mamba env create --strict-channel-priority --override-channels -f environment.yml
```

Then activate the tested project paths and state directories with:

```bash
source activate_env.sh
```

Run the complete smoke test at any time with:

```bash
python scripts/smoke_test.py
```

Exact installed versions and verification notes are in
`installation_report.txt`. The fully resolved package URLs are in
`environment.osx-64.lock.txt`.
