#!/usr/bin/env bash
# Source this file from the project root: source activate_env.sh
set -eu

if [ -n "${BASH_VERSION:-}" ]; then
    SCRIPT_PATH="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION:-}" ]; then
    SCRIPT_PATH="${(%):-%N}"
else
    SCRIPT_PATH="$0"
fi
PROJECT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
ENV_PREFIX="$PROJECT_DIR/.conda/envs/aml-variant-demo"

export CONDA_PREFIX="$ENV_PREFIX"
export CONDA_DEFAULT_ENV="aml-variant-demo"
export CONDA_SUBDIR="osx-64"
export CONDA_ENVS_PATH="$PROJECT_DIR/.conda/envs"
export CONDA_PKGS_DIRS="$PROJECT_DIR/.conda/pkgs_x86"
export JAVA_HOME="$ENV_PREFIX/lib/jvm"
export NXF_HOME="$PROJECT_DIR/.nextflow"
export CLOUDSDK_CONFIG="$PROJECT_DIR/.gcloud"
export PATH="$ENV_PREFIX/bin:$JAVA_HOME/bin:$PATH"

echo "Activated aml-variant-demo at $ENV_PREFIX"
