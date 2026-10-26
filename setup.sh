#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

if ! command -v conda >/dev/null; then
    echo "ERROR: conda no está disponible. Instala Miniforge: https://conda-forge.org/download/" >&2
    exit 1
fi

# Per-user env and package cache: on shared conda installs, the common pkgs/ can hold
# folders other users own (Permission denied), and a shared envs/ would let one user's
# setup delete another's env.
export CONDA_PKGS_DIRS="$HOME/.conda/pkgs"
ENV_PREFIX="$HOME/.conda/envs/usher-env"
# Older conda exits non-zero if it can't fetch channel notices after a successful
# install, which aborts the script under set -e.
export CONDA_NUMBER_CHANNEL_NOTICES=0

# The lock pins exact builds so every server gets the same env; no solver step.
# Override the user's .condarc channels: the lock downloads from repo.prefix.dev, and
# conda's ToS and notices checks would otherwise query conda.anaconda.org, which some
# servers can't reach.
conda create -y -p "$ENV_PREFIX" --file envs/usher-env.linux-64.lock \
    --override-channels -c https://repo.prefix.dev/conda-forge -c https://repo.prefix.dev/bioconda
# taxoniumtools is only on PyPI; this lock pins it and all its dependencies by hash.
conda run -p "$ENV_PREFIX" python -m pip install --no-deps --require-hashes -r envs/taxonium-requirements.lock

echo
echo "Env installed in $ENV_PREFIX. For running the pipeline:"
echo "    conda activate $ENV_PREFIX"
echo "    python pipeline.py"
