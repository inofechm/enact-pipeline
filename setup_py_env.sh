#!/bin/bash
eval "$(conda shell.bash hook)"

set -e

PY_ENV_PATH=$1

# Create Python environment
if ! conda info --envs | grep -q "$PY_ENV_PATH"; then
    echo "Environment $PY_ENV_PATH does not exist. Creating..."
    conda create --prefix $PY_ENV_PATH python=3.9
    conda activate $PY_ENV_PATH
    pip install -r requirements.txt
fi
