#!/bin/bash
eval "$(conda shell.bash hook)"

set -e

PY_ENV_PATH=$1

# Run ENACT pipeline
conda activate $PY_ENV_PATH
python -m src.eval.cell_annotation_eval

