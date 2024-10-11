#!/bin/bash
eval "$(conda shell.bash hook)"

set -e

PY_ENV_PATH=$1

conda activate $PY_ENV_PATH

FILE_URL="https://zenodo.org/records/13887921/files/ENACT_supporting_files.zip"
OUTPUT_FILE="ENACT_supporting_files.zip"

# Download ENACT supporting files if they are not present
if [ -f "$OUTPUT_FILE" ]; then
    echo "$OUTPUT_FILE already exists. Skipping download."
else
    echo "$OUTPUT_FILE is downloading."
    wget -O $OUTPUT_FILE $FILE_URL
    unzip $OUTPUT_FILE
fi


# Need to add step to download files from Zenodo to ENACT_supporting_files (in repo home directory)
# Run ENACT pipeline to test all combinations of Bin-to-cell assignment and cell annotation methods - Order of experiments matters, don't change!
python -m src.enact.pipeline --configs_path ENACT_supporting_files/public_data/human_colorectal/config_files/naive-celltypist.yaml
python -m src.enact.pipeline --configs_path ENACT_supporting_files/public_data/human_colorectal/config_files/naive-cellassign.yaml
python -m src.enact.pipeline --configs_path ENACT_supporting_files/public_data/human_colorectal/config_files/weighted_by_area-celltypist.yaml
python -m src.enact.pipeline --configs_path ENACT_supporting_files/public_data/human_colorectal/config_files/weighted_by_area-cellassign.yaml
python -m src.enact.pipeline --configs_path ENACT_supporting_files/public_data/human_colorectal/config_files/weighted_by_transcript-celltypist.yaml
python -m src.enact.pipeline --configs_path ENACT_supporting_files/public_data/human_colorectal/config_files/weighted_by_transcript-cellassign.yaml
python -m src.enact.pipeline --configs_path ENACT_supporting_files/public_data/human_colorectal/config_files/weighted_by_cluster-celltypist.yaml
python -m src.enact.pipeline --configs_path ENACT_supporting_files/public_data/human_colorectal/config_files/weighted_by_cluster-cellassign.yaml
