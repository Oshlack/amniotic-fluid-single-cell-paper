#!/bin/bash

source /Users/grimacalandra/miniconda3/etc/profile.d/conda.sh

conda activate celltypist_py310

# Project specific variables ---------------------------------------------------

PROJECT="your_project_path"

INDATA="${PROJECT}/data/Gerli_et_al_2024/Gerli_et_al_2024.Counts.csv"
OUTDIR="${PROJECT}/celltypist/output_Gerli_et_al_2024_pan_fetal"
mkdir ${OUTDIR}

# Run celltypist ---------------------------------------------------------------

echo "Running CellTypist..."
celltypist --indata ${INDATA} --outdir ${OUTDIR} --majority-voting --transpose-input --model Pan_Fetal_Human.pkl

echo "All done! :)"

conda deactivate