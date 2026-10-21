#!/bin/bash

source /Users/grimacalandra/miniconda3/etc/profile.d/conda.sh

conda activate celltypist_py310

# Project specific variables ---------------------------------------------------

PROJECT="your_project_path"

INDATA="${PROJECT}/data/workflow/AmnioSC_ALL_optimised_ref.Preprocessed_merged.Counts.csv"
OUTDIR="${PROJECT}/celltypist/output_ALL_optimised_ref_pan_fetal"
mkdir ${OUTDIR}

# Run celltypist ---------------------------------------------------------------

echo "Running CellTypist..."
celltypist --indata ${INDATA} --outdir ${OUTDIR} --majority-voting --transpose-input --model Pan_Fetal_Human.pkl

echo "All done! :)"

conda deactivate