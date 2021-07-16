#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/corces/allan/variantapp/scripts/job_output
#$ -cwd
#$ -j y
#$ -l mem_free=10G
#$ -l h_rt=24:00:00

BASE_DIR=/wynton/home/corces/allan/variantapp
REFERENCE_DIR=/wynton/home/corces/allan/variantapp/reference

PREDICTIONS_DIR=$BASE_DIR/predictions
SHAP_DIR=$BASE_DIR/shap
DATA_DIR=$BASE_DIR/data
MODEL_DIR=$BASE_DIR/models

INPUT_DATA=$DATA_DIR/input_data.json
CV_SPLITS=$DATA_DIR/splits.json

INPUT_BW=$DATA_DIR/tracks/C24_track.bw
# PEAKS_F=$DATA_DIR/peaks/C24_peaks.bed
PEAKS_F=$DATA_DIR/peaks/app.bed

CHROM_SIZES=$REFERENCE_DIR/hg38.chrom.sizes
REFERENCE_GENOME=$REFERENCE_DIR/hg38.genome.fa

mkdir -p $PREDICTIONS_DIR

INPUT_SEQ_LEN=2114
OUTPUT_LEN=1000

predict \
    --model $MODEL_DIR/bpnet-hintv1/model.h5 \
    --chrom-sizes $CHROM_SIZES \
    --reference-genome $REFERENCE_GENOME \
    --exponentiate-counts \
    --output-dir $PREDICTIONS_DIR \
    --input-data $INPUT_DATA \
    --input-seq-len $INPUT_SEQ_LEN \
    --output-len $OUTPUT_LEN \
    --output-window-size $OUTPUT_LEN \
    --write-buffer-size 4000 \
    --batch-size 1  \
    --predict-peaks

LOGITS_FILE=$PREDICTIONS_DIR/'bpnet_task0.bw'
COUNTS_FILE=$PREDICTIONS_DIR/'bpnet_task0_exponentiated_counts.bw'

logits2profile \
        --logits-file $LOGITS_FILE \
        --counts-file $COUNTS_FILE \
        --output-directory $PREDICTIONS_DIR \
        --output-filename predicted_profile \
        --peaks $PEAKS_F \
        --chrom-sizes $CHROM_SIZES