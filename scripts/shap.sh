#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/corces/allan/variantapp/scripts/job_output
#$ -cwd
#$ -j y
#$ -l mem_free=10G
#$ -l h_rt=24:00:00

BASE_DIR=/wynton/home/corces/allan/variantapp/ATAC
REFERENCE_DIR=/wynton/home/corces/allan/variantapp/reference
PREDICTIONS_DIR=$BASE_DIR/predictions
DATA_DIR=$BASE_DIR/data
MODEL_DIR=$BASE_DIR/model

INPUT_DATA=$BASE_DIR/data/input_data.json
CV_SPLITS=$BASE_DIR/data/splits.json

INPUT_BW=$DATA_DIR/Cluster24.bpnet.unstranded.bw
PEAKS_F=$DATA_DIR/Cluster24.peaks.bed

CHROM_SIZES=$REFERENCE_DIR/hg38.chrom.sizes
REFERENCE_GENOME=$REFERENCE_DIR/hg38.genome.fa

SHAP_DIR=$BASE_DIR/shap
mkdir -p $SHAP_DIR

shap_scores \
    --reference-genome $REFERENCE_GENOME \
    --model $MODEL_DIR/bpnet.256.15.001_split000.h5 \
    --bed-file $PEAKS_F \
    --chroms chr1 \
    --output-dir $SHAP_DIR \
    --input-seq-len 2114 \
    --control-len 1000 \
    --task-id 0 \
    --control-info $INPUT_DATA

peaks_valid_scores=$SHAP_DIR/peaks_valid_scores.bed
counts_scores_h5=$SHAP_DIR/counts_scores.h5
counts_out_pfx=$SHAP_DIR/count_scores

python hdf5_to_bw.py \
    -h5 $counts_scores_h5 \
    -r $peaks_valid_scores \
    -c $CHROM_SIZES \
    -o $counts_out_pfx.bw \
    -s $counts_out_pfx.stats.txt \
    -t 1

profile_scores_h5=$SHAP_DIR/profile_scores.h5
profile_out_pfx=$SHAP_DIR/profile_scores

python hdf5_to_bw.py \
    -h5 $profile_scores_h5 \
    -r $peaks_valid_scores \
    -c $CHROM_SIZES \
    -o $profile_out_pfx.bw \
    -s $profile_out_pfx.stats.txt \
    -t 1