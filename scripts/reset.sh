#!/bin/bash

BASE_DIR=/wynton/home/corces/allan/BPMVF/ATAC
REFERENCE_DIR=/wynton/home/corces/allan/BPMVF/reference
DEPLOY_DIR=/wynton/home/corces/allan/BPMVF/deployment

DATA_DIR=$BASE_DIR/data
MODEL_DIR=$BASE_DIR/models
CHROM_SIZES=$REFERENCE_DIR/hg38.chrom.sizes
REFERENCE_GENOME=$REFERENCE_DIR/hg38.genome.fa
CV_SPLITS=$BASE_DIR/splits.json
INPUT_DATA=$BASE_DIR/input_data.json
PREDICTIONS_DIR=$BASE_DIR/predictions
BOUNDS_DIR=$BASE_DIR/bounds
METRICS_DIR=$BASE_DIR/metrics
SHAP_DIR=$BASE_DIR/shap

rm -r $BOUNDS_DIR
rm -r $METRICS_DIR
rm -r $PREDICTIONS_DIR

mkdir $BASE_DIR/predictions
mkdir $BASE_DIR/bounds
mkdir $BASE_DIR/metrics