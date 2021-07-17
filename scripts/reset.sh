BASE_DIR=/wynton/home/corces/allan/variantapp
REFERENCE_DIR=/wynton/home/corces/allan/variantapp/reference

PREDICTIONS_DIR=$BASE_DIR/predictions
SHAP_DIR=$BASE_DIR/shap
DATA_DIR=$BASE_DIR/data
MODEL_DIR=$BASE_DIR/models

INPUT_DATA=$DATA_DIR/input_data.json
CV_SPLITS=$DATA_DIR/splits.json

INPUT_BW=$DATA_DIR/tracks/C24_track.bw
PEAKS_F=$DATA_DIR/peaks/app.bed

CHROM_SIZES=$REFERENCE_DIR/hg38.chrom.sizes
REFERENCE_GENOME=$REFERENCE_DIR/hg38.genome.fa

rm -r /Users/alexlan/Desktop/Research/variantapp/static/images/app
mkdir -p /Users/alexlan/Desktop/Research/variantapp/static/images/app