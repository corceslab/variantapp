#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/corces/allan/variantapp
#$ -cwd
#$ -j y
#$ -l mem_free=10G
#$ -l h_rt=24:00:00

python3 predict.py