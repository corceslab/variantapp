#!/bin/bash

tabix https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.all_motifs.v1.0.bed.gz $1 > static/text/motifs.txt