#!/bin/bash

while getopts c:l: flag
    do
        case "${flag}" in
                c) chrom=${OPTARG}
                        ;;
                l) loc=${OPTARG}
                         ;;
                *) echo "Invalid option: -$flag" ;;
        esac
    done

printf "chrom:"
echo $chrom
printf "loc:"
echo $loc
QUERY="${chrom}:${loc}-${loc}"
echo $QUERY
tabix https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bed.gz $QUERY > static/text/motifs.txt