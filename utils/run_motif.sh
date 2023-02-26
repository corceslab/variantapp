#!/bin/bash

while getopts c:p:f: flag
    do
        case "${flag}" in
                c) chrom=${OPTARG}
                        ;;
                p) loc=${OPTARG}
                        ;;
                f) out=${OPTARG}
                        ;;
                *) echo "Invalid option: -$flag" ;;
        esac
    done

printf "chrom:"
echo $chrom
printf "input loc:"
echo $loc
printf "output loc:"
echo $out
QUERY="${chrom}:${loc}-${loc}"
OUTPUT="static/GIDBcache/motifs/${out}.txt"
echo $QUERY
tabix motifs/hg38.archetype_motifs.v1.0.bed.gz $QUERY > $OUTPUT