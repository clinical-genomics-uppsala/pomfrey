#!/bin/bash
module load snakemake/6.6.1
snakefile=$1
config=$2
outSvg=$3

snakemake --rulegraph  -s ${snakefile} --configfile ${config} | dot -Tsvg > ${outSvg}
