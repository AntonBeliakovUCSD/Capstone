#!/bin/bash

# Check if the correct number of arguments was passed
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <sumstats_path> <weights_path> <weights_dir_path> <name> <output_directory>"
  exit 1
fi

# Assign command-line arguments to variables
SUMSTATS=$1
WEIGHTS=$2
WEIGHTS_DIR=$3
NAME=$4
OUTPUT_DIR=$5

# Ensure the output directory exists
mkdir -p $OUTPUT_DIR

# Loop through chromosomes 1 to 22
for CHR in {1..22}; do
    OUT="${OUTPUT_DIR}/${NAME}.${CHR}.dat"
    
    Rscript src/FUSION.assoc_test.R \
    --sumstats $SUMSTATS \
    --weights $WEIGHTS \
    --weights_dir $WEIGHTS_DIR \
    --ref_ld_chr data/LDREF/1000G.EUR. \
    --chr $CHR \
    --out $OUT
    
    echo "Processed chromosome $CHR"
done
