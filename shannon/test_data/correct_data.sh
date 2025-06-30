#!/bin/bash

CONTEXTS=("CpG" "CHG" "CHH")
CHROMOSOMES=("1" "2" "3" "4" "5" "Mt" "Pt")

for CONTEXT in "${CONTEXTS[@]}"; do
  for CHR in "${CHROMOSOMES[@]}"; do
    CSV="${CONTEXT}_chr${CHR}.csv"
    echo "label,url" > "$CSV"

    for SAMPLE_DIR in *_se; do
      SAMPLE_NAME="$SAMPLE_DIR"
      FILE_PATH="${SAMPLE_DIR}/${CONTEXT}_${SAMPLE_NAME}.bismark_chr_${CHR}.sorted.cov.gz"
      
      if [[ -f "$FILE_PATH" ]]; then
        echo "${SAMPLE_NAME},test_data/${FILE_PATH}" >> "$CSV"
      fi
    done
  done
done
