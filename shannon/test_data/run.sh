#!/bin/bash

BASE_DIR="test_data"

CONTEXTS=("CpG" "CHG" "CHH")
CHROMOSOMES=("1" "2" "3" "4" "5" "Mt" "Pt")

for CONTEXT in "${CONTEXTS[@]}"; do
  for CHR in "${CHROMOSOMES[@]}"; do
    OUT_CSV="${CONTEXT}_chr${CHR}.csv"
    echo "label,url" > "$OUT_CSV"
    
    for SAMPLE_PATH in "$BASE_DIR"/*_se; do
      SAMPLE_NAME=$(basename "$SAMPLE_PATH")
      FILE_PATH="$SAMPLE_PATH/${CONTEXT}_${SAMPLE_NAME}.bismark_chr_${CHR}.cov.gz"
      
      if [ -f "$FILE_PATH" ]; then
        echo "$SAMPLE_NAME,$FILE_PATH" >> "$OUT_CSV"
      fi
    done
  done
done
