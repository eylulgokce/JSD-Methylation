#!/bin/bash

for SAMPLE_DIR in *_se; do
  echo "Entering $SAMPLE_DIR"
  for FILE in "$SAMPLE_DIR"/*.cov.gz; do
    echo "Processing $FILE"

    gunzip -f "$FILE"

    RAW_FILE="${FILE%.gz}"
    SORTED_FILE="${RAW_FILE%.cov}.sorted.cov"
    sort -k1,1 -k2,2n "$RAW_FILE" > "$SORTED_FILE"

    bgzip -f "$SORTED_FILE"

    tabix -s 1 -b 2 -e 2 "${SORTED_FILE}.gz"
  done
done
