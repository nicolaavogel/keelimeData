#!/usr/bin/env bash
set -u  # remove -e so script doesn't exit on failure

# --- CONFIGURATION ---
INPUT_DIR="/projects/wintherpedersen/people/bfj994/keelime/data"
SPADES_CMD="spades.py"

# --- MAIN LOOP ---
echo "[$(date)] Starting SPAdes assemblies..."

for FILEPATH in "$INPUT_DIR"/*.fq.gz; do
  FILENAME=$(basename "$FILEPATH")
  BASENAME="${FILENAME%.fq.gz}"
  OUTPUT_DIR="${BASENAME}_spades2"

  echo "[$(date)] Running SPAdes on $FILENAME -> $OUTPUT_DIR"
  
  if ! $SPADES_CMD -s "$FILEPATH" -o "$OUTPUT_DIR" -k 21 ; then
    echo "[$(date)] ERROR: SPAdes failed on $FILENAME. Skipping." >&2
    continue
  fi
done

echo "[$(date)] All SPAdes jobs attempted."
