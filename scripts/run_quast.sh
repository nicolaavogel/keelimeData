#!/usr/bin/env bash
set -eu

REF_DIR="/projects/wintherpedersen/people/bfj994/keelime/sim"  # Update if needed
THREADS=10

echo "[$(date)] Starting QUAST runs..."

for spades_dir in *_spades2; do
    contigs="$spades_dir/contigs.fasta"
    if [[ ! -f "$contigs" ]]; then
        echo "[$(date)] WARNING: $contigs not found, skipping $spades_dir" >&2
        continue
    fi

    base=${spades_dir%_spades}
    quast_out="quast_${base}"
    
    # Determine reference
    if [[ "$base" == PBHigh* || "$base" == PBNone* ]]; then
        ref="PB.fa"
    else
        # Extract prefix before "High"
        prefix="${base%%High*}"
        ref="${prefix}.fa"
    fi

    # Check if reference file exists
    if [[ ! -f "$REF_DIR/$ref" ]]; then
        echo "[$(date)] ERROR: Reference $ref not found for $spades_dir" >&2
        continue
    fi

    echo "[$(date)] Running QUAST on $spades_dir using $ref..."
    quast.py -o "$quast_out" -r "$REF_DIR/$ref" \
        --threads "$THREADS" --fragmented --min-identity 89.99 "$contigs"
done

echo "[$(date)] All QUAST jobs completed."

