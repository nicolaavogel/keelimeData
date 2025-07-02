#!/bin/bash
#SBATCH --job-name=edit_distance
#SBATCH --output=logs/edit_distance_%A_%a.out
#SBATCH --error=logs/edit_distance_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-90  # Adjust based on number of MSA files - 1

# Define your list of MSA files
MSA_LIST=(
PB100kHigh_fin_25000normConsensus_prank.best.fas
PB100kHigh_fin_25000reckConsensus_prank.best.fas
PB100kHigh_fin_25000strictConsensus_prank.best.fas
PB100kHigh_fin_3000normConsensus_prank.best.fas
PB100kHigh_fin_3000reckConsensus_prank.best.fas
PB100kHigh_fin_3000strictConsensus_prank.best.fas
PB10kHigh_fin_25000normConsensus_prank.best.fas
PB10kHigh_fin_25000reckConsensus_prank.best.fas
PB10kHigh_fin_25000strictConsensus_prank.best.fas
PB10kHigh_fin_3000normConsensus_prank.best.fas
PB10kHigh_fin_3000reckConsensus_prank.best.fas
PB10kHigh_fin_3000strictConsensus_prank.best.fas
PB10myaHigh_fin_25000normConsensus_prank.best.fas
PB10myaHigh_fin_25000reckConsensus_prank.best.fas
PB10myaHigh_fin_25000strictConsensus_prank.best.fas
PB10myaHigh_fin_3000normConsensus_prank.best.fas
PB10myaHigh_fin_3000reckConsensus_prank.best.fas
PB10myaHigh_fin_3000strictConsensus_prank.best.fas
PB1myaHigh_fin_25000normConsensus_prank.best.fas
PB1myaHigh_fin_25000reckConsensus_prank.best.fas
PB1myaHigh_fin_25000strictConsensus_prank.best.fas
PB1myaHigh_fin_3000normConsensus_prank.best.fas
PB1myaHigh_fin_3000reckConsensus_prank.best.fas
PB1myaHigh_fin_3000strictConsensus_prank.best.fas
PB20myaHigh_fin_25000normConsensus_prank.best.fas
PB20myaHigh_fin_25000reckConsensus_prank.best.fas
PB20myaHigh_fin_25000strictConsensus_prank.best.fas
PB20myaHigh_fin_3000normConsensus_prank.best.fas
PB20myaHigh_fin_3000reckConsensus_prank.best.fas
PB20myaHigh_fin_3000strictConsensus_prank.best.fas
PB3myaHigh_fin_25000normConsensus_prank.best.fas
PB3myaHigh_fin_25000reckConsensus_prank.best.fas
PB3myaHigh_fin_25000strictConsensus_prank.best.fas
PB3myaHigh_fin_3000normConsensus_prank.best.fas
PB3myaHigh_fin_3000reckConsensus_prank.best.fas
PB3myaHigh_fin_3000strictConsensus_prank.best.fas
PB5myaHigh_fin_25000normConsensus_prank.best.fas
PB5myaHigh_fin_25000reckConsensus_prank.best.fas
PB5myaHigh_fin_25000strictConsensus_prank.best.fas
PB5myaHigh_fin_3000normConsensus_prank.best.fas
PB5myaHigh_fin_3000reckConsensus_prank.best.fas
PB5myaHigh_fin_3000strictConsensus_prank.best.fas
PB500kHigh_fin_25000normConsensus_prank.best.fas
PB500kHigh_fin_25000reckConsensus_prank.best.fas
PB500kHigh_fin_25000strictConsensus_prank.best.fas
PB500kHigh_fin_3000normConsensus_prank.best.fas
PB500kHigh_fin_3000reckConsensus_prank.best.fas
PB500kHigh_fin_3000strictConsensus_prank.best.fas)  # Add all file paths here

# Get current MSA file based on SLURM_ARRAY_TASK_ID
MSA_FILE=${MSA_LIST[$SLURM_ARRAY_TASK_ID]}
BASELINE_INDEX=0  # or set dynamically if needed
OUTPUT_FILE="results/$(basename "${MSA_FILE%.fas}")_new_edit_distance.txt"

# Make sure output directories exist
mkdir -p logs results

# Run your script
python editDist.py "$MSA_FILE" "$BASELINE_INDEX" "$OUTPUT_FILE"
