import sys
from Bio import AlignIO

def calculate_edit_distance_from_msa(msa_file, baseline_index, output_file):
    # Read the MSA file
    alignment = AlignIO.read(msa_file, "fasta")

    baseline_seq = alignment[baseline_index]
    edit_distances = []

    # Initialize genome coverage counters
    total_positions = len(baseline_seq.seq.replace('-', ''))  # Exclude gaps in the baseline from total positions

    for i, record in enumerate(alignment):
        if i == baseline_index:
            continue

        aligned_seq1 = str(baseline_seq.seq)
        aligned_seq2 = str(record.seq)

        # Initialize counters
        insertions = 0
        correct_matches = 0

        # Count the number of valid insertions and correct matches
        for res1, res2 in zip(aligned_seq1, aligned_seq2):
            if res1 == '-' and res2 in 'ACGT':  # Insertions into baseline
                insertions += 1
            elif res1 in 'ACGT' and res1 == res2:  # Correct matches
                correct_matches += 1

        # Store the results for each sequence comparison
        edit_distances.append((record.id, insertions, correct_matches))

    # Prepare the output
    result = f"Baseline Sequence: {baseline_seq.id}\n"
    for record_id, insertions, correct_matches in edit_distances:
        result += f"\nComparison with {record_id}:\n"
        result += f"Insertions: {insertions}\n"
        result += f"Correct Matches: {correct_matches}\n"

    # Write results to the output file
    with open(output_file, 'w') as file:
        file.write(result)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <msa_file> <baseline_index> <output_file>")
        sys.exit(1)

    msa_file = sys.argv[1]
    baseline_index = int(sys.argv[2])
    output_file = sys.argv[3]

    calculate_edit_distance_from_msa(msa_file, baseline_index, output_file)
