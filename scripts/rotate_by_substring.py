import gzip
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to rotate a sequence to start with the given substring
def rotate_sequence(seq, substring):
    pos = seq.find(substring)
    if pos == -1:
        return None
    return seq[pos:] + seq[:pos]

# Function to reverse complement a sequence
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

# Function to find a suitable substring and rotate the first sequence
def find_and_rotate(seq_list):
    second_seq = seq_list[1]
    for start in range(len(second_seq) - 19):
        substring = second_seq[start:start + 20]
        print(f"Trying substring from second sequence: {substring}")
        
        rotated_seq = rotate_sequence(seq_list[0], substring)
        if rotated_seq:
            print(f"Found matching substring in first sequence: {substring}")
            # Rotate second sequence to start with the same substring
            rotated_second_seq = rotate_sequence(second_seq, substring)
            return rotated_seq, rotated_second_seq

        reverse_comp_substring = reverse_complement(substring)
        rotated_seq = rotate_sequence(seq_list[0], reverse_comp_substring)
        if rotated_seq:
            reverse_comp_second_seq = reverse_complement(second_seq)  # Reverse complement the entire second sequence
            print(f"Found reverse complemented substring in first sequence: {reverse_comp_substring}")
            print("Reverse complementing the entire second sequence and rotating both sequences...")
            print(f"Original second sequence: {second_seq}")
            print(f"Reverse complemented second sequence: {reverse_comp_second_seq}")
            # Rotate second sequence to start with the same reverse complemented substring
            rotated_second_seq = rotate_sequence(reverse_comp_second_seq, reverse_comp_substring)
            return rotated_seq, rotated_second_seq

    raise ValueError("No suitable substring found in the first sequence, even after trying reverse complements.")

# Main function
def main(fasta_file, output_file):
    seq_list = []
    headers = []

    # Read the sequences from the input FASTA file
    with gzip.open(fasta_file, "rt") as file:
        for record in SeqIO.parse(file, "fasta"):
            seq_list.append(str(record.seq))
            headers.append(record.id)

    # Ensure there are at least two sequences
    if len(seq_list) < 2:
        raise ValueError("The input file must contain at least two sequences.")

    try:
        rotated_seq, rotated_second_seq = find_and_rotate(seq_list)
        rotated_records = [
            SeqRecord(Seq(rotated_seq), id=headers[0], description=""),
            SeqRecord(Seq(rotated_second_seq), id=headers[1], description="")
        ]
    except ValueError as e:
        print(e)
        return

    # Add any additional sequences to the output
    for i in range(2, len(seq_list)):
        rotated_records.append(SeqRecord(Seq(seq_list[i]), id=headers[i], description=""))

    # Write the rotated sequences to the output FASTA file
    with gzip.open(output_file, "wt") as file:
        SeqIO.write(rotated_records, file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rotate the sequences in a FASTA file to start with a suitable 20 bases substring from the second sequence, or its reverse complement if necessary.")
    parser.add_argument('fasta_file', help="Input FASTA file (gzipped)")
    parser.add_argument('output_file', help="Output FASTA file (gzipped)")
    args = parser.parse_args()

    main(args.fasta_file, args.output_file)
