import random
import sys
import argparse

# Define mutation rates
coding_region_mutation_rate = 0.0001338  # Mutation rate for the coding region
d_loop_mutation_rate = 0.00135  # Mutation rate for the D-loop region
titv = 1/22

def transition(b):
    if b == 'A':
        return 'G'
    if b == 'G':
        return 'A'
    if b == 'C':
        return 'T'
    if b == 'T':
        return 'C'
    sys.stderr.write("transition error.\n")
    sys.exit(1)

def transversion(b):
    if b in ['A', 'G']:
        return ['C', 'T'][random.randint(0, 1)]
    if b in ['C', 'T']:
        return ['A', 'G'][random.randint(0, 1)]
    sys.stderr.write("transversion error.\n")
    sys.exit(1)

def mutate(seq, mutation_rate, mut):
    newseq = ""
    for i in seq:
        pmut = random.random()
        if pmut < mutation_rate:
            mut += 1
            print(mut)
            ptitv = random.random()
            if ptitv < (float(titv) / (float(titv) + float(1.0))):  # is a transition
                newseq += transition(i)
            else:
                newseq += transversion(i)
        else:
            newseq += i
    return newseq

def reverse_time_simulation(genome, generations):
    coding_region_end = 13300
    coding_region_start = 14300

    # Simulate reverse mutations for each generation
    m = 0
    for g in range(1, generations + 1):
        sys.stderr.write(f"Generation #{g}\n")
        mutation_rate = d_loop_mutation_rate if coding_region_start <= g <= coding_region_end else coding_region_mutation_rate
        genome = mutate(genome, mutation_rate, m)
    return genome

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate reverse time genetic mutations.")
    parser.add_argument("reffile", help="Path to the reference FASTA file.")
    parser.add_argument("gen", type=int, help="Number of generations to simulate.")
    args = parser.parse_args()

    # Reading the reference genome from a specified file
    with open(args.reffile, "r") as filefa:
        header = filefa.readline().strip()
        genome = filefa.read().replace('\n', '').upper()

    # Reverse time simulation
    resulting_genome = reverse_time_simulation(genome, args.gen)

    # Output
    with open("out_reverse.fa", "w") as outfastafp:
        outfastafp.write(f">{header}_reverse\n")
        outfastafp.write(resulting_genome + "\n")

    sys.stderr.write("Done\n")
