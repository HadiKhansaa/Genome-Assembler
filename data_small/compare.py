def read_sequence(file_path):
    with open(file_path, 'r') as file:
        return file.read().replace('\n', '')

def compare_sequences(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    similarity = (matches / len(seq1)) * 100 if len(seq1) > 0 else 0
    return matches, mismatches, similarity

def main():
    contig_file = 'ddd.txt'  # Path to your contig file
    reference_genome_file = 'genome.txt'  # Path to your reference genome file

    contigs = read_sequence(contig_file)
    reference_genome = read_sequence(reference_genome_file)

    # Compare the first contig with the reference genome
    # Note: This is a simple comparison. In real cases, you might need to align sequences.
    matches, mismatches, similarity = compare_sequences(contigs, reference_genome)

    print(f"Matches: {matches}")
    print(f"Mismatches: {mismatches}")
    print(f"Similarity: {similarity:.2f}%")

if __name__ == "__main__":
    main()
