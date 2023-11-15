def parse_sam_file(sam_file):
    with open(sam_file, 'r') as file:
        reads = []
        for line in file:
            if line.startswith('@'):
                continue
            parts = line.split('\t')
            position = int(parts[3])
            sequence = parts[9]
            reads.append((position, sequence))
        return sorted(reads, key=lambda x: x[0])

def merge_reads_and_identify_contigs(reads, gap_threshold=1000):
    contigs = []
    current_contig = ''
    current_position = 1
    for position, sequence in reads:
        if position - current_position > gap_threshold:
            if current_contig:
                contigs.append(current_contig)
            current_contig = sequence
            current_position = position + len(sequence)
        else:
            overlap = current_position - position
            current_contig += sequence[overlap:]
            current_position += len(sequence) - overlap
    if current_contig:
        contigs.append(current_contig)
    return contigs

def assemble_contigs_from_sam(sam_file, output_file, gap_threshold=1000):
    reads = parse_sam_file(sam_file)
    contigs = merge_reads_and_identify_contigs(reads, gap_threshold)
    with open(output_file, 'w') as file:
        for contig in contigs:
            file.write(contig + '\n')


# Use the function with your SAM file
assemble_contigs_from_sam('our_reads.sam', 'ddd.txt')
