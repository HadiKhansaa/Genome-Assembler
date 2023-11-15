def fastq_to_text(input_fastq, output_text):
    with open(input_fastq, 'r') as fastq_file, open(output_text, 'w') as text_file:
        # Iterate through the FASTQ file
        for line_number, line in enumerate(fastq_file):
            # Check if the line is a DNA sequence line (skip every 4th line)
            if line_number % 4 == 1:
                # Write the DNA sequence to the text file
                text_file.write(line.strip() + '\n')

# Replace 'output_tiny_30xCov1.fq' and 'output.txt' with your actual file names
fastq_to_text('output_tiny_30xCov1.fq', 'output.txt')
