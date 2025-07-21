from Bio import SeqIO

def trim_genomes(fasta_file, positions_file, output_file):
    positions = {}
    with open(positions_file, 'r') as pos_file:
        for line in pos_file:
            genome_name, start, end = line.strip().split()
            positions[genome_name] = (int(start), int(end))

    trimmed_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_name = record.id
        if genome_name in positions:
            start, end = positions[genome_name]
            trimmed_seq = record.seq[start-1:end]
            trimmed_record = record[start-1:end]
            trimmed_records.append(trimmed_record)

    SeqIO.write(trimmed_records, output_file, "fasta")

# Example usage
fasta_file = "linear_chromosomes.fasta"
positions_file = "trimming_position.txt"
output_file = "trimmed_chromosomes.fasta"
trim_genomes(fasta_file, positions_file, output_file)
