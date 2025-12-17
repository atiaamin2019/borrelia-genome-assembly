import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def create_blast_db(reference_file):
    print(f"Creating BLAST database for {reference_file}")
    db_command = [
        "makeblastdb",
        "-in", reference_file,
        "-dbtype", "nucl"
    ]
    subprocess.run(db_command, check=True)

def run_blast(query_file, db_file, output_file):
    print(f"Running BLAST: query={query_file}, db={db_file}, output={output_file}")
    blast_command = [
        "blastn",
        "-query", query_file,
        "-db", db_file,
        "-out", output_file,
        "-outfmt", "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"
    ]
    subprocess.run(blast_command, check=True)

def parse_blast_output(blast_output_file, is_small_contig1):
    print(f"Parsing BLAST output: {blast_output_file}")
    df = pd.read_csv(blast_output_file, header=None)
    df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sstrand"]

    if is_small_contig1:
        # For small_contig1 (Left_end_blast.csv)
        plus_strand = df[df["sstrand"] == "plus"]
        if plus_strand.empty:
            raise ValueError(f"No 'plus' strand alignment found in {blast_output_file}")

        candidate = plus_strand[plus_strand["sstart"] == 1]
        if candidate.empty:
            raise ValueError(f"No suitable alignment found in {blast_output_file} for small_contig1")

        overlap_start = candidate["qstart"].iloc[0]
        overlap_end = candidate["qend"].iloc[0]

    else:
        # For small_contig2 (Right_end_blast.csv), use the existing logic
        plus_strand = df[df["sstrand"] == "plus"]
        if plus_strand.empty:
            raise ValueError(f"No 'plus' strand alignment found in {blast_output_file}")

        candidate = plus_strand[plus_strand["send"] == plus_strand["slen"]]
        if candidate.empty:
            raise ValueError(f"No suitable alignment found where subject end equals subject length in {blast_output_file}")

        overlap_start = candidate["qstart"].iloc[0]
        overlap_end = candidate["qend"].iloc[0]

    print(f"Determined overlap start: {overlap_start}, end: {overlap_end}")

    return overlap_start, overlap_end

def merge_contigs(large_contig, small_contig_1, small_contig_2, output_dir):
    print(f"Merging contigs for: {large_contig}, {small_contig_1}, {small_contig_2}")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    left_blast_output = os.path.join(output_dir, "Left_end_blast.csv")
    right_blast_output = os.path.join(output_dir, "Right_end_blast.csv")

    # Create BLAST database for the large contig
    create_blast_db(large_contig)

    # Run BLAST for both small contigs against the large contig
    run_blast(small_contig_1, large_contig, left_blast_output)
    run_blast(small_contig_2, large_contig, right_blast_output)

    # Parse BLAST output to get overlap coordinates
    overlap_start_1, overlap_end_1 = parse_blast_output(left_blast_output, is_small_contig1=True)
    overlap_start_2, overlap_end_2 = parse_blast_output(right_blast_output, is_small_contig1=False)

    # Read sequences
    print("Reading sequences...")
    large_seq = SeqIO.read(large_contig, "fasta").seq
    small_seq_1 = SeqIO.read(small_contig_1, "fasta").seq
    small_seq_2 = SeqIO.read(small_contig_2, "fasta").seq

    print("Constructing merged sequence...")
    # Construct the merged sequence
    merged_seq = small_seq_1[:overlap_start_1] + large_seq + small_seq_2[overlap_end_2:]

    # Create a new SeqRecord
    merged_record = SeqIO.SeqRecord(Seq(merged_seq), id="merged_contig", description="Merged contig")

    # Write to file
    output_file = os.path.join(output_dir, "merged_contig.fasta")
    SeqIO.write(merged_record, output_file, "fasta")
    print(f"Written merged sequence to {output_file}")

if __name__ == "__main__":
    # Expecting command-line arguments: large_contig, small_contig_1, small_contig_2, output_dir
    if len(sys.argv) != 5:
        print("Usage: python merge_contigs_blast_modified.py <large_contig> <small_contig_1> <small_contig_2> <output_dir>")
        sys.exit(1)

    large_contig = sys.argv[1]
    small_contig_1 = sys.argv[2]
    small_contig_2 = sys.argv[3]
    output_dir = sys.argv[4]

    # Call merge_contigs function with provided arguments
    merge_contigs(large_contig, small_contig_1, small_contig_2, output_dir)
