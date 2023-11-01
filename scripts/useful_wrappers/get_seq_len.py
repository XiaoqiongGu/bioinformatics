import sys
from Bio import SeqIO


def fasta_to_csv(fasta_file, csv_file):
    with open(csv_file, 'w') as output:
        # Write the header to the CSV file
        output.write("Sequence Name,Sequence Length\n")

        # Parse the FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_name = record.id
            seq_length = len(record.seq)
            output.write(f"{seq_name},{seq_length}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py path_to_fasta_file path_to_output_csv")
        sys.exit(1)

    fasta_file_path = sys.argv[1]
    csv_file_path = sys.argv[2]

    fasta_to_csv(fasta_file_path, csv_file_path)
