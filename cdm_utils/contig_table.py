import hashlib
import gzip
from Bio import SeqIO

class ContigTable:
    def __init__(self, assembly_paths_file):
        self.assembly_paths_file = assembly_paths_file
        self.contig_stats = []

    def compute_md5(self, sequence):
        """
        Compute the MD5 checksum of the contig sequence.
        """
        md5_hash = hashlib.md5()
        md5_hash.update(sequence.encode('utf-8'))
        return md5_hash.hexdigest()

    def calculate_contig_stats(self, fasta_file, assembly_id):
        """
        Calculate statistics for each contig in the assembly file.
        Handles both compressed (.gz) and uncompressed files.
        """
        # Use gzip.open if the file is compressed, otherwise use open
        open_func = gzip.open if fasta_file.endswith('.gz') else open

        with open_func(fasta_file, 'rt') as handle:
            sequences = SeqIO.parse(handle, "fasta")
            for seq_record in sequences:
                contig_name = seq_record.id
                sequence = str(seq_record.seq).upper()
                length = len(sequence)
                gc_content = (sequence.count('G') + sequence.count('C')) / length if length > 0 else 0
                contig_id = self.compute_md5(sequence)  # Compute MD5 based on the sequence content

                self.contig_stats.append({
                    'id': contig_id,
                    'contig_name': contig_name,
                    'length': length,
                    'gc_content': gc_content,
                    'assembly_id': assembly_id,
                    'fasta_file': fasta_file
                })

    def process_assemblies(self):
        """
        Process each assembly file path from the assembly_paths_file, calculate contig stats,
        and store the results in contig_stats.
        """
        with open(self.assembly_paths_file, 'r') as f:
            assembly_paths = f.read().splitlines()

        for assembly_file in assembly_paths:
            try:
                # Compute MD5 checksum for assembly file content for consistent ID
                assembly_id = self.compute_md5_from_file(assembly_file)

                # Calculate the contig statistics
                self.calculate_contig_stats(assembly_file, assembly_id)
                print(f"Processed contigs for assembly file {assembly_file}.")
            except Exception as e:
                print(f"Error processing assembly file {assembly_file}: {e}")

    def compute_md5_from_file(self, assembly_file):
        """
        Compute the MD5 checksum of the entire assembly file content,
        ensuring the same result for both compressed and uncompressed versions.
        """
        md5_hash = hashlib.md5()
        open_func = gzip.open if assembly_file.endswith('.gz') else open

        with open_func(assembly_file, 'rt') as handle:
            for line in handle:
                md5_hash.update(line.encode('utf-8'))

        return md5_hash.hexdigest()

    def write_to_tsv(self, output_file):
        """
        Write the contig statistics to a TSV file for database loading.
        """
        with open(output_file, 'w') as tsvfile:
            # Write the header row
            header = ['id', 'contig_name', 'length', 'gc_content', 'assembly_id', 'fasta_file']
            tsvfile.write('\t'.join(header) + '\n')

            # Write the contig data
            for stat in self.contig_stats:
                row = [
                    stat['id'],
                    stat['contig_name'],
                    str(stat['length']),
                    f"{stat['gc_content']:.4f}",  # Format GC content to four decimal places
                    stat['assembly_id'],
                    stat['fasta_file']
                ]
                tsvfile.write('\t'.join(row) + '\n')


# Example usage
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python contig_table.py <assembly_paths_file>")
        sys.exit(1)

    assembly_paths_file = sys.argv[1]  # Read the input file path from command line argument
    contigs_output_tsv_path = "contigs_output.tsv"  # Define the output file path

    # Create an instance of ContigTable
    contig_table = ContigTable(assembly_paths_file)

    # Process the assemblies and calculate the contig statistics
    contig_table.process_assemblies()

    # Write results to TSV file
    contig_table.write_to_tsv(contigs_output_tsv_path)

    print(f"Contig statistics have been written to {contigs_output_tsv_path}.")

