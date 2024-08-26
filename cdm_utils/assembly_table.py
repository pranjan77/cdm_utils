import hashlib
import sys
import gzip
from .bbmap_assembly_stats import BBMapAssemblyStats

class AssemblyTable:
    def __init__(self, assembly_paths_file):
        self.assembly_paths_file = assembly_paths_file
        self.assemblies = []  # Store assembly statistics here
        self.bbmap_parser = BBMapAssemblyStats()

    def compute_md5(self, assembly_file):
        """
        Compute the MD5 checksum of the assembly file content,
        ensuring the same result for both compressed and uncompressed versions.
        """
        md5_hash = hashlib.md5()

        # Check if the file is compressed (ends with .gz)
        open_func = gzip.open if assembly_file.endswith('.gz') else open

        with open_func(assembly_file, 'rt') as f:  # Open in text mode
            for line in f:
                md5_hash.update(line.encode('utf-8'))  # Update hash with line content

        return md5_hash.hexdigest()

    def run_bbmap_and_parse(self, assembly_file):
        """
        Run BBMap assembly stats and parse the output.
        """
        bbmap_output = self.bbmap_parser.run_bbmap_stats(assembly_file)
        if bbmap_output:
            # Parse the BBMap stats.sh output
            self.bbmap_parser.parse_bbmap_output(bbmap_output)
            # Get and return the parsed stats
            return self.bbmap_parser.get_stats()
        else:
            print(f"BBMap stats.sh output is empty or invalid for {assembly_file}.")
            return None

    def add_assembly(self, md5sum, assembly_file, parsed_data):
        """
        Add an assembly record with stats to the assemblies list.
        """
        if not parsed_data:
            print(f"No parsed data for {assembly_file}, skipping.")
            return

        # Extract relevant fields for the assembly table
        assembly_record = {
            'id': md5sum,
            'assembly_file': assembly_file,
            'A_content': parsed_data.get('A_content', None),
            'C_content': parsed_data.get('C_content', None),
            'G_content': parsed_data.get('G_content', None),
            'T_content': parsed_data.get('T_content', None),
            'N_content': parsed_data.get('N_content', None),
            'IUPAC_content': parsed_data.get('IUPAC_content', None),
            'Other_content': parsed_data.get('Other_content', None),
            'GC_content': parsed_data.get('GC_content', None),
            'GC_stdev': parsed_data.get('GC_stdev', None),
            'scaffold_total': parsed_data.get('scaffold_total', None),
            'contig_total': parsed_data.get('contig_total', None),
            'scaffold_sequence_total': parsed_data.get('scaffold_sequence_total', None),
            'contig_sequence_total': parsed_data.get('contig_sequence_total', None),
            'contig_gap_percentage': parsed_data.get('contig_gap_percentage', None),
            'scaffold_N50': parsed_data.get('scaffold_N50', None),
            'scaffold_L50': parsed_data.get('scaffold_L50', None),
            'contig_N50': parsed_data.get('contig_N50', None),
            'contig_L50': parsed_data.get('contig_L50', None),
            'scaffold_N90': parsed_data.get('scaffold_N90', None),
            'scaffold_L90': parsed_data.get('scaffold_L90', None),
            'max_scaffold_length': parsed_data.get('max_scaffold_length', None),
            'max_contig_length': parsed_data.get('max_contig_length', None),
            'large_scaffold_count_gt_50kb': parsed_data.get('large_scaffold_count_gt_50kb', None),
            'percent_genome_in_large_scaffolds_gt_50kb': parsed_data.get('percent_genome_in_large_scaffolds_gt_50kb', None)
        }

        self.assemblies.append(assembly_record)

    def write_to_tsv(self, output_file):
        """
        Write the assembly statistics to a TSV file for database loading.
        """
        with open(output_file, 'w') as tsvfile:
            # Write the header row
            header = [
                "id", "assembly_file", "A_content", "C_content", "G_content", "T_content", "N_content",
                "IUPAC_content", "Other_content", "GC_content", "GC_stdev",
                "scaffold_total", "contig_total", "scaffold_sequence_total",
                "contig_sequence_total", "contig_gap_percentage", "scaffold_N50",
                "scaffold_L50", "contig_N50", "contig_L50", "scaffold_N90",
                "scaffold_L90", "max_scaffold_length", "max_contig_length",
                "large_scaffold_count_gt_50kb", "percent_genome_in_large_scaffolds_gt_50kb"
            ]
            tsvfile.write('\t'.join(header) + '\n')

            # Write the assembly data
            for assembly in self.assemblies:
                row = [str(assembly[col] if assembly[col] is not None else '') for col in header]
                tsvfile.write('\t'.join(row) + '\n')

    def process_assemblies(self, output_file):
        """
        Process each assembly file path from the assembly_paths_file, run BBMap stats, and save the results to a TSV file.

        Parameters:
        - output_file: Path to the output TSV file.
        """
        with open(self.assembly_paths_file, 'r') as f:
            assembly_paths = f.read().splitlines()

        for assembly_file in assembly_paths:
            try:
                # Compute MD5 checksum
                md5sum = self.compute_md5(assembly_file)

                # Run BBMap and parse output
                parsed_data = self.run_bbmap_and_parse(assembly_file)

                if parsed_data:
                    # Add assembly data to the table
                    self.add_assembly(md5sum, assembly_file, parsed_data)
                    print(f"Processed and added data for assembly file {assembly_file}.")
                else:
                    print(f"Failed to process data for assembly file {assembly_file}.")
            except Exception as e:
                print(f"Error processing assembly file {assembly_file}: {e}")

        # Write results to TSV file
        self.write_to_tsv(output_file)
        print(f"TSV data successfully saved to {output_file}")


# Example usage with sys.argv[1] for input file
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python assembly_table.py <assembly_paths_file>")
        sys.exit(1)

    assembly_paths_file = sys.argv[1]  # Read the input file path from command line argument
    output_tsv_path = "assembly_output.tsv"  # Define the output file path

    assembly_table = AssemblyTable(assembly_paths_file)
    assembly_table.process_assemblies(output_tsv_path)

