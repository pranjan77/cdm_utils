import subprocess
import sys
import os
import uuid
import re

class ProdigalAnnotation:
    def __init__(self, assembly_file, prefix, output_dir):
        self.assembly_file = assembly_file
        self.prefix = prefix
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.decompressed_file = os.path.join(self.output_dir, f"{prefix}_decompressed.fna")
        self.gff_output = os.path.join(self.output_dir, f"{self.prefix}_original.gff")
        self.faa_output = os.path.join(self.output_dir, f"{self.prefix}_prodigal.faa")
        self.updated_gff_output = os.path.join(self.output_dir, f"{self.prefix}_prodigal.gff")

    def run_command(self, command):
        """Helper method to run a shell command."""
        try:
            subprocess.run(command, shell=True, check=True)
            print(f"Executed: {command}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing command: {command}")
            print(e)  # Print the error message for debugging
            sys.exit(1)

    def prepare_assembly_file(self):
        """Prepare the assembly file by decompressing if gzipped, or use the uncompressed file directly."""
        if not os.path.isfile(self.assembly_file):
            print(f"Error: The file {self.assembly_file} does not exist.")
            sys.exit(1)

        if self.assembly_file.endswith('.gz'):
            command = f"gzip -dc {self.assembly_file} > {self.decompressed_file}"
            self.run_command(command)
            # Use the decompressed file
            self.assembly_file = self.decompressed_file
        else:
            # Use the provided uncompressed file directly
            self.decompressed_file = self.assembly_file

    def run_prodigal(self):
        """Run Prodigal with a specified prefix, outputting to a UUID directory."""
        if not os.path.isfile(self.decompressed_file):
            print(f"Error: The decompressed file {self.decompressed_file} does not exist.")
            sys.exit(1)
        command = f"prodigal -i {self.decompressed_file} -o {self.gff_output} -a {self.faa_output} -f gff"
        self.run_command(command)

        print(f"GFF output: {self.gff_output}")
        print(f"Protein FASTA output: {self.faa_output}")

    def update_gff_ids(self):
        """Update the GFF file with gene entries and corresponding CDS entries with Parent attributes."""
        with open(self.gff_output, 'r') as infile, open(self.updated_gff_output, 'w') as outfile:
            cds_counter = 1  # Initialize counter for CDS

            for line in infile:
                if line.startswith("#"):
                    outfile.write(line)  # Write headers and comments as-is
                else:
                    fields = line.strip().split("\t")
                    seq_id = fields[0]  # Get the sequence ID
                    start = fields[3]
                    end = fields[4]
                    strand = fields[6]

                    # Extract the existing ID for CDS from the attributes
                    attributes = fields[-1]
                    match = re.search(r'ID=([^;]+)', attributes)
                    if match:
                        original_cds_id = match.group(1)

                        # Create a new CDS ID by prefixing with the sequence ID
                        new_cds_id = f"{seq_id}_{cds_counter}"

                        # Create a gene ID by appending '_gene' to the new CDS ID
                        gene_id = f"{new_cds_id}_gene"

                        # Create a new gene entry
                        gene_entry = "\t".join([
                            seq_id, "Prodigal_v2.6.3", "gene", start, end, ".", strand, ".",
                            f"ID={gene_id};Name={gene_id}"
                        ]) + "\n"
                        outfile.write(gene_entry)

                        # Construct the protein_id by appending "_prot" to the new CDS ID
                        protein_id = f"{new_cds_id}_prot"

                        # Update the CDS entry with the new ID, Parent attribute, and protein_id
                        attributes = re.sub(r'ID=[^;]+', f'ID={new_cds_id};Parent={gene_id};protein_id={protein_id}', attributes)
                        fields[-1] = attributes
                        outfile.write("\t".join(fields) + "\n")

                        # Increment the CDS counter for unique ID assignment
                        cds_counter += 1

        print(f"Updated GFF file with gene and CDS entries saved as {self.updated_gff_output}")

    def update_faa_file(self):
        """Update the FAA file to match the updated protein IDs in the GFF file."""
        modified_faa_lines = []
        with open(self.faa_output, 'r') as faa_file:
            for line in faa_file:
                if line.startswith('>'):
                    protein_id = line.split()[0][1:]  # Extract protein ID from fasta header
                    new_protein_id = f"{protein_id}_prot"  # Update protein ID with _prot suffix
                    modified_faa_lines.append(f">{new_protein_id} {line[len(protein_id)+1:]}\n")  # Replace old ID with new ID
                else:
                    modified_faa_lines.append(line)

        # Write modified FAA back to file
        with open(self.faa_output, 'w') as faa_file:
            faa_file.writelines(modified_faa_lines)
        print(f"Modified FAA file saved with updated protein IDs.")

    def clean_up(self):
        """Clean up the decompressed file if it was originally gzipped."""
        if self.assembly_file.endswith('.gz') and os.path.exists(self.decompressed_file):
            try:
                os.remove(self.decompressed_file)
                print(f"Removed file: {self.decompressed_file}")
            except OSError as e:
                print(f"Error removing file {self.decompressed_file}: {e}")

    def run(self):
        """Main method to execute the full workflow."""
        self.prepare_assembly_file()
        self.run_prodigal()
        self.update_gff_ids()
        self.update_faa_file()  # Update the FAA file to match the updated protein IDs
        self.clean_up()
        print(f"Final outputs are stored in: {self.output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python prodigal_annotation.py <path_to_assembly_file> <prefix> <output_dir>")
        sys.exit(1)

    assembly_file = sys.argv[1]
    prefix = sys.argv[2]
    output_dir = sys.argv[3]

    annotation = ProdigalAnnotation(assembly_file, prefix, output_dir)
    annotation.run()

