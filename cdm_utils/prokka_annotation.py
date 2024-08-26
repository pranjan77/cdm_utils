import subprocess
import sys
import os
import uuid
import shutil

class ProkkaAnnotation:
    def __init__(self, assembly_file, prefix, output_dir):
        self.assembly_file = assembly_file
        self.prefix = prefix
        self.output_dir = output_dir
        self.temp_fna_dir = os.path.join(os.path.dirname(self.output_dir), str(uuid.uuid4()))  # Temporary directory for .fna file
        self.temp_prokka_dir = os.path.join(os.path.dirname(self.output_dir), str(uuid.uuid4()))  # Temporary directory for Prokka output
        self.temp_fna_path = os.path.join(self.temp_fna_dir, "temp.fna")  # Decompressed file path in temp_fna_dir

    def run_command(self, command):
        """Helper method to run a shell command and print output."""
        try:
            subprocess.run(command, shell=True, check=True, executable='/bin/bash')
            print(f"Executed: {command}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing command: {command}")
            print(e)
            sys.exit(1)

    def prepare_directories(self):
        """Check if directories exist and create necessary directories."""
        if os.path.exists(self.output_dir):
            print(f"Error: Output directory {self.output_dir} already exists. Please provide a non-existing directory.")
            sys.exit(1)
        os.makedirs(self.output_dir, exist_ok=False)  # Create output directory
        os.makedirs(self.temp_fna_dir, exist_ok=False)  # Create temporary directory for .fna file
        os.makedirs(self.temp_prokka_dir, exist_ok=False)  # Create temporary directory for Prokka output

    def prepare_assembly_file(self):
        """Prepare the assembly file for Prokka, handling both gzipped and uncompressed formats."""
        if not os.path.isfile(self.assembly_file):
            print(f"Error: The file {self.assembly_file} does not exist.")
            sys.exit(1)

        if self.assembly_file.endswith('.gz'):
            command = f"gzip -dc {self.assembly_file} > {self.temp_fna_path}"
            self.run_command(command)
            # Check if the decompressed file was created and is not empty
            if not os.path.exists(self.temp_fna_path) or os.stat(self.temp_fna_path).st_size == 0:
                print(f"Error: The file {self.temp_fna_path} was not created or is empty.")
                sys.exit(1)
        else:
            self.temp_fna_path = self.assembly_file  # Use the uncompressed file directly

    def run_prokka(self):
        """Run Prokka in a temporary directory and copy results to the output directory."""
        prokka_command = f"prokka --outdir {self.temp_prokka_dir} --prefix {self.prefix} --addgenes --force {self.temp_fna_path}"
        self.run_command(prokka_command)

        # Copy GFF and FAA files to output directory
        gff_file = os.path.join(self.temp_prokka_dir, f"{self.prefix}.gff")
        faa_file = os.path.join(self.temp_prokka_dir, f"{self.prefix}.faa")
        if os.path.exists(gff_file) and os.path.exists(faa_file):
            shutil.copy(gff_file, self.output_dir)
            shutil.copy(faa_file, self.output_dir)
            print(f"Copied {gff_file} and {faa_file} to {self.output_dir}.")
        else:
            print(f"Error: Prokka did not produce the expected output files.")

    def modify_gff_file(self):
        """Modify the GFF file to add protein IDs to CDS entries."""
        gff_file_path = os.path.join(self.output_dir, f"{self.prefix}.gff")

        if not os.path.exists(gff_file_path):
            print(f"Error: GFF file not found in output directory.")
            sys.exit(1)

        modified_gff_lines = []
        with open(gff_file_path, 'r') as gff_file:
            for line in gff_file:
                if line.startswith('#') or '\tCDS\t' not in line:
                    modified_gff_lines.append(line)
                    continue

                fields = line.strip().split('\t')
                attributes = fields[-1]
                cds_id = None
                if 'ID=' in attributes:
                    cds_id = attributes.split('ID=')[1].split(';')[0]

                if cds_id:
                    protein_id = f"{cds_id}_prot"  # Construct protein_id from CDS ID
                    attributes += f";protein_id={protein_id}"
                    fields[-1] = attributes
                    modified_gff_lines.append('\t'.join(fields) + '\n')
                else:
                    modified_gff_lines.append(line)

        # Write modified GFF back to file
        with open(gff_file_path, 'w') as gff_file:
            gff_file.writelines(modified_gff_lines)
        print(f"Modified GFF file saved with protein IDs added.")

    def modify_faa_file(self):
        """Modify the FAA file to match protein IDs with the updated GFF file."""
        faa_file_path = os.path.join(self.output_dir, f"{self.prefix}.faa")

        if not os.path.exists(faa_file_path):
            print(f"Error: FAA file not found in output directory.")
            sys.exit(1)

        modified_faa_lines = []
        with open(faa_file_path, 'r') as faa_file:
            for line in faa_file:
                if line.startswith('>'):
                    protein_id = line.split()[0][1:]  # Extract protein ID from fasta header
                    new_protein_id = f"{protein_id}_prot"  # Update protein ID with _prot suffix
                    modified_faa_lines.append(f">{new_protein_id} {line[len(protein_id)+1:]}")  # Replace old ID with new ID
                else:
                    modified_faa_lines.append(line)

        # Write modified FAA back to file
        with open(faa_file_path, 'w') as faa_file:
            faa_file.writelines(modified_faa_lines)
        print(f"Modified FAA file saved with updated protein IDs.")

    def clean_up(self):
        """Remove the temporary directories after Prokka runs."""
        if os.path.exists(self.temp_fna_dir) and self.assembly_file.endswith('.gz'):
            shutil.rmtree(self.temp_fna_dir)  # Remove the entire temporary .fna directory
            print(f"Removed temporary directory: {self.temp_fna_dir}")
        if os.path.exists(self.temp_prokka_dir):
            shutil.rmtree(self.temp_prokka_dir)  # Remove the entire temporary Prokka directory
            print(f"Removed temporary directory: {self.temp_prokka_dir}")

    def run(self):
        """Main method to execute the full workflow."""
        self.prepare_directories()
        self.prepare_assembly_file()
        self.run_prokka()
        self.modify_gff_file()  # Modify the GFF file to add protein IDs
        self.modify_faa_file()  # Modify the FAA file to match the updated protein IDs
        self.clean_up()
        print(f"All outputs are stored in: {self.output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python run_prokka.py <path_to_assembly_fasta_file> <prefix> <output_dir>")
        sys.exit(1)

    assembly_file = sys.argv[1]
    prefix = sys.argv[2]
    output_dir = sys.argv[3]

    annotation = ProkkaAnnotation(assembly_file, prefix, output_dir)
    annotation.run()

