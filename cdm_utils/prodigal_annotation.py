import subprocess
import sys
import os
import uuid

class ProdigalAnnotation:
    def __init__(self, gz_file, prefix, output_dir):
        self.gz_file = gz_file
        self.prefix = prefix
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.decompressed_file = os.path.join(self.output_dir, f"{prefix}_decompressed.fna")
        self.gff_output = os.path.join(self.output_dir, f"{self.prefix}_prodigal.gff")
        self.faa_output = os.path.join(self.output_dir, f"{self.prefix}_prodigal.faa")

    def run_command(self, command):
        """Helper method to run a shell command."""
        try:
            subprocess.run(command, shell=True, check=True)
            print(f"Executed: {command}")
        except subprocess.CalledProcessError as e:
            print(f"Error executing command: {command}")
            sys.exit(1)

    def decompress_gz(self):
        """Decompress the gzipped file."""
        command = f"zcat {self.gz_file} > {self.decompressed_file}"
        self.run_command(command)

    def run_prodigal(self):
        """Run Prodigal with a specified prefix, outputting to a UUID directory."""
        command = f"prodigal -i {self.decompressed_file} -o {self.gff_output} -a {self.faa_output} -f gff"
        self.run_command(command)

        print(f"GFF output: {self.gff_output}")
        print(f"Protein FASTA output: {self.faa_output}")

    def clean_up(self):
        """Clean up all files in the output directory except the GFF and FAA files."""
        os.remove(self.decompressed_file)
        #print(f"Removed file: {self.decompressed_file}")

    def run(self):
        """Main method to execute the full workflow."""
        self.decompress_gz()
        self.run_prodigal()
        self.clean_up()
        print(f"Final outputs are stored in: {self.output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python prodigal_annotation.py <path_to_gz_file> <prefix> <output_dir>")
        sys.exit(1)

    gz_file = sys.argv[1]
    prefix = sys.argv[2]
    output_dir = sys.argv[3]

    annotation = ProdigalAnnotation(gz_file, prefix, output_dir)
    annotation.run()

