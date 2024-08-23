import subprocess
import configparser
import sys
import os
import tempfile
import shutil
class BBMapAssemblyStats:
    def __init__(self, config_file='config.ini'):
        # Initialize an empty dictionary to store the parsed data

        self.stats = {}
        stats_path = shutil.which('stats.sh')
        if stats_path:
            self.stats_path = stats_path
            print(f"Path to stats.sh: {stats_path}")
        else:
            print("stats.sh not found in PATH")
           # Get the path to the stats.sh binary

    def run_bbmap_stats(self, assembly_file):
        """
        Run BBMap's stats.sh on the provided assembly file and return the output.
        """
        try:
            # Create a temporary file to store the output
            args = [self.stats_sh_path, f'in={assembly_file}']
            # Run the stats.sh command and redirect output to the temp file
            result = subprocess.run(args,
                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
            )

           # Check the output from both stdout and stderr
            stdout_output = result.stdout.decode('utf-8') if result.stdout else "No stdout output"
           #stderr_output = result.stderr.decode('utf-8') if result.stderr else "No stderr output"

            return stdout_output

        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running BBMap stats.sh: {e}")
            return None



    def parse_bbmap_output(self, output):
        lines = output.splitlines()

        # Remove empty lines and trim each line
        lines = [line.strip() for line in lines if line.strip()]

        for i, line in enumerate(lines):
            if line.startswith("A\tC\tG\tT\tN\tIUPAC\tOther\tGC\tGC_stdev"):
                # The next line contains the values
                values = lines[i + 1].split("\t")
                self.stats['A_content'] = float(values[0])
                self.stats['C_content'] = float(values[1])
                self.stats['G_content'] = float(values[2])
                self.stats['T_content'] = float(values[3])
                self.stats['N_content'] = float(values[4])
                self.stats['IUPAC_content'] = float(values[5])
                self.stats['Other_content'] = float(values[6])
                self.stats['GC_content'] = float(values[7])
                self.stats['GC_stdev'] = float(values[8])

            elif line.startswith("Main genome scaffold total:"):
                self.stats['scaffold_total'] = int(line.split("\t")[1])
            elif line.startswith("Main genome contig total:"):
                self.stats['contig_total'] = int(line.split("\t")[1])
            elif line.startswith("Main genome scaffold sequence total:"):
                self.stats['scaffold_sequence_total'] = line.split("\t")[1]
            elif line.startswith("Main genome contig sequence total:"):
                parts = line.split("\t")
                self.stats['contig_sequence_total'] = parts[1]
                self.stats['contig_gap_percentage'] = parts[2].strip()
            elif line.startswith("Main genome scaffold N/L50:"):
                nl50 = line.split("\t")[1].split("/")
                self.stats['scaffold_N50'] = f"{nl50[0]} {nl50[1].split()[1]}"
                self.stats['scaffold_L50'] = nl50[1]
            elif line.startswith("Main genome contig N/L50:"):
                nl50 = line.split("\t")[1].split("/")
                self.stats['contig_N50'] = f"{nl50[0]} {nl50[1].split()[1]}"
                self.stats['contig_L50'] = nl50[1]
            elif line.startswith("Main genome scaffold N/L90:"):
                nl90 = line.split("\t")[1].split("/")
                self.stats['scaffold_N90'] = f"{nl90[0]} {nl90[1].split()[1]}"
                self.stats['scaffold_L90'] = nl90[1]
            elif line.startswith("Main genome contig N/L90:"):
                nl90 = line.split("\t")[1].split("/")
                self.stats['contig_N90'] = f"{nl90[0]} {nl90[1].split()[1]}"
                self.stats['contig_L90'] = nl90[1]
            elif line.startswith("Max scaffold length:"):
                self.stats['max_scaffold_length'] = line.split("\t")[1]
            elif line.startswith("Max contig length:"):
                self.stats['max_contig_length'] = line.split("\t")[1]
            elif line.startswith("Number of scaffolds > 50 KB:"):
                self.stats['large_scaffold_count_gt_50kb'] = int(line.split("\t")[1])
            elif line.startswith("% main genome in scaffolds > 50 KB:"):
                self.stats['percent_genome_in_large_scaffolds_gt_50kb'] = float(line.split("\t")[1].strip('%'))







    def get_stats(self):
        # Return the parsed genome stats as a dictionary
        return self.stats

if __name__ == "__main__":
    # Example usage: running BBMap stats and parsing the output
    assembly_file = sys.argv[1]  
    bbmap_parser = BBMapAssemblyStats(config_file="../config.ini")
    
    # Run BBMap stats.sh
    bbmap_output = bbmap_parser.run_bbmap_stats(assembly_file)
    
    if bbmap_output:
        # Parse the BBMap stats.sh output
        bbmap_parser.parse_bbmap_output(bbmap_output)
        
        # Get and print the parsed stats
        parsed_data = bbmap_parser.get_stats()
        print(parsed_data)

