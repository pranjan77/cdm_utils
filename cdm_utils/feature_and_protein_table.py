import csv
import hashlib
import argparse
import os
import gzip

# Define SO terms mapping
so_terms = {
    "gene": "SO:0000704",
    "pseudogene": "SO:0000336",
    "ncRNA_gene": "SO:0001263",
    "mRNA": "SO:0000234",
    "CDS": "SO:0000316",
    "exon": "SO:0000147",
    "five_prime_UTR": "SO:0000204",
    "three_prime_UTR": "SO:0000205",
    "ncRNA": "SO:0000655",
    "rRNA": "SO:0000252",
    "tRNA": "SO:0000253",
    "SRP_RNA": "SO:0000590",
    "RNase_P_RNA": "SO:0000386",
    "riboswitch": "SO:0000035",
    "direct_repeat": "SO:0000319",
    "origin_of_replication": "SO:0000296",
    "CRISPR": "SO:0001459",
    "mobile_genetic_element": "SO:0001037",
    "region": "SO:0000001",
    "sequence_feature": "SO:0000110"
}

class GFFParser:
    def __init__(self, assembly_file, gff_file, protein_file):
        self.assembly_file = assembly_file
        self.gff_file = gff_file
        self.protein_file = protein_file
        self.features = []
        self.feature_associations = []
        self.feature_protein_associations = []
        self.assembly_md5 = None
        self.contig_md5s = {}
        self.protein_ids = set()

    @staticmethod
    def generate_file_md5(filepath, blocksize=65536):
        """Generate the MD5 checksum of a file's decompressed content."""
        md5 = hashlib.md5()
        open_func = gzip.open if filepath.endswith('.gz') else open
        try:
            with open_func(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
                for block in iter(lambda: f.read(blocksize), ''):
                    md5.update(block.encode('utf-8'))
            return md5.hexdigest()
        except Exception as e:
            print(f"Error generating MD5 for {filepath}: {e}")
            return None

    @staticmethod
    def generate_hash_id(seq_id, start, end, feature_type, file_md5, attribute_value=None):
        """Generate a hash-based ID for each feature, including filename, file MD5 checksum, and an attribute."""
        unique_string = f"{seq_id}_{start}_{end}_{feature_type}_{file_md5}"
        if attribute_value:
            unique_string += f"_{attribute_value}"
        return hashlib.md5(unique_string.encode()).hexdigest()

    @staticmethod
    def parse_attributes(attributes_str):
        """Parse the GFF3 attributes field into a dictionary."""
        attributes = {}
        for attribute in attributes_str.split(';'):
            if '=' in attribute:
                key, value = attribute.split('=', 1)
                attributes[key] = value
        return attributes

    def calculate_md5_checksums(self):
        """Calculate MD5 checksums for the assembly and its contigs."""
        print(f"Calculating MD5 for assembly: {self.assembly_file}")
        self.assembly_md5 = self.generate_file_md5(self.assembly_file)
        if not self.assembly_md5:
            print(f"Error calculating MD5 for assembly file {self.assembly_file}")
            return

        open_func = gzip.open if self.assembly_file.endswith('.gz') else open

        try:
            with open_func(self.assembly_file, 'rt', encoding='utf-8', errors='ignore') as file:
                current_contig = None
                contig_data = []

                for line in file:
                    if line.startswith('>'):
                        if current_contig and contig_data:
                            # Calculate MD5 for the previous contig
                            contig_md5 = hashlib.md5(''.join(contig_data).encode()).hexdigest()
                            self.contig_md5s[current_contig] = contig_md5
                            print(f"Calculated MD5 for contig {current_contig}: {contig_md5}")
                        current_contig = line[1:].strip().split()[0]  # Get the contig name without '>'
                        contig_data = []  # Reset contig data
                    else:
                        contig_data.append(line.strip())

                # Calculate MD5 for the last contig
                if current_contig and contig_data:
                    contig_md5 = hashlib.md5(''.join(contig_data).encode()).hexdigest()
                    self.contig_md5s[current_contig] = contig_md5
                    print(f"Calculated MD5 for contig {current_contig}: {contig_md5}")
        except Exception as e:
            print(f"Error reading assembly file {self.assembly_file}: {e}")

    def prepare_gff3_data(self):
        """Prepare data for insertion into the database."""
        print(f"Preparing GFF3 data from: {self.gff_file}")
        file_md5 = self.generate_file_md5(self.gff_file)
        if not file_md5:
            print(f"Error calculating MD5 for GFF file {self.gff_file}")
            return

        open_func = gzip.open if self.gff_file.endswith('.gz') else open

        try:
            with open_func(self.gff_file, 'rt', encoding='utf-8', errors='ignore') as file:
                reader = csv.reader(file, delimiter='\t')

                for row in reader:
                    if row[0].startswith('#') or len(row) < 9:
                        continue

                    seq_id = row[0]
                    source = row[1]
                    feature_type = row[2]
                    start = int(row[3])
                    end = int(row[4])
                    score = row[5] if row[5] != '.' else None
                    strand = row[6] if row[6] in ['+', '-'] else None
                    phase = row[7] if row[7] in ['0', '1', '2'] else None
                    attributes_str = row[8]

                    # Parse attributes
                    attributes = self.parse_attributes(attributes_str)
                    feature_id_value = attributes.get('ID', None)
                    parent_value = attributes.get('Parent', None)
                    protein_id = attributes.get('protein_id', None)

                    if feature_type == "CDS" and protein_id:
                        self.protein_ids.add(protein_id)

                    # Generate a unique hash ID for each feature
                    feature_id = self.generate_hash_id(seq_id, start, end, feature_type, file_md5, feature_id_value)

                    # Prepare feature data including MD5 of the assembly and contig
                    feature_data = {
                        'feature_uid': feature_id,
                        'seq_id': seq_id,
                        'feature_type': feature_type,
                        'feature_ontology': so_terms.get(feature_type, ""),
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'score': score,
                        'phase': phase,
                        'original_id': feature_id_value,
                        'parent': parent_value,
                        'assembly_md5': self.assembly_md5,
                        'contig_md5': self.contig_md5s.get(seq_id, ""),  # Use the contig MD5 if available
                        'protein_id': protein_id  # Add protein_id to feature data
                    }
                    self.features.append(feature_data)

                    # Add all other attributes to the associations
                    for key, value in attributes.items():
                        self.feature_associations.append({
                            'feature_id': feature_id,
                            'key': key,
                            'value': value
                        })

            print(f"Finished preparing GFF3 data. Total features: {len(self.features)}")
        except Exception as e:
            print(f"Error reading GFF file {self.gff_file}: {e}")

    def prepare_protein_associations(self):
        """Prepare protein associations data from the protein file."""
        print(f"Preparing protein associations from: {self.protein_file}")
        open_func = gzip.open if self.protein_file.endswith('.gz') else open

        try:
            with open_func(self.protein_file, 'rt', encoding='utf-8', errors='ignore') as file:
                current_protein_id = None
                protein_sequence = []

                for line in file:
                    if line.startswith('>'):
                        if current_protein_id and protein_sequence:
                            protein_md5 = hashlib.md5(''.join(protein_sequence).encode()).hexdigest()
                            if current_protein_id not in self.protein_ids:
                                raise ValueError(f"Protein ID {current_protein_id} in FAA file does not match any protein_id in GFF file.")
                            self.feature_protein_associations.append({
                                'protein_id': current_protein_id,
                                'protein_md5': protein_md5
                            })
                            print(f"Protein ID {current_protein_id} with MD5 {protein_md5}")
                        current_protein_id = line[1:].strip().split()[0]
                        protein_sequence = []
                    else:
                        protein_sequence.append(line.strip())

                # Capture the last protein sequence
                if current_protein_id and protein_sequence:
                    protein_md5 = hashlib.md5(''.join(protein_sequence).encode()).hexdigest()
                    if current_protein_id not in self.protein_ids:
                        raise ValueError(f"Protein ID {current_protein_id} in FAA file does not match any protein_id in GFF file.")
                    self.feature_protein_associations.append({
                        'protein_id': current_protein_id,
                        'protein_md5': protein_md5
                    })
                    print(f"Protein ID {current_protein_id} with MD5 {protein_md5}")

            print(f"Finished preparing protein associations. Total proteins: {len(self.feature_protein_associations)}")
        except Exception as e:
            print(f"Error reading protein file {self.protein_file}: {e}")

    def match_proteins_to_features(self):
        """Match proteins to features based on the protein ID in the GFF attributes."""
        matched_proteins = []
        print("Matching proteins to features...")
        for feature in self.features:
            protein_id = feature['protein_id']
            matching_proteins = [assoc for assoc in self.feature_protein_associations if assoc['protein_id'] == protein_id]
            if matching_proteins:
                for match in matching_proteins:
                    matched_proteins.append({
                        'feature_id': feature['feature_uid'],
                        'protein_id': match['protein_id'],
                        'protein_md5': match['protein_md5']
                    })
                print(f"Feature ID {feature['feature_uid']} matched with protein ID {protein_id}")
            else:
                print(f"No matching protein for feature ID {feature['feature_uid']}")

        self.feature_protein_associations = matched_proteins  # Update the list with matched data
        print(f"Finished matching proteins to features. Total matches: {len(self.feature_protein_associations)}")

    def save_as_tsv(self, features_tsv, associations_tsv, protein_associations_tsv):
        """Save data as TSV files."""
        try:
            # Save features data
            with open(features_tsv, 'w', newline='') as f_out:
                writer = csv.DictWriter(f_out, fieldnames=['feature_uid', 'seq_id', 'feature_type', 'feature_ontology', 'start', 'end', 'strand', 'score', 'phase', 'original_id', 'parent', 'assembly_md5', 'contig_md5', 'protein_id'], delimiter='\t')
                writer.writeheader()
                writer.writerows(self.features)
            print(f"Features saved to {features_tsv}")

            # Save feature associations data
            with open(associations_tsv, 'w', newline='') as f_out:
                writer = csv.DictWriter(f_out, fieldnames=['feature_id', 'key', 'value'], delimiter='\t')
                writer.writeheader()
                writer.writerows(self.feature_associations)
            print(f"Feature associations saved to {associations_tsv}")

            # Save feature-protein associations data
            with open(protein_associations_tsv, 'w', newline='') as f_out:
                writer = csv.DictWriter(f_out, fieldnames=['feature_id', 'protein_id', 'protein_md5'], delimiter='\t')
                writer.writeheader()
                writer.writerows(self.feature_protein_associations)
            print(f"Feature-protein associations saved to {protein_associations_tsv}")

        except Exception as e:
            print(f"Error saving TSV files: {e}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Parse GFF3, assembly, and protein files and generate TSV outputs.")
    parser.add_argument('input_tsv', type=str, help='Input TSV file containing file paths for assembly, GFF, and protein files')
    parser.add_argument('--delimiter', type=str, choices=['space', 'tab'], default='tab', help='Delimiter used in the input TSV file (space or tab)')
    parser.add_argument('--features_output', type=str, default='features.tsv', help='Output TSV file for features')
    parser.add_argument('--associations_output', type=str, default='feature_associations.tsv', help='Output TSV file for feature associations')
    parser.add_argument('--protein_associations_output', type=str, default='feature_protein_associations.tsv', help='Output TSV file for feature-protein associations')

    args = parser.parse_args()

    # Set the delimiter based on user input
    delimiter = ' ' if args.delimiter == 'space' else '\t'

    # Read the input TSV file to get file paths
    with open(args.input_tsv, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=delimiter)
        for row in reader:
            if len(row) < 3:
                print("Skipping row due to missing file paths.")
                continue  # Skip rows that don't have all three file paths
            assembly_file, gff_file, protein_file = row
            print(f"Processing: Assembly: {assembly_file}, GFF: {gff_file}, Protein: {protein_file}")
            parser = GFFParser(assembly_file, gff_file, protein_file)
            parser.calculate_md5_checksums()
            parser.prepare_gff3_data()
            parser.prepare_protein_associations()
            parser.match_proteins_to_features()
            parser.save_as_tsv(args.features_output, args.associations_output, args.protein_associations_output)

if __name__ == "__main__":
    main()

