import pandas as pd
import hashlib
import argparse

# Helper function to generate unique MD5 hash IDs
def generate_md5(*args):
    hash_input = ''.join(args)
    return hashlib.md5(hash_input.encode('utf-8')).hexdigest()

# Define the ObservationAndAssembly class
class ObservationAndAssembly:
    def __init__(self):
        self.assemblies = []
        self.observations = []

    def add_assembly(self, assembly_accession, assembly_name, assembly_level, sample_id, observation_id):
        assembly_entry = {
            'assembly_accession': assembly_accession,
            'assembly_name': assembly_name,
            'assembly_level': assembly_level,
            'sample_id': sample_id,
            'measurement_id': observation_id  # Renamed from measurement_id to observation_id
        }
        self.assemblies.append(assembly_entry)

    def add_observation(self, observation_id, protocol_id, sample_id, value, data_file):
        observation_entry = {
            'observation_id': observation_id,
            'protocol_id': protocol_id,
            'sample_id': sample_id,
            'value': value,
            'data_file': data_file
        }
        self.observations.append(observation_entry)

    def to_dataframe(self):
        assembly_df = pd.DataFrame(self.assemblies).drop_duplicates()
        observation_df = pd.DataFrame(self.observations)
        return assembly_df, observation_df

# Main function to handle command-line arguments and process files
def main(input_file_path, assembly_output_path, observation_output_path, protocol_output_path):
    # Load the TSV file
    df = pd.read_csv(input_file_path, sep='\t')

    # Initialize list to store protocol data
    protocol_data = []

    # Create an instance of ObservationAndAssembly
    observation_and_assembly = ObservationAndAssembly()

    # Iterate over the dataframe to fill in the assembly, observation, and protocol table data
    for index, row in df.iterrows():
        # Extract relevant assembly and sample information
        sample_id = row['id']
        assembly_accession = row['assembly_accession']
        assembly_name = row['assembly_name']
        assembly_level = row['assembly_level']

        protocol_id = row.get('protocol_id', None)  # Assuming protocol_id is a column in the TSV file
        protocol_name = row.get('protocol_name', None)  # Assuming protocol_name is a column in the TSV file
        protocol_description = row.get('protocol_description', None)  # Assuming protocol_description is a column in the TSV file
        data_file = row.get('data_file', None)  # Assuming data_file is a column in the TSV file

        # Create a unique observation_id
        observation_id = generate_md5(sample_id, assembly_accession)

        # Add assembly entry using the ObservationAndAssembly class
        observation_and_assembly.add_assembly(
            assembly_accession=assembly_accession,
            assembly_name=assembly_name,
            assembly_level=assembly_level,
            sample_id=sample_id,
            observation_id=observation_id
        )

        # Add observation entry using the ObservationAndAssembly class
        observation_and_assembly.add_observation(
            observation_id=observation_id,
            protocol_id=protocol_id,
            sample_id=sample_id,
            value=assembly_accession,  # Here, 'value' is set to assembly_accession; adjust based on your needs
            data_file=data_file
        )

        # Create protocol table entry if protocol_id is not None
        if protocol_id:
            protocol_entry = {
                'protocol_id': protocol_id,
                'protocol_name': protocol_name,
                'protocol_description': protocol_description
            }
            protocol_data.append(protocol_entry)

    # Convert lists to DataFrames
    assembly_df, observation_df = observation_and_assembly.to_dataframe()
    protocol_df = pd.DataFrame(protocol_data).drop_duplicates(subset=['protocol_id'])  # Remove duplicates by protocol_id

    # Save the DataFrames to TSV files
    assembly_df.to_csv(assembly_output_path, sep='\t', index=False)
    observation_df.to_csv(observation_output_path, sep='\t', index=False)
    protocol_df.to_csv(protocol_output_path, sep='\t', index=False)

    # Output the paths to the TSV files created
    print(f"Assembly data saved to: {assembly_output_path}")
    print(f"Observation data saved to: {observation_output_path}")
    print(f"Protocol data saved to: {protocol_output_path}")

if __name__ == "__main__":
    # Set up argument parser with default values
    parser = argparse.ArgumentParser(description="Process input and output file paths for observation and assembly data.")
    parser.add_argument(
        "input_file",
        help="Path to the input TSV file.",
        nargs="?",
        default="sample_details.tsv"  # Default input file path
    )
    parser.add_argument(
        "assembly_output",
        help="Path to the output TSV file for assembly data.",
        nargs="?",
        default="assembly.tsv"  # Default output file path for assembly data
    )
    parser.add_argument(
        "observation_output",
        help="Path to the output TSV file for observation data.",
        nargs="?",
        default="observation.tsv"  # Default output file path for observation data
    )
    parser.add_argument(
        "protocol_output",
        help="Path to the output TSV file for protocol data.",
        nargs="?",
        default="protocol.tsv"  # Default output file path for protocol data
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call main function with the parsed arguments
    main(args.input_file, args.assembly_output, args.observation_output, args.protocol_output)

