import json
import hashlib
import pandas as pd
import sys
import re

class NCBIJSONLParser:
    def __init__(self, input_file='fastgenomics.jsonl', sample_details_file='sample_details.tsv',
                 sample_attributes_file='sample_attributes.tsv', sample_xref_file='sample_xref.tsv'):
        self.input_file = input_file
        self.sample_details_file = sample_details_file
        self.sample_attributes_file = sample_attributes_file
        self.sample_xref_file = sample_xref_file

        # Initialize lists to store extracted data
        self.sample_details_data = []
        self.sample_attributes_data = []
        self.sample_xref_data = []

    @staticmethod
    def generate_md5(*args):
        hash_input = ''.join(args)
        return hashlib.md5(hash_input.encode('utf-8')).hexdigest()

    @staticmethod
    def safe_float_conversion(value, default=None):
        try:
            return float(value)
        except (ValueError, TypeError):
            return default

    @staticmethod
    def parse_lat_lon(lat_lon):
        """Parse latitude and longitude from a string like '22.4932 N 113.8762 E'."""
        if not lat_lon:
            return None, None

        # Regular expression to match the latitude and longitude values
        match = re.match(r"([+-]?\d+\.\d+)\s*([NS])\s*([+-]?\d+\.\d+)\s*([EW])", lat_lon)
        if not match:
            return None, None

        lat, lat_dir, lon, lon_dir = match.groups()
        lat = float(lat) * (1 if lat_dir == 'N' else -1)
        lon = float(lon) * (1 if lon_dir == 'E' else -1)
        return lat, lon

    def parse(self):
        # Read the JSONL file and extract relevant information
        with open(self.input_file, 'r') as file:
            for line in file:
                record = json.loads(line.strip())

                # Extracting Source, Geolocation, and Biosample Info
                project_accession = record.get('assemblyInfo', {}).get('bioprojectAccession', '')
                project_title = next((bp.get('title', '') for lineage in record.get('assemblyInfo', {}).get('bioprojectLineage', []) for bp in lineage.get('bioprojects', [])), '')
                project_source = "fastgenomics"
                submitter = record.get('submitter', 'Unknown')

                biosample = record.get('assemblyInfo', {}).get('biosample', {})
                geo_loc_name = biosample.get('geoLocName', None) if biosample.get('geoLocName', '').strip().lower() not in ['missing', 'not determined'] else None
                lat_lon = biosample.get('latLon', None) if biosample.get('latLon', '').strip().lower() not in ['missing', 'not determined'] else None
                
                # Split latitude and longitude using the parse_lat_lon method
                latitude, longitude = self.parse_lat_lon(lat_lon)
                
                elevation = self.safe_float_conversion(next((attr.get('value') for attr in biosample.get('attributes', []) if attr.get('name', '').lower() == 'elevation'), None))
                depth = self.safe_float_conversion(next((attr.get('value') for attr in biosample.get('attributes', []) if attr.get('name', '').lower() == 'depth'), None))
                collection_date = biosample.get('collectionDate', '')
                host = biosample.get('host', 'Unknown')
                sample_accession = biosample.get('accession', '')
                sample_parent_accession = biosample.get('parent_accession', None)

                # Extracting strain attribute correctly from a list of attributes
                strain_attribute = next((attr.get('value') for attr in biosample.get('attributes', []) if attr.get('name', '').lower() == 'strain'), 'Unknown')

                # Extracting environment package and model information
                environment_package = biosample.get('package', 'Unknown')
                models = '; '.join(biosample.get('models', []))

                # Assembly-related details
                assembly_accession = record.get('accession', '')
                assembly_name = record.get('assemblyInfo', {}).get('assemblyName', '')
                assembly_level = record.get('assemblyInfo', {}).get('assemblyLevel', '')

                # Generate a unique ID for each sample
                sample_id = self.generate_md5(sample_accession)

                # Create a sample details entry
                sample_details_entry = {
                    'id': sample_id,
                    'source_project_source': project_source,
                    'source_project_accession': project_accession,
                    'source_project_title': project_title,
                    'source_submitter': submitter,
                    'geolocation_geo_loc_name': geo_loc_name,
                    'geolocation_latitude': latitude,
                    'geolocation_longitude': longitude,
                    'geolocation_elevation': elevation,
                    'geolocation_depth': depth,
                    'biosample_collection_date': collection_date,
                    'biosample_host': host,
                    'biosample_accession': sample_accession,
                    'biosample_isolate_strain': strain_attribute,
                    'biosample_parent_accession': sample_parent_accession,
                    'biosample_environment_package': environment_package,
                    'biosample_models': models,
                    'assembly_accession': assembly_accession,
                    'assembly_name': assembly_name,
                    'assembly_level': assembly_level,
                    'sample_parent_id': None  # Placeholder, adjust if parent info is available
                }
                self.sample_details_data.append(sample_details_entry)

                # Extract data for Sample Attributes table
                for attribute in biosample.get('attributes', []):
                    sample_attribute_id = self.generate_md5(sample_id, attribute.get('name', ''))
                    sample_attributes_entry = {
                        'id': sample_attribute_id,
                        'sample_id': sample_id,
                        'metadata_key': attribute.get('name', ''),
                        'metadata_key_ontology': '',
                        'metadata_value': attribute.get('value', ''),
                        'metadata_value_ontology': '',
                        'unit': '',
                        'unit_ontology': ''
                    }
                    self.sample_attributes_data.append(sample_attributes_entry)

                # Extract cross-reference data for sample_xref table
                for sample_id_info in biosample.get('sampleIds', []):
                    db = sample_id_info.get('db', '').strip()
                    accession = sample_id_info.get('value', '')
                    if db and db.lower() != 'unknown' and accession:
                        sample_xref_entry = {
                            'sample_id': sample_id,
                            'db': db,
                            'accession': accession
                        }
                        self.sample_xref_data.append(sample_xref_entry)

    def save_to_tsv(self):
        # Convert the extracted data to dataframes and save to TSV files
        pd.DataFrame(self.sample_details_data).to_csv(self.sample_details_file, sep='\t', index=False)
        pd.DataFrame(self.sample_attributes_data).to_csv(self.sample_attributes_file, sep='\t', index=False)
        pd.DataFrame(self.sample_xref_data).to_csv(self.sample_xref_file, sep='\t', index=False)
        print(f"Data saved to {self.sample_details_file}, {self.sample_attributes_file}, {self.sample_xref_file}")

if __name__ == '__main__':
    # Default paths
    input_file_path = 'fastgenomics.jsonl'
    sample_details_path = 'sample_details.tsv'
    sample_attributes_path = 'sample_attributes.tsv'
    sample_xref_path = 'sample_xref.tsv'

    # If user provides paths via command line, use those instead
    if len(sys.argv) > 1:
        input_file_path = sys.argv[1]
    if len(sys.argv) > 2:
        sample_details_path = sys.argv[2]
    if len(sys.argv) > 3:
        sample_attributes_path = sys.argv[3]
    if len(sys.argv) > 4:
        sample_xref_path = sys.argv[4]

    # Create an instance of NCBIJSONLParser
    parser = NCBIJSONLParser(input_file_path, sample_details_path, sample_attributes_path, sample_xref_path)
    parser.parse()
    parser.save_to_tsv()

