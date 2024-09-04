import json
import hashlib
import pandas as pd
import sys
import re
import argparse


class NCBIJSONLParser:
    def __init__(self, input_file='fastgenomics.jsonl', sample_details_file='sample_details.tsv',
                 sample_attributes_file='sample_attributes.tsv', source_details_file='source_details.tsv',observation_details_file='observation_details.tsv' ):
        self.input_file = input_file
        self.sample_details_file = sample_details_file
        self.sample_attributes_file = sample_attributes_file
        self.source_details_file = source_details_file
        self.observation_details_file = observation_details_file


        # Initialize lists to store extracted data
        self.sample_details_data = []
        self.sample_attributes_data = []
        self.source_details_data = []
        self.observation_details_data = []


    @staticmethod
    def generate_md5(*args):
        """Generate an MD5 hash for a set of input arguments."""
        hash_input = ''.join(args)
        return hashlib.md5(hash_input.encode('utf-8')).hexdigest()

    @staticmethod
    def safe_float_conversion(value, default=None):
        """Safely convert a value to a float."""
        try:
            return float(value)
        except (ValueError, TypeError):
            return default

    @staticmethod
    def parse_lat_lon(lat_lon):
        """Parse latitude and longitude from a string like '22.4932 N 113.8762 E'."""
        if not lat_lon or lat_lon.lower() in ['missing', 'not determined']:
            return None, None

        # Regular expression to match latitude and longitude values
        match = re.match(r"([+-]?\d+\.\d+)\s*([NS])\s*([+-]?\d+\.\d+)\s*([EW])", lat_lon)
        if not match:
            return None, None

        lat, lat_dir, lon, lon_dir = match.groups()
        lat = float(lat) * (1 if lat_dir == 'N' else -1)
        lon = float(lon) * (1 if lon_dir == 'E' else -1)
        return lat, lon

    def parse(self):
        """Parse the JSONL file and extract relevant information."""
        with open(self.input_file, 'r') as file:
            for line in file:
                record = json.loads(line.strip())

                # Extract Source, Geolocation, and Biosample Info
                project_accession = record.get('assemblyInfo', {}).get('bioprojectAccession', '')
                project_id = self.generate_md5(project_accession)

                project_title = next((bp.get('title', '') for lineage in record.get('assemblyInfo', {}).get('bioprojectLineage', []) for bp in lineage.get('bioprojects', [])), '')
                submitter = record.get('assemblyInfo', {}).get('submitter', 'Unknown')

                biosample = record.get('assemblyInfo', {}).get('biosample', {})
                #geo_loc_name = biosample.get('geoLocName', None) if biosample.get('geoLocName', '').strip().lower() not in ['missing', 'not determined'] else None
                
                # Parse latitude and longitude using the parse_lat_lon method
                #lat_lon = biosample.get('latLon', None)

                lat_lon = self.safe_float_conversion(next((attr.get('value') for attr in biosample.get('attributes', []) if attr.get('name', '').lower() == 'latLon'), None))
                latitude, longitude = self.parse_lat_lon(lat_lon)

                #elevation = self.safe_float_conversion(next((attr.get('value') for attr in biosample.get('attributes', []) if attr.get('name', '').lower() == 'elevation'), None))
                #depth = self.safe_float_conversion(next((attr.get('value') for attr in biosample.get('attributes', []) if attr.get('name', '').lower() == 'depth'), None))
                #collection_date = biosample.get('collectionDate', '')
                #host = biosample.get('host', 'Unknown')
                sample_accession = biosample.get('accession', '')
                sample_parent_accession = biosample.get('parent_accession', None)
                sample_id = self.generate_md5(sample_accession)

                # Extract environment package and model information
                environment_package = biosample.get('package', 'Unknown')
                models = '; '.join(biosample.get('models', []))

                # Assembly-related details
                assembly_accession = record.get('accession', '')
                assembly_name = record.get('assemblyInfo', {}).get('assemblyName', '')
                assembly_level = record.get('assemblyInfo', {}).get('assemblyLevel', '')

                # Generate a unique ID for each sample

                # Initialize additional fields with default values
                annotations = json.dumps(biosample.get('annotations', []))
                add_date = biosample.get('submissionDate')
                mod_date = biosample.get('lastUpdated')
                #emsl_biosample_identifiers = json.dumps(biosample.get('sampleIds', []))
                # Initialize all variables with None
                latitude = longitude = depth = elevation = env_broad_scale_id = None
                env_local_scale_id = env_medium_id = collection_date = ecosystem = None
                ecosystem_category = ecosystem_type = ecosystem_subtype = specific_ecosystem = None
                geo_loc_name = host = is_metagenomic= derived_from= None
    
    # Extract cross-reference data for sample_xref table
                xref = list()
                for sample_id_info in biosample.get('sampleIds', []):
                    db = sample_id_info.get('db', '').strip()
                    accession = sample_id_info.get('value', '')
                    if db and db.lower() != 'unknown' and accession:
                        xref.append(db + ":" + accession)

                alternate_identifiers = ', '.join(xref)
                



                # Extract attributes relevant to the biosample
                attributes = biosample.get('attributes', [])
                other_attributes = dict()
                for attribute in attributes:
                    name = attribute.get('name', '').lower()
                    value = attribute.get('value', '')
                    if name == 'lat_lon' and value.lower() != "missing":
                        latitude, longitude = self.parse_lat_lon(value)
                    elif name == 'depth':
                        depth = self.safe_float_conversion(value)
                    elif name == 'elevation':
                        elevation = self.safe_float_conversion(value)
                    elif name == 'env_broad_scale':
                        env_broad_scale_id = value if value.lower() != "missing" else None
                    elif name == 'env_local_scale':
                        env_local_scale_id = value if value.lower() != "missing" else None
                    elif name == 'env_medium':
                        env_medium_id = value if value.lower() != "missing" else None
                    elif name == 'collection_date':
                        collection_date = value if value.lower() != "missing" else None
                    elif name == 'ecosystem':
                        ecosystem = value if value.lower() != "missing" else None
                    elif name == 'ecosystem_category':
                        ecosystem_category = value if value.lower() != "missing" else None
                    elif name == 'ecosystem_type':
                        ecosystem_type = value if value.lower() != "missing" else None
                    elif name == 'ecosystem_subtype':
                        ecosystem_subtype = value if value.lower() != "missing" else None
                    elif name == 'specific_ecosystem':
                        specific_ecosystem = value if value.lower() != "missing" else None
                    elif name == 'geo_loc_name':
                        geo_loc_name = value if value.lower() != "missing" else None

                    elif name == 'geographic location (elevation)':
                        elevation = value if value.lower() != "missing" else None

                    elif name == 'geographic location (depth)':
                        depth = value if value.lower() != "missing" else None

                    elif name == 'geographic location (latitude)':
                        latitude = value if value.lower() != "missing" else None

                    elif name == 'geographic location (longitude)':
                        longitude = value if value.lower() != "missing" else None

                    elif name == 'geographic location (region and locality)':
                        geo_loc_name = value if value.lower() != "missing" else None

                    elif name == 'host':
                        host = value if value.lower() != "missing" else None

                    elif name == 'derived_from':
                        v = value if value.lower() != "missing" else None
                        derived_from = ", ".join(re.findall(r"SAMN\d{8}", v))


                    elif name == 'metagenomic':
                        is_metagenomic = value if value.lower() != "missing" else None


                    elif attribute.get('value') != 'missing':
                        sample_attributes_entry = {
                            'sample_id': sample_id,
                            'metadata_key': attribute.get('name', ''),
                            'metadata_key_ontology': '',
                            'metadata_value': attribute.get('value', ''),
                            'metadata_value_ontology': '',
                            'unit': '',
                            'unit_ontology': ''
                        }
                        self.sample_attributes_data.append(sample_attributes_entry)
                # Create a sample details entry
                sample_details_entry = {
                    'id': sample_id,
                    'name': biosample.get('description', {}).get('title', ''),
                    'description': biosample.get('description', {}).get('comment', ''),
                    'accession': sample_accession,
                    'source_id': project_id,
                    'alternate_identifiers': alternate_identifiers,
                    'annotations': annotations,
                    'add_date': add_date,
                    'mod_date': mod_date,
                    'collection_date': collection_date,
                    'depth': depth,
                    'env_broad_scale_id': env_broad_scale_id, 
                    'env_local_scale_id': env_local_scale_id, 
                    'env_medium_id': env_medium_id, 
                    'latitude': latitude, 
                    'longitude': longitude,
                    'study_id': project_accession,
                    'ecosystem': ecosystem,
                    'ecosystem_category': ecosystem_category,
                    'ecosystem_type': ecosystem_type,
                    'ecosystem_subtype': ecosystem_subtype,
                    'specific_ecosystem': specific_ecosystem,
                    'is_metagenomic': is_metagenomic,
                    'location': geo_loc_name,
                    'elevation': elevation,
                    'host': host,
                    'derived_from': derived_from,
                    'environment_package': environment_package,
                    'models': models,
                    'sample_parent_id': None  # Placeholder, adjust if parent info is available
                }
                self.sample_details_data.append(sample_details_entry)

                source_details_entry = {
                    'id': project_id,
                    'accession': project_accession,
                    'title': project_title,
                    'submitter': submitter
                }  
                self.source_details_data.append(source_details_entry) 

             

                observation_details_entry = {
                    'id': sample_id,
                    'assembly_accession': assembly_accession,
                    'assembly_name': assembly_name,
                    'assembly_level': assembly_level,
                }  
                self.observation_details_data.append(observation_details_entry) 

    def save_to_tsv(self):
        """Convert the extracted data to dataframes and save to TSV files."""
        pd.DataFrame(self.sample_details_data).to_csv(self.sample_details_file, sep='\t', index=False)
        pd.DataFrame(self.sample_attributes_data).to_csv(self.sample_attributes_file, sep='\t', index=False)
        pd.DataFrame(self.source_details_data).to_csv(self.source_details_file, sep='\t', index=False)
        pd.DataFrame(self.observation_details_data).to_csv(self.observation_details_file, sep='\t', index=False)
        print(f"Data saved to {self.sample_details_file}, {self.sample_attributes_file}, {self.source_details_file}, {self.observation_details_file}")

if __name__ == '__main__':
 
    # Initialize the parser
    parser = argparse.ArgumentParser(description='Process file paths for the NCBI JSONL parser.')

    # Add arguments
    parser.add_argument('--input_file_path', type=str, default='fastgenomics.jsonl', help='Path to the input JSONL file')
    parser.add_argument('--sample_details_path', type=str, default='sample_details.tsv', help='Path to the output sample details TSV file')
    parser.add_argument('--sample_attributes_path', type=str, default='sample_attributes.tsv', help='Path to the output sample attributes TSV file')
    parser.add_argument('--source_details_path', type=str, default='source_details.tsv', help='Path to the output source details TSV file')
    parser.add_argument('--observation_details_path', type=str, default='observation_details.tsv', help='Path to the output observation details TSV file')

    # Parse the arguments
    args = parser.parse_args()

    # Assign the parsed arguments to variables
    input_file_path = args.input_file_path
    sample_details_path = args.sample_details_path
    sample_attributes_path = args.sample_attributes_path
    source_details_path = args.source_details_path
    observation_details_path = args.observation_details_path

    # Create an instance of NCBIJSONLParser
    parser_instance = NCBIJSONLParser(input_file_path, sample_details_path, sample_attributes_path, source_details_path, observation_details_path)
    parser_instance.parse()
    parser_instance.save_to_tsv()
