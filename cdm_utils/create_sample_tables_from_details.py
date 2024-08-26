import pandas as pd

class SampleTable:
    def __init__(self, input_file='sample_details.tsv'):
        self.input_file = input_file
        self.df = pd.read_csv(self.input_file, sep='\t')
        self.project_df = None
        self.sample_df = None
        self.isolate_df = None
        self.cultivation_df = None

    def create_project_table(self):
        """Create Project Table from the input dataframe."""
        project_columns = [
            'source_project_source', 'source_project_accession', 'source_project_title', 'source_submitter'
        ]
        self.project_df = self.df[project_columns].drop_duplicates().rename(columns={
            'source_project_source': 'project_source',
            'source_project_accession': 'project_accession',
            'source_project_title': 'project_title',
            'source_submitter': 'submitter'
        })

        # Generate project_id based on project_accession
        self.project_df['project_id'] = self.project_df['project_accession']

    def create_sample_table(self):
        """Create Sample Table from the input dataframe."""
        sample_columns = [
            'id', 'biosample_accession', 'biosample_collection_date', 'biosample_host',
            'geolocation_geo_loc_name', 'geolocation_latitude', 'geolocation_longitude', 'geolocation_elevation',
            'geolocation_depth', 'biosample_environment_package', 'biosample_models',
            'biosample_parent_accession', 'source_project_accession'
        ]
        self.sample_df = self.df[sample_columns].rename(columns={
            'id': 'sample_id',
            'biosample_accession': 'accession',
            'biosample_collection_date': 'collection_date',
            'biosample_host': 'host',
            'geolocation_geo_loc_name': 'geo_loc_name',
            'geolocation_latitude': 'latitude',
            'geolocation_longitude': 'longitude',
            'geolocation_elevation': 'elevation',
            'geolocation_depth': 'depth',
            'biosample_environment_package': 'environment_package',
            'biosample_models': 'models',
            'biosample_parent_accession': 'parent_accession',
            'source_project_accession': 'project_id'
        })

    def create_isolate_table(self):
        """Create Isolate Table from the input dataframe."""
        isolate_columns = ['id', 'biosample_isolate_strain']
        self.isolate_df = self.df[isolate_columns].rename(columns={
            'id': 'sample_id',
            'biosample_isolate_strain': 'strain'
        })

        # Add isolate_id as a unique identifier for each isolate entry, typically could be the same as sample_id or a new generated ID
        self.isolate_df['isolate_id'] = self.isolate_df['sample_id']  # Use sample_id as isolate_id here

    def create_cultivation_table(self):
        """Create Cultivation Table if 'biosample_cultivation' exists in the input dataframe."""
        if 'biosample_cultivation' in self.df.columns:
            cultivation_columns = ['id', 'biosample_cultivation']
            self.cultivation_df = self.df[cultivation_columns].rename(columns={
                'id': 'sample_id',
                'biosample_cultivation': 'cultivation_details'
            })
        else:
            self.cultivation_df = pd.DataFrame()  # Create an empty DataFrame if no cultivation data

    def save_tables(self, sample_output='sample.tsv', project_output='project.tsv',
                    isolate_output='isolate.tsv', cultivation_output='cultivation.tsv'):
        """Save the tables to TSV files."""
        if self.sample_df is not None:
            self.sample_df.to_csv(sample_output, sep='\t', index=False)
        if self.project_df is not None:
            self.project_df.to_csv(project_output, sep='\t', index=False)
        if self.isolate_df is not None:
            self.isolate_df.to_csv(isolate_output, sep='\t', index=False)
        if self.cultivation_df is not None and not self.cultivation_df.empty:
            self.cultivation_df.to_csv(cultivation_output, sep='\t', index=False)

    def process_and_save(self):
        """Run the complete process of creating tables and saving them to TSV files."""
        self.create_project_table()
        self.create_sample_table()
        self.create_isolate_table()
        self.create_cultivation_table()
        self.save_tables()


if __name__ == '__main__':
    # Define input and output file paths
    input_file_path = 'sample_details.tsv'
    sample_output_path = 'sample.tsv'
    project_output_path = 'project.tsv'
    isolate_output_path = 'isolate.tsv'
    cultivation_output_path = 'cultivation.tsv'

    # Initialize SampleTable with the input file
    sample_table = SampleTable(input_file=input_file_path)

    # Process data and save to output files
    sample_table.process_and_save()

    # Output the paths to the TSV files created
    print(f"Files saved to: {sample_output_path}, {project_output_path}, {isolate_output_path}, "
          f"{cultivation_output_path if not sample_table.cultivation_df.empty else 'No cultivation data found'}")

