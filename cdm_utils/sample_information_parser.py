import json
import argparse

def parse_metadata(json_data):
    """
    Parses relevant metadata from the provided JSON data.
    
    Args:
        json_data (dict): A dictionary containing the JSON data.
    
    Returns:
        dict: A dictionary with the extracted metadata.
    """
    metadata = {}

    # Accession
    metadata['accession'] = json_data.get('accession')

    # Annotation Info
    annotation_info = json_data.get('annotationInfo', {})
    metadata['annotation_method'] = annotation_info.get('method')
    metadata['annotation_pipeline'] = annotation_info.get('pipeline')
    metadata['annotation_provider'] = annotation_info.get('provider')
    metadata['annotation_release_date'] = annotation_info.get('releaseDate')
    metadata['annotation_software_version'] = annotation_info.get('softwareVersion')

    # Assembly Info
    assembly_info = json_data.get('assemblyInfo', {})
    metadata['assembly_level'] = assembly_info.get('assemblyLevel')
    metadata['assembly_method'] = assembly_info.get('assemblyMethod')
    metadata['assembly_status'] = assembly_info.get('assemblyStatus')
    metadata['assembly_type'] = assembly_info.get('assemblyType')
    metadata['sequencing_tech'] = assembly_info.get('sequencingTech')
    metadata['bioproject_accession'] = assembly_info.get('bioprojectAccession')
    
    # Biosample Info
    biosample = assembly_info.get('biosample', {})
    metadata['biosample_accession'] = biosample.get('accession')
    metadata['strain'] = biosample.get('strain')
    metadata['host'] = biosample.get('host')
    metadata['geo_loc_name'] = biosample.get('geoLocName')
    metadata['isolation_source'] = biosample.get('isolationSource')

    # Assembly Stats
    assembly_stats = json_data.get('assemblyStats', {})
    metadata['contig_n50'] = assembly_stats.get('contigN50')
    metadata['scaffold_n50'] = assembly_stats.get('scaffoldN50')
    metadata['gc_percent'] = assembly_stats.get('gcPercent')
    metadata['genome_coverage'] = assembly_stats.get('genomeCoverage')

    # CheckM Info
    checkm_info = json_data.get('checkmInfo', {})
    metadata['completeness'] = checkm_info.get('completeness')
    metadata['contamination'] = checkm_info.get('contamination')

    # Organism Info
    organism = json_data.get('organism', {})
    metadata['organism_name'] = organism.get('organismName')
    metadata['tax_id'] = organism.get('taxId')

    # Type Material
    type_material = json_data.get('typeMaterial', {})
    metadata['type_material_label'] = type_material.get('typeLabel')
    metadata['type_material_display_text'] = type_material.get('typeDisplayText')

    # WGS Info
    wgs_info = json_data.get('wgsInfo', {})
    metadata['wgs_project_accession'] = wgs_info.get('wgsProjectAccession')
    metadata['master_wgs_url'] = wgs_info.get('masterWgsUrl')
    metadata['wgs_contigs_url'] = wgs_info.get('wgsContigsUrl')

    return metadata

def main(input_file, output_file):
    # Read JSON data from a file
    with open(input_file, 'r') as file:
        json_data = json.load(file)

    # Parse the metadata
    metadata = parse_metadata(json_data)

    # Output the metadata as JSON
    with open(output_file, 'w') as out_file:
        json.dump(metadata, out_file, indent=4)

    print(f"Extracted metadata has been written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sample metadata from JSON and output as JSON.")
    parser.add_argument('input_file', type=str, help='Input file containing JSON data.')
    parser.add_argument('output_file', type=str, help='Output file to write the extracted metadata as JSON.')
    
    args = parser.parse_args()

    main(args.input_file, args.output_file)

