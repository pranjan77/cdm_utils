import unittest
import os
from unittest.mock import patch, MagicMock
from io import StringIO

from cdm_utils.bbmap_assembly_stats import BBMapAssemblyStats

class TestBBMapAssemblyStatsWithFiles(unittest.TestCase):

    def setUp(self):
        # Path to the assembly file inside the tests/data directory
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.bbmap_output_file = os.path.join(self.data_dir, 'bbmap_output.txt')
        # Load the sample output from the file
        with open(self.bbmap_output_file, 'r') as file:
            self.bbmap_output = file.read()

    def test_run_bbmap_stats_with_real_files(self):
        # This method is a placeholder for running tests with real files.
        # You can implement this as needed, or if not needed, you can remove it.
        pass


    def test_parse_bbmap_output(self):
        # Create instance of BBMapAssemblyStats
        bbmap_parser = BBMapAssemblyStats()

        # Parse the output (loaded in setUp)
        bbmap_parser.parse_bbmap_output(self.bbmap_output)

        # Get the parsed stats
        stats = bbmap_parser.get_stats()

        expected_result = {
            'A_content': 0.1862, 'C_content': 0.3191, 'G_content': 0.3149, 'T_content': 0.1798, 'N_content': 0.0004,
            'IUPAC_content': 0.0, 'Other_content': 0.0, 'GC_content': 0.634, 'GC_stdev': 0.0805, 'scaffold_total': 13,
            'contig_total': 20, 'scaffold_sequence_total': '1.879 MB', 'contig_sequence_total': '1.878 MB  ',
            'contig_gap_percentage': '0.042% gap', 'scaffold_N50': '2 KB', 'scaffold_L50': '247.538 KB', 
            'contig_N50': '5 KB', 'contig_L50': '169.915 KB', 'scaffold_N90': '6 KB', 'scaffold_L90': '136.211 KB', 
            'contig_N90': '10 KB', 'contig_L90': '95.394 KB', 'max_scaffold_length': '859.216 KB', 
            'max_contig_length': '314.009 KB', 'large_scaffold_count_gt_50kb': 6, 
            'percent_genome_in_large_scaffolds_gt_50kb': 96.55
        }

        # Check if the parsed data matches the expected result
        self.assertEqual(stats, expected_result)

if __name__ == '__main__':
    unittest.main()

