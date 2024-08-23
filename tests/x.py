import unittest
from unittest.mock import patch, MagicMock
from io import StringIO

class TestBBMapAssemblyStats(unittest.TestCase):

    def setUp(self):
        # Load the sample output from the file
        with open('tests/data/bbmap_output.txt', 'r') as file:
            self.bbmap_output = file.read()

    @patch('subprocess.run')
    def test_run_bbmap_stats(self, mock_run):
        # Mock the subprocess.run method using the content from bbmap_output.txt
        mock_run.return_value = MagicMock(stdout=self.bbmap_output.encode('utf-8'), stderr=b'', returncode=0)

        # Create instance of BBMapAssemblyStats
        bbmap_parser = BBMapAssemblyStats(config_file='config.ini')
        
        # Create a temporary assembly file path (mocked)
        assembly_file = 'test_assembly.fasta'

        # Run the run_bbmap_stats method
        output = bbmap_parser.run_bbmap_stats(assembly_file)

        # Check if the output matches the mocked content
        self.assertEqual(output, self.bbmap_output)

        # Ensure subprocess.run was called with correct parameters
        mock_run.assert_called_once_with([bbmap_parser.stats_sh_path, f'in={assembly_file}'],
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    def test_parse_bbmap_output(self):
        # Create instance of BBMapAssemblyStats
        bbmap_parser = BBMapAssemblyStats()

        # Parse the output (loaded in setUp)
        bbmap_parser.parse_bbmap_output(self.bbmap_output)

        # Get the parsed stats
        stats = bbmap_parser.get_stats()

        # Assertions to check if the stats dictionary is correct
        self.assertEqual(stats['A_content'], 0.1862)
        self.assertEqual(stats['C_content'], 0.3191)
        self.assertEqual(stats['G_content'], 0.3149)
        self.assertEqual(stats['T_content'], 0.1798)
        self.assertEqual(stats['N_content'], 0.0004)
        self.assertEqual(stats['GC_content'], 0.6340)
        self.assertEqual(stats['GC_stdev'], 0.0805)
        self.assertEqual(stats['scaffold_total'], 13)
        self.assertEqual(stats['contig_total'], 20)
        self.assertEqual(stats['scaffold_sequence_total'], "1.879 MB")
        self.assertEqual(stats['contig_sequence_total'], "1.878 MB")
        self.assertEqual(stats['contig_gap_percentage'], "0.042% gap")
        self.assertEqual(stats['scaffold_N50'], "2 KB")
        self.assertEqual(stats['scaffold_L50'], "247.538 KB")
        self.assertEqual(stats['contig_N50'], "5 KB")
        self.assertEqual(stats['contig_L50'], "169.915 KB")
        self.assertEqual(stats['scaffold_N90'], "6 KB")
        self.assertEqual(stats['scaffold_L90'], "136.211 KB")
        self.assertEqual(stats['contig_N90'], "10 KB")
        self.assertEqual(stats['contig_L90'], "95.394 KB")
        self.assertEqual(stats['max_scaffold_length'], "859.216 KB")
        self.assertEqual(stats['max_contig_length'], "314.009 KB")
        self.assertEqual(stats['large_scaffold_count_gt_50kb'], 6)
        self.assertEqual(stats['percent_genome_in_large_scaffolds_gt_50kb'], 96.55)

if __name__ == '__main__':
    unittest.main()

