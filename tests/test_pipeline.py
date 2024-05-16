# tests/test_pipeline.py

import unittest
from unittest.mock import patch, MagicMock
from superseeker.pipeline import run_pipeline
from superseeker.data_manipulation import (
    vcf_to_pyclone_input,
    pyclone_to_vcf,
    extract_sseeker_stats,
    make_dot_files
)

class TestSuperSeekerPipeline(unittest.TestCase):
    @patch('subprocess.run')
    def test_run_pipeline(self, mock_subprocess_run):
        # Mock the subprocess.run to always return success
        mock_subprocess_run.return_value = MagicMock(returncode=0)

        # Call the run_pipeline function with test parameters
        run_pipeline(
            patient='test_patient',
            vcf_file='test.vcf',
            patient_sex='M',
            restarts=10,
            clusters=5,
            cn_neutral=True,
            cn_override=False,
            germfilter=True
        )

        # Assert subprocess.run was called (you can check specific calls if needed)
        self.assertTrue(mock_subprocess_run.called)

    def test_vcf_to_pyclone_input(self):
        # Example test for vcf_to_pyclone_input
        vcf_to_pyclone_input(
            patient='test_patient',
            vcf_file='test.vcf',
            facets_dir='test_dir',
            output_file='output.tsv',
            patient_sex='M',
            cn_neutral=True,
            cn_override=False,
            germfilter=True
        )
        # Add assertions as needed based on the expected behavior of vcf_to_pyclone_input

    def test_pyclone_to_vcf(self):
        # Example test for pyclone_to_vcf
        pyclone_to_vcf(
            patient='test_patient',
            vcf_file='test.vcf',
            clustered_file='clustered.tsv',
            output_file='output.vcf'
        )
        # Add assertions as needed based on the expected behavior of pyclone_to_vcf

    def test_extract_sseeker_stats(self):
        # Example test for extract_sseeker_stats
        extract_sseeker_stats(
            stats_file='stats.txt',
            output_file='output.txt'
        )
        # Add assertions as needed based on the expected behavior of extract_sseeker_stats

    def test_make_dot_files(self):
        # Example test for make_dot_files
        make_dot_files(
            subclones_vcf='subclones.vcf',
            tmp_graph_files='tmp_graph_files'
        )
        # Add assertions as needed based on the expected behavior of make_dot_files

if __name__ == '__main__':
    unittest.main()
