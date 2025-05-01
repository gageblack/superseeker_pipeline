# tests/test_pipeline.py

import unittest
import os
import tempfile
from unittest.mock import patch, MagicMock
from superseeker.pipeline import run_pipeline, cluster_variants, run_superseeker

class TestPipeline(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.patient = "test_patient"
        self.vcf_file = os.path.join(self.test_dir, "test.vcf")
        self.facets_dir = os.path.join(self.test_dir, "facets")
        os.makedirs(self.facets_dir, exist_ok=True)
        
        # Create a minimal VCF file
        with open(self.vcf_file, "w") as f:
            f.write("""##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="Annotation">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
chr1\t100\t.\tA\tT\t.\tPASS\tANN=MODIFIER\tGT:AO:RO:DP\t0/1:10:90:100\t0/1:20:80:100
""")
            
        # Create a minimal FACETS output file
        with open(os.path.join(self.facets_dir, "SAMPLE1.cncf.txt"), "w") as f:
            f.write("""chrom\tstart\tend\tcell_fraction\ttotal_cn\tminor_cn
1\t1\t1000\t0.5\t2\t1
""")
            
        with open(os.path.join(self.facets_dir, "SAMPLE2.cncf.txt"), "w") as f:
            f.write("""chrom\tstart\tend\tcell_fraction\ttotal_cn\tminor_cn
1\t1\t1000\t0.5\t2\t1
""")

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.test_dir)

    @patch('subprocess.run')
    def test_run_pipeline(self, mock_subprocess):
        """Test the complete pipeline run."""
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        run_pipeline(
            patient=self.patient,
            vcf_file=self.vcf_file,
            facets_dir=self.facets_dir,
            patient_sex='F',
            restarts=10,
            clusters=5
        )
        
        # Check if output directory was created
        output_dir = f"{self.patient}_superseeker_results"
        self.assertTrue(os.path.exists(output_dir))
        
        # Check if expected files were created
        expected_files = [
            f"{self.patient}.pyclone_input.tsv",
            f"{self.patient}.pyclone.clustered.tsv",
            f"{self.patient}.somatic.clustered.vcf",
            f"{self.patient}.cluster_lines.png",
            f"{self.patient}.subclones.vcf",
            f"{self.patient}.solutions.pdf",
            f"{self.patient}.evolution.txt"
        ]
        
        for file in expected_files:
            self.assertTrue(os.path.exists(os.path.join(output_dir, file)))

    @patch('subprocess.run')
    def test_cluster_variants(self, mock_subprocess):
        """Test the variant clustering step."""
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        cluster_variants(
            patient=self.patient,
            vcf_file=self.vcf_file,
            facets_dir=self.facets_dir,
            patient_sex='F'
        )
        
        # Check if output directory was created
        output_dir = f"{self.patient}_superseeker_results"
        self.assertTrue(os.path.exists(output_dir))
        
        # Check if expected files were created
        expected_files = [
            f"{self.patient}.pyclone_input.tsv",
            f"{self.patient}.pyclone.clustered.tsv",
            f"{self.patient}.somatic.clustered.vcf",
            f"{self.patient}.cluster_lines.png"
        ]
        
        for file in expected_files:
            self.assertTrue(os.path.exists(os.path.join(output_dir, file)))

    @patch('subprocess.run')
    def test_run_superseeker(self, mock_subprocess):
        """Test the SuperSeeker step."""
        mock_subprocess.return_value = MagicMock(returncode=0)
        
        # First run clustering to create required input
        cluster_variants(
            patient=self.patient,
            vcf_file=self.vcf_file,
            facets_dir=self.facets_dir,
            patient_sex='F'
        )
        
        # Then run SuperSeeker
        run_superseeker(
            patient=self.patient,
            vcf_file=self.vcf_file
        )
        
        # Check if expected files were created
        output_dir = f"{self.patient}_superseeker_results"
        expected_files = [
            f"{self.patient}.subclones.vcf",
            f"{self.patient}.solutions.pdf",
            f"{self.patient}.evolution.txt"
        ]
        
        for file in expected_files:
            self.assertTrue(os.path.exists(os.path.join(output_dir, file)))

    def test_invalid_vcf_file(self):
        """Test handling of invalid VCF file."""
        with self.assertRaises(FileNotFoundError):
            run_pipeline(
                patient=self.patient,
                vcf_file="nonexistent.vcf",
                facets_dir=self.facets_dir
            )

    def test_invalid_facets_dir(self):
        """Test handling of invalid FACETS directory."""
        with self.assertRaises(FileNotFoundError):
            run_pipeline(
                patient=self.patient,
                vcf_file=self.vcf_file,
                facets_dir="nonexistent_dir"
            )

    def test_invalid_patient_sex(self):
        """Test handling of invalid patient sex."""
        with self.assertRaises(ValueError):
            run_pipeline(
                patient=self.patient,
                vcf_file=self.vcf_file,
                facets_dir=self.facets_dir,
                patient_sex='X'  # Invalid sex
            )

if __name__ == '__main__':
    unittest.main()
