import unittest
import os
import tempfile
from superseeker.data_processing import (
    get_CN_info,
    get_copy_numbers,
    vcf_to_pyclone_input,
    pyclone_to_vcf,
    identify_evolution,
    make_dot_files,
    make_line_plot,
    Subclone
)

class TestDataProcessing(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.samples = ["SAMPLE1", "SAMPLE2"]
        self.facets_dir = os.path.join(self.test_dir, "facets")
        os.makedirs(self.facets_dir, exist_ok=True)
        
        # Create test FACETS output files
        for sample in self.samples:
            with open(os.path.join(self.facets_dir, f"{sample}.cncf.txt"), "w") as f:
                f.write("""chrom\tstart\tend\tcell_fraction\ttotal_cn\tminor_cn
1\t1\t1000\t0.5\t2\t1
X\t1\t1000\t0.5\t2\t1
""")
        
        # Create test VCF file
        self.vcf_file = os.path.join(self.test_dir, "test.vcf")
        with open(self.vcf_file, "w") as f:
            f.write("""##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="Annotation">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
chr1\t100\t.\tA\tT\t.\tPASS\tANN=MODIFIER\tGT:AO:RO:DP\t0/1:10:90:100\t0/1:20:80:100
chrX\t200\t.\tG\tC\t.\tPASS\tANN=HIGH\tGT:AO:RO:DP\t0/1:15:85:100\t0/1:25:75:100
""")

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.test_dir)

    def test_get_CN_info(self):
        """Test copy number information extraction."""
        output_list, X_cn = get_CN_info(self.samples, self.facets_dir, "F")
        
        # Check output structure
        self.assertEqual(len(output_list), len(self.samples))
        self.assertEqual(X_cn, "2")  # Female X chromosome copy number
        
        # Check content of first sample's CN data
        first_sample_data = output_list[0]
        self.assertEqual(len(first_sample_data), 2)  # Two regions
        self.assertEqual(first_sample_data[0][0], "1")  # Chromosome
        self.assertEqual(first_sample_data[0][4], "2")  # Total CN
        self.assertEqual(first_sample_data[0][5], "1")  # Major CN
        self.assertEqual(first_sample_data[0][6], "1")  # Minor CN

    def test_get_copy_numbers(self):
        """Test copy number retrieval for specific positions."""
        sample_cn_data = [
            ["1", "1", "1000", "0.5", "2", "1", "1"],
            ["X", "1", "1000", "0.5", "2", "1", "1"]
        ]
        
        # Test autosomal position
        major, minor = get_copy_numbers(sample_cn_data, "1", "500")
        self.assertEqual(major, "1")
        self.assertEqual(minor, "1")
        
        # Test X chromosome position
        major, minor = get_copy_numbers(sample_cn_data, "X", "500", X_normal_cn=2)
        self.assertEqual(major, "1")
        self.assertEqual(minor, "1")

    def test_vcf_to_pyclone_input(self):
        """Test VCF to PyClone input conversion."""
        output_file = os.path.join(self.test_dir, "pyclone_input.tsv")
        
        vcf_to_pyclone_input(
            self.vcf_file,
            self.facets_dir,
            output_file,
            "F",
            cn_neutral=False,
            cn_override=False,
            germfilter=False
        )
        
        # Check if output file exists and has correct format
        self.assertTrue(os.path.exists(output_file))
        
        # Read and validate output
        with open(output_file, "r") as f:
            header = f.readline().strip()
            self.assertEqual(header, "mutation_id\tsample_id\tref_counts\talt_counts\tmajor_cn\tminor_cn\tnormal_cn")
            
            # Check first line of data
            first_line = f.readline().strip()
            fields = first_line.split("\t")
            self.assertEqual(len(fields), 7)
            self.assertTrue(fields[0].startswith("chr1:100:T"))

    def test_pyclone_to_vcf(self):
        """Test PyClone to VCF conversion."""
        # Create test PyClone output
        pyclone_file = os.path.join(self.test_dir, "pyclone_output.tsv")
        with open(pyclone_file, "w") as f:
            f.write("""mutation_id\tsample_id\tcluster_id
chr1:100:T\tSAMPLE1\t1
chr1:100:T\tSAMPLE2\t1
chrX:200:C\tSAMPLE1\t2
chrX:200:C\tSAMPLE2\t2
""")
        
        output_vcf = os.path.join(self.test_dir, "output.vcf")
        pyclone_to_vcf(self.vcf_file, pyclone_file, output_vcf)
        
        # Check if output file exists
        self.assertTrue(os.path.exists(output_vcf))
        
        # Validate VCF format and cluster assignments
        with open(output_vcf, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                info = fields[7]
                self.assertIn("AFCLU=", info)

    def test_identify_evolution(self):
        """Test evolution pattern identification."""
        # Create test stats file
        stats_file = os.path.join(self.test_dir, "stats.txt")
        with open(stats_file, "w") as f:
            f.write("""Cluster1\tCluster2\tCluster3
0.1\t0.2\t0.3
0.2\t0.1\t0.4
0.3\t0.1\t0.5
""")
        
        output_file = os.path.join(self.test_dir, "evolution.txt")
        identify_evolution(stats_file, output_file)
        
        # Check if output file exists and contains expected content
        self.assertTrue(os.path.exists(output_file))
        with open(output_file, "r") as f:
            content = f.read()
            self.assertIn("Evolution type:", content)

    def test_make_dot_files(self):
        """Test DOT file generation for tree visualization."""
        # Create test SuperSeeker VCF
        subclones_vcf = os.path.join(self.test_dir, "subclones.vcf")
        with open(subclones_vcf, "w") as f:
            f.write("""##subclone="1->2, 2->3"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\tPASS\tAFCLU=1
""")
        
        tmp_graph_files = os.path.join(self.test_dir, "tmp_graph_files")
        os.makedirs(tmp_graph_files, exist_ok=True)
        
        make_dot_files(subclones_vcf, tmp_graph_files)
        
        # Check if DOT file was created
        dot_file = os.path.join(tmp_graph_files, "solution1.gv")
        self.assertTrue(os.path.exists(dot_file))
        
        # Validate DOT file content
        with open(dot_file, "r") as f:
            content = f.read()
            self.assertIn("digraph", content)
            self.assertIn("1->2", content)
            self.assertIn("2->3", content)

    def test_make_line_plot(self):
        """Test VAF line plot generation."""
        output_plot = os.path.join(self.test_dir, "vaf_plot.png")
        
        make_line_plot(
            self.vcf_file,
            "Test VAF Plot",
            show=False,
            save=True,
            plot_file_name=output_plot
        )
        
        # Check if plot was created
        self.assertTrue(os.path.exists(output_plot))

    def test_subclone_class(self):
        """Test Subclone class functionality."""
        subclone = Subclone("1")
        
        # Test VAF addition
        subclone.add_vaf(0.1)
        subclone.add_vaf(0.2)
        self.assertEqual(len(subclone.vafs), 2)
        self.assertEqual(subclone.vafs[0], 0.1)
        self.assertEqual(subclone.vafs[1], 0.2)
        
        # Test evolution flags
        self.assertFalse(subclone.selection)
        self.assertFalse(subclone.emergence)
        self.assertFalse(subclone.replacement)

if __name__ == '__main__':
    unittest.main() 