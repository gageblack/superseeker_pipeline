# superseeker/pipeline.py
import subprocess
import os
import logging
from pathlib import Path
from typing import Optional
from .data_processing import (
    vcf_to_pyclone_input,
    pyclone_to_vcf,
    identify_evolution,
    make_dot_files,
    make_line_plot
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def validate_inputs(patient: str, vcf_file: str, facets_dir: Optional[str] = None) -> None:
    """Validate input parameters and files.
    
    Args:
        patient: Patient identifier
        vcf_file: Path to VCF file
        facets_dir: Path to FACETS output directory
        
    Raises:
        ValueError: If any input validation fails
        FileNotFoundError: If required files don't exist
    """
    if not patient or not isinstance(patient, str):
        raise ValueError("Patient identifier must be a non-empty string")
    
    vcf_path = Path(vcf_file)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")
    if not vcf_path.suffix == '.vcf':
        raise ValueError("VCF file must have .vcf extension")
    
    if facets_dir:
        facets_path = Path(facets_dir)
        if not facets_path.exists():
            raise FileNotFoundError(f"FACETS directory not found: {facets_dir}")

def run_pipeline(
    patient: str,
    vcf_file: str,
    facets_dir: str = '',
    patient_sex: str = '',
    restarts: int = 100,
    clusters: int = 10,
    cn_neutral: bool = False,
    cn_override: bool = False,
    germfilter: bool = True
) -> None:
    """Run the complete SuperSeeker pipeline.
    
    Args:
        patient: Patient identifier
        vcf_file: Path to VCF file
        facets_dir: Path to FACETS output directory
        patient_sex: Patient sex ('F' or 'M')
        restarts: Number of PyClone-VI restarts
        clusters: Number of clusters to use
        cn_neutral: Whether to exclude variants in CNV regions
        cn_override: Whether to override copy number information
        germfilter: Whether to filter germline variants
        
    Raises:
        ValueError: If input validation fails
        FileNotFoundError: If required files don't exist
        RuntimeError: If any pipeline step fails
    """
    try:
        validate_inputs(patient, vcf_file, facets_dir)
        
        output_dir = f'{patient}_superseeker_results'
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info(f"Starting pipeline for patient {patient}")
        
        # Convert VCF to PyClone input
        pyclone_input_file = f'{output_dir}/{patient}.pyclone_input.tsv'
        logger.info("Converting VCF to PyClone input format")
        vcf_to_pyclone_input(vcf_file, facets_dir, pyclone_input_file, patient_sex, cn_neutral, cn_override, germfilter)
        
        # Run PyClone
        clustered_file = f'{output_dir}/{patient}.pyclone.clustered.tsv'
        logger.info("Running PyClone-VI clustering")
        try:
            subprocess.run([
                'pyclone-vi', 'fit', '-i', pyclone_input_file, '-o', f'{output_dir}/{patient}.h5',
                '-d', 'beta-binomial', '-r', str(restarts), '-c', str(clusters)
            ], check=True)
            subprocess.run([
                'pyclone-vi', 'write-results-file', '-i', f'{output_dir}/{patient}.h5', '-o', clustered_file
            ], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"PyClone-VI failed: {str(e)}")
        
        # Add PyClone clusters to VCF
        clustered_vcf_file = f'{output_dir}/{patient}.somatic.clustered.vcf'
        logger.info("Adding PyClone clusters to VCF")
        pyclone_to_vcf(vcf_file, clustered_file, clustered_vcf_file)
        
        # Make visual representation of clustering
        logger.info("Generating VAF line plot")
        cluster_lines_pdf = f'{output_dir}/{patient}.cluster_lines.png'
        make_line_plot(clustered_vcf_file, "VAFs", show=False, save=True, plot_file_name=cluster_lines_pdf)
        
        # Run SuperSeeker
        subclones_vcf = f'{output_dir}/{patient}.subclones.vcf'
        logger.info("Running SuperSeeker")
        try:
            subprocess.run(['superseeker', clustered_vcf_file, '>', subclones_vcf], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SuperSeeker failed: {str(e)}")
        
        # Make a visual representation of the trees found
        tmp_graph_files = f'{output_dir}/tmp_graph_files'
        os.makedirs(tmp_graph_files, exist_ok=True)
        make_dot_files(subclones_vcf, tmp_graph_files)
        try:
            subprocess.run(['dot', '-Tpdf', tmp_graph_files, '-o', f'{output_dir}/{patient}.solutions.pdf'], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Graphviz failed: {str(e)}")
        
        # Extract SuperSeeker stats
        stats_file = f'{output_dir}/{patient}.stats.txt'
        evolution_file = f'{output_dir}/{patient}.evolution.txt'
        identify_evolution(stats_file, evolution_file)
        
        logger.info("Pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        raise

def cluster_variants(
    patient: str,
    vcf_file: str,
    facets_dir: str = '',
    patient_sex: str = '',
    restarts: int = 100,
    clusters: int = 10,
    cn_neutral: bool = False,
    cn_override: bool = False,
    germfilter: bool = True
) -> None:
    """Run only the variant clustering step using PyClone-VI.
    
    Args:
        patient: Patient identifier
        vcf_file: Path to VCF file
        facets_dir: Path to FACETS output directory
        patient_sex: Patient sex ('F' or 'M')
        restarts: Number of PyClone-VI restarts
        clusters: Number of clusters to use
        cn_neutral: Whether to exclude variants in CNV regions
        cn_override: Whether to override copy number information
        germfilter: Whether to filter germline variants
        
    Raises:
        ValueError: If input validation fails
        FileNotFoundError: If required files don't exist
        RuntimeError: If clustering step fails
    """
    try:
        validate_inputs(patient, vcf_file, facets_dir)
        
        output_dir = f'{patient}_superseeker_results'
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info(f"Starting variant clustering for patient {patient}")
        
        # Convert VCF to PyClone input
        pyclone_input_file = f'{output_dir}/{patient}.pyclone_input.tsv'
        logger.info("Converting VCF to PyClone input format")
        vcf_to_pyclone_input(vcf_file, facets_dir, pyclone_input_file, patient_sex, cn_neutral, cn_override, germfilter)
        
        # Run PyClone
        clustered_file = f'{output_dir}/{patient}.pyclone.clustered.tsv'
        logger.info("Running PyClone-VI clustering")
        try:
            subprocess.run([
                'pyclone-vi', 'fit', '-i', pyclone_input_file, '-o', f'{output_dir}/{patient}.h5',
                '-d', 'beta-binomial', '-r', str(restarts), '-c', str(clusters)
            ], check=True)
            subprocess.run([
                'pyclone-vi', 'write-results-file', '-i', f'{output_dir}/{patient}.h5', '-o', clustered_file
            ], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"PyClone-VI failed: {str(e)}")
        
        # Add PyClone clusters to VCF
        clustered_vcf_file = f'{output_dir}/{patient}.somatic.clustered.vcf'
        logger.info("Adding PyClone clusters to VCF")
        pyclone_to_vcf(vcf_file, clustered_file, clustered_vcf_file)
        
        # Make visual representation of clustering
        logger.info("Generating VAF line plot")
        cluster_lines_pdf = f'{output_dir}/{patient}.cluster_lines.png'
        make_line_plot(clustered_vcf_file, "VAFs", show=False, save=True, plot_file_name=cluster_lines_pdf)
        
        logger.info("Variant clustering completed successfully")
        
    except Exception as e:
        logger.error(f"Variant clustering failed: {str(e)}")
        raise

def run_superseeker(
    patient: str,
    vcf_file: str,
    facets_dir: str = '',
    patient_sex: str = '',
    restarts: int = 100,
    clusters: int = 10,
    cn_neutral: bool = False,
    cn_override: bool = False,
    germfilter: bool = True
) -> None:
    """Run only the SuperSeeker step on pre-clustered variants.
    
    Args:
        patient: Patient identifier
        vcf_file: Path to VCF file
        facets_dir: Path to FACETS output directory
        patient_sex: Patient sex ('F' or 'M')
        restarts: Number of PyClone-VI restarts
        clusters: Number of clusters to use
        cn_neutral: Whether to exclude variants in CNV regions
        cn_override: Whether to override copy number information
        germfilter: Whether to filter germline variants
        
    Raises:
        ValueError: If input validation fails
        FileNotFoundError: If required files don't exist
        RuntimeError: If SuperSeeker step fails
    """
    try:
        validate_inputs(patient, vcf_file, facets_dir)
        
        output_dir = f'{patient}_superseeker_results'
        clustered_vcf_file = f'{output_dir}/{patient}.somatic.clustered.vcf'
        
        if not os.path.isfile(clustered_vcf_file):
            raise FileNotFoundError(f"Clustered VCF file not found: {clustered_vcf_file}. Please run cluster_variants first.")
        
        logger.info(f"Starting SuperSeeker analysis for patient {patient}")
        
        # Run SuperSeeker
        subclones_vcf = f'{output_dir}/{patient}.subclones.vcf'
        logger.info("Running SuperSeeker")
        try:
            subprocess.run(['superseeker', clustered_vcf_file, '>', subclones_vcf], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SuperSeeker failed: {str(e)}")
        
        # Make a visual representation of the trees found
        tmp_graph_files = f'{output_dir}/tmp_graph_files'
        os.makedirs(tmp_graph_files, exist_ok=True)
        make_dot_files(subclones_vcf, tmp_graph_files)
        try:
            subprocess.run(['dot', '-Tpdf', tmp_graph_files, '-o', f'{output_dir}/{patient}.solutions.pdf'], check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Graphviz failed: {str(e)}")
        
        # Extract SuperSeeker stats
        stats_file = f'{output_dir}/{patient}.stats.txt'
        evolution_file = f'{output_dir}/{patient}.evolution.txt'
        identify_evolution(stats_file, evolution_file)
        
        logger.info("SuperSeeker analysis completed successfully")
        
    except Exception as e:
        logger.error(f"SuperSeeker analysis failed: {str(e)}")
        raise
