# superseeker/data_processing.py

import pandas as pd
import matplotlib.pyplot as plt
import logging
from pathlib import Path
from typing import Dict, List, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_CN_info(samples: List[str], facets_dir: str, patient_sex: str) -> Tuple[List[List[str]], str]:
    """Extract copy number information from FACETS output files.
    
    Args:
        samples: List of sample identifiers
        facets_dir: Path to FACETS output directory
        patient_sex: Patient sex ('F' or 'M')
        
    Returns:
        Tuple containing:
            - List of copy number entries for each sample
            - Normal copy number for X chromosome
            
    Raises:
        FileNotFoundError: If FACETS output files are not found
        ValueError: If patient sex is invalid
    """
    if patient_sex not in ['F', 'M']:
        raise ValueError("Patient sex must be 'F' or 'M'")
        
    output_list = []
    x_cn_list: Dict[str, int] = {}
    
    for samp in samples:
        cn_file_path = Path(facets_dir) / f"{samp}.cncf.txt"
        if not cn_file_path.exists():
            raise FileNotFoundError(f"FACETS output file not found: {cn_file_path}")
            
        logger.info(f"Processing FACETS output for sample {samp}")
        entry = []
        
        with open(cn_file_path, "r") as cn_file:
            for line in cn_file:
                info = line.strip().split("\t")
                if info[0] == "chrom":
                    continue
                    
                chrom = info[0]
                start = info[9]
                end = info[10]
                cell_fraction = info[11]
                total_cn = info[12]
                minor_cn = info[13]
                
                # Handle X chromosome copy number
                if chrom == "23" and patient_sex not in ["M", "Y"]:
                    dist = int(end) - int(start)
                    x_cn_list[total_cn] = x_cn_list.get(total_cn, 0) + dist
                    
                if minor_cn == "NA":
                    continue
                    
                major_cn = str(int(total_cn) - int(minor_cn))
                if major_cn == "0":
                    major_cn = "1"
                    minor_cn = "0"
                    
                entry.append([chrom, start, end, cell_fraction, total_cn, major_cn, minor_cn])
                
        output_list.append(entry)
        
    # Determine X chromosome copy number
    if patient_sex == "M":
        X_cn = "1"
    elif patient_sex == "F" or not x_cn_list:
        X_cn = "2"
    else:
        X_cn = max(x_cn_list.items(), key=lambda x: x[1])[0]
        
    return output_list, X_cn

def get_copy_numbers(
    sample_cn_data: List[List[str]],
    chr: str,
    position: str,
    X_normal_cn: int = 2
) -> Tuple[str, str]:
    """Get copy numbers for a specific genomic position.
    
    Args:
        sample_cn_data: Copy number data for a sample
        chr: Chromosome
        position: Genomic position
        X_normal_cn: Normal copy number for X chromosome
        
    Returns:
        Tuple of (major_cn, minor_cn)
    """
    major = "1"
    minor = "1"
    
    if len(chr) > 3 and chr[:3] == "chr":
        chr = chr[3:]
        
    for val in sample_cn_data:
        if val[0] != chr:
            continue
        if int(val[1]) < int(position) and int(val[2]) > int(position):
            major = val[5]
            minor = val[6]
            return major, minor
            
    if chr == "X":
        minor = str(int(X_normal_cn) - 1)
        
    return major, minor

def get_variant_output_lines(
    variant_line: str,
    samples: List[str],
    output_file,
    HIGH_IMPACT: bool = False,
    germfilter: bool = False,
    cn_override: bool = False,
    cn_neutral: bool = False,
    patient_sex: str = "F",
    X_normal_cn: int = 2,
    samples_cn_lists: List[List[List[str]]] = []
) -> None:
    """Process a variant line and write output to file.
    
    Args:
        variant_line: Line from VCF file
        samples: List of sample identifiers
        output_file: File object to write output to
        HIGH_IMPACT: Whether to filter for high impact variants
        germfilter: Whether to filter germline variants
        cn_override: Whether to override copy number information
        cn_neutral: Whether to exclude variants in CNV regions
        patient_sex: Patient sex ('F' or 'M')
        X_normal_cn: Normal copy number for X chromosome
        samples_cn_lists: Copy number data for each sample
    """
    fields = variant_line.strip().split("\t")
    mutation_id = f"{fields[0]}:{fields[1]}:{fields[4]}"
    info = fields[7].split(";")
    
    if HIGH_IMPACT:
        for i in reversed(info):
            if i[:4] == "ANN=":
                ann = i.split("|")
                for a in ann:
                    if a == "HIGH" or a == "MODERATE":
                        break
                    if a == "MODIFIER" or a == "LOW":
                        return
                break
                
    Format = fields[8].split(':')
    AO_index = Format.index("AO")
    RO_index = Format.index("RO")
    
    if germfilter:
        germline_alt = fields[9].split(":")[AO_index]
        germline_ref = fields[9].split(":")[RO_index]
        Germline_AF = int(germline_alt)/(int(germline_alt)+int(germline_ref))
        if Germline_AF >= 0.01:
            pass
            
    output = ""
    for i, sample_id in enumerate(samples):
        if germfilter:
            sample_info = fields[10+i].split(":")
        else:
            sample_info = fields[9+i].split(":")
            
        ref_counts = sample_info[RO_index]
        alt_counts = sample_info[AO_index]
        
        if cn_override:
            major_cn, minor_cn = "1", "1"
            normal_cn = "1" if fields[0] == "chrX" and patient_sex == "M" else "2"
        else:
            major_cn, minor_cn = get_copy_numbers(samples_cn_lists[i], fields[0], fields[1])
            normal_cn = X_normal_cn if fields[0] == "chrX" else "2"
            
            if cn_neutral and int(major_cn) + int(minor_cn) != int(normal_cn):
                logger.debug(f"Skipping variant in CNV region: {mutation_id}")
                return
                
        output += f"{mutation_id}\t{sample_id}\t{ref_counts}\t{alt_counts}\t{major_cn}\t{minor_cn}\t{normal_cn}\n"
        
    output_file.write(output)

def vcf_to_pyclone_input(
    vcf_file_name: str,
    facets_dir: str,
    output_file_name: str,
    patient_sex: str,
    cn_neutral: bool,
    cn_override: bool,
    germfilter: bool
) -> None:
    """Convert VCF file to PyClone-VI input format.
    
    Args:
        vcf_file_name: Path to VCF file
        facets_dir: Path to FACETS output directory
        output_file_name: Path to output file
        patient_sex: Patient sex ('F' or 'M')
        cn_neutral: Whether to exclude variants in CNV regions
        cn_override: Whether to override copy number information
        germfilter: Whether to filter germline variants
        
    Raises:
        FileNotFoundError: If input files are not found
        ValueError: If VCF file format is invalid
    """
    vcf_path = Path(vcf_file_name)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_file_name}")
        
    if vcf_path.suffix == ".vcf.gz":
        raise ValueError("Compressed VCF files are not supported. Please decompress first.")
    elif vcf_path.suffix != ".vcf":
        raise ValueError("Input file must be a VCF file")
        
    logger.info(f"Converting VCF file to PyClone input format: {vcf_file_name}")
    
    with open(vcf_file_name, "r") as vcf_file, open(output_file_name, "w") as output_file:
        variant_lines = []
        header_line = ""
        
        for line in vcf_file:
            if line[:2] == "##":
                continue
            if line[0] == "#":
                header_line = line.strip()
            else:
                variant_lines.append(line.strip())
                
        columns = header_line.split("\t")
        if germfilter:
            germline_sample = columns[9]
            samples = columns[10:]
        else:
            samples = columns[9:]
            
        logger.info(f"Processing {len(samples)} samples")
        
        samples_cn_lists = []
        X_normal_cn = 2
        
        if not cn_override:
            samples_cn_lists, X_normal_cn = get_CN_info(samples, facets_dir, patient_sex)
            
        output_file.write("mutation_id\tsample_id\tref_counts\talt_counts\tmajor_cn\tminor_cn\tnormal_cn\n")
        
        for line in variant_lines:
            get_variant_output_lines(
                line, samples, output_file, False, germfilter,
                cn_override, cn_neutral, patient_sex, X_normal_cn, samples_cn_lists
            )

def pyclone_to_vcf(
    vcf_file: str,
    clustered_file: str,
    output_file_name: str,
    HIGH_IMPACT: bool = False
) -> None:
    """Add PyClone cluster assignments to VCF file.
    
    Args:
        vcf_file: Path to input VCF file
        clustered_file: Path to PyClone clustering results
        output_file_name: Path to output VCF file
        HIGH_IMPACT: Whether to filter for high impact variants
        
    Raises:
        FileNotFoundError: If input files are not found
    """
    vcf_path = Path(vcf_file)
    clustered_path = Path(clustered_file)
    
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")
    if not clustered_path.exists():
        raise FileNotFoundError(f"Clustered file not found: {clustered_file}")
        
    logger.info("Adding PyClone clusters to VCF file")
    
    cluster_assignments = {}
    with open(clustered_file, 'r') as pyclone_file:
        for line in pyclone_file:
            fields = line.strip().split('\t')
            if fields[0] == "mutation_id":
                continue
            elif fields[0] not in cluster_assignments:
                cluster_assignments[fields[0]] = [fields[2]]
                
    info_line = "##INFO=<ID=AFCLU,Number=1,Type=String,Description=\"The allele frequency cluster this variant belongs to\">"
    
    with open(vcf_file, 'r') as vcf_file, open(output_file_name, "w") as output_file:
        for line in vcf_file:
            if line.strip()[:2] == "##":
                output_file.write(line.strip() + '\n')
            elif line.strip()[0] == "#":
                output_file.write(info_line + '\n')
                output_file.write(line.strip() + '\n')
            else:
                fields = line.split('\t')
                mut_id = f"{fields[0]}:{fields[1]}:{fields[4]}"
                
                if mut_id not in cluster_assignments:
                    logger.debug(f"Skipping variant not found in clustering results: {mut_id}")
                    continue
                    
                skip = False
                info = fields[7].split(";")
                
                if HIGH_IMPACT:
                    for i in reversed(info):
                        if i[:4] == "ANN=":
                            ann = i.split("|")
                            for a in ann:
                                if a == "HIGH" or a == "MODERATE":
                                    break
                                if a == "MODIFIER" or a == "LOW":
                                    skip = True
                            break
                            
                if not skip:
                    info_fields = fields[7].split(";")
                    cluster = cluster_assignments[mut_id]
                    
                    if info_fields[-1].split("=")[0] == "AFCLU":
                        info_fields[-1] = f"AFCLU={cluster[0]}"
                    else:
                        info_fields.append(f"AFCLU={cluster[0]}")
                        
                    fields[7] = ";".join(info_fields)
                    output_file.write("\t".join(fields))

class Subclone:
    """Class representing a subclone with its evolutionary properties."""
    
    def __init__(self, ID: str):
        """Initialize a new Subclone.
        
        Args:
            ID: Subclone identifier
        """
        self.ID = ID
        self.vafs: List[float] = []
        self.selection = False
        self.emergence = False
        self.replacement = False
        
    def add_vaf(self, vaf: float) -> None:
        """Add a variant allele frequency to the subclone.
        
        Args:
            vaf: Variant allele frequency
        """
        self.vafs.append(vaf)

def find_replacement(subclones: Dict[str, Subclone]) -> None:
    """Identify subclones that show replacement pattern.
    
    Args:
        subclones: Dictionary of Subclone objects
    """
    for ID in subclones:
        highest_at_end = True
        not_highest_at_start = False
        
        if subclones[ID].emergence:
            for check in subclones:
                if check == ID:
                    continue
                if subclones[check].vafs[-1] >= subclones[ID].vafs[-1]:
                    highest_at_end = False
                if subclones[check].vafs[1] > subclones[ID].vafs[1]:
                    not_highest_at_start = True
                    
            if not_highest_at_start and highest_at_end:
                subclones[ID].replacement = True

def find_evolution(subclones: Dict[str, Subclone]) -> None:
    """Identify evolutionary patterns in subclones.
    
    Args:
        subclones: Dictionary of Subclone objects
    """
    for ID in subclones:
        i = 1  # Start after germline
        change = 0.0
        
        while i < len(subclones[ID].vafs) - 1:
            change += float(subclones[ID].vafs[i+1]) - float(subclones[ID].vafs[i])
            
            if change >= 0.1:
                subclones[ID].emergence = True
            if change <= -0.1:
                subclones[ID].selection = True
                
            i += 1

def identify_evolution(stats_file: str, output_file: str) -> None:
    """Identify evolutionary patterns from SuperSeeker stats.
    
    Args:
        stats_file: Path to SuperSeeker stats file
        output_file: Path to output file
        
    Raises:
        FileNotFoundError: If input file is not found
    """
    stats_path = Path(stats_file)
    if not stats_path.exists():
        raise FileNotFoundError(f"Stats file not found: {stats_file}")
        
    logger.info("Identifying evolutionary patterns")
    
    with open(stats_file, 'r') as input_file, open(output_file, 'w') as output_file:
        subclones = {}
        clusters = input_file.readline().strip().split("\t")
        
        for cluster in clusters:
            subclones[cluster] = Subclone(cluster)
            
        for line in input_file:
            if line[:4] == "Move":
                break
                
            sample = line.strip().split("\t")
            if len(sample) < 4:
                logger.warning("Not enough samples to identify evolution")
                break
                
            for i, val in enumerate(sample[1:], 1):
                subclones[str(i)].add_vaf(float(val))
                
        find_evolution(subclones)
        find_replacement(subclones)
        
        selection = False
        emergence = False
        replacement = False
        
        selection_list = []
        emergence_list = []
        replacement_list = []
        
        for ID in subclones:
            if subclones[ID].selection:
                selection = True
                selection_list.append(ID)
            if subclones[ID].emergence:
                emergence = True
                emergence_list.append(ID)
            if subclones[ID].replacement:
                replacement = True
                replacement_list.append(ID)
                
        if replacement:
            evolution_type = "Replacement"
        elif emergence:
            evolution_type = "Positive Selection"
        elif selection:
            evolution_type = "Negative Selection"
        else:
            evolution_type = "No Evolution"
            
        logger.info(f"Evolution type: {evolution_type}")
        output_file.write(f"{evolution_type}\n")
        
        output_file.write(f"Subclones with Negative Selection: {','.join(selection_list)}\n")
        output_file.write(f"Subclones with Positive Selection: {','.join(emergence_list)}\n")
        output_file.write(f"Subclones with Replacement: {','.join(replacement_list)}\n")

def make_dot_files(subclones_vcf: str, tmp_graph_files: str) -> None:
    """Create DOT files for evolutionary tree visualization.
    
    Args:
        subclones_vcf: Path to SuperSeeker VCF file
        tmp_graph_files: Path to output directory for DOT files
        
    Raises:
        FileNotFoundError: If input file is not found
    """
    vcf_path = Path(subclones_vcf)
    if not vcf_path.exists():
        raise FileNotFoundError(f"SuperSeeker VCF file not found: {subclones_vcf}")
        
    logger.info("Creating DOT files for tree visualization")
    
    with open(subclones_vcf, "r") as infile:
        i = 0
        for line in infile:
            if line[0] != "#":
                break
            if line[0:10] == "##subclone":
                i += 1
                edges = line.split("\"")[1]
                outfile_path = Path(tmp_graph_files) / f"solution{i}.gv"
                
                with open(outfile_path, "w") as outfile:
                    outfile.write("digraph D{\n")
                    for edge in edges.split(", "):
                        outfile.write(f"{edge}\n")
                    outfile.write(f"label=\"Solution {i}\"\nlabelloc=\"t\"\n}}\n")

def make_line_plot(
    vcf_file_name: str,
    title: str,
    show: bool = False,
    save: bool = True,
    plot_file_name: str = "vaf.png"
) -> None:
    """Create line plot of variant allele frequencies.
    
    Args:
        vcf_file_name: Path to VCF file
        title: Plot title
        show: Whether to display the plot
        save: Whether to save the plot
        plot_file_name: Path to save plot
        
    Raises:
        FileNotFoundError: If input file is not found
    """
    vcf_path = Path(vcf_file_name)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_file_name}")
        
    logger.info("Creating VAF line plot")
    
    # Create DataFrame of variant allele frequencies
    variant_df = pd.DataFrame(columns=["Sample", "Position", "BP_Change", "Cluster", "Allele_Frequency"])
    
    with open(vcf_file_name, "r") as infile:
        header = None
        for line in infile:
            if line[0] == "#":
                header = line.strip().split()
            else:
                fields = line.strip().split()
                clu = fields[7].split(";")[-1].split("=")[1]
                Format = fields[8].split(":")
                AO_index = Format.index("AO")
                DP_index = Format.index("DP")
                
                for samp_col_num in range(9, len(fields)):
                    sample_name = header[samp_col_num]
                    dp = fields[samp_col_num].split(":")[DP_index]
                    ao = fields[samp_col_num].split(":")[AO_index]
                    af = float(ao)/float(dp)
                    variant_df.loc[len(variant_df.index)] = [
                        sample_name,
                        f"{fields[0]}:{fields[1]}",
                        f"{fields[3]}->{fields[4]}",
                        clu,
                        af
                    ]
    
    # Create line plot
    sample_order = variant_df['Sample'].unique()
    allele_freq_matrix = variant_df.pivot(
        index='Sample',
        columns='Position',
        values='Allele_Frequency'
    ).reindex(sample_order)
    
    unique_clusters = variant_df['Cluster'].unique()
    cluster_colors = {
        cluster: plt.get_cmap('gist_rainbow')(i / len(unique_clusters))
        for i, cluster in enumerate(unique_clusters)
    }
    
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 20,
        'axes.labelsize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14
    })
    
    plt.figure(figsize=(16, 10))
    
    for position in allele_freq_matrix.columns:
        mutation_cluster = variant_df[variant_df['Position'] == position]['Cluster'].iloc[0]
        color = cluster_colors[mutation_cluster]
        plt.plot(
            allele_freq_matrix.index,
            allele_freq_matrix[position],
            label=position,
            linestyle='-',
            linewidth=2,
            alpha=0.7,
            color=color
        )
    
    plt.xticks(rotation=90)
    plt.title(title)
    plt.ylabel('Cellular Prevalence')
    plt.xlabel('Sample ID')
    plt.legend(
        title='Mutation ID',
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        fontsize='small'
    )
    plt.tight_layout()
    
    if show:
        plt.show()
    if save:
        plt.savefig(plot_file_name, bbox_inches='tight')
        logger.info(f"Plot saved to {plot_file_name}")