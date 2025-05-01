# SuperSeeker Pipeline

The SuperSeeker pipeline is a Python library for identifying subclonal evolution in cancer. It integrates various tools and scripts to process VCF files, cluster variants using PyClone-VI, and reconstruct subclonal evolutionary trees using SuperSeeker.

## Features

- Convert VCF files to PyClone-VI input format
- Cluster genetic variants using PyClone-VI
- Reconstruct subclonal evolutionary trees using SuperSeeker
- Visualize variant allele frequencies and evolutionary trees
- Identify patterns of evolution (selection, emergence, replacement)

## Installation

### Prerequisites

- Python 3.6 or higher
- Conda (for installing PyClone-VI)

### Step 1: Install External Tools

#### Installing PyClone-VI

PyClone-VI is used for clustering genetic variants. You can install it using Conda:

```bash
conda install -c bioconda pyclone-vi
```

#### Installing SuperSeeker

SuperSeeker is used for reconstructing the evolutionary history of subclones. You can install it from its source repository:

```bash
git clone https://github.com/yiq/SuperSeeker.git
cd SuperSeeker
pip install -r requirements.txt
```

### Step 2: Install the SuperSeeker Pipeline Package

Install the SuperSeeker Pipeline Python package:

```bash
pip install superseeker-pipeline
```

## Usage

### Basic Usage

The simplest way to use the pipeline is through the `run_pipeline` function:

```python
from superseeker.pipeline import run_pipeline

run_pipeline(
    patient='patient1',
    vcf_file='input_file.vcf',
    facets_dir='facets_output_directory',
    patient_sex='F',  # 'F' for female, 'M' for male
    restarts=100,     # Number of PyClone-VI restarts
    clusters=10,      # Number of clusters to use
    cn_neutral=False, # Whether to exclude variants in CNV regions
    cn_override=False,# Whether to override copy number information
    germfilter=True   # Whether to filter germline variants
)
```

### Advanced Usage

You can also run the pipeline in stages:

```python
from superseeker.pipeline import cluster_variants, run_superseeker

# First, cluster the variants
cluster_variants(
    patient='patient1',
    vcf_file='input_file.vcf',
    facets_dir='facets_output_directory',
    patient_sex='F'
)

# Then, run SuperSeeker on the clustered variants
run_superseeker(
    patient='patient1',
    vcf_file='input_file.vcf'
)
```

### Output Files

The pipeline generates several output files in a directory named `{patient}_superseeker_results`:

- `{patient}.pyclone_input.tsv`: Input file for PyClone-VI
- `{patient}.pyclone.clustered.tsv`: PyClone-VI clustering results
- `{patient}.somatic.clustered.vcf`: VCF file with cluster assignments
- `{patient}.cluster_lines.png`: Visualization of variant allele frequencies
- `{patient}.subclones.vcf`: SuperSeeker output with subclone assignments
- `{patient}.solutions.pdf`: Visual representation of evolutionary trees
- `{patient}.evolution.txt`: Summary of evolutionary patterns

## API Documentation

### Main Functions

#### `run_pipeline(patient, vcf_file, facets_dir='', patient_sex='', restarts=100, clusters=10, cn_neutral=False, cn_override=False, germfilter=True)`

Run the complete pipeline from VCF to evolutionary tree reconstruction.

#### `cluster_variants(patient, vcf_file, facets_dir='', patient_sex='', restarts=100, clusters=10, cn_neutral=False, cn_override=False, germfilter=True)`

Run only the variant clustering step using PyClone-VI.

#### `run_superseeker(patient, vcf_file, facets_dir='', patient_sex='', restarts=100, clusters=10, cn_neutral=False, cn_override=False, germfilter=True)`

Run only the SuperSeeker step on pre-clustered variants.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details.
