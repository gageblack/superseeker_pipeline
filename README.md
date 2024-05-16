
# SuperSeeker

SuperSeeker is a Python library for identifying subclonal evolution in cancer. It integrates various tools and scripts to process VCF files, cluster variants using PyClone-VI, and reconstruct subclonal evolutionary trees using SuperSeeker.

## Installation

### Step 1: Install External Tools

#### Installing PyClone-VI

PyClone-VI is used for clustering genetic variants. You can install it using Conda:

```bash
conda install -c bioconda pyclone-vi
```

#### Installing SuperSeeker

SuperSeeker is used for reconstructing the evolutionary history of subclones. You can install it from its source repository:

```bash
git clone https://github.com/yourusername/superseeker.git
cd superseeker
make install
```

### Step 2: Install Python Dependencies

Clone the SuperSeeker repository and navigate to the project directory:

```bash
git clone https://github.com/yourusername/SuperSeeker.git
cd SuperSeeker
```

Use the `requirements.txt` file to install the necessary Python packages:

```bash
pip install -r requirements.txt
```

### Step 3: Install the SuperSeeker Python Package

Finally, install the SuperSeeker Python package:

```bash
pip install -e .
```

This will install the package in "editable" mode, meaning changes to the code will be reflected without needing to reinstall the package. This is useful for development.

## Usage

To use the SuperSeeker pipeline, you can run the `run_pipeline` function from the `superseeker.pipeline` module. Here's an example:

```python
from superseeker.pipeline import run_pipeline

run_pipeline(
    patient='patient1',
    vcf_file='input_file.vcf',
    patient_sex='M',
    restarts=100,
    clusters=10,
    cn_neutral=True
)
```

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details.
