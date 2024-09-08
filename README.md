# identify-shared-SVs

The `identify-shared-SVs.py` script identifies structural variant (SV) sites shared between two VCF files. It outputs sites from the first file that are deemed to be shared with the second file.

## Requirements

- [Python 3](https://www.python.org/) version 3.8 or higher
- [pysam](https://pysam.readthedocs.io/en/latest/) version 0.22.0

## Installation

Clone the repository:

```bash
git clone git@github.com:paulstothard/identify-shared-SVs.git
```

Create a virtual environment and install the required packages:

```bash
cd identify-shared-SVs
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Test the script:

```bash
python identify-shared-SVs.py --vcf-file1 sample-input/manta.test.vcf.gz \
--vcf-file2 sample-input/smoove.test.vcf.gz \
--outfile sample-output/shared-SVs.vcf.gz \
--shared-variants-file sample-output/shared-variants.txt \
--genotype-overlap-percent 90 \
--position-overlap-percent 90 \
--not-shared-if-opposing-homozygotes \
--progress-count 1000
```

### Command-line Options

- `--vcf-file1` (required): Path to the first input VCF file (e.g., `vcf1.vcf.gz`).
- `--vcf-file2` (required): Path to the second input VCF file (e.g., `vcf2.vcf.gz`).
- `--outfile`: Path to the output VCF file (e.g., `shared_sites.vcf.gz`).
- `--shared-variants-file`: Path to the optional tab-delimited output file providing identifiers of shared sites (e.g., `shared_variants.txt`).
- `--genotype-overlap-percent`: Minimum genotype overlap percentage to be classified as a shared site (default: 90).
- `--position-overlap-percent`: Minimum position overlap percentage to be classified as a shared site (default: 90).
- `--not-shared-if-opposing-homozygotes`: Sites aren't classified as shared if opposing homozygotes are detected (default: False). Use this option to enable this behavior.
- `--progress-count`: Frequency of progress updates (every N variants) (default: 1000).
- `--version`: Show the program version and exit.
