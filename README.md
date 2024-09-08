# identify-shared-SVs

The `identify-shared-SVs.py` script identifies structural variant (SV) sites shared between two VCF files. It outputs sites from the first file that are deemed to be shared with the second file.

## Requirements

- [Python 3](https://www.python.org/) version 3.8 or higher
- [pysam](https://pysam.readthedocs.io/en/latest/) version 0.22.0

## Usage

First create a virtual environment and install the required packages:

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Then run the script, for example:

```bash
python identify-shared-SVs.py --vcf-file1 sample-input/manta.test.vcf.gz --vcf-file2 sample-input/smoove.test.vcf.gz --outfile sample-output/shared-SVs.vcf.gz --overwrite
```

To see available options, run:

```bash
python identify-shared-SVs.py --help
```
