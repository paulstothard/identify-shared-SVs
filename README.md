# identify-shared-SVs

## Summary

The `identify-shared-SVs.py` script identifies structural variant (SV) sites shared between two VCF files. It outputs sites from the first file that are deemed to be shared with the second file. Sharing status is based on SV type, position overlap, and the percentage of matching genotypes. An optional `--not-shared-if-opposing-homozygotes` flag can be used to require that the sites have no opposing homozygote genotypes to be considered shared. The script accepts VCF and compressed VCF files (with or without an index) as input. The script can output compressed and uncompressed VCF files (based on the extension of the output file specified) and an optional tab-delimited file listing the identifiers (or identifying information if IDs are missing) of shared sites.

## Determining shared sites

The SV comparison is done by comparing each SV in `--vcf-file1` to those in `--vcf-file2`. For a site in `--vcf-file1` to be deemed shared with a site in `--vcf-file2`, the following conditions must be met:

- The `SVTYPE` values have to be the same.

- The position overlap has to be equal to or greater than the `--position-overlap-percent` option setting (default: 90%). Position overlap is determined by comparing the start and stop positions of the two variants under consideration. The overlap is calculated as the length where the two variants intersect. This overlap length is then divided by the total length that spans both variants to calculate the overlap percentage.

- The genotype overlap has to be equal to or greater than the `--genotype-overlap-percent` option setting (default: 90%). Genotype overlap is determined by comparing the genotypes of the two variants under consideration. The overlap is calculated as the number of matching genotypes divided by the total number of genotypes.

- If the `--not-shared-if-opposing-homozygotes` option is used, the two variants under consideration must not have opposing homozygote genotypes. For example, if one variant has a `0/0` genotype and the other has a `1/1` genotype, the variants are not considered shared.

## Output

The script can write the sites from `--vcf-file1` that are shared with `--vcf-file2` to a VCF file specified by the `--outfile` option. Note that a single site in the output VCF file can be shared with multiple sites in `--vcf-file2`. The `--shared-variants-file` option can be used to write a tab-delimited file showing the identifiers or identifying information of all `--vcf-file2` sites that are shared with `--vcf-file1` sites.

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

### Command-line options

- `--vcf-file1` (required): Path to the first input VCF file (e.g., `vcf1.vcf.gz`).
- `--vcf-file2` (required): Path to the second input VCF file (e.g., `vcf2.vcf.gz`).
- `--outfile`: Path to the output VCF file (e.g., `shared_sites.vcf.gz`).
- `--shared-variants-file`: Path to the optional tab-delimited output file providing identifiers of shared sites (e.g., `shared_variants.txt`).
- `--genotype-overlap-percent`: Minimum genotype overlap percentage to be classified as a shared site (default: 90).
- `--position-overlap-percent`: Minimum position overlap percentage to be classified as a shared site (default: 90).
- `--not-shared-if-opposing-homozygotes`: Sites aren't classified as shared if opposing homozygotes are detected (default: False). Use this option to enable this behavior.
- `--progress-count`: Frequency of progress updates (every N variants) (default: 1000).
- `--version`: Show the program version and exit.
