#!/usr/bin/env bash
python identify-shared-SVs.py --vcf-file1 sample-input/manta.test.vcf.gz \
--vcf-file2 sample-input/smoove.test.vcf.gz \
--outfile sample-output/shared-SVs.vcf.gz \
--shared-variants-file sample-output/shared-variants.txt \
--genotype-overlap-percent 90 \
--position-overlap-percent 90 \
--not-shared-if-opposing-homozygotes \
--progress-count 1000