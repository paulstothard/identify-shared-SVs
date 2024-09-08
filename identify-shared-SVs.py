#!/usr/bin/env python3
"""
Author:
    Paul Stothard
"""

import pysam
import os
from pathlib import Path
import argparse


def get_1_based_start_and_end_positions(pysam_record):
    return pysam_record.pos, pysam_record.stop


def check_input_file(vcf_file):
    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")


def check_samples_order(vcf1, vcf2):
    samples1 = list(vcf1.header.samples)
    samples2 = list(vcf2.header.samples)

    if samples1 != samples2:
        print(f"Samples in first VCF file: {samples1}")
        print(f"Samples in second VCF file: {samples2}")
        raise ValueError(
            "Sample mismatch: VCF files have different samples or the samples are in a different order."
        )


def is_position_overlap(
    variant1_start,
    variant1_stop,
    variant2_start,
    variant2_stop,
    position_overlap_percent,
):
    overlap_start = max(variant1_start, variant2_start)
    overlap_end = min(variant1_stop, variant2_stop)
    overlap = max(0, overlap_end - overlap_start + 1)
    total_length = max(variant1_stop, variant2_stop) - min(
        variant1_start, variant2_start
    )
    return overlap / total_length >= position_overlap_percent / 100


def is_genotype_overlap(variant1, variant2, genotype_overlap_percent):
    matching_samples = [
        sample
        for sample in variant1.samples
        if variant1.samples[sample]["GT"] == variant2.samples[sample]["GT"]
    ]
    return (
        len(matching_samples) / len(variant1.samples) >= genotype_overlap_percent / 100,
        matching_samples,
    )


def is_homozygous(genotype):
    return (
        genotype is not None
        and genotype[0] is not None
        and genotype[1] is not None
        and genotype[0] == genotype[1]
    )


def has_opposing_homozygotes(variant1, variant2):
    return any(
        is_homozygous(variant1.samples[sample]["GT"])
        and is_homozygous(variant2.samples[sample]["GT"])
        and variant1.samples[sample]["GT"] != variant2.samples[sample]["GT"]
        for sample in variant1.samples
    )


def ensure_vcf_index(vcf_file):
    vcf_file = Path(vcf_file)
    if vcf_file.suffix == ".gz":
        compressed_vcf_file = vcf_file
    else:
        compressed_vcf_file = vcf_file.with_suffix(vcf_file.suffix + ".gz")

    index_file_tbi = Path(str(compressed_vcf_file) + ".tbi")
    index_file_csi = Path(str(compressed_vcf_file) + ".csi")

    def is_file_older(file1, file2):
        return (
            file1.exists()
            and file2.exists()
            and os.path.getmtime(file1) < os.path.getmtime(file2)
        )

    if vcf_file.suffix != ".gz":
        if not compressed_vcf_file.exists():
            pysam.tabix_index(
                str(vcf_file), preset="vcf", force=True, keep_original=True
            )
        else:
            if is_file_older(compressed_vcf_file, vcf_file):
                raise ValueError(
                    "The compressed VCF file is older than the uncompressed VCF file."
                )

    index_needs_creation = False
    if not index_file_tbi.exists() and not index_file_csi.exists():
        index_needs_creation = True
    elif index_file_tbi.exists() and is_file_older(index_file_tbi, compressed_vcf_file):
        index_needs_creation = True
    elif index_file_csi.exists() and is_file_older(index_file_csi, compressed_vcf_file):
        index_needs_creation = True

    if index_needs_creation:
        pysam.tabix_index(
            str(compressed_vcf_file), preset="vcf", force=True, keep_original=True
        )

    return compressed_vcf_file


def get_shared_SV_sites(
    vcf_file1,
    vcf_file2,
    genotype_overlap_percent=90,
    position_overlap_percent=90,
    outfile=None,
    not_shared_if_opposing_homozygotes=True,
    progress_count=1000,
    overwrite=True,
):
    check_input_file(vcf_file1)
    check_input_file(vcf_file2)

    if not overwrite and outfile and os.path.exists(outfile):
        raise FileExistsError(
            f"Output file '{outfile}' already exists. Use '--overwrite' to overwrite it."
        )

    if overwrite and outfile:
        if os.path.exists(outfile):
            os.remove(outfile)
        if os.path.exists(str(outfile) + ".tbi"):
            os.remove(str(outfile) + ".tbi")

    vcf_file1 = ensure_vcf_index(vcf_file1)
    vcf_file2 = ensure_vcf_index(vcf_file2)

    vcf1 = pysam.VariantFile(vcf_file1)
    vcf2 = pysam.VariantFile(vcf_file2)

    check_samples_order(vcf1, vcf2)

    compress_and_index = str(outfile).endswith(".vcf.gz")

    out_vcf = None
    if outfile:
        mode = "w" if not compress_and_index else "wz"
        out_vcf = pysam.VariantFile(str(outfile), mode, header=vcf1.header)

    variant_count = 0
    unique_variants = set()
    unique_variants_by_type = {}

    for idx, variant1 in enumerate(vcf1):
        if progress_count is not None and (idx + 1) % progress_count == 0:
            print(f"Processed {idx + 1} sites from {vcf_file1}")

        variant1_start, variant1_stop = get_1_based_start_and_end_positions(variant1)
        variant_key = (
            variant1.chrom,
            variant1_start,
            variant1.ref,
            tuple(variant1.alts),
        )

        variant_processed = False

        # Calculate variant size
        variant_size = variant1_stop - variant1_start + 1
        # Calculate scan distance dynamically based on position overlap percent
        scan_distance = int(variant_size * (100 - position_overlap_percent) / 100)

        # Adjust start and end positions based on the calculated scan distance
        start = max(0, variant1_start - scan_distance - 1)
        end = variant1_stop + scan_distance

        for variant2 in vcf2.fetch(variant1.chrom, start, end):
            variant2_start, variant2_stop = get_1_based_start_and_end_positions(
                variant2
            )

            if variant1.info["SVTYPE"] == variant2.info["SVTYPE"]:
                if is_position_overlap(
                    variant1_start,
                    variant1_stop,
                    variant2_start,
                    variant2_stop,
                    position_overlap_percent,
                ):
                    genotype_match, matching_samples = is_genotype_overlap(
                        variant1, variant2, genotype_overlap_percent
                    )
                    if genotype_match:
                        if (
                            not_shared_if_opposing_homozygotes
                            and has_opposing_homozygotes(variant1, variant2)
                        ):
                            continue

                        if variant_key not in unique_variants:
                            unique_variants.add(variant_key)
                            unique_variants_by_type[variant1.info["SVTYPE"]] = (
                                unique_variants_by_type.get(variant1.info["SVTYPE"], 0)
                                + 1
                            )

                        if not variant_processed:
                            if out_vcf:
                                out_vcf.write(variant1)
                            variant_count += 1
                            variant_processed = True

        vcf2.reset()

    if out_vcf:
        out_vcf.close()

        if compress_and_index:
            pysam.tabix_index(str(outfile), preset="vcf", force=True)

        print(
            f"{variant_count} variants from {vcf_file1} shared with {vcf_file2} written out to {outfile}."
        )

    print(f"Variants in {vcf_file1} shared with {vcf_file2} by SVTYPE:")
    for svtype, unique_count in unique_variants_by_type.items():
        print(f"{svtype}: Count = {unique_count}")

    print(f"Total sites: {len(unique_variants)}")


def main():
    parser = argparse.ArgumentParser(
        description="Writes sites from the first file that are deemed to be shared with the second file."
    )
    parser.add_argument(
        "--vcf-file1",
        required=True,
        type=str,
        help="Path to the first input VCF file (e.g., vcf1.vcf.gz)",
    )
    parser.add_argument(
        "--vcf-file2",
        required=True,
        type=str,
        help="Path to the second input VCF file (e.g., vcf2.vcf.gz)",
    )
    parser.add_argument(
        "--genotype-overlap-percent",
        type=int,
        default=90,
        help="Minimum genotype overlap percentage to be classified as a shared site (default: 90)",
    )
    parser.add_argument(
        "--position-overlap-percent",
        type=int,
        default=90,
        help="Minimum position overlap percentage to be classified as a shared site (default: 90)",
    )
    parser.add_argument(
        "--outfile",
        type=str,
        help="Path to the output VCF file (e.g., shared_sites.vcf.gz)",
    )
    parser.add_argument(
        "--not-shared-if-opposing-homozygotes",
        action="store_true",
        default=True,
        help="Sites aren't classified as shared if opposing homozygotes are detected (default: True)",
    )
    parser.add_argument(
        "--progress-count",
        type=int,
        default=1000,
        help="Frequency of progress updates (every N variants) (default: 1000)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite the output file if it exists (default: False)",
    )

    args = parser.parse_args()

    get_shared_SV_sites(
        vcf_file1=args.vcf_file1,
        vcf_file2=args.vcf_file2,
        genotype_overlap_percent=args.genotype_overlap_percent,
        position_overlap_percent=args.position_overlap_percent,
        outfile=args.outfile,
        not_shared_if_opposing_homozygotes=args.not_shared_if_opposing_homozygotes,
        progress_count=args.progress_count,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()