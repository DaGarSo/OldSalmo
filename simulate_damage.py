#!/usr/bin/env python3

"""
Simulate ancient DNA damage patterns in a VCF file.

This script introduces C>T and G>A substitutions into reference alleles
at a user-defined rate, mimicking common ancient DNA deamination patterns.

Input and output files should be provided as arguments.
A random seed can be specified to ensure reproducibility.
"""

import gzip
import random
import argparse


def damage_gt(gt, ref, alt, damage_rate):
    """
    Introduce simulated deamination damage into diploid genotypes.

    C>T and G>A substitutions are simulated by randomly converting reference
    alleles (0) into alternative alleles (1) according to the specified
    damage rate.

    Non-standard genotypes are returned unchanged.
    """

    # Only process simple diploid genotypes
    if gt not in ["0/0", "0/1", "1/0", "1/1"]:
        return gt

    # C>T damage pattern
    if ref == "C" and alt == "T":
        target = "1"

    # G>A damage pattern
    elif ref == "G" and alt == "A":
        target = "1"

    else:
        return gt

    alleles = gt.replace("|", "/").split("/")

    damaged_alleles = []

    for allele in alleles:

        # Convert reference alleles into damaged alleles at the specified rate
        if allele == "0" and random.random() < damage_rate:
            damaged_alleles.append(target)
        else:
            damaged_alleles.append(allele)

    return "/".join(damaged_alleles)


def main():

    parser = argparse.ArgumentParser(
        description="Simulate ancient DNA damage patterns in a VCF file."
    )

    parser.add_argument(
        "input_vcf",
        help="Input VCF file (gzip compressed)"
    )

    parser.add_argument(
        "output_vcf",
        help="Output damaged VCF file (gzip compressed)"
    )

    parser.add_argument(
        "damage_rate",
        type=float,
        help="Probability of introducing damage (e.g. 0.05)"
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Random seed for reproducibility (default: 12345)"
    )

    args = parser.parse_args()

    # Set random seed for reproducibility
    random.seed(args.seed)

    with gzip.open(args.input_vcf, "rt") as fin, \
         gzip.open(args.output_vcf, "wt") as fout:

        for line in fin:

            # Preserve VCF header lines
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip().split("\t")

            chrom, pos, vid, ref, alt = fields[:5]

            samples = fields[9:]
            samples_new = []

            for sample in samples:

                genotype_fields = sample.split(":")
                gt = genotype_fields[0]

                new_gt = damage_gt(
                    gt,
                    ref,
                    alt,
                    args.damage_rate
                )

                # Preserve additional FORMAT fields
                remaining_fields = genotype_fields[1:]

                if remaining_fields:
                    new_sample = new_gt + ":" + ":".join(remaining_fields)
                else:
                    new_sample = new_gt

                samples_new.append(new_sample)

            fields[9:] = samples_new

            fout.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
