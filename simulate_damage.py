#!/usr/bin/env python3

import gzip
import random
import argparse
import sys

############################################################
# ARGUMENTS
############################################################

parser = argparse.ArgumentParser(
    description="Simulate ancient DNA damage in haploid VCF genotypes"
)

parser.add_argument(
    "input_vcf",
    help="Input VCF file (gzipped or uncompressed)"
)

parser.add_argument(
    "output_vcf",
    help="Output VCF file (use .gz for gzip compression)"
)

parser.add_argument(
    "damage_rate",
    type=float,
    help="Damage probability (e.g. 0.05)"
)

parser.add_argument(
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility"
)

args = parser.parse_args()


input_vcf = args.input_vcf
output_vcf = args.output_vcf
damage_rate = args.damage_rate


if args.seed is not None:
    random.seed(args.seed)


############################################################
# FUNCTION
############################################################

def damage_gt(gt, ref, alt):

    # Only haploid genotypes
    if gt not in ["0", "1"]:
        return gt

    ########################################################
    # C -> T and T -> C
    ########################################################

    if ref == "C" and alt == "T":

        # Reference C damaged into T
        if gt == "0" and random.random() < damage_rate:
            return "1"

        return gt


    elif ref == "T" and alt == "C":

        # ALT C reverts to T
        if gt == "1" and random.random() < damage_rate:
            return "0"

        return gt


    ########################################################
    # G -> A and A -> G
    ########################################################

    elif ref == "G" and alt == "A":

        # Reference G damaged into A
        if gt == "0" and random.random() < damage_rate:
            return "1"

        return gt


    elif ref == "A" and alt == "G":

        # ALT G reverts to A
        if gt == "1" and random.random() < damage_rate:
            return "0"

        return gt


    ########################################################
    # Other SNPs unchanged
    ########################################################

    return gt



############################################################
# OPEN INPUT / OUTPUT
############################################################

def open_input(filename):

    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


def open_output(filename):

    if filename.endswith(".gz"):
        return gzip.open(filename, "wt")
    else:
        return open(filename, "w")



############################################################
# PROCESS VCF
############################################################

with open_input(input_vcf) as fin, open_output(output_vcf) as fout:

    for line in fin:

        # Keep headers
        if line.startswith("#"):
            fout.write(line)
            continue


        fields = line.rstrip("\n").split("\t")

        ref = fields[3]
        alt = fields[4]


        # Skip multiallelic variants
        if "," in alt:
            fout.write(line)
            continue


        samples_new = []


        for sample in fields[9:]:

            parts = sample.split(":")

            gt = parts[0]

            new_gt = damage_gt(
                gt,
                ref,
                alt
            )

            parts[0] = new_gt

            samples_new.append(":".join(parts))


        fields[9:] = samples_new

        fout.write("\t".join(fields) + "\n")


print(
    f"Finished. Damage rate={damage_rate}, seed={args.seed}",
    file=sys.stderr
)
