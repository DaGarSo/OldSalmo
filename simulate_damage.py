#!/usr/bin/env python3

import gzip
import random
import sys

############################################################
# ARGUMENTS
############################################################

input_vcf = sys.argv[1]
output_vcf = sys.argv[2]
damage_rate = float(sys.argv[3])

# Optional random seed
if len(sys.argv) > 4:
    random.seed(int(sys.argv[4]))

############################################################
# FUNCTION
############################################################

def damage_gt(gt, ref, alt):

    # Only haploid genotypes
    if gt not in ["0", "1"]:
        return gt

    ############################################################
    # C <-> T
    ############################################################

    if ref == "C" and alt == "T":

        # Reference C becomes ALT T
        if gt == "0" and random.random() < damage_rate:
            return "1"

        return gt

    elif ref == "T" and alt == "C":

        # ALT C becomes reference T
        if gt == "1" and random.random() < damage_rate:
            return "0"

        return gt

    ############################################################
    # G <-> A
    ############################################################

    elif ref == "G" and alt == "A":

        # Reference G becomes ALT A
        if gt == "0" and random.random() < damage_rate:
            return "1"

        return gt

    elif ref == "A" and alt == "G":

        # ALT G becomes reference A
        if gt == "1" and random.random() < damage_rate:
            return "0"

        return gt

    ############################################################
    # Other SNPs
    ############################################################

    return gt

############################################################
# PROCESS VCF
############################################################

# Input can be gzipped (BGZF is also readable by gzip)
with gzip.open(input_vcf, "rt") as fin, open(output_vcf, "w") as fout:

    for line in fin:

        if line.startswith("#"):
            fout.write(line)
            continue

        fields = line.rstrip().split("\t")

        ref = fields[3]
        alt = fields[4]

        # Skip multiallelic sites
        if "," in alt:
            fout.write(line)
            continue

        samples_new = []

        for sample in fields[9:]:

            parts = sample.split(":")

            gt = parts[0]

            new_gt = damage_gt(gt, ref, alt)

            parts[0] = new_gt

            samples_new.append(":".join(parts))

        fields[9:] = samples_new

        fout.write("\t".join(fields) + "\n")
