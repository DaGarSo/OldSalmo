#!/usr/bin/env python3

import gzip
import random
import argparse


############################################################
# ARGUMENTS
############################################################

parser = argparse.ArgumentParser(
    description="Convert diploid VCF genotypes into pseudo-haploid genotypes"
)

parser.add_argument(
    "input_vcf",
    help="Input VCF file (gzipped)"
)

parser.add_argument(
    "output_vcf",
    help="Output pseudo-haploid VCF file"
)

parser.add_argument(
    "--seed",
    type=int,
    default=12345,
    help="Random seed for reproducibility (default: 12345)"
)

args = parser.parse_args()


input_vcf = args.input_vcf
output_vcf = args.output_vcf


############################################################
# RANDOM SEED
############################################################

random.seed(args.seed)



############################################################
# PSEUDO-HAPLOIDIZATION FUNCTION
############################################################

def pseudo_gt(gt):

    """
    Convert a diploid genotype into a pseudo-haploid genotype by randomly
    selecting one allele.

    Missing or non-diploid genotypes are converted to missing data.
    """

    alleles = gt.replace('|', '/').split('/')

    if len(alleles) == 2 and all(a != '.' for a in alleles):
        return random.choice(alleles)

    return "."



############################################################
# READ AND WRITE VCF
############################################################

with gzip.open(input_vcf, 'rt') as fin, open(output_vcf, 'w') as fout:

    for line in fin:

        if line.startswith('#'):

            fout.write(line)

        else:

            fields = line.rstrip().split('\t')

            format_fields = fields[8].split(':')
            gt_index = format_fields.index('GT')

            samples = fields[9:]

            new_samples = []

            for sample in samples:

                sample_fields = sample.split(':')

                original_gt = sample_fields[gt_index]

                hap_gt = pseudo_gt(original_gt)

                new_samples.append(hap_gt)


            fout.write(
                '\t'.join(fields[:8]) +
                '\tGT\t' +
                '\t'.join(new_samples) +
                '\n'
            )
