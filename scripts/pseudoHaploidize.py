import random
import gzip

# Input and output VCF files
# Modify these paths according to your input VCF and desired output filename
input_vcf = "pruned.filtered.vcf.gz"
output_vcf = "pruned.pseudoHapl.vcf"

# Random seed for reproducibility
# Change this value if a different random sampling is desired
random_seed = 12345
random.seed(random_seed)


def pseudo_gt(gt):
    """
    Convert a diploid genotype into a pseudo-haploid genotype by randomly
    selecting one of the two alleles.

    Missing or non-diploid genotypes are converted to missing data ('.').
    """
    alleles = gt.replace('|', '/').split('/')

    if len(alleles) == 2 and all(a != '.' for a in alleles):
        return random.choice(alleles)

    return '.'


# Read compressed input VCF and write pseudo-haploid VCF
with gzip.open(input_vcf, 'rt') as fin, open(output_vcf, 'w') as fout:

    for line in fin:

        # Preserve VCF header lines
        if line.startswith('#'):
            fout.write(line)

        else:
            fields = line.strip().split('\t')

            # Identify the position of the genotype (GT) field
            format_fields = fields[8].split(':')
            gt_index = format_fields.index('GT')

            samples = fields[9:]
            new_samples = []

            for sample in samples:
                sample_fields = sample.split(':')
                original_gt = sample_fields[gt_index]

                # Randomly select one allele to generate pseudo-haploid genotype
                hap_gt = pseudo_gt(original_gt)

                new_samples.append(hap_gt)

            # Write VCF line retaining only the GT field
            fout.write(
                '\t'.join(fields[:8]) +
                '\tGT\t' +
                '\t'.join(new_samples) +
                '\n'
            )
