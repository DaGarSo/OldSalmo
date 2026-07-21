# OldSalmo

Custom scripts used for the manuscript:

**"What Modern Genomes Overlook: Ancient DNA Rewrites the Phylogeography of European Atlantic Salmon (*Salmo salar*)"**  
García-Souto et al.

This repository contains custom Python scripts developed for processing and simulating genomic data in the context of ancient DNA analyses of Atlantic salmon.

---

# Scripts

## 1. Pseudo-haploidization of VCF files

### `pseudoHaploidize.py`

This script converts diploid genotypes from a VCF file into pseudo-haploid genotypes by randomly sampling one allele from each heterozygous genotype.

Pseudo-haploidization is commonly applied to low-coverage ancient DNA datasets to reduce biases associated with genotype uncertainty and differences in sequencing depth between ancient and modern samples.

The script:

- Retains the original VCF header information.
- Converts heterozygous genotypes into randomly sampled homozygous states.
- Keeps only the genotype (`GT`) field in the output VCF.

---

## Input and output

Before running the script, modify the input and output filenames in the script:

```python
input_vcf = "pruned.filtered.vcf.gz"
output_vcf = "pruned.pseudoHapl.vcf"
Reproducibility

Pseudo-haploidization involves random allele sampling. To ensure reproducibility, a fixed random seed is defined:

random_seed = 12345
random.seed(random_seed)

Changing the seed value will generate a different pseudo-haploid representation.

Usage

Run the script with:

python pseudoHaploidize.py
2. Ancient DNA damage simulation
simulate_damage.py

## 2. Simulation of ancient DNA damage

### `simulate_damage.py`

This script introduces simulated ancient DNA damage patterns into a VCF file.

The simulation models common post-mortem deamination patterns observed in ancient DNA datasets:

C→T substitutions
G→A substitutions

Reference alleles are randomly converted into alternative alleles according to a user-defined damage probability. This allows testing the potential impact of ancient DNA damage on downstream population genomic analyses.

Only simple diploid genotypes are modified:

0/0
0/1
1/0
1/1

Other genotype configurations are retained unchanged.

Input and output

The script requires three mandatory arguments:

Argument	Description
Input VCF	Input VCF file (gzip compressed)
Output VCF	Output damaged VCF file (gzip compressed)
Damage rate	Probability of introducing simulated damage at susceptible sites
Usage example
python simulate_damage.py \
    input.vcf.gz \
    damaged.vcf.gz \
    0.05

where:

0.05 = 5% probability of introducing damage

at susceptible sites.

Reproducibility

Because damage introduction is based on random sampling, a random seed can be specified using the --seed option:

python simulate_damage.py \
    input.vcf.gz \
    damaged.vcf.gz \
    0.05 \
    --seed 12345

Using the same seed will generate identical simulated damage patterns.

```markdown
# Requirements

- Python ≥ 3.8
- Standard Python libraries:
  - `random`
  - `gzip`
  - `argparse`
