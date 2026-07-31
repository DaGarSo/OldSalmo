# OldSalmo

Custom scripts used for the manuscript:

**"What Modern Genomes Overlook: Ancient DNA Rewrites the Phylogeography of European Atlantic Salmon (Salmo salar)"**  
García-Souto et al.

---

## Overview

This repository contains custom scripts developed for the processing, simulation, quality control, and visualization of genomic data in the context of ancient DNA analyses of Atlantic salmon (*Salmo salar*).

The repository includes:

- Custom Python scripts for ancient DNA data processing and simulation.
- Shell scripts for genomic quality-control analyses.
- R scripts to reproduce all main and supplementary figures presented in the manuscript.
- Code used for sensitivity analyses evaluating the impact of missing data and ancient DNA damage on population genomic inference.

All scripts are provided to facilitate reproducibility of the analyses described in the manuscript.

---

# Custom Scripts

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
```

## Reproducibility

Pseudo-haploidization involves random allele sampling. To ensure reproducibility, a fixed random seed is defined:

```python
random_seed = 12345
random.seed(random_seed)
```

Changing the seed value will generate a different pseudo-haploid representation.

## Usage

Run the script with:

```bash
python pseudoHaploidize.py
```

---

## 2. Simulation of ancient DNA damage

### `simulate_damage.py`

This script introduces simulated ancient DNA damage patterns into a VCF file.

The simulation models common post-mortem deamination patterns observed in ancient DNA datasets:

- C→T substitutions
- G→A substitutions

Reference alleles are randomly converted into alternative alleles according to a user-defined damage probability. This allows testing the potential impact of ancient DNA damage on downstream population genomic analyses.

Only simple diploid genotypes are modified:

```text
0/0
0/1
1/0
1/1
```

Other genotype configurations are retained unchanged.

---

## Input and output

The script requires three mandatory arguments:

| Argument | Description |
|---|---|
| Input VCF | Input VCF file (gzip compressed) |
| Output VCF | Output damaged VCF file (gzip compressed) |
| Damage rate | Probability of introducing simulated damage at susceptible sites |

---

## Usage example

```bash
python simulate_damage.py \
    input.vcf.gz \
    damaged.vcf.gz \
    0.05
```

where:

```text
0.05 = 5% probability of introducing damage
```

at susceptible sites.

---

## Reproducibility

Because damage introduction is based on random sampling, a random seed can be specified using the `--seed` option:

```bash
python simulate_damage.py \
    input.vcf.gz \
    damaged.vcf.gz \
    0.05 \
    --seed 12345
```

## 3. Ti/Tv calculation

### `scripts/calculate_ti_tv_by_sample.sh`

Script used to calculate transition/transversion ratios for each individual from a VCF file.

The script extracts:

REF and ALT alleles.
Individual genotypes.
Counts of nucleotide substitutions.
Transition (Ti) and transversion (Tv) rates.
Ti/Tv ratio per sample.

---

## Usage example

```bash scripts/calculate_ti_tv_by_sample.sh input.vcf.gz output_prefix```

---

## Outputs:

output_prefix.snp_gt_matrix.tsv
output_prefix.ti_tv_by_sample.tsv
