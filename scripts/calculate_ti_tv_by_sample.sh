#!/bin/bash

###############################################################
# OldSalmo project
#
# Calculate transition/transversion (Ti/Tv) ratios
# by individual from a VCF file
#
# Input:
#   VCF file
#
# Output:
#   <prefix>.snp_gt_matrix.tsv
#   <prefix>.ti_tv_by_sample.tsv
#
# Usage:
#   bash calculate_ti_tv_by_sample.sh input.vcf.gz output_prefix
#
###############################################################


set -euo pipefail


###############################################################
# Arguments
###############################################################


if [ $# -lt 1 ]; then

    echo "Usage:"
    echo "bash calculate_ti_tv_by_sample.sh input.vcf.gz [output_prefix]"
    exit 1

fi


VCF=$1

OUTPREFIX=${2:-ti_tv_results}



###############################################################
# Extract REF, ALT and genotypes
###############################################################


bcftools query \
-f '%REF\t%ALT[\t%GT]\n' \
"$VCF" \
> "${OUTPREFIX}.snp_gt_matrix.tsv"



###############################################################
# Count substitutions per sample
###############################################################


awk '

BEGIN {

split(
"AG,GA,CT,TC,AC,CA,AT,TA,GC,CG,GT,TG",
mut_list,
","
);

}


NR == 1 {

n_samples = NF - 2;

for (i in mut_list) {

m = mut_list[i];

for (s = 1; s <= n_samples; s++) {

count[m][s]=0;

}

}

next;

}



{

ref=$1;

alt=$2;

mut=ref alt;


for (i=3; i<=NF; i++) {


gt=$i;

idx=i-2;


if (

gt!="./." &&
gt!="." &&
gt!="0/0" &&
gt!="0|0"

)

{

count[mut][idx]++;

}

}

}



END {


print "Sample\tAG\tGA\tCT\tTC\tAC\tCA\tAT\tTA\tGC\tCG\tGT\tTG\tTi\tTv\tTi/Tv";


for (s=1; s<=n_samples; s++) {


ti =
count["AG"][s] +
count["GA"][s] +
count["CT"][s] +
count["TC"][s];


tv=0;


for (i in mut_list) {


m=mut_list[i];


if (

m!="AG" &&
m!="GA" &&
m!="CT" &&
m!="TC"

)

tv += count[m][s];


}


ratio=(tv==0) ? "Inf" : ti/tv;


printf "Sample%d",s;


for (i in mut_list)

printf "\t%d",count[mut_list[i]][s];


printf "\t%d\t%d\t%.2f\n",ti,tv,ratio;


}


}

' "${OUTPREFIX}.snp_gt_matrix.tsv" \
> "${OUTPREFIX}.ti_tv_by_sample.tsv"



echo "Ti/Tv calculation completed"
echo ""
echo "Generated files:"
echo "- ${OUTPREFIX}.snp_gt_matrix.tsv"
echo "- ${OUTPREFIX}.ti_tv_by_sample.tsv"