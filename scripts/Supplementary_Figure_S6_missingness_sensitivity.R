###############################################################
# Supplementary Figure S6
#
# Missingness sensitivity analysis:
# Ancient vs modern Spain-Cant PCA displacement
#
# OldSalmo project
#
# Input:
#   data/sensitivity_missingness/
#
# Output:
#   figures/Supplementary_Fig_S6.pdf
#
###############################################################

##############################
# Libraries
##############################

library(tidyverse)
library(patchwork)
library(here)

##############################
# Paths
##############################

input_dir <- here(
  "data",
  "sensitivity_missingness"
)

output_file <- here(
  "figures",
  "Supplementary_Fig_S6.pdf"
)


############################################################
# PARAMETERS
############################################################

missing_levels <- c(
  "M0.40",
  "M0.50",
  "M0.75",
  "M0.80",
  "M0.90"
)


############################################################
# REGION ASSIGNMENT FUNCTION
############################################################

assign_region <- function(x){

case_when(

str_detect(x,"Teno") ~ "North-Norway",

str_detect(x,"Tornio") ~ "Finland",

str_detect(x,"^1_") ~ "Spain-Atl",
str_detect(x,"^2_") ~ "Spain-Atl",
str_detect(x,"^3_") ~ "Spain-Atl",

str_detect(x,"^4_") ~ "Spain-Cant",
str_detect(x,"^5_") ~ "Spain-Cant",
str_detect(x,"^6_") ~ "Spain-Cant",

str_detect(x,"^7_") ~ "Spain-Atl",
str_detect(x,"^8_") ~ "Spain-Atl",
str_detect(x,"^9_") ~ "Spain-Atl",
str_detect(x,"^10_") ~ "Spain-Atl",

str_detect(x,"^11_") ~ "UK-Scotland",
str_detect(x,"^12_") ~ "UK-Scotland",

str_detect(x,"^13_") ~ "Spain-Atl",
str_detect(x,"^14_") ~ "Spain-Atl",
str_detect(x,"^15_") ~ "Spain-Atl",
str_detect(x,"^16_") ~ "Spain-Atl",

str_detect(
x,
"SRR5585853|SRR5585854|SRR5585855|
SRR5585856|SRR5585857|SRR5585858|
SRR5585859|SRR5585860|SRR5585861"
) ~ "UK-Scotland",


str_detect(
x,
"SRR12588279|SRR12588280|SRR12588286|
SRR12588289|SRR12588290|SRR12588291|
SRR12588292|SRR2170984|SRR2174328|
SRR2174329|SRR2174330|SRR2174331|
SRR2174332|SRR2174333|SRR2174334|
SRR2174357|SRR2174358|SRR2174359|
SRR2174360|SRR2174361|SRR2174362|
SRR2174363|SRR2174364|SRR2174366|
SRR2174367|SRR2174368|SRR2174369|
SRR2174449|SRR2174450|SRR2174451|
SRR2174452|SRR3669078|SRR3669756|
SRR3669927"
) ~ "South-Norway",

TRUE ~ NA_character_

)

}


############################################################
# LOAD REFERENCE PCA
############################################################

ref <- read.table(
  file.path(
    input_dir,
    "M0.40_reference_coordinates.tsv"
  ),
  header = TRUE,
  stringsAsFactors = FALSE
)


ref <- ref %>%
  mutate(
    region = assign_region(sample)
  )


############################################################
# FUNCTION TO LOAD SIMULATIONS
############################################################

process_missingness <- function(level){

message("Processing ", level)


proj <- read.table(
  file.path(
    input_dir,
    paste0(level,"_projected_coordinates.tsv")
  ),
  header = TRUE,
  stringsAsFactors = FALSE
)


proj <- proj %>%
mutate(

OriginalSample =
str_remove(
sample,
paste0("_",level,"_R[0-9]+$")
),

replicate =
str_extract(
sample,
"(?<=_R)[0-9]+$"
) %>%
as.numeric(),

missingness = level

)


proj <- proj %>%
mutate(
region =
assign_region(OriginalSample)
)


proj <- proj %>%
left_join(
ref,
by=c(
"OriginalSample"="sample"
),
suffix=c(
"_proj",
"_ref"
)
)


return(proj)

}


############################################################
# LOAD ALL MISSINGNESS LEVELS
############################################################

all_dat <- map_dfr(
missing_levels,
process_missingness
)


all_dat <- all_dat %>%
mutate(
PC1_proj = -PC1_proj,
PC2_proj = -PC2_proj,
PC1_ref  = -PC1_ref,
PC2_ref  = -PC2_ref
)


############################################################
# Colours
############################################################

missing_cols <- c(
"M0.40" = "#56B4E9",
"M0.50" = "#009E73",
"M0.75" = "#E69F00",
"M0.80" = "#D55E00",
"M0.90" = "#CC79A7"
)


############################################################
# SIMULATED SPAIN-CANT CENTROIDS
############################################################

sim_centroids <- all_dat %>%

filter(
region_ref=="Spain-Cant"
) %>%

group_by(
missingness,
replicate
) %>%

summarise(

PC1 =
mean(PC1_proj,na.rm=TRUE),

PC2 =
mean(PC2_proj,na.rm=TRUE),

.groups="drop"

)


############################################################
# REFERENCE CENTROID
############################################################

ref_centroid <- all_dat %>%

filter(
region_ref=="Spain-Cant"
) %>%

summarise(

PC1 =
mean(PC1_ref,na.rm=TRUE),

PC2 =
mean(PC2_ref,na.rm=TRUE)

)


############################################################
# DISTANCES
############################################################

sim_distances <- sim_centroids %>%

mutate(

distance =
sqrt(
(PC1-ref_centroid$PC1)^2 +
(PC2-ref_centroid$PC2)^2
)

)


############################################################
# Ancient samples
############################################################

ancient <- read.table(
file.path(
input_dir,
"pca.evec"
),
header=FALSE,
skip=1,
stringsAsFactors=FALSE
)


colnames(ancient) <- c(
"sample",
paste0("PC",1:10),
"group"
)


ancient <- ancient %>%
filter(
group=="Ancient"
)


ancient_centroid <- ancient %>%

summarise(

PC1=mean(PC1),
PC2=mean(PC2)

)


observed_distance <- sqrt(

(ancient_centroid$PC1 -
ref_centroid$PC1)^2 +

(ancient_centroid$PC2 -
ref_centroid$PC2)^2

)


############################################################
# Plot
############################################################

all_dat <- all_dat %>%

mutate(
region_ref = case_when(

str_detect(OriginalSample,"^SRR217") &
is.na(region_ref) ~
"UK-Scotland",

TRUE ~ region_ref

)

) %>%

filter(
!is.na(region_ref)
)


best_rep <- sim_distances %>%

filter(
missingness=="M0.90"
) %>%

mutate(
dev=abs(distance-median(distance))
) %>%

slice_min(
dev,
n=1
) %>%

pull(replicate)



proj_plot <- all_dat %>%

filter(
replicate==best_rep
) %>%

mutate(
missingness=factor(
missingness,
levels=c(
"M0.90",
"M0.80",
"M0.75",
"M0.50",
"M0.40"
)
)

)


p_left <- ggplot() +

geom_segment(
data=proj_plot,
aes(
x=PC1_ref,
y=PC2_ref,
xend=PC1_proj,
yend=PC2_proj,
colour=missingness
),
alpha=0.30,
linewidth=0.35
)+

facet_wrap(~region_ref)+

scale_colour_manual(
values=missing_cols,
guide="none"
)+

theme_bw(
base_size=14
)+

labs(
x="PC1",
y="PC2"
)



p_right <- ggplot(
sim_distances,
aes(
x=distance,
fill=missingness
)
)+

geom_density(
aes(colour=missingness),
alpha=0.30,
linewidth=0.8
)+

geom_vline(
xintercept=observed_distance,
colour="red",
linewidth=1,
linetype=2
)+

scale_fill_manual(
values=missing_cols,
name="Missingness"
)+

scale_colour_manual(
values=missing_cols,
guide="none"
)+

theme_bw(
base_size=14
)+

theme(
legend.position="bottom"
)+

labs(
x="Centroid displacement",
y="Density"
)


final_plot <-

(p_left | p_right) +

plot_layout(
guides="collect"
)+

plot_annotation(
tag_levels="A"
)&

theme(
legend.position="bottom",
plot.tag=element_text(
size=16,
face="bold"
)
)


############################################################
# Save figure
############################################################

ggsave(
filename=output_file,
plot=final_plot,
width=12,
height=6,
units="in",
device=cairo_pdf
)


############################################################
# Summary table
############################################################

summary_table <- sim_distances %>%

group_by(missingness) %>%

summarise(

n=n(),

Mean_Euclidean_Distance=mean(distance),

SD=sd(distance),

Median=median(distance),

q95=quantile(distance,0.95),

Maximum=max(distance),

.groups="drop"

) %>%

mutate(

Observed_Ancient_Distance=observed_distance,

Observed_vs_Mean=
observed_distance/
Mean_Euclidean_Distance

)


summary_table