###############################################################
# Supplementary Figure S8
#
# aDNA damage sensitivity analysis
# Ancient vs modern Spain-Cant PCA displacement
#
# OldSalmo project
#
# Input:
#   data/sensitivity_damage/
#
# Output:
#   figures/Supplementary_Fig_S8.pdf
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
  "sensitivity_damage"
)


output_file <- here(
  "figures",
  "Supplementary_Fig_S8.pdf"
)



############################################################
# LOAD PCA PROJECTIONS
############################################################


dat <- read.table(
  file.path(
    input_dir,
    "all_damage_projections.evec.txt"
  ),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)



############################################################
# ASSIGN REGIONS
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
SRR5585856|SRR5585861"
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



dat <- dat %>%
mutate(
  region = assign_region(sample)
)



############################################################
# SPLIT BASE AND PROJECTIONS
############################################################


base <- dat %>%

filter(
  type == "base"
) %>%

select(
  sample,
  PC1,
  PC2,
  region
)



projection <- dat %>%

filter(
  type == "projection"
) %>%

select(
  sample,
  PC1,
  PC2,
  damage,
  replicate,
  region
) %>%

mutate(
  damage = sprintf("%.3f", damage)
)



############################################################
# CHECK DAMAGE LEVELS
############################################################


print(
unique(
projection$damage
)
)



############################################################
# SPAIN-CANT BASE CENTROID
############################################################


ref_centroid <- base %>%

filter(
  region == "Spain-Cant"
) %>%

summarise(

PC1 = mean(PC1, na.rm = TRUE),

PC2 = mean(PC2, na.rm = TRUE)

)



############################################################
# SIMULATED SPAIN-CANT CENTROIDS
############################################################


sim_centroids <- projection %>%

filter(
  region == "Spain-Cant"
) %>%

group_by(
  damage,
  replicate
) %>%

summarise(

PC1 = mean(PC1, na.rm = TRUE),

PC2 = mean(PC2, na.rm = TRUE),

.groups = "drop"

)



############################################################
# DISTANCES CAUSED BY DAMAGE
############################################################


sim_distances <- sim_centroids %>%

mutate(

distance =
sqrt(

(PC1 - ref_centroid$PC1)^2 +

(PC2 - ref_centroid$PC2)^2

)

)



############################################################
# LOAD ANCIENT PCA
############################################################


ancient <- read.table(

file.path(
  input_dir,
  "pca.evec"
),

header = FALSE,

skip = 1,

stringsAsFactors = FALSE

)



colnames(ancient) <- c(

"sample",

paste0("PC",1:10),

"group"

)



ancient <- ancient %>%

filter(
  group == "Ancient"
)



############################################################
# ANCIENT CENTROID
############################################################


ancient_centroid <- ancient %>%

summarise(

PC1 = mean(PC1, na.rm = TRUE),

PC2 = mean(PC2, na.rm = TRUE)

)



############################################################
# OBSERVED DISTANCE
############################################################


observed_distance <- sqrt(

(ancient_centroid$PC1 -
 ref_centroid$PC1)^2 +

(ancient_centroid$PC2 -
 ref_centroid$PC2)^2

)


print(observed_distance)

############################################################
# SUMMARY OF SIMULATED DAMAGE EFFECTS
############################################################


sim_summary <- sim_distances %>%

group_by(
  damage
) %>%

summarise(

mean_distance = mean(distance),

q95 = quantile(distance,0.95),

p_value = mean(distance >= observed_distance),

.groups = "drop"

)


print(sim_summary)



############################################################
# DAMAGE LEVELS
############################################################


damage_levels <- c(
"0.001",
"0.003",
"0.005",
"0.010"
)


damage_labels <- c(
"0.001"="1%",
"0.003"="3%",
"0.005"="5%",
"0.010"="10%"
)



sim_distances <- sim_distances %>%

mutate(

damage = factor(
  damage,
  levels = damage_levels
)

)



############################################################
# SELECT BEST REPLICATE
# Median displacement at 10% damage
############################################################


best_rep <- sim_distances %>%

filter(
  damage == "0.010"
) %>%

mutate(

dev =
abs(
distance -
median(distance, na.rm = TRUE)
)

) %>%

slice_min(
  dev,
  n = 1
) %>%

pull(replicate)



message(
"Best replicate (10% damage): ",
best_rep
)



############################################################
# PREPARE VECTORS
############################################################


base_vectors <- base %>%

mutate(

OriginalSample =
str_remove(
sample,
"_base$"
)

) %>%

select(
OriginalSample,
PC1,
PC2,
region
)



projection_vectors <- projection %>%

mutate(

OriginalSample =
str_remove(
sample,
"_d[0-9]+_r[0-9]+$"
)

)



############################################################
# FIX REGIONS
############################################################


projection_vectors <- projection_vectors %>%

mutate(

region =
case_when(

str_detect(OriginalSample,"^SRR217") &
is.na(region) ~
"UK-Scotland",

TRUE ~ region

)

)



############################################################
# JOIN PROJECTIONS WITH BASE
############################################################


proj_plot_all <- projection_vectors %>%

left_join(

base_vectors,

by = "OriginalSample",

suffix = c(
"_proj",
"_ref"
)

)



print(
table(
is.na(proj_plot_all$PC1_ref)
)
)



############################################################
# REMOVE UNKNOWN REGIONS
############################################################


proj_plot_all <- proj_plot_all %>%

filter(
!is.na(region_ref)
)



############################################################
# SELECT BEST REPLICATE DATA
############################################################


proj_plot <- proj_plot_all %>%

filter(
replicate == best_rep
)



############################################################
# INVERT PC AXES
############################################################


proj_plot <- proj_plot %>%

mutate(

PC1_proj = -PC1_proj,

PC2_proj = -PC2_proj,

PC1_ref = -PC1_ref,

PC2_ref = -PC2_ref

)



############################################################
# DAMAGE ORDER FOR PLOTTING
############################################################


damage_levels_plot <- c(
"0.010",
"0.005",
"0.003",
"0.001"
)



proj_plot <- proj_plot %>%

mutate(

damage = factor(
damage,
levels = damage_levels_plot
)

) %>%

arrange(
damage
)



############################################################
# COLOURS
############################################################


damage_cols <- c(

"0.001"="#56B4E9",

"0.003"="#009E73",

"0.005"="#E69F00",

"0.010"="#D55E00"

)



############################################################
# PANEL A
############################################################


p_left <- ggplot() +


geom_segment(

data = proj_plot,

aes(

x = PC1_ref,

y = PC2_ref,

xend = PC1_proj,

yend = PC2_proj,

colour = damage

),

alpha = 0.35,

linewidth = 0.35

) +


facet_wrap(
~region_ref
) +


scale_x_continuous(
breaks = c(-0.05,0)
) +


scale_y_continuous(
breaks = c(-0.05,0)
) +


scale_colour_manual(

values = damage_cols,

labels = damage_labels,

breaks = damage_levels,

name = "DNA damage"

) +


theme_bw(
base_size = 14
) +


theme(

strip.background =
element_rect(
fill = "grey90",
colour = "black",
linewidth = 0.5
),

panel.border =
element_rect(
colour = "black",
fill = NA,
linewidth = 0.5
),

legend.position = "bottom",

legend.justification = "center"

) +


labs(

x = "PC1",

y = "PC2"

)



############################################################
# PANEL B
############################################################


sim_distances <- sim_distances %>%

mutate(

damage = factor(

as.character(damage),

levels = damage_levels

)

)



p_right <- ggplot(

sim_distances,

aes(

x = distance,

fill = damage

)

) +


geom_density(

aes(
colour = damage
),

alpha = 0.35,

linewidth = 0.8

) +


geom_vline(

xintercept = observed_distance,

colour = "red",

linewidth = 1,

linetype = 2

) +


scale_fill_manual(

values = damage_cols,

labels = damage_labels,

guide = "none"

) +


scale_colour_manual(

values = damage_cols,

labels = damage_labels,

guide = "none"

) +


theme_bw(
base_size = 14
) +


theme(
legend.position = "none"
) +


labs(

x = "Centroid displacement",

y = "Density"

)



############################################################
# FINAL FIGURE
############################################################


final_plot <-

(p_left | p_right) +

plot_layout(

widths = c(1,1),

guides = "collect"

) +

plot_annotation(

tag_levels = "A"

) &

theme(

legend.position = "bottom",

legend.justification = "center",

plot.tag =
element_text(
size = 16,
face = "bold"
)

)



final_plot



############################################################
# SAVE FIGURE
############################################################


ggsave(

filename = output_file,

plot = final_plot,

width = 12,

height = 6,

units = "in",

device = cairo_pdf

)



############################################################
# SUMMARY TABLE
############################################################


summary_table <- sim_distances %>%

group_by(
damage
) %>%

summarise(

n = n(),

Mean_Euclidean_Distance =
mean(distance),

SD =
sd(distance),

Median =
median(distance),

q95 =
quantile(distance,0.95),

Maximum =
max(distance),

Empirical_P =
mean(distance >= observed_distance),

.groups = "drop"

) %>%

mutate(

Observed_Ancient_Distance =
observed_distance,

Observed_vs_Mean =
observed_distance /
Mean_Euclidean_Distance

)



summary_table



############################################################
# PUBLICATION TABLE
############################################################


summary_table_pub <- summary_table %>%

mutate(

Damage = c(
"1%",
"3%",
"5%",
"10%"
),


Mean_Euclidean_Distance =
sprintf(
"%.6f",
Mean_Euclidean_Distance
),


SD =
sprintf(
"%.6f",
SD
),


Median =
sprintf(
"%.6f",
Median
),


q95 =
sprintf(
"%.6f",
q95
),


Maximum =
sprintf(
"%.6f",
Maximum
),


Observed_Ancient_Distance =
sprintf(
"%.6f",
Observed_Ancient_Distance
),


Observed_vs_Mean =
sprintf(
"%.1f×",
Observed_vs_Mean
),


Empirical_P =
ifelse(
Empirical_P == 0,
"<0.01",
sprintf("%.3f",Empirical_P)
)

) %>%


select(

Damage,

n,

Mean_Euclidean_Distance,

SD,

Median,

q95,

Maximum,

Observed_Ancient_Distance,

Observed_vs_Mean,

Empirical_P

)



summary_table_pub



############################################################
# END
############################################################