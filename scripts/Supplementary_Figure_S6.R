###############################################################
# Supplementary Figure S6
#
# ADMIXTURE ancestry proportions (K = 3)
#
# OldSalmo project
#
# Input:
#   data/admixture/
#
# Output:
#   figures/Supplementary_Fig_S6.pdf
#
###############################################################


##############################
# Libraries
##############################

library(tidyverse)
library(here)



##############################
# Paths
##############################

q_file <- here(
  "data",
  "admixture",
  "MERGED.noAmerica.clean.3.Q"
)


fam_file <- here(
  "data",
  "admixture",
  "MERGED.noAmerica.clean.fam"
)


metadata_file <- here(
  "data",
  "admixture",
  "metadata.txt"
)


output_file <- here(
  "figures",
  "Supplementary_Fig_S6.pdf"
)



############################################################
# LOAD DATA
############################################################


Q <- read.table(
  q_file,
  header = FALSE
)


FAM <- read.table(
  fam_file,
  header = FALSE
)


Q <- cbind(
  FAM,
  Q
)


Q$V1 <- sub(
  " .*",
  "",
  Q$V1
)


Q <- Q[,c(
  1,
  7:ncol(Q)
)]


colnames(Q) <- c(
  "Sample",
  "PC1",
  "PC2",
  "PC3"
)



############################################################
# METADATA
############################################################


metadata <- read.table(

  metadata_file,

  header = TRUE,

  sep = "\t",

  stringsAsFactors = FALSE,

  fill = TRUE

)


Q <- Q %>%

left_join(
  metadata,
  by = "Sample"
)



############################################################
# FILTERING AND POPULATION ASSIGNMENT
############################################################


Q$Area[
  grep(
    "Tornio",
    Q$Sample
  )
] <- "Finland"



Q <- Q %>%

mutate(

Area = if_else(
  is.na(Area),
  "Oldsalmo",
  Area
)

) %>%

filter(

!Area %in%
c(
"North-America",
"Russia",
"Sweden"
)

) %>%

filter(

!is.na(Area),

!is.na(PC1),

!is.na(PC2),

!is.na(PC3)

)



############################################################
# SUBSAMPLING
############################################################


max_samples <- 100


Q_sub <- Q %>%

group_by(
  Area
) %>%

group_modify(

~slice_sample(
.x,
n=min(
max_samples,
nrow(.x)
)
)

) %>%

ungroup()



Q_sub <- Q_sub %>%

filter(
Area != "Norway"
)



############################################################
# LONG FORMAT
############################################################


Q_long <- Q_sub %>%

pivot_longer(

cols = c(
PC1,
PC2,
PC3
),

names_to = "Ancestry",

values_to = "Proportion"

)



Q_long <- Q_long %>%

mutate(

Sample = factor(

Sample,

levels =
Q_long %>%

distinct(
Sample,
Area
) %>%

arrange(
Area
) %>%

pull(
Sample
)

)

)



Q_long$Area <- factor(

Q_long$Area,

levels = c(

"Finland",

"North-Norway",

"South-Norway",

"UK-Scotland",

"Oldsalmo",

"Spain-Cant",

"Spain-Atl"

)

)



############################################################
# LABELS
############################################################


area_labels <- c(

"Finland" =
"Finland",

"North-Norway" =
"North Norway",

"South-Norway" =
"South Norway",

"UK-Scotland" =
"UK Scotland",

"Oldsalmo" =
"Ancient Asón",

"Spain-Cant" =
"Spain Cantabrian",

"Spain-Atl" =
"Spain Atlantic"

)



############################################################
# PLOT
############################################################


q3plotsub <- ggplot(

Q_long,

aes(

x = Sample,

y = Proportion,

fill = Ancestry

)

) +

geom_bar(

stat = "identity",

width = 1

) +


facet_grid(

. ~ Area,

scales = "free_x",

space = "free_x",

labeller =
labeller(
Area = area_labels
)

) +


scale_fill_brewer(

palette = "Set2"

) +


scale_y_continuous(

expand = c(0,0),

breaks = seq(
0,
1,
0.25
)

) +


labs(

x = "Individuals",

y = "Ancestry proportion (K = 3)",

fill = "Component"

) +


theme_classic(

base_size = 11

) +


theme(

axis.text.x =
element_blank(),

axis.ticks.x =
element_blank(),

axis.title.x =
element_text(
face="bold",
margin=margin(t=10)
),

axis.title.y =
element_text(
face="bold",
margin=margin(r=10)
),

strip.text.x =
element_text(
face="bold",
size=10
),

strip.background =
element_blank(),

panel.spacing =
unit(
0.15,
"lines"
),

legend.position =
"right",

legend.title =
element_text(
face="bold"
)

)



q3plotsub



############################################################
# EXPORT
############################################################


ggsave(

filename = output_file,

plot = q3plotsub,

width = 10,

height = 2.5,

dpi = 300

)