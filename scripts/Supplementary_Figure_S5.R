###############################################################
# Supplementary Figure S5
#
# ADMIXTURE ancestry proportions (K = 4)
# Including American populations
#
# OldSalmo project
#
# Input:
#   data/admixture/MERGED.America.clean.4.Q
#
# Output:
#   figures/Supplementary_Fig_S5.pdf
#
###############################################################


##############################
# Libraries
##############################

library(tidyr)
library(dplyr)
library(ggplot2)
library(here)



##############################
# Paths
##############################

input_file <- here(
  "data",
  "admixture",
  "MERGED.America.clean.4.Q"
)


output_file <- here(
  "figures",
  "Supplementary_Fig_S5.pdf"
)



############################################################
# LOAD DATA
############################################################


df <- read.table(

  input_file,

  header = FALSE,

  sep = "",

  fill = TRUE,

  stringsAsFactors = FALSE

)



colnames(df) <- c(

"ID",

"Population",

paste0(
"Q",
1:(ncol(df)-2)
)

)



############################################################
# FILTER POPULATIONS
############################################################


df <- df %>%

filter(

!Population %in%
c(
"Russia",
"Sweden"
)

)



############################################################
# RENAME POPULATIONS
############################################################


df <- df %>%

mutate(

Population = case_when(

Population == "Norway" &
grepl("^S", ID)
~
"NW-South",


Population == "Norway"
~
"NW-North",


TRUE
~
Population

)

)



############################################################
# POPULATION ORDER
############################################################


pop_order <- c(

"America",

"Finland",

"NW-North",

"NW-South",

"UK-Scotland",

"Oldsalmo",

"Spain-Cant",

"Spain-Atl"

)



df$Population <- factor(

df$Population,

levels = pop_order

)



############################################################
# BALANCED SUBSAMPLING
############################################################


max_samples <- 100


Q_sub <- df %>%

group_by(
Population
) %>%

group_modify(

~slice_sample(
.x,
n = min(
max_samples,
nrow(.x)
)
)

) %>%

ungroup()



############################################################
# LONG FORMAT
############################################################


Q_long <- Q_sub %>%

pivot_longer(

cols = starts_with("Q"),

names_to = "Ancestry",

values_to = "Proportion"

)



############################################################
# SAMPLE ORDER
############################################################


Q_long <- Q_long %>%

mutate(

ID = factor(

ID,

levels =
Q_long %>%

distinct(
ID,
Population
) %>%

arrange(
Population
) %>%

pull(ID)

)

)



############################################################
# REMOVE UNKNOWN POPULATIONS
############################################################


Q_long <- Q_long %>%

filter(

!Population %in%
c(
"Russia",
"Sweden"
)

) %>%

filter(
complete.cases(.)
)



############################################################
# PLOT
############################################################


q4plotsub <- ggplot(

Q_long,

aes(

x = ID,

y = Proportion,

fill = Ancestry

)

) +


geom_bar(

stat = "identity",

width = 1

) +


facet_grid(

. ~ Population,

scales = "free_x",

space = "free_x"

) +


scale_fill_brewer(

palette = "Set2"

) +


labs(

x = "Samples",

y = "Ancestry proportion (K = 4)",

fill = "Component"

) +


theme_minimal(

base_size = 12

) +


theme(

axis.text.x =
element_blank(),

axis.ticks.x =
element_blank(),

panel.spacing =
unit(
0.1,
"lines"
),

panel.grid =
element_blank(),

strip.text =
element_text(
size = 11,
face = "bold"
),

legend.title =
element_blank(),

legend.text =
element_text(
size = 10
)

)



q4plotsub



############################################################
# EXPORT
############################################################


ggsave(

filename = output_file,

plot = q4plotsub,

width = 14,

height = 5,

units = "in",

device = cairo_pdf

)