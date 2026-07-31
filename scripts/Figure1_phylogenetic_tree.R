###############################################################
# Figure 1. Mitochondrial phylogenetic tree
#
# OldSalmo project
#
# Description:
# Generates the mitochondrial phylogenetic tree shown
# in Figure 1 from the RAxML tree with bootstrap support.
#
# Input:
#   data/phylogeny/salmon_mt_rooted.raxml.support
#
# Output:
#   figures/Figure1_phylogenetic_tree.pdf
#
###############################################################

##############################
# Libraries
##############################

library(ape)
library(ggtree)
library(tidyverse)
library(here)
library(grid)

##############################
# Input / output paths
##############################

tree_file <- here(
  "data",
  "phylogeny",
  "salmon_mt_rooted.raxml.support"
)

output_file <- here(
  "figures",
  "Figure1_phylogenetic_tree.pdf"
)

##############################
# Read tree
##############################

tree <- read.tree(tree_file)

##############################
# Clean bootstrap labels
##############################

tree$node.label[tree$node.label == ""] <- NA

##############################
# Haplotypes from modern Asón
##############################

iberian_haps <- c(
  "Hap016",
  "Hap023",
  "Hap092",
  "Hap095",
  "Hap096"
)

##############################
# Metadata
##############################

tip_df <- tibble(
  label = tree$tip.label
) %>%
  mutate(
    Category = case_when(

      label %in% iberian_haps ~
        "Verspoor 2012 - Modern Asón",

      str_detect(label, "^Hap") ~
        "Verspoor 2012 - Europe",

      str_detect(label, "^S1|^S2") ~
        "Ancient Asón",

      str_detect(
        label,
        regex("trutta|trout", ignore_case = TRUE)
      ) ~
        "Salmo trutta",

      str_detect(
        label,
        regex("SRR2821|JQ390055", ignore_case = TRUE)
      ) ~
        "American",

      TRUE ~
        "European"
    )
  )

##############################
# Colour palette
##############################

tip_colors <- c(
  "Verspoor 2012 - Europe"      = "grey55",
  "Verspoor 2012 - Modern Asón" = "#E69F00",
  "European"                    = "#4DAF4A",
  "American"                    = "#377EB8",
  "Ancient Asón"                = "#D73027",
  "Salmo trutta"                = "black"
)

##############################
# Build tree
##############################

p <- ggtree(
  tree,
  layout = "circular",
  branch.length = "none"
) %<+% tip_df

##############################
# Bootstrap labels
##############################

boot_df <- p$data %>%
  filter(
    !isTip,
    !is.na(label)
  ) %>%
  mutate(
    bootstrap = suppressWarnings(
      as.numeric(label)
    )
  ) %>%
  filter(
    !is.na(bootstrap),
    bootstrap >= 50
  )

##############################
# Plot
##############################

figure1 <-

  p +

  geom_text(
    data = boot_df,
    aes(
      x = x,
      y = y,
      label = bootstrap
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 3.2
  ) +

  geom_tippoint(
    aes(colour = Category),
    size = 5,
    show.legend = TRUE
  ) +

  geom_tiplab(
    aes(colour = Category),
    size = 4,
    offset = 1.0,
    show.legend = FALSE
  ) +

  scale_colour_manual(
    values = tip_colors,
    name = NULL
  ) +

  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "cm"),
    plot.margin = margin(15, 15, 15, 15)
  )

##############################
# Save figure
##############################

ggsave(
  filename = output_file,
  plot = figure1,
  width = 14,
  height = 14,
  units = "in",
  dpi = 600,
  device = cairo_pdf,
  bg = "white"
)

##############################
# Display
##############################

figure1
