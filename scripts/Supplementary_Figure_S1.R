###############################################################
# Supplementary Figure S1
#
# Map Europe only
#
# OldSalmo project
#
# Input:
#   data/sampling_map/Loc.txt
#
# Output:
#   figures/Supplementary_Fig_S1.pdf
#
###############################################################


##############################
# Libraries
##############################

library(sf)
library(ggplot2)
library(rnaturalearth)
library(ggrepel)
library(ggspatial)
library(grid)
library(here)



##############################
# Paths
##############################

input_file <- here(
  "data",
  "sampling_map",
  "Loc.txt"
)


output_file <- here(
  "figures",
  "Supplementary_Fig_S1.pdf"
)



############################################################
# DATA
############################################################


sites <- read.table(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE
)


sites$River <- gsub(
  "_",
  " ",
  sites$River
)



############################################################
# REMOVE NON-EUROPEAN SITES
############################################################


sites <- subset(
  sites,
  !Country %in% c(
    "USA",
    "Canada"
  )
)



############################################################
# RENUMBER SITES
############################################################


sites$ID <- seq_len(
  nrow(sites)
)



sites_sf <- st_as_sf(
  sites,
  coords = c(
    "Longitude",
    "Latitude"
  ),
  crs = 4326
)



############################################################
# MAP DATA
############################################################


world <- ne_countries(
  scale = 10,
  returnclass = "sf"
)



############################################################
# COLOURS
############################################################


land.col   <- "#F2F2F2"
water.col  <- "white"
border.col <- "#E8E6E1"
grid.col   <- "#E8E6E1"
point.col  <- "#8B0000"



############################################################
# GRATICULE
############################################################


graticule <- st_graticule(
  lat = seq(
    30,
    75,
    by = 10
  ),
  lon = seq(
    -30,
    60,
    by = 10
  ),
  crs = 4326
)



############################################################
# EUROPE MAP
############################################################


map <- ggplot() +


theme_void() +


geom_sf(
  data = graticule,
  colour = grid.col,
  linewidth = 0.15,
  alpha = 0.8
) +


geom_sf(
  data = world,
  fill = land.col,
  colour = border.col,
  linewidth = 0.25
) +


geom_sf(
  data = sites_sf,
  shape = 21,
  size = 3.2,
  fill = point.col,
  colour = "white",
  stroke = 0.5
) +


geom_text_repel(
  data = sites,
  aes(
    Longitude,
    Latitude,
    label = ID
  ),

  seed = 42,

  size = 3.3,

  fontface = "bold",

  colour = "black",

  box.padding = 0.35,

  point.padding = 0.25,

  force = 2,

  min.segment.length = 0,

  segment.size = 0.25,

  segment.colour = "grey40",

  max.overlaps = Inf
) +


annotation_scale(
  location = "bl",
  width_hint = 0.18,
  text_cex = 0.7,
  line_width = 0.4
) +


annotation_north_arrow(
  location = "tr",
  which_north = "true",
  height = unit(
    0.7,
    "cm"
  ),
  width = unit(
    0.7,
    "cm"
  ),
  style = north_arrow_minimal()
) +


coord_sf(
  xlim = c(
    -25,
    60
  ),
  ylim = c(
    34,
    72
  ),
  expand = FALSE
) +


theme_minimal(
  base_family = "Arial"
) +


theme(

panel.background =
element_rect(
  fill = "white",
  colour = NA
),

panel.grid =
element_blank(),

axis.text =
element_text(
  colour = "black",
  size = 8
),

plot.margin =
unit(
  c(
    0.2,
    0.2,
    0.2,
    0.2
  ),
  "cm"
)

)



map



############################################################
# EXPORT
############################################################


ggsave(

filename = output_file,

plot = map,

width = 11,

height = 8,

dpi = 600,

device = cairo_pdf

)