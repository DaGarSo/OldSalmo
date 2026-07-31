############################################################
# Figure 3
#
# Mitochondrial nucleotide diversity (π) and effective
# population size (Ne) estimates for Atlantic salmon
#
# OldSalmo project
#
# Calculates:
# - Corrected nucleotide diversity (π)
# - Bootstrap confidence intervals
# - Effective population size (Ne)
#   assuming generation times of 3 and 4 years
#
# Input:
# data/Nef/Filtered_mtDNA_alignment.fasta
#
# Output:
# figures/Figure3_pi_Nef.pdf
############################################################


############################################################
# Libraries
############################################################

library(ape)
library(pegas)
library(dplyr)
library(tidyr)
library(ggplot2)


############################################################
# Paths
############################################################

fasta_file <- "data/Nef/Filtered_mtDNA_alignment.fasta"

output_file <- "figures/Figure3_pi_Nef.pdf"



############################################################
# Load mitochondrial sequences
############################################################

seqs <- read.FASTA(
  fasta_file
)

names(seqs) <- gsub(
  " ",
  "_",
  names(seqs)
)


# Remove problematic ancient sample
seqs <- seqs[
  !names(seqs) %in%
    "S22643"
]


############################################################
# Ancient samples
############################################################

muestras_OLD <- c(
  "S16680",
  "S16681",
  "S16682",
  "S22641",
  "S22642",
  "S22644",
  "S22645"
)



############################################################
# Population assignment
############################################################

pop_name_map <- c(
  "AS" = "AS",
  "AG" = "AG",
  "EL" = "EL",
  "TE" = "TE",
  "TA" = "TA",
  "BV" = "BV",
  "CW" = "CW",
  "DB" = "DB",
  "ST" = "ST",
  "F"  = "F",
  "OR" = "OR",
  "LX" = "LX",
  "NU" = "NU",
  "OY" = "OY",
  "UF" = "UF",
  "NE" = "NE",
  "TW" = "TW",
  "HF" = "HF",
  "SG" = "SG",
  "BJ" = "BJ",
  "ER" = "ER",
  "KO" = "KO",
  "NA" = "NA",
  "RY" = "RY",
  "PG" = "PG",
  "PE" = "PE",
  "TN" = "TN"
)


nombres <- names(seqs)


poblaciones <- sapply(
  nombres,
  function(x){

    if(x %in% muestras_OLD){
      return("OLD")
    }

    prefix <- sub(
      "^([A-Z]+).*",
      "\\1",
      x
    )

    if(prefix %in% names(pop_name_map)){
      return(pop_name_map[[prefix]])
    }

    return(NA)

  }
)


valid_idx <- !is.na(poblaciones)


sub_seqs <- seqs[valid_idx]

poblaciones <- poblaciones[valid_idx]


attr(
  sub_seqs,
  "pop"
) <- factor(poblaciones)


print(
  table(attr(sub_seqs,"pop"))
)



############################################################
# Parameters
############################################################

mu <- 1e-7

gen3 <- 3

gen4 <- 4

nboot <- 10000

set.seed(123)



############################################################
# Diversity and Ne estimation
############################################################


poblaciones_unicas <- unique(
  attr(sub_seqs,"pop")
)


resultados <- data.frame()


boot_long <- data.frame()



for(pop in poblaciones_unicas){

  cat(
    "Processing:",
    pop,
    "\n"
  )


  seq_pop <- sub_seqs[
    attr(sub_seqs,"pop")==pop
  ]


  N <- length(seq_pop)


  dna_pop <- as.DNAbin(seq_pop)


  pi_raw <- nuc.div(
    dna_pop
  )


  pi_corr <- pi_raw *
    (N/(N-1))


  boot_vals_corr <- numeric(
    nboot
  )


  for(i in 1:nboot){

    idx <- sample(
      seq_len(N),
      replace=TRUE
    )


    dna_boot <- dna_pop[idx]


    pi_b <- nuc.div(
      dna_boot
    )


    boot_vals_corr[i] <-
      pi_b *
      (length(idx)/(length(idx)-1))

  }


  nef_gen3 <- pi_corr /
    (2*mu*gen3)


  nef_gen4 <- pi_corr /
    (2*mu*gen4)


  boot_nef_gen3 <-
    boot_vals_corr /
    (2*mu*gen3)


  boot_nef_gen4 <-
    boot_vals_corr /
    (2*mu*gen4)

###############################################################
# GEOGRAPHICAL ORDER AND METADATA
###############################################################

library(tibble)
library(tidyr)
library(scales)
library(here)


# Geographic order south -> north

orden_sur_norte <- c(
  "AS",
  "OLD",
  "AG",
  "EL",
  "TE",
  "TA",
  "BV",
  "CW",
  "DB",
  "ST",
  "F",
  "OR",
  "LX",
  "NU",
  "OY",
  "UF",
  "NE",
  "TW",
  "HF",
  "SG",
  "BJ",
  "ER",
  "KO",
  "NA",
  "RY",
  "PG",
  "PE",
  "TN"
)



###############################################################
# RIVER AND COUNTRY METADATA
###############################################################


rio_pais <- tibble(

Poblacion = c(
"AS","OLD","AG","EL","TE","TA",
"BV","CW","DB","ST",
"F","OR","LX","NU","OY",
"UF","NE","TW",
"HF","SG",
"BJ","ER","KO","NA",
"RY","PG","PE","TN"
),


Rio = c(

"Ason - Modern",
"Ason - Ancient",
"Allier",
"Elorn",
"Teign",
"Taw",
"Blackwater",
"Conwy",
"Dacre Beck",
"Stinchar",
"Feochan",
"Orchy",
"Laxford",
"North Uist",
"Oykel",
"Ugie",
"North Esk",
"Tweed",
"Hofsa",
"Sog",
"Bjerkreimselva",
"Eirvassdraget",
"Komagelva",
"Namsen",
"Rynda",
"Pongoma",
"Pechora",
"Teno"

),


Pais = c(

"Spain",
"Spain",
"France",
"France",
"UK",
"UK",
"UK",
"UK",
"UK",
"UK",
"Scotland",
"Scotland",
"Scotland",
"Scotland",
"Scotland",
"Scotland",
"Scotland",
"UK",
"Iceland",
"Norway",
"Norway",
"Norway",
"Norway",
"Norway",
"Russia",
"Russia",
"Russia",
"Norway"

)

)



# Combine UK + Scotland

rio_pais <- rio_pais %>%
mutate(
Pais =
ifelse(
Pais %in% c("UK","Scotland"),
"Great Britain",
Pais
)
)



###############################################################
# PREPARE DATA FOR PLOTTING
###############################################################


resultados_long <- resultados %>%


left_join(
rio_pais,
by="Poblacion"
) %>%


mutate(

Poblacion =
factor(
Poblacion,
levels=orden_sur_norte
),

Rio =
factor(
Rio,
levels =
rio_pais$Rio[
match(
orden_sur_norte,
rio_pais$Poblacion
)
]

)

) %>%


rename(

Nef_gen3_value = Nef_gen3,
Nef_gen4_value = Nef_gen4

) %>%


pivot_longer(

cols =
starts_with("Nef_gen"),

names_to =
c(
"Generacion",
".value"
),

names_pattern =
"Nef_gen(3|4)_(value|lower|upper)"

) %>%


rename(

Nef=value

) %>%


mutate(

Generacion =
case_when(

Generacion=="3" ~ "3 years",

Generacion=="4" ~ "4 years"

)

)



###############################################################
# EXPORT RESULTS
###############################################################


dir.create(
here("data","mtDNA"),
recursive=TRUE,
showWarnings=FALSE
)


write.csv(

resultados,

here(
"data",
"mtDNA",
"results_pi_Nef.csv"
),

row.names=FALSE

)


write.csv(

boot_long,

here(
"data",
"mtDNA",
"bootstrap_pi_Nef_long.csv"
),

row.names=FALSE

)



###############################################################
# FIGURE 3 STYLE
###############################################################


palette_pub <- c(

"3 years"="darkorange",

"4 years"="steelblue"

)



###############################################################
# FIGURE 3 STYLE
###############################################################


library(ggplot2)

# Tema común unificado
theme_pub <- theme_bw(base_size = 20) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(face = "bold", size = 30),
    plot.title = element_text(face = "bold", size = 40),
    plot.subtitle = element_text(size = 30),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 28),
    legend.text = element_text(size = 20)
  )

palette_pub <- c("3 years" = "darkorange", "4 years" = "steelblue")

p <- ggplot(resultados_long, aes(x = Rio, y = Nef, color = Generacion)) +
  geom_point(size = 4, alpha = 0.85,
             position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.18, size = 0.9, alpha = 0.6,
                position = position_dodge(width = 0.6)) +
  scale_color_manual(values = palette_pub, name = "Generation time") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "River",
    y = expression(N[ef]),
    color = "Generation time"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major.y = element_line(size = 0.3, color = "grey80"),
    panel.grid.minor = element_blank()
  )


###############################################################
# EXPORT FIGURE
###############################################################


dir.create(

here(
"figures"
),

showWarnings=FALSE

)



ggsave(

filename = here(
"figures",
"Figure3_pi_Nef"
),

plot=p,

device=cairo_pdf,

width=25,

height=15,

units="in",

dpi=600

)
