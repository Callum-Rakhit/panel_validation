# TODO(Callum)
#  - Take input directly from snappy
#  - Automate file import and naming (no hardcoding)

##### Load/Install relevant packages #####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("ggplot2", "reshape2", "wesanderson", "tidyverse", "scales", "doParallel", 
              "devtools", "dplyr", "gtable", "grid", "gridExtra", "data.table"))

# Load relevant data amplicon/coverage/sample data
amplicon_coverage_filenames <- Sys.glob(paths = "/mnt/shared_data/work/three_runs_together/*amplicon_coverage")
amplicon_coverage_sampleIDs <- paste0(basename(amplicon_coverage_filenames))
amplicon_coverage_list <- lapply(amplicon_coverage_filenames, function(i){read.table(file = i, header = T)})
amplicon_coverage_melted <- do.call(rbind, amplicon_coverage_list)
amplicon_coverage_melted$id <- factor(rep(amplicon_coverage_sampleIDs, each = sapply(amplicon_coverage_list, nrow)))

amplicon_coverage_average_per_run <- setDT(amplicon_coverage_melted)[ , .(mean_coverage = mean(BARCODE_COUNT)), by = id]
amplicon_coverage_melted$mean_coverage <- amplicon_coverage_average_per_run$mean_coverage[match(
  amplicon_coverage_melted$id, amplicon_coverage_average_per_run$id)]


# Get list for specific genes/categories
CEBPA_amplicon_coverage_melted <- amplicon_coverage_melted[grep("^CEBPA", amplicon_coverage_melted$PRIMER),]

Adult_Solid_Gene_List <- c(
  "^BRAF",
  "^BRCA1",
  "^BRCA2",
  "^CDKN2A",
  "^EGFR",
  "^HRAS",
  "^KIT",
  "^KRAS",
  "^MET",
  "^MLH1",
  "^NRAS",
  "^PDGFRA",
  "^SMARCA4",
  "^TP53"
)

Adult_Solid_amplicon_coverage_melted <- amplicon_coverage_melted[(
  grep(paste(Adult_Solid_Gene_List, collapse="|"), amplicon_coverage_melted$PRIMER)
  )]

# Pick colours
colour_palette <- rep(x = wesanderson::wes_palettes$Darjeeling1, times = 10)

# Create the plotting function
RunCov.vs.AvCov <- function(dataframe, id, BARCODE_COUNT, mean_coverage){
  ggplot(dataframe) + 
  geom_boxplot(aes(reorder(x = id, X = BARCODE_COUNT),
                 y = BARCODE_COUNT, color = "#000000"), coef = 6) +
  geom_line(aes(reorder(x = id, X = BARCODE_COUNT),
                y = mean_coverage, group = 1)) +
  # scale_fill_manual(values = colour_palette) +
  # scale_y_log10(limits = c(1, 10000)) +
  xlab("Various different runs") +
  ylab("Coverage (average coverage for all amplicons is the line)") + 
  ggtitle("Coverage per run vs average coverage") +
  theme(
    # Lengends to the top
    legend.position = "none",
    # Remove the y-axis
    # axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank())
}

Deviation.from.Av <- function(dataframe, PRIMER, BARCODE_COUNT, mean_coverage){
  ggplot(dataframe) +
  geom_boxplot(aes(x = PRIMER, y = (BARCODE_COUNT - mean_coverage), color = PRIMER), coef = 6) +
  scale_fill_manual(values = colour_palette) +
  xlab("Sample ID") +
  theme(
    # Lengends to the top
    legend.position = "none",
    # Remove the y-axis
    axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank())
}

# Create the plot
output <- RunCov.vs.AvCov(CEBPA_amplicon_coverage_melted, id, BARCODE_COUNT, mean_coverage)
ggsave("~/Desktop/CEPRA_1.pdf", output, width = 16*1.25, height = 9*1.25)
output <- Deviation.from.Av(CEBPA_amplicon_coverage_melted, PRIMER, BARCODE_COUNT, mean_coverage)
ggsave("~/Desktop/CEPRA_2.pdf", output, width = 16*1.25, height = 9*1.25)

output <- RunCov.vs.AvCov(Adult_Solid_amplicon_coverage_melted, id, BARCODE_COUNT, mean_coverage)
ggsave("~/Desktop/AS_1.pdf", output, width = 16*1.25, height = 9*1.25)
output <- Deviation.from.Av(Adult_Solid_amplicon_coverage_melted, PRIMER, BARCODE_COUNT, mean_coverage)
ggsave("~/Desktop/AS_2.pdf", output, width = 16*1.25, height = 9*1.25)

Adult_Solid_amplicon_coverage_melted$diff <- (Adult_Solid_amplicon_coverage_melted$BARCODE_COUNT/Adult_Solid_amplicon_coverage_melted$mean_coverage)
CEBPA_amplicon_coverage_melted$diff <- (CEBPA_amplicon_coverage_melted$BARCODE_COUNT/CEBPA_amplicon_coverage_melted$mean_coverage)

boxplot(Adult_Solid_amplicon_coverage_melted$diff)
mean(CEBPA_amplicon_coverage_melted$diff)
mean(Adult_Solid_amplicon_coverage_melted$diff)
mean(CEBPA_amplicon_coverage_melted$diff)

# Read in low coverage (<200) regions from snappy
low_coverage_filenames <- Sys.glob(paths = "/mnt/shared_data/work/three_runs_together/*low_coverage")

# Extract all low coverage samples
low_coverage_regions <- lapply(
  low_coverage_filenames, function(i) {
    read.table(i, header = T)
  })

# Identify low coverage regions found in all samples
common_low_cov_primers <- Reduce(
  f = intersect, 
  x = lapply(low_coverage_regions, '[[', 4)  # 4th element is a list of names of the low coverage amplicon regions
)

common_low_cov_primers<- as.data.frame(common_low_cov_primers)
common_low_cov_primers <- common_low_cov_primers[!grepl("^chr", common_low_cov_primers$common_low_cov_primers),]
View(common_low_cov_primers)

