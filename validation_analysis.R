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

install_github("kassambara/easyGgplot2")  # Need devtools to use this function
library(easyGgplot2)

##### Load relevant data #####

# Amplicon/coverage/sample data
# amplicon_coverage_filenames <- Sys.glob(
#   paths = "/mnt/shared_data/work/metrics_extraction_for_validation/*amplicon_coverage")
amplicon_coverage_filenames <- Sys.glob(
  paths = "/mnt/shared_data/work/metrics_extraction_for_validation_12_samples/*amplicon_coverage")
amplicon_coverage_sampleIDs <- paste0(basename(amplicon_coverage_filenames))
amplicon_coverage_list <- lapply(amplicon_coverage_filenames, function(i){read.table(file = i, header = T)})
amplicon_coverage_melted <- do.call(rbind, amplicon_coverage_list)
amplicon_coverage_melted$id <- factor(rep(amplicon_coverage_sampleIDs, each = sapply(amplicon_coverage_list, nrow)))

# ROI coverage percentage information
# coverage_percentages_filenames <- Sys.glob(
#   paths = "/mnt/shared_data/work/metrics_extraction_for_validation/*coverage_percentages*")
coverage_percentages_filenames <- Sys.glob(
  paths = "/mnt/shared_data/work/metrics_extraction_for_validation_12_samples/*coverage_percentages*")
coverage_percentage_sampleIDs <- paste0(basename(coverage_percentages_filenames))
coverage_percentage_list <- lapply(coverage_percentages_filenames, function(i){read.table(file = i, header = T)})
coverage_percentage_melted <- do.call(rbind, coverage_percentage_list)
coverage_percentage_melted$id <- factor(rep(
  coverage_percentage_sampleIDs, each = sapply(coverage_percentage_list, nrow)))
colnames(coverage_percentage_melted) <- c("Amplicon", "1x", "10x", "20x", "40x", "80x", 
                                          "160x", "200x", "320x", "640x", "sampleID")

# Crude coverage summary
total_reads <- t(as.data.frame(lapply(unique(amplicon_coverage_melted$id), function(i) {
  amplicon_coverage_melted[amplicon_coverage_melted$id == i,] %>%
    select(BARCODE_COUNT) %>% sum })))
total_reads <- as.data.frame(round((total_reads/1000000), digits = 1))
total_reads$sampleID <- (unique(coverage_percentage_melted$sampleID))
colnames(total_reads) <- c("totalreads", "sampleID")
percentages_totalreads_merged <- merge(coverage_percentage_melted, total_reads, by = "sampleID")

# Get the observed VAFs for the horizon controls
variant_frequency_filenames <- Sys.glob(
  paths = "/mnt/shared_data/work/metrics_extraction_for_validation_12_samples/*VAF_frequencies_bare")
variant_frequency_list <- list()
invisible(lapply(variant_frequency_filenames, function(i){read.table(file = i, header = F) %>%
    {. ->> variant_frequency_list[[paste0(basename(i))]]} }))
variant_frequency_melted <- do.call(rbind, variant_frequency_list)
variant_frequency_melted <- setDT(variant_frequency_melted, keep.rownames = T)[]
variant_frequency_melted$rn <- gsub("_default.*", "", variant_frequency_melted$rn)
colnames(variant_frequency_melted) <- c("sampleID", "location", "VAF")

# Get the expected VAFs for the controls samples run
known_VAFs <- read.csv("~/panel_validation/Horizon_Control_Locations_AF.csv")
known_VAFs$Exact_Location <- gsub(" ", "", x = known_VAFs$Exact_Location)

# Filter the observations based on those found in the horizon controls 
VAF <- list()
for(sample.filename in names(variant_frequency_list)){
  sample.name <- (gsub('_S.*', '', sample.filename))
  sample.name.inc.conc <- (gsub('_S.*', '', sample.filename))
  sample.name <- (gsub('-.*', '', sample.name))
  sample.name <- (gsub('18F', '18F-', sample.name))
  temp_df <- as.data.frame(variant_frequency_list[[sample.filename]])
  df_name <- list()
  for(i in known_VAFs[known_VAFs$ControlID == sample.name, ]$Exact_Location){
    df_name[[i]] <- temp_df[temp_df$V1 == i, ]}
  VAFs_HorizonID <- do.call(rbind, df_name)
  colnames(VAFs_HorizonID) <- c("Exact_Location", sample.name)
  VAFs_HorizonID <- merge(known_VAFs, VAFs_HorizonID, by = "Exact_Location", all = T)
  VAFs_HorizonID <- VAFs_HorizonID[VAFs_HorizonID$ControlID == sample.name, ]
  VAFs_HorizonID$AF <- gsub(pattern = "%", replacement = "", x = VAFs_HorizonID$AF)
  VAFs_HorizonID$AF <- as.numeric(VAFs_HorizonID$AF)/100
  print(VAFs_HorizonID)
  VAF[[sample.name.inc.conc]] <- VAFs_HorizonID
}

##### Pick colours #####

colour_palette <- rep(x = wesanderson::wes_palettes$Darjeeling1, times = 10)

##### Generate the plotting functions #####

# Plot percentage of target regions achieving various levels of coverage
CoverageDepth <- function(dataframe, coverage, coverage.depth, sample){
  ggplot(dataframe) +
    geom_boxplot(aes(x = coverage.depth, y = coverage, fill = coverage.depth), outlier.shape = NA, notch = F) +
    xlab(sample) +
    scale_fill_manual(values = colour_palette) +
    theme(# Lengends to the top
      legend.position = "none",
      # Remove the y-axis
      axis.title.y = element_blank(),
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      panel.grid.minor = element_blank(),
      # Remove panel background
      panel.background = element_blank(),
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 0, hjust = 0))
}

# Plot distribution of barcode reads for all samples
SampleCoverageDistrubtion <- function(dataframe, coverage, sample){
  ggplot(dataframe) +
    geom_boxplot(aes(x = reorder(x = sample, X = coverage), y = coverage, fill = sample), 
                 outlier.shape = NA, notch = T) +
    ggtitle("Disribution of reads per amplicon for each sample\n") +
    # scale_fill_manual(values = colour_palette) +
    coord_cartesian(ylim = quantile(coverage, c(0.01, 0.99), na.rm = T)) +
    scale_fill_manual(values = colour_palette) +
    ylab("IQR of barcode-adjusted reads for all amplicons\n") +
    xlab("Sample") +
    labs(fill = "") + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          # Lengends to the top
          legend.position = "none",
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Rotate the x-axis labels 90 degrees
          axis.text.x = element_text(angle = 90, hjust = 0)
    )
}

# Plot barcode reads for all samples by amplicon
AmpliconCoverageDistrubtion <- function(dataframe, coverage, amplicon, sample){
  ggplot(dataframe) + 
    geom_point(aes(reorder(x = amplicon, X = coverage), y = coverage, color = sample)) +
    scale_fill_manual(values = colour_palette) +
    scale_y_log10(limits = c(1, 10000)) +
    xlab(sample) +
    theme(
      # Lengends to the top
      legend.position = "none",
      # Remove the y-axis
      axis.title.y = element_blank(),
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

# Plot percentage of amplicons achieving 320x coverage 
PercentageAt320 <- function(dataframe, coverage_number, read_depth){
  ggplot(dataframe) +
    geom_boxplot(aes(x = as.factor(read_depth), y = coverage_number, fill = as.factor(read_depth)),
                 outlier.shape = NA, notch = F) +
    xlab("Total Read Depth (in Millions)") +
    ylab("Percentage of Amplicons Achieving 320x Coverage") +
    scale_fill_manual(values = colour_palette) +
    ggtitle("Total Number of Reads and Percentage of Amplicons Achieving 320x Coverage") +
    theme(# Lengends to the top
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      # Remove the y-axis
      axis.title.y = element_blank(),
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      panel.grid.minor = element_blank(),
      # Remove panel background
      panel.background = element_blank(),
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 0, hjust = 0))
}

# Plot the VAFs
VAFPlot <- function(dataframe, HorizonID, HorizonID_AF){
  ggplot(dataframe, aes(x = AF, y = HorizonID_AF, 
                        color = ifelse(is.na(HorizonID_AF), "Missing", "Present"), 
                        shape = ifelse(is.na(HorizonID_AF), "Missing", "Present"),
                        size = ifelse(is.na(HorizonID_AF), "Missing", "Present"))) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_point(shape = 1, size = 2) +
    geom_point(aes(y = AF)) +
    geom_segment(aes(xend = AF, yend = AF), size = 0.5, linetype = 2, colour = "black", show.legend = F) +
    scale_shape_manual(name = NULL, values = c(Missing = 4, Present = 19)) +
    scale_color_manual(name = NULL, values = c(Missing = "red", Present = "black")) +
    scale_size_manual(name = NULL, values = c(Missing = 4, Present = 2)) +
    xlab(label = "True allelic frequency (black dots)") +
    ylab(label = "Observed allelic frequency (white dots)") +
    ggtitle(label = paste("Observed variant allelic frequency versus expected - ", HorizonID, sep = "")) +
    xlim(0, .3) +
    ylim(0, .3) +
    theme(
      # Centre title
      plot.title = element_text(hjust = 0.5),
      # Lengends to the top
      legend.position = "top",
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      panel.grid.minor = element_blank(),
      # Remove panel background
      panel.background = element_blank())
}

##### Generate the plots #####

# Percentage coverage per sample plot
p <- list()

lapply(unique(coverage_percentage_melted$sampleID), function(i) {
  coverage_percentage_melted[coverage_percentage_melted$sampleID == i,] %>% 
    melt %>% CoverageDepth(., .$value, .$variable, i) %>% {. ->> p[[i]]} })  

output <- grid.arrange(grobs = p, top = textGrob(
  "Percent of amplicons achieving 'x' level of coverage", vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)),
  left = textGrob("IQR of amplicon coverage (%)", rot = 90, vjust = 1))

ggsave("~/Desktop/Rplot.pdf", output, width = 16*1.25, height = 9*1.25)

# Coverage per amplicon plot
pdf("~/Desktop/Rplot.pdf", width = 16*1.5, height = 9*1.5)

SampleCoverageDistrubtion(
  amplicon_coverage_melted,
  amplicon_coverage_melted$BARCODE_COUNT,
  amplicon_coverage_melted$id)

graphics.off()

# Coverage per amplicon log adjusted plot
p <- list()

lapply(unique(amplicon_coverage_melted$id), function(i) {
  amplicon_coverage_melted[amplicon_coverage_melted$id == i,] %>% 
    AmpliconCoverageDistrubtion(., .$BARCODE_COUNT, .$PRIMER, i) %>% {. ->> p[[i]] } })  

output <- grid.arrange(grobs = p, top = textGrob(
  "Vertical read depth for each amplicon", vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)),
  left = textGrob("Barcode-adjusted read depth", rot = 90, vjust = 1))

ggsave("~/Desktop/Rplot.pdf", output, width = 16*1.25, height = 9*1.25)

# Plot 4a - Horizon control HD200
output <- VAFPlot(VAF$HD200, "HD200", VAF$HD200$HD200) 
ggsave("~/Desktop/Rplot.pdf", output, width = 16*1, height = 9*1)

# Plot 4b - Horizon control HD701
output <- VAFPlot(VAF$HD701, "HD701", VAF$HD701$HD701)
ggsave("~/Desktop/Rplot.pdf", output, width = 16*1, height = 9*1)

# Plot 4c - Horizon control HD798
output <- VAFPlot(VAF$HD798, "HD798", VAF$HD798$HD798)
ggsave("~/Desktop/Rplot.pdf", output, width = 16*1, height = 9*1)

# Plot the 320x percentage plot
output <- PercentageAt320(
  percentages_totalreads_merged, percentages_totalreads_merged$`320x`, percentages_totalreads_merged$totalreads)

ggsave("~/Desktop/Rplot.pdf", output, width = 16*1, height = 9*1)

##### Identify common zero coverage primers #####

# Read in zero coverage regions from snappy
zero_coverage_regions_s1_L1 <- read.table(
  "/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L001_default_zero_coverage", header = T)
zero_coverage_regions_s1_L2 <- read.table(
  "/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L002_default_zero_coverage", header = T)
zero_coverage_regions_s1_L3 <- read.table(
  "/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L003_default_zero_coverage", header = T)
zero_coverage_regions_s1_L4 <- read.table(
  "/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L004_default_zero_coverage", header = T)

# Export zero coverage regions found in all samples
common_zero_cov_primers <- base::Reduce(
  base::intersect, base::list(
    zero_coverage_regions_s1_L1$name, 
    zero_coverage_regions_s1_L2$name, 
    zero_coverage_regions_s1_L3$name,
    zero_coverage_regions_s1_L4$name))

write.table(common_zero_cov_primers, file = "~/Desktop/zero_coverage_amplicons")

##### Identify common low coverage primers #####

# Read in low coverage (<200) regions from snappy
below_200_coverage_regions_s1_L1 <- read.table(
  "/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L001_default_low_coverage", header = T)
below_200_coverage_regions_s1_L2 <- read.table(
  "/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L002_default_low_coverage", header = T)
below_200_coverage_regions_s1_L3 <- read.table(
  "/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L003_default_low_coverage", header = T)
below_200_coverage_regions_s1_L4 <- read.table(
  "/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L004_default_low_coverage", header = T)

# Export zero coverage regions found in all samples
common_zero_cov_primers <-base::Reduce(
  base::intersect, base::list(
    below_200_coverage_regions_s1_L1$name, below_200_coverage_regions_s1_L2$name, 
    below_200_coverage_regions_s1_L3$name, below_200_coverage_regions_s1_L4$name))

write.table(common_zero_cov_primers, file = "~/Desktop/zero_coverage_amplicons")

# Extract percentage of target regions achieving various levels of coverage
percentage_coverage_regions_s1 <- colMeans(percentage_coverage_regions_s1[,2:10])
percentage_coverage_regions_s2 <- colMeans(percentage_coverage_regions_s2[,2:10])
percentage_coverage_regions_s3 <- colMeans(percentage_coverage_regions_s3[,2:10])
percentage_coverage_regions_s4 <- colMeans(percentage_coverage_regions_s4[,2:10])
coverage <- as.numeric(c("1", "10", "20", "40", "80", "160", "200", "320", "640"))

df <- data_frame(coverage, percentage_coverage_regions_s1, percentage_coverage_regions_s2, 
                 percentage_coverage_regions_s3, percentage_coverage_regions_s4)
sample_coverage <- melt(df, id.vars = "coverage")
colnames(sample_coverage) <- c("read_depth", "sample", "coverage_percentage")

# Pull out all the VAFs and plot them
s1_VAF <- read.table("~/S1_VAF", header = F, col.names = c("Location", "VAF"))
s1_VAF$Sample <- 1
s2_VAF <- read.table("~/S2_VAF", header = F, col.names = c("Location", "VAF"))
s2_VAF$Sample <- 2
s3_VAF <- read.table("~/S3_VAF", header = F, col.names = c("Location", "VAF"))
s3_VAF$Sample <- 3
s4_VAF <- read.table("~/S4_VAF", header = F, col.names = c("Location", "VAF"))
s4_VAF$Sample <- 4

VAF <- rbind(s1_VAF, s2_VAF, s3_VAF, s4_VAF)

ggplot2.histogram(data = VAF, xName = 'VAF', groupName = 'Sample', legendPosition = "top",
                  faceting = T, facetingVarNames = 'Sample', alpha = 0.5, position = "identity")
