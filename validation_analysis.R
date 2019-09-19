# TODO(Callum)
#   Streamline graph generation scripts
#     - Automate file import and naming (no hardcoding)

##### Load/Install relevant packages #####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04

GetPackages <- function(required.packages) {
  packages.not.installed <- 
    required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){
    install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("ggplot2", "reshape2", "wesanderson", "tidyverse", "scales", "doParallel", 
              "devtools", "dplyr", "gtable", "grid", "gridExtra"))
install_github("kassambara/easyGgplot2")  # Need devtools to use this function
library(easyGgplot2)

##### Load relevant data #####

# ROI coverage percentage information
coverage_percentages_filenames <- Sys.glob(paths = "/mnt/shared_data/work/metrics_extraction_for_validation/*coverage_percentages*")
coverage_percentage_sample_names <- paste0(basename(coverage_percentages_filenames))
coverage_percentage_list <- lapply(coverage_percentages_filenames, function(i){read.table(file = i, header = T)})
coverage_percentage_melted <- do.call(rbind, coverage_percentage_list)
coverage_percentage_melted$id <- factor(rep(coverage_percentage_sample_names, each = sapply(coverage_percentage_list, nrow)))
colnames(coverage_percentage_melted) <- c("Amplicon", "1x", "10x", "20x", "40x", "80x", "160x", "200x", "320x", "640x", "sampleID")

# Amplicon/coverage/sample data
amplicon_coverage_filenames <- Sys.glob(paths = "/mnt/shared_data/work/metrics_extraction_for_validation/*amplicon_coverage")
amplicon_coverage_sample_names <- paste0(basename(amplicon_coverage_filenames))
amplicon_coverage_list <- lapply(amplicon_coverage_filenames, function(i){read.table(file = i, header = T)})
amplicon_coverage_melted <- do.call(rbind, amplicon_coverage_list)
amplicon_coverage_melted$id <- factor(rep(amplicon_coverage_sample_names, each = sapply(amplicon_coverage_list, nrow)))

##### Pick colours #####

colour_palette <- wesanderson::wes_palettes$Darjeeling1

##### Generate the plotting functions #####

# Plot percentage of target regions achieving various levels of coverage
CoverageDepth <- function(dataframe, coverage, coverage.depth, sample){
  ggplot(dataframe) +
    geom_boxplot(aes(x = coverage.depth, y = coverage, fill = coverage.depth), outlier.shape = NA, notch = F) +
    xlab(sample) +
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
    geom_boxplot(aes(x = reorder(x = sample, X = coverage), y = coverage, fill = sample), outlier.shape = NA, notch = T) +
    ggtitle("Distribution of total barcode reads per amplicon across each sample\n") +
    # scale_fill_manual(values = colour_palette) +
    # coord_cartesian(ylim = quantile(coverage, c(0.1, 0.9), na.rm = T)) +
    scale_y_continuous(breaks = seq(0, 4000, 500)) +
    ylab("Number of barcode reads for a specific amplicon\n") +
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
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
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


##### Generate the plots #####

# Plot1
p <- list()
lapply(unique(coverage_percentage_melted$sampleID), function(i) {
  coverage_percentage_melted[coverage_percentage_melted$sampleID == i,] %>% 
  melt %>%
  CoverageDepth(., .$value, .$variable, i) %>%
  {. ->> p[[i]] }
})  
output <- grid.arrange(grobs = p,
             top = textGrob("Percentage of amplicons achieving various levels of coverage", vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)),
             left = textGrob("IQR of amplicon coverage", rot = 90, vjust = 1))
ggsave("~/Desktop/Rplot.pdf", output, width = 16*1.25, height = 9*1.25)

# Plot 2
pdf("~/Desktop/Rplot.pdf", width = 16*1.5, height = 9*1.5)
SampleCoverageDistrubtion(
  amplicon_coverage_melted,
  amplicon_coverage_melted$BARCODE_COUNT,
  amplicon_coverage_melted$id)
graphics.off()

# Plot 3
p <- list()
lapply(unique(amplicon_coverage_melted$id), function(i) {
  amplicon_coverage_melted[amplicon_coverage_melted$id == i,] %>%
    AmpliconCoverageDistrubtion(., .$BARCODE_COUNT, .$PRIMER, i) %>%
    {. ->> p[[i]] }
})  
output <- grid.arrange(
  grobs = p,
  top = textGrob("Vertical read depth for each amplicon", vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)),
  left = textGrob("Barcode-adjusted read depth", rot = 90, vjust = 1))
ggsave("~/Desktop/Rplot.png", device = "png", output, width = 16*1.25, height = 9*1.25)





placeholder






















##### Unfinished attempts to plot VAF graphs #####

# Load the datavariant_frequency <- Sys.glob(paths = "/mnt/shared_data/work/metrics_extraction_for_validation/*VAF_frequencies_bare*")
variant_frequency <- lapply(variant_frequency, read.table)
variant_frequency <- do.call(rbind.data.frame, variant_frequency)
variant_frequency <- variant_frequency[-1,]
hist(as.numeric(variant_frequency$V2), breaks = 100, xlim = c(0, 1))

##### Identify common zero coverage primers #####

zero_coverage_regions_s1_L1 <- read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L001_default_zero_coverage", header = T)
zero_coverage_regions_s1_L2 <- read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L002_default_zero_coverage", header = T)
zero_coverage_regions_s1_L3 <- read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L003_default_zero_coverage", header = T)
zero_coverage_regions_s1_L4 <- read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L004_default_zero_coverage", header = T)

# Export zero coverage regions found in all samples
common_zero_cov_primers <- base::Reduce(
  base::intersect, base::list(
    zero_coverage_regions_s1_L1$name, 
    zero_coverage_regions_s1_L2$name, 
    zero_coverage_regions_s1_L3$name,
    zero_coverage_regions_s1_L4$name)
)

write.table(common_zero_cov_primers, file = "~/Desktop/zero_coverage_amplicons")

##### Identify common low coverage primers #####

below_200_coverage_regions_s1_L1 <- read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L001_default_low_coverage", header = T)
below_200_coverage_regions_s1_L2 <- read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L002_default_low_coverage", header = T)
below_200_coverage_regions_s1_L3 <- read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L003_default_low_coverage", header = T)
below_200_coverage_regions_s1_L4 <- read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L004_default_low_coverage", header = T)

# Export zero coverage regions found in all samples
common_zero_cov_primers <-base::Reduce(
  base::intersect, base::list(below_200_coverage_regions_s1_L1$name, 
                              below_200_coverage_regions_s1_L2$name, 
                              below_200_coverage_regions_s1_L3$name, 
                              below_200_coverage_regions_s1_L4$name))

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

# Import the per primer pair coverage metrics
barcode_coverage_s1_per_amplicon <- 
  read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-80_S1_L001_default_amplicon_coverage", 
             header = T, col.names = c("Amplicon", "Sample_1_L1"))
barcode_coverage_s2_per_amplicon <- 
  read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-20_S4_L002_default_amplicon_coverage",
             header = T, col.names = c("Amplicon", "Sample_1_L2"))
barcode_coverage_s3_per_amplicon <-
  read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-20_S4_L003_default_amplicon_coverage",
             header = T, col.names = c("Amplicon", "Sample_1_L3"))
barcode_coverage_s4_per_amplicon <-
  read.table("/mnt/shared_data/work/metrics_extraction_for_validation/18F199-20_S4_L004_default_amplicon_coverage", skip = 7,
             header = T, col.names = c("Amplicon", "Sample_1_L4"))

barcode_coverage_df_per_amplicon <- 
  list(barcode_coverage_s1_per_amplicon, 
       barcode_coverage_s2_per_amplicon,
       barcode_coverage_s3_per_amplicon, 
       barcode_coverage_s4_per_amplicon
  ) %>% 
  reduce(full_join, by = "Amplicon")

# Plot all barcode reads for all amplicon by sample
amplicon_sample_coverage <- melt(barcode_coverage_df_per_amplicon, measure.vars = c(
  "Sample_1_L1", "Sample_1_L2", "Sample_1_L3", "Sample_1_L4"
))

colnames(amplicon_sample_coverage) <- c("amplicon", "sample", "coverage")

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

