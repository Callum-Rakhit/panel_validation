# TODO(Callum)
#   Automate the graph generation scripts
#     - Automate file import and naming

##### Load/Install relevant packages #####

GetPackages <- function(required.packages) {
  packages.not.installed <- 
    required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){
    install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("ggplot2", "reshape2", "wesanderson", "tidyverse", "scales", "doParallel", "devtools"))

install_github("kassambara/easyGgplot2")
library(easyGgplot2)

##### Load relevant data #####

coverage_metrics <- read.table(
  "analysis/20190727_183445/18F-199_S1/default/metrics/*.metrics.targetcoverage")

# Identify common zero coverage primers
zero_coverage_regions_s1 <- read.table("~/outfile_S104", header = T)
zero_coverage_regions_s2 <- read.table("~/outfile_S204", header = T)
zero_coverage_regions_s3 <- read.table("~/outfile_S304", header = T)
zero_coverage_regions_s4 <- read.table("~/outfile_S404", header = T)

# Export zero coverage regions
zero_summary <- rbind(zero_coverage_regions_s1,
                      zero_coverage_regions_s2,
                      zero_coverage_regions_s3,
                      zero_coverage_regions_s4)
write.table(zero_summary, file = "~/Desktop/zero_coverage_amplicons")

common_zero_cov_primers <- base::Reduce(
  base::intersect, base::list(
    zero_coverage_regions_s1$name, 
    zero_coverage_regions_s2$name, 
    zero_coverage_regions_s3$name,
    zero_coverage_regions_s4$name)
)

# Identify common low coverage primers
below_200_coverage_regions_s1 <- read.table("~/outfile_S105", header = T)
below_200_coverage_regions_s2 <- read.table("~/outfile_S205", header = T)
below_200_coverage_regions_s3 <- read.table("~/outfile_S305", header = T)
below_200_coverage_regions_s4 <- read.table("~/outfile_S405", header = T)

common_zero_cov_primers <-base::Reduce(
  base::intersect, base::list(below_200_coverage_regions_s1$name, 
                              below_200_coverage_regions_s2$name, 
                              below_200_coverage_regions_s3$name, 
                              below_200_coverage_regions_s4$name))


# Identify common low coverage primers
percentage_coverage_regions_s1 <- read.table("~/outfile_S106", header = T)
percentage_coverage_regions_s2 <- read.table("~/outfile_S206", header = T)
percentage_coverage_regions_s3 <- read.table("~/outfile_S306", header = T)
percentage_coverage_regions_s4 <- read.table("~/outfile_S406", header = T)

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

##### Produce functions which generate the plots #####

# Plot percentage of target regions achieving various levels of coverage
CoverageDepth <- function(dataframe, read_depth, coverage_percentage, sample){
  ggplot(dataframe, aes(read_depth, coverage_percentage, col = sample)) +
    xlab("Average coverage") +
    ylab("Fraction of the region of interest included") +
    geom_path() +
    ylim(0.75, 1) +
    theme_minimal()
}

# Import the per primer pair coverage metrics
barcode_coverage_s1_per_amplicon <- 
  read.table("~/sample_1_primer_barcode_coverage_inc_amplicon", skip = 7, header = F, col.names = c("Amplicon", "Sample_1"))
barcode_coverage_s2_per_amplicon <- 
  read.table("~/sample_2_primer_barcode_coverage_inc_amplicon", skip = 7, header = F, col.names = c("Amplicon", "Sample_2"))
barcode_coverage_s3_per_amplicon <- 
  read.table("~/sample_3_primer_barcode_coverage_inc_amplicon", skip = 7, header = F, col.names = c("Amplicon", "Sample_3"))
barcode_coverage_s4_per_amplicon <- 
  read.table("~/sample_4_primer_barcode_coverage_inc_amplicon", skip = 7, header = F, col.names = c("Amplicon", "Sample_4"))

barcode_coverage_df_per_amplicon <- 
  list(barcode_coverage_s1_per_amplicon, 
       barcode_coverage_s2_per_amplicon,
       barcode_coverage_s3_per_amplicon,
       barcode_coverage_s4_per_amplicon) %>% 
  reduce(full_join, by = "Amplicon")

# Plot all barcode reads for all amplicon by sample
amplicon_sample_coverage <- melt(barcode_coverage_df_per_amplicon, measure.vars = c(
  "Sample_1", "Sample_2", "Sample_3", "Sample_4"))
colnames(amplicon_sample_coverage) <- c("amplicon", "sample", "coverage")

SampleCoverageDistrubtion <- function(dataframe, amplicon, coverage, sample){
  ggplot(dataframe) +
    geom_boxplot(aes(x = sample, y = coverage, fill = sample), outlier.shape = NA, notch = T) +
    ggtitle("Distribution of total barcode reads per amplicon across each sample\n") +
    scale_fill_manual(values = colour_palette) +
    coord_cartesian(ylim = quantile(dataframe$coverage, c(0.1, 0.9))) +
    scale_y_continuous(breaks = seq(0, 4000, 500)) +
    ylab("Number of barcode reads for a specific amplicon\n") +
    xlab("\nSample") +
    labs(fill = "") + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          # Lengends to the top
          legend.position = "right",
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          axis.text.x = element_blank()
    )
}

# Plot barcode reads for all samples by amplicon
AmpliconCoverageDistrubtion <- function(amplicon_sample_coverage, amplicon, value){
  ggplot(amplicon_sample_coverage, aes(x = reorder(Amplicon, value, FUN = mean), y = value)) +
    geom_boxplot(outlier.size = 0.1) +
    geom_jitter(position = position_jitter(0.2)) +
    ggtitle("Distribution of total barcode reads per amplicon across each sample\n") +
    ylab("Average barcode reads across all samples\n") +
    xlab("\nAmplicon") + 
    labs(fill = "") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          # Lengends to the top
          legend.position = "right",
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          axis.text.x = element_blank()
    )
}

##### Pick colours #####

colour_palette <- wesanderson::wes_palettes$Darjeeling1

##### Generate the plots #####

CoverageDepth(sample_coverage, sample_coverage$read_depth, 
              sample_coverage$coverage_percentage, sample_coverage$sample)

SampleCoverageDistrubtion(
  amplicon_sample_coverage, amplicon_sample_coverage$amplicon, 
  amplicon_sample_coverage$sample, amplicon_sample_coverage$coverage)

AmpliconCoverage(amplicon_sample_coverage, amplicon_sample_coverage$amplicon, amplicon_sample_coverage$coverage)

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

