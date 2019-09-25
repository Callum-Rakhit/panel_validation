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
    geom_boxplot(aes(x = reorder(x = sample, X = coverage), y = coverage, fill = sample), outlier.shape = NA, notch = T) +
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


##### Generate the plots #####

# Crude coverage summary
# lapply(unique(amplicon_coverage_melted$id), function(i) {
#   amplicon_coverage_melted[amplicon_coverage_melted$id == i,] %>%
#     select(BARCODE_COUNT) %>%
#     sum
# })

# Plot1
p <- list()
lapply(unique(coverage_percentage_melted$sampleID), function(i) {
  coverage_percentage_melted[coverage_percentage_melted$sampleID == i,] %>% 
  melt %>%
  CoverageDepth(., .$value, .$variable, i) %>%
  {. ->> p[[i]] }
})  
output <- grid.arrange(grobs = p,
             top = textGrob("Percent of amplicons achieving 'x' level of coverage", vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)),
             left = textGrob("IQR of amplicon coverage (%)", rot = 90, vjust = 1))
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
ggsave("~/Desktop/Rplot.pdf", output, width = 16*1.25, height = 9*1.25)

placeholder































##### Unfinished attempts to plot VAF graphs #####

VAFs <- read.csv("~/panel_validation/Horizon_Control_Locations_AF.csv")

HD200 <- as.data.frame(variant_frequency_list[["HD200_S11_default_VAF_frequencies_bare"]])
HD701 <- as.data.frame(variant_frequency_list[["HD701_S10_default_VAF_frequencies_bare"]])
HD798 <- as.data.frame(variant_frequency_list[["HD798_S9_default_VAF_frequencies_bare"]])

HD200_VAFs <- list()

for(i in VAFs[VAFs$ControlID == "HD200", ]$Exact_Location){ 
  HD200_VAFs[[i]] <- HD200[HD200$V1 == i, ] }

HD200_VAFs <- do.call(rbind, HD200_VAFs)
colnames(HD200_VAFs) <- c("Exact_Location", "HD200_AF")

VAFs_HD200 <- merge(VAFs, HD200_VAFs, by = 'Exact_Location', all = T)
VAFs_HD200$AF <- gsub(pattern = "%", replacement = "", x = VAFs_HD200$AF)
VAFs_HD200$AF <- as.numeric(VAFs_HD200$AF)/100

# Plot residuals
# VAF.lm <- lm(AF ~ HD200_AF, data = VAFs_HD200[VAFs_HD200$ControlID == "HD200", ], na.action = na.exclude)
VAF.lm <- lm(AF ~ HD200_AF, data = VAFs_HD200[VAFs_HD200$ControlID == "HD200", ], na.action = na.exclude)
VAF.pred <- predict(VAF.lm)
VAF.actual <- VAFs_HD200[VAFs_HD200$ControlID == "HD200", ]$AF
VAF.res <- resid(VAF.lm)

View(VAF.pred)
View(VAFs_HD200[VAFs_HD200$ControlID == "HD200", ])

ggplot(VAFs_HD200[VAFs_HD200$ControlID == "HD200", ],
       aes(x = AF, y = HD200_AF, 
           color = ifelse(is.na(HD200_AF), "Missing", "Present"), 
           shape = ifelse(is.na(HD200_AF), "Missing", "Present"),
           size = ifelse(is.na(HD200_AF), "Missing", "Present"))) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.2) +
  geom_point(shape = 1, size = 2) +
  geom_point(aes(y = AF)) +
  geom_segment(aes(xend = AF, yend = AF), alpha = .2, size = 1, linetype = 2) +
  scale_shape_manual(name = "", values = c(Missing = 4, Present = 19)) +
  scale_color_manual(name = "", values = c(Missing = "red", Present = "black")) +
  scale_size_manual(name = "", values = c(Missing = 3, Present = 2)) +
  xlab(label = "True AF (black dots)") +
  ylab(label = "Observed AF (white dots)") +
  theme(
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

ggsave("~/Desktop/Rplot.pdf", output, width = 16*1.25, height = 9*1.25)
  
hold

ggpubr::show_point_shapes()

ggplot(d, aes(x = hp, y = mpg)) +
  geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
  geom_segment(aes(xend = hp, yend = predicted), alpha = .2) +  # alpha to fade lines
  geom_point() +
  geom_point(aes(y = predicted), shape = 1) +
  theme_bw()  # Add theme for cleaner look


a[a$V1 == "1:100292538", ]
a[a$V1 == "4:55599333", ]
a[a$V1 == "7:55241707", ]
a[a$V1 == "4:55599333", ]

# Load the datavariant_frequency <- Sys.glob(paths = "/mnt/shared_data/work/metrics_extraction_for_validation/*VAF_frequencies_bare*")
variant_frequency_filenames <- Sys.glob(paths = "/mnt/shared_data/work/metrics_extraction_for_validation_12_samples/*VAF_frequencies_bare")
variant_frequency_list <- list()
lapply(variant_frequency_filenames, function(i){read.table(file = i, header = F) %>%
    { . ->> variant_frequency_list[[paste0(basename(i))]] }
    })
variant_frequency_melted <- do.call(rbind, variant_frequency_list)



p <- list()
lapply(unique(amplicon_coverage_melted$id), function(i) {
  amplicon_coverage_melted[amplicon_coverage_melted$id == i,] %>%
    AmpliconCoverageDistrubtion(., .$BARCODE_COUNT, .$PRIMER, i) %>%
    {. ->> p[[i]] }
})  
variant_frequency_melted$id <- factor(rep(variant_frequency_sample_names, each = sapply(variant_frequency_filenames, nrow)))
hist(variant_frequency_melted$V2)

##### Identify common zero coverage primers #####

# Read in data files from snappy
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

