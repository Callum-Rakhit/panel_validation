# TODO(Callum)
#  - Remove labels
#  - Add bounds with "grey zone"
#  - Convert divergance in "ratio"

##### Load/Install relevant packages #####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("ggplot2", "nlme", "devtools", "tidyverse", "reshape2", "wesanderson", "gmodels"))

# Load in the amplicon/coverage/sample data
amplicon_coverage_filenames <- Sys.glob(paths = "/mnt/shared_data/work/three_runs_together/*default_amplicon_coverage")
amplicon_coverage_filenames <- amplicon_coverage_filenames[1:length(amplicon_coverage_filenames)-1]
amplicon_coverage_sampleIDs <- paste0(basename(amplicon_coverage_filenames))
amplicon_coverage_list <- lapply(amplicon_coverage_filenames, function(i){read.table(file = i, header = T)})
amplicon_coverage_melted <- do.call(rbind, amplicon_coverage_list)
amplicon_coverage_melted$id <- factor(rep(amplicon_coverage_sampleIDs, each = sapply(amplicon_coverage_list, nrow)))
colnames(amplicon_coverage_melted) <- c("PRIMER", "COVERAGE", "ID")
rm(list = c("amplicon_coverage_filenames", "amplicon_coverage_sampleIDs", "amplicon_coverage_list"))

# Load in the chromosomal location data
bed_file <- read.table(file = "/mnt/shared_data/snappy/snappy/tas/ccp1/ccp1_primer.bed", header = F, )
colnames(bed_file) <- c("CHR", "START", "STOP", "PRIMER", "SCORE", "STRAND")

#######################################################################################################################

# Formatting data for mahattan plot
manhattan_data <- merge(x = amplicon_coverage_melted, y = bed_file, by = "PRIMER")
manhattan_data$CHR <- as.numeric(manhattan_data$CHR)
manhattan_data <- manhattan_data[order(manhattan_data$CHR),]
manhattan_data$CHR <- ifelse(manhattan_data$CHR == 23, "X", manhattan_data$CHR)

nCHR <- length(unique(manhattan_data$CHR))
manhattan_data$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(manhattan_data$CHR)){
  nbp[i] <- max(manhattan_data[manhattan_data$CHR == i,]$START)
  manhattan_data[manhattan_data$CHR == i, "BPcum"] <- manhattan_data[manhattan_data$CHR == i, "START"] + s
  s <- s + nbp[i]
}

axis.set <- manhattan_data %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- max(manhattan_data$COVERAGE) + 2 

ggplot(subset_manhattan, 
                   aes(x = subset_manhattan$BPcum, y = subset_manhattan$COVERAGE, 
                       color = as.factor(subset_manhattan$CHR), size = subset_manhattan$COVERAGE)) +
  geom_smooth() +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2500)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = "Chromsome", y = "Coverage") + 
  theme_minimal() +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

##### Proportionality plot #####

amplicon_proportions_all_samples <- as.data.frame(amplicon_coverage_melted$PRIMER)
for(i in sort(unique(amplicon_coverage_melted$ID))){
  df <- amplicon_coverage_melted[(amplicon_coverage_melted$ID == i), ]
  options(scipen = 999)
  df$i <- df$COVERAGE/sum(df$COVERAGE)
  amplicon_proportions_all_samples <- cbind(amplicon_proportions_all_samples, df$i)
}

example_plot <- amplicon_coverage_melted[(amplicon_coverage_melted$ID == "18F-199_S1_default_amplicon_coverage"), ]
options(scipen = 999)
example_plot$divergance199 <- example_plot$COVERAGE/sum(example_plot$COVERAGE)

amplicon_median_proportion <- melt(amplicon_proportions_all_samples)
amplicon_median_proportion <- amplicon_median_proportion %>% 
  group_by(amplicon_median_proportion$`amplicon_coverage_melted$PRIMER`) %>% 
  summarise(average = mean(value))
colnames(amplicon_median_proportion) <- c("PRIMER", "MEDIAN_DIVERGANCE")

example_plot <- merge(x = example_plot, y = amplicon_median_proportion, by = "PRIMER")
example_plot$ratio <- example_plot$divergance199/example_plot$MEDIAN_DIVERGANCE
colnames(example_plot) <- c("PRIMER", "COVERAGE", "ID", "DIVERGANCE", "MEDIAN_DIVERGANCE", "RATIO")
example_plot <- merge(x = example_plot, y = bed_file, by = "PRIMER")

# sample$BP <- factor(manhattan_data$START, levels = manhattan_data[ !duplicated(
#   manhattan_data[,"START"]), "START"][order(manhattan_data[!duplicated(manhattan_data[ , "START"]),  "CHR"] )])
# 
# manhattan_data[manhattan_data$ID == "18F-199_S1_default_amplicon_coverage",]

# Single_Sample_Proportional_Read_Depth
# 13.333 inch by 7.5

# Make proportionality Manhatten

example_plot$CHR <- as.numeric(example_plot$CHR)
example_plot <- example_plot[order(example_plot$CHR),]
example_plot$CHR <- ifelse(example_plot$CHR == 23, "X", example_plot$CHR)
nCHR <- length(unique(example_plot$CHR))
example_plot$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(example_plot$CHR)){
  nbp[i] <- max(example_plot[example_plot$CHR == i,]$START)
  example_plot[example_plot$CHR == i, "BPcum"] <- example_plot[example_plot$CHR == i, "START"] + s
  s <- s + nbp[i]
}
axis.set <- example_plot %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

# Remove 0 and infinate values
example_plot <- example_plot[(example_plot$RATIO != 0), ]
example_plot <- example_plot[(example_plot$RATIO != Inf), ]
example_plot <- example_plot[!(is.na(example_plot$RATIO)), ]

# Calculate a 95% CI for the ratio value
hist(example_plot$RATIO, breaks = 1000)
ci(example_plot$RATIO)
mean(example_plot$RATIO)
median(example_plot$RATIO)

# Create the final plot
ggplot(example_plot, aes(x = example_plot$BPcum, y = example_plot$RATIO, 
                         color = as.factor(example_plot$CHR))) +
  geom_point(size = 2) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_color_manual(values = rep(wes_palettes$Darjeeling1, 5)) +
  labs(x = "Chromsome", y = "Proportional increase/decrease in barcoded reads") + 
  # geom_text(aes(label = ifelse(RATIO > 20 | RATIO <0.01, as.character(PRIMER), '')), hjust = -0.05, vjust = 0) +
  theme_minimal() +
  ggtitle("The proportion of molecularly barcoded reads for each primer in a single sample vs the median expected values") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
  )

##### LOH plot #####
START <- data.frame(do.call('rbind', strsplit(as.character(variant_frequency_melted$location), ':', fixed = T)))
# LOH_data <- within(variant_frequency_melted, START <- data.frame(do.call('rbind', strsplit(as.character(location), ':', fixed = T))))
LOH_data <- cbind(START$X1, START$X2, variant_frequency_melted$VAF, variant_frequency_melted$coverage)
LOH_data <- as.data.frame(LOH_data)
colnames(LOH_data) <- c("CHR", "START", "AF", "COVERAGE")

LOH_data$CHR <- as.numeric(LOH_data$CHR)
LOH_data <- LOH_data[order(LOH_data$CHR),]
LOH_data$CHR <- ifelse(LOH_data$CHR == 23, "X", LOH_data$CHR)
nCHR <- length(unique(LOH_data$CHR))
LOH_data$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(LOH_data$CHR)){
  nbp[i] <- max(LOH_data[LOH_data$CHR == i,]$START)
  LOH_data[LOH_data$CHR == i, "BPcum"] <- LOH_data[LOH_data$CHR == i, "START"] + s
  s <- s + nbp[i]
}
axis.set <- LOH_data %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

# Remove 0 and infinate values
LOH_data <- LOH_data[!(is.na(LOH_data$AF)), ]
LOH_data <- LOH_data[(LOH_data$AF <= 1.1), ]

hist(LOH_data$AF, breaks = 100)

# Create the final plot
install.packages("plotly")
library(plotly)

ggplot(LOH_data, aes(x = LOH_data$BPcum, y = LOH_data$AF, 
                         color = as.factor(LOH_data$CHR))) +
  geom_point(size = 1) +
  # geom_density_2d() + 
  # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_color_manual(values = rep(wes_palettes$Darjeeling1, 5)) +
  labs(x = "Chromsome", y = "Allelic frequency") + 
  # geom_text(aes(label = ifelse(RATIO > 20 | RATIO <0.01, as.character(PRIMER), '')), hjust = -0.05, vjust = 0) +
  theme_minimal() +
  ggtitle("Allelic frequencies across the genome for a sample run through the CCP") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
  )

placeholder




