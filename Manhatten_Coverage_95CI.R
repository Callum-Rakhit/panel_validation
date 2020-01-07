# TODO(Callum)
#  - Make the graph

##### Load/Install relevant packages #####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("ggplot2", "nlme", "devtools", "tidyverse", "reshape2", "wesanderson"))

# Amplicon/coverage/sample data
amplicon_coverage_filenames <- Sys.glob(paths = "/mnt/shared_data/work/three_runs_together/*default_amplicon_coverage")
amplicon_coverage_filenames <- Sys.glob(paths = "/mnt/shared_data/work/metrics_extraction_for_validation_4_samples/*default_amplicon_coverage")

amplicon_coverage_filenames <- amplicon_coverage_filenames[1:length(
  amplicon_coverage_filenames)-1]  # Remove last element in list (undetermined)
amplicon_coverage_sampleIDs <- paste0(basename(amplicon_coverage_filenames))
amplicon_coverage_list <- lapply(amplicon_coverage_filenames, function(i){read.table(file = i, header = T)})
amplicon_coverage_melted <- do.call(rbind, amplicon_coverage_list)
amplicon_coverage_melted$id <- factor(rep(amplicon_coverage_sampleIDs, each = sapply(amplicon_coverage_list, nrow)))
rm(list = c("amplicon_coverage_filenames", "amplicon_coverage_sampleIDs", "amplicon_coverage_list"))

bed_file <- read.table(file = "/mnt/shared_data/snappy/snappy/tas/ccp1/ccp1_primer.bed", header = F, )

colnames(bed_file) <- c("CHR", "START", "STOP", "PRIMER", "SCORE", "STRAND")
colnames(amplicon_coverage_melted) <- c("PRIMER", "COVERAGE", "ID")

manhattan_data <- merge(x = amplicon_coverage_melted, y = bed_file, by = "PRIMER")

manhattan.plot(as.numeric(manhattan_data$CHR), manhattan_data$START, manhattan_data$COVERAGE)

sort(unique(as.numeric(manhattan_data$CHR)))
class(manhattan_data$CHR)

sample$BP <- factor(manhattan_data$START, levels = manhattan_data[ !duplicated(
  manhattan_data[,"START"]), "START"][order(manhattan_data[!duplicated(manhattan_data[ , "START"]),  "CHR"] )])

#######################################################################################################################

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

subset_manhattan <- manhattan_data[manhattan_data$ID == "18F-199_S1_default_amplicon_coverage",]

manhplot <- ggplot(manhattan_data, 
                   aes(x = manhattan_data$BPcum, y = manhattan_data$COVERAGE, 
                       color = as.factor(manhattan_data$CHR), size = manhattan_data$COVERAGE)) +
  geom_smooth() +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(250, 750)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = "Chromsome", y = "Coverage") + 
  theme_minimal() +
  
  # geom_smooth(data = manhattan_data[manhattan_data$ID == "18F-199_S1_default_amplicon_coverage",], 
  #             aes(x = manhattan_data[manhattan_data$ID == "18F-199_S1_default_amplicon_coverage",]$BPcum, 
  #                 y = manhattan_data[manhattan_data$ID == "18F-199_S1_default_amplicon_coverage",]$COVERAGE, 
  #                 color = as.factor(manhattan_data[manhattan_data$ID == "18F-199_S1_default_amplicon_coverage",]$CHR),
  #                 # color = "one big line",
  #                 size = manhattan_data[manhattan_data$ID == "18F-199_S1_default_amplicon_coverage",]$COVERAGE)) +
  
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

print(manhplot)

manhplot <- ggplot(subset_manhattan, 
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

print(manhplot)

##### Proportionality plot #####

bed_file <- read.table(file = "/mnt/shared_data/snappy/snappy/tas/ccp1/ccp1_primer.bed", header = F, )
colnames(bed_file) <- c("CHR", "START", "STOP", "PRIMER", "SCORE", "STRAND")
colnames(example) <- c("PRIMER", "blah", "bleh", "divergance")
example <- merge(x = example, y = bed_file, by = "PRIMER")
sample$BP <- factor(manhattan_data$START, levels = manhattan_data[ !duplicated(
  manhattan_data[,"START"]), "START"][order(manhattan_data[!duplicated(manhattan_data[ , "START"]),  "CHR"] )])




manhattan_data[manhattan_data$ID == "18F-199_S1_default_amplicon_coverage",]



amplicon_all_proportions <- as.data.frame(amplicon_coverage_melted$PRIMER)

for(i in sort(unique(amplicon_coverage_melted$id))){
  df <- amplicon_coverage_melted[(amplicon_coverage_melted$id == i), ]
  # options(scipen = 999)
  df$i <- df$BARCODE_COUNT/sum(df$BARCODE_COUNT)
  amplicon_all_proportions <- cbind(amplicon_all_proportions, df$i)
}

amplicon_all_proportions_melted <- melt(amplicon_all_proportions)

class(amplicon_all_proportions_melted$`amplicon_coverage_melted$PRIMER`)

options(scipen = 999)

amplicon_median_proportion <- amplicon_all_proportions_melted %>% 
  group_by(amplicon_all_proportions_melted$`amplicon_coverage_melted$PRIMER`) %>% 
  summarise(average = median(value))

example <- cbind(amplicon_median_proportion, df$i)
example$divergance <- example$average - example$`df$i`
colnames(example) <- c("PRIMER", "blah", "bleh", "divergance")







# Single_Sample_Proportional_Read_Depth
# 13.333 inch by 7.5

example$divergance2

ifelse(example$divergance>0.025, as.character(example$primer), '')
  
ggplot(example, aes(x = primer, y = divergance, color = divergance, label = primer)) +
  geom_point() +
  labs(x = "Primer", y = "Proportional increase/decrease in barcoded reads") + 
  theme_minimal() +
  geom_text(aes(label = ifelse(divergance>0.025, as.character(primer), '')), hjust = 0, vjust = 0) +
  scale_y_continuous(limits = c(-0.001, 0.004)) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_color_manual(breaks = c("-0.001", "0", "0.001"), values = wes_palettes$Darjeeling1) +
  # geom_tile("The proportion of molecularly barcoded reads for each primer in a single sample vs the median expected values") +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, size = 0.01, vjust = 0.5)
  )

ggf

























# Make proportionality Manhatten

example$CHR <- as.numeric(example$CHR)
example <- example[order(example$CHR),]
example$CHR <- ifelse(example$CHR == 23, "X", example$CHR)

nCHR <- length(unique(example$CHR))
example$BPcum <- NA
s <- 0
nbp <- c()

for (i in unique(example$CHR)){
  nbp[i] <- max(example[example$CHR == i,]$START)
  example[example$CHR == i, "BPcum"] <- example[example$CHR == i, "START"] + s
  s <- s + nbp[i]
}

axis.set <- example %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
# ylim <- max(example$divergance) + 2 

ggplot(example, aes(x = example$BPcum, y = example$divergance, color = as.factor(example$CHR), label = example$PRIMER)) +
  geom_point(size = 2) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(limits = c(-0.0005, 0.004)) +
  scale_color_manual(values = rep(wes_palettes$Darjeeling1, 5)) +
  labs(x = "Chromsome", y = "Proportional increase/decrease in barcoded reads") + 
  geom_text(aes(label = ifelse(divergance > 0.00025, as.character(PRIMER), '')), hjust = -0.05, vjust = 0) +
  theme_minimal() +
  ggtitle("The proportion of molecularly barcoded reads for each primer in a single sample vs the median expected values") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
  )

sss





