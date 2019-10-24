##### Load/Install relevant packages #####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("ggplot2", "reshape2", "wesanderson", "tidyverse", "scales", "doParallel", 
              "devtools", "dplyr", "gtable", "grid", "gridExtra", "data.table", "rlist"))

##### Load Data #####

# Get the observed VAFs for the horizon controls
variant_frequency_filenames <- Sys.glob(paths = "/mnt/shared_data/work/three_runs_together/*VAF_frequencies_bare")
variant_frequency_list <- list()
invisible(lapply(variant_frequency_filenames, function(i){read.table(file = i, header = F) %>%
    {. ->> variant_frequency_list[[paste0(basename(i))]]} }))

##### Process Data #####

# Function to plot the VAFs
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
    # xlim(0, .3) +
    # ylim(0, .3) +
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

# Get the expected VAFs for the controls samples run
known_VAFs <- read.csv("~/panel_validation/Horizon_Control_Locations_AF.csv")
known_VAFs$Exact_Location <- gsub(" ", "", x = known_VAFs$Exact_Location)

# Filter the observations based on those found in the horizon controls 
VAF <- list()

for(sample.filename in names(variant_frequency_list)){
  sample.name.inc.conc <- (gsub('_r.*', '', sample.filename))
  sample.name.inc.conc <- (gsub('_d.*', '', sample.name.inc.conc))
  sample.name.inc.conc <- (gsub("_S.*_", "_", sample.name.inc.conc))
  sample.name <- (gsub('_S.*', '', sample.filename))
  sample.name <- (gsub("18F1", "18F-1", sample.name))
  sample.name <- (gsub("18F2", "18F-2", sample.name))
  sample.name <- (gsub("18F3", "18F-3", sample.name))
  sample.name <- (gsub("18F-18F", "18F", sample.name))
  sample.name <- (gsub('9-20', '9', sample.name))
  sample.name <- (gsub('9-40', '9', sample.name))
  sample.name <- (gsub('9-60', '9', sample.name))
  sample.name <- (gsub('9-80', '9', sample.name))
  sample.name <- (gsub('3-20', '3', sample.name))
  sample.name <- (gsub('3-40', '3', sample.name))
  sample.name <- (gsub('3-60', '3', sample.name))
  sample.name <- (gsub('3-80', '3', sample.name))
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
  VAF[[sample.name.inc.conc]] <- VAFs_HorizonID
}

known_VAFs$Exact_Location

# Plot 4a - Horizon control HD200
VAFPlot(VAF$HD200_S10, "HD200", VAF$HD200_S10$HD200)
VAFPlot(VAF$HD200_S11, "HD200", VAF$HD200_S11$HD200)
VAFPlot(VAF$HD200_5M, "HD200", VAF$HD200_5M$HD200)
VAFPlot(VAF$HD200_2.5M, "HD200", VAF$HD200_2.5M$HD200)
VAFPlot(VAF$HD200_1.25M, "HD200", VAF$HD200_1.25M$HD200)

# Plot 4b - Horizon control HD701
VAFPlot(VAF$HD701_S10, "HD701", VAF$HD701_S10$HD701)
VAFPlot(VAF$HD701_S11, "HD701", VAF$HD701_S11$HD701)
VAF$HD701_S11
VAFPlot(VAF$HD701_5M, "HD701", VAF$HD701_5M$HD701)
VAFPlot(VAF$HD701_2.5M, "HD701", VAF$HD701_2.5M$HD701)
VAFPlot(VAF$HD701_1.25M, "HD701", VAF$HD701_1.25M$HD701)

# Plot 4c - Horizon control HD798
VAFPlot(VAF$HD798_S9, "HD798", VAF$HD798_S9$HD798)
VAFPlot(VAF$HD798_5M, "HD798", VAF$HD798_5M$HD798)
VAFPlot(VAF$HD798_2.5M, "HD798", VAF$HD798_2.5M$HD798)
VAFPlot(VAF$HD798_1.25M, "HD798", VAF$HD798_1.25M$HD798)
VAF$HD798_1.25M

HD798_comparison <- merge(variant_frequency_list$HD701_S10_default_VAF_frequencies_bare,
                          variant_frequency_list$HD701_S10_1.25M_resample_VAF_frequencies_bare,
                          by = "V1", all = T)
colnames(HD798_comparison) <- c("Location", "AF", "Downsampled_AF")

HD798_downsampled_rejects <- HD798_comparison[(is.na(HD798_comparison$AF)),]
hist(HD798_downsampled_rejects$Downsampled_AF, breaks = 20)

HD798_comparison <- HD798_comparison[!(is.na(HD798_comparison$AF)),]

VAFPlot(HD798_comparison, "HD798 10M total reads vs 1.25M", HD798_comparison$Downsampled_AF)

sum(is.na(HD798_comparison$Downsampled_AF))
sum(is.na(HD798_comparison$AF))

summary(HD798_comparison$Downsampled_AF)
hist(HD798_comparison$Downsampled_AF, breaks = 40)
hist(HD798_comparison$AF, breaks = 40)

##### Total reads versus VAF #####





