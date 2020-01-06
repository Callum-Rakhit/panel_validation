# TODO(Callum)
#  - Make the graph

##### Load/Install relevant packages #####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("lattice", "devtools"))

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

manhattan.plot <- function(chr, pos, pvalue, sig.level = NA, annotate = NULL, ann.default = list(), 
                           should.thin = T, thin.pos.places = 2, thin.logp.places = 2, xlab = "Chromosome",
                           # ylab = expression(-log[10](p-value)),
                           # ylab = expression(p-value),
                           ylab = expression(Coverage), 
                           # col = c("gray", "darkgray"), 
                           col = c("red", "blue"),
                           panel.extra = NULL, pch = 20, cex = 0.8, ...) {
  
  if (length(chr) == 0) stop("chromosome vector is empty")
  if (length(pos) == 0) stop("position vector is empty")
  if (length(pvalue) == 0) stop("pvalue vector is empty")
  
  # Make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[, drop = T]
  }
  
  # Make sure positions are in kbp
  if (any(pos > 1e6)) pos <- pos/1e6;
  
  # Calculate absolute genomic position from relative chromosomal positions
  posmin <- tapply(pos, chr, min);
  posmax <- tapply(pos, chr, max);
  posshift <- head(c(0, cumsum(posmax)), -1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos <- function(cchr, cpos) {
    p <- posshift[as.character(cchr)] + cpos
    return(p)
  }
  
  # Parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default <- list(x = "peak", y = "peak", adj = NULL, pos = 3, offset = 0.5, 
                      col = NULL, fontface = NULL, fontsize = NULL, show = F)
  parse.label <- function(rawval, groupname) {
    r <- list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval) >= 1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times = length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]] <- list(pch = pch, col = col, cex = cex, fill = col, label = label.default)
  
  if (length(ann.settings) > 1) { 
    lcols <- trellis.par.get("superpose.symbol")$col 
    lfills <- trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch = pch, 
                              col = lcols[(i-2) %% length(lcols) + 1 ], 
                              fill = lfills[(i-2) %% length(lfills) + 1 ], 
                              cex = cex, label = label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings) <- levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate) > 1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1]) != "")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols <- user.cols - 1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label <- parse.label(annotate[[user.cols[i]]]$label, levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]] <- modifyList(ann.settings[[ann.cols[i]]], annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  # Reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      # logp=round(-log10(pvalue),thin.logp.places),
      logp = round(pvalue, thin.logp.places),
      pos = round(genpos, thin.pos.places), 
      chr = chr,
      grp = grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    # logp <- -log10(pvalue)
    logp <- pvalue
  }
  rm(pos, pvalue)
  gc()
  
  # Custom axis to print chromosome names
  axis.chr <- function(side, ...) {
    if(side == "bottom") {
      panel.axis(side = side, outside = T,
                 at = ((posmax+posmin)/2 + posshift),
                 labels=levels(chr), 
                 ticks = F, rot = 0,
                 check.overlap = F
      )
    } else if (side == "top" || side == "right") {
      panel.axis(side = side, draw.labels = F, ticks = F);
    }
    else {
      axis.default(side = side, ...);
    }
  }
  
  # Make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr <- function(x, y , ...) { 
    A <- list();
    # maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    maxy <- ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0))) + 0.5;
    A$ylim = c(0, maxy);
    A;
  }
  
  # xyplot(logp ~ genpos, chr = chr, groups = grp,
  #        axis = axis.chr, ann.settings = ann.settings,
  #        prepanel = prepanel.chr, scales = list(axs = "i"),
  #        panel = function(x, y, ..., getgenpos));
  # 
  # ggplot(logp ~ genpos)
  
  xyplot(logp ~ genpos, chr = chr, groups = grp,
         axis = axis.chr, ann.settings = ann.settings,
         prepanel = prepanel.chr, scales = list(axs="i"),
         panel = function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             # Add significance line (if requested)
             # panel.abline(h=-log10(sig.level), lty=2);
             panel.abline(h = sig.level, lty = 2);
           }
           panel.superpose(x, y, ..., getgenpos = getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x, y, getgenpos, ...)
           }
         },
         panel.groups = function(x, y, ..., subscripts, group.number) {
           A <- list(...)

           # Allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)

           # Draw labels (if requested)
           if(gs$label$show) {
             gt <- gs$label
             names(gt)[which(names(gt) == "text")] <- "labels"
             gt$show <- NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x == "peak") {gt$x <- x[peak]}
                 if(gt$x == "center") {gt$x <- center}
               }
               if (is.character(gt$y)) {
                 if(gt$y == "peak") {gt$y <- y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x <- A$getgenpos(gt$x[[1]], gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab = xlab, ylab = ylab,
         panel.extra = panel.extra, getgenpos = getGenPos, ...
  );
}

manhattan.plot(as.numeric(manhattan_data$CHR), manhattan_data$START, manhattan_data$COVERAGE)

sort(unique(as.numeric(manhattan_data$CHR)))
class(manhattan_data$CHR)

require(ggplot2)
require(nlme)

sample$BP <- factor(manhattan_data$START,
                    levels = manhattan_data[ !duplicated(manhattan_data[,"START"]), "START"][order(
                      manhattan_data[!duplicated(manhattan_data[ , "START"]),  "CHR"] )]
)

logp ~ genpos

ggplot(manhattan_data, aes(START, COVERA GE)) + geom_point() + geom_smooth()

ggplot(manhattan_data, aes(START, COVERAGE)) + geom_smooth()

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

































