#!/usr/bin/env bash
LOCATION=$1  # Specify location as first argument
FILENAMESUFFIX1="_R1.fastq.gz"
FILENAMESUFFIX2="_R2.fastq.gz"

# list unique sample names alongside R1/R2 status
for i in *R1_001.fastq.gz; do cat $i >> "${i%_L00**}_R1.fastq.gz"; done
for i in *R2_001.fastq.gz; do cat $i >> "${i%_L00**}_R2.fastq.gz"; done

# Make files executable/owned by all
chmod 777 *

# Make a directory for the sample
for i in *$FILENAMESUFFIX1; do mkdir "${i%_R1.fastq.gz}"; done

# Place the relevant R1 files in the relevant directories
for i in *$FILENAMESUFFIX1; do mv "$i" "${i%_R1.fastq.gz}"; done

# Place the relevant R2 files in the relevant directories
for i in *$FILENAMESUFFIX2; do mv "$i" "${i%_R2.fastq.gz}"; done
