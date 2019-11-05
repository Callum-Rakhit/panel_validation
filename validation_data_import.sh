# TODO(Callum)
#   Pipe instead of copying, altering and removing files
#   Remove hardcoded references to a particular run

# Create a folder to dump the output into
mkdir -p $SNAPPYWORK/metrics_extraction_for_validation

##### Get the coverage metrics #####

# Copy, rename and split the relevant metrics file for each downsampling analysis performed for each sample
function aggregate_coverage_files {
  for i in $(ls $1); do for x in $(ls $1/$i); 
  do cp $1/$i/$x/metrics/*.metrics.targetcoverage $2 && 
  for y in $1/$i/$x/metrics/*.metrics.targetcoverage; 
  do z=$(basename $y) && a=$(echo $z | cut -c1-8) && 
  mv $2/$a* $2/${i}_${x}_metrics && 
  sed -n '/# NOT COVERED/,/# LOW COVERAGE (<200X)/ p' $2/${i}_${x}_metrics | head -n -2 | sed '1d' > $2/${i}_${x}_zero_coverage && 
  sed -n '/# LOW COVERAGE (<200X)/,/# PERCENT COVERED/ p' $2/${i}_${x}_metrics | head -n -2  | sed '1d' > $2/${i}_${x}_low_coverage && 
  sed -e '1,/# PERCENT COVERED/ d' $2/${i}_${x}_metrics | head -n -1 > $2/${i}_${x}_coverage_percentages;
  done; 
  done; 
  done
}


##### Get the coverage per primer pair metrics #####

# Copy, rename and split the amplicon coverage metrics file
function aggregate_amplion_coverage {
  for i in $(ls $1); do for x in $(ls $1/$i); 
  do cp $1/$i/$x/metrics/*.counts.PRIMER $2 &&
  for y in $1/$i/$x/metrics/*.counts.PRIMER;
  do z=$(basename $y) && a=$(echo $z | cut -c1-8) &&
  mv $2/$a* $2/${i}_${x}_amplicon_coverage_temp &&
  awk '{print $1, $3}' $2/${i}_${x}_amplicon_coverage_temp | sed '1,6d' | head -n -1 > $2/${i}_${x}_amplicon_coverage;
  done;
  done;
  done
}


##### Get the VAF/variant information  #####

# Copy, rename and split the amplicon coverage metrics file
function aggregate_VAF {
  for i in $(ls $1); do for x in $(ls $1/$i); 
  do cp $1/$i/$x/variants/*.VD.filtered.vcf.gz $2 &&
  for y in $1/$i/$x/variants/*.VD.filtered.vcf.gz;
  do z=$(basename $y) && a=$(echo $z | cut -c1-8) &&
  mv $2/$a* $2/${i}_${x}_VAF_frequencies_verbose &&
  bcftools query -f '%CHROM:%POS %AF %DP\n' $2/${i}_${x}_VAF_frequencies_verbose > $2/${i}_${x}_VAF_frequencies_bare;
  done;
  done;
  done
}

##### Run functions #####

# Final function to call all previous functions by iterating for the sample folders
function make_plots {
  aggregate_coverage_files $1 $2
  aggregate_amplion_coverage  $1 $2
  aggregate_VAF $1 $2
  # rm *temp* && rm *verbose* && rm *metrics*
}

make_plots $SNAPPYWORK/analysis/20191003_190447/ $SNAPPYWORK/metrics_extraction_for_validation/