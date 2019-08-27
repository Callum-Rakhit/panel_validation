##### Get the coverage metrics #####

# Copy, rename and split the relevant metrics file for each downsampling analysis performed for each sample
for i in $(ls analysis/20190727_183445/); 
  do for x in $(ls analysis/20190727_183445/$i); 
    do cp analysis/20190727_183445/$i/$x/metrics/*.metrics.targ* ./metrics_extraction_for_validation/
      && for y in analysis/20190727_183445/$i/$x/metrics/*.metrics.targ*;
        do z=$(basename $y) && a=$(echo $z | cut -c1-8) 
        && mv ./metrics_extraction_for_validation/$a* ./metrics_extraction_for_validation/${i}_${x}_metrics
        && csplit --digits=2 --quiet --prefix=./metrics_extraction_for_validation/${i}_${x}_metrics_split ./metrics_extraction_for_validation/${i}_${x}_metrics "/#/+1" "{*}";
    done;
  done;
done

##### Get the coverage per primer pair metrics #####

# awk '{print $1, $3}' analysis/20190727_183445/18F-199_S1/default/metrics/63f8c42f.counts.PRIMER > sample_1_primer_barcode_coverage_inc_amplicon
# awk '{print $1, $3}' analysis/20190727_183445/18F-208_S2/default/metrics/912d7e9c.counts.PRIMER > sample_2_primer_barcode_coverage_inc_amplicon
# awk '{print $1, $3}' analysis/20190727_183445/18F-181_S3/default/metrics/c5763424.counts.PRIMER > sample_3_primer_barcode_coverage_inc_amplicon
# awk '{print $1, $3}' analysis/20190727_183445/18F-381_S4/default/metrics/2c132364.counts.PRIMER > sample_4_primer_barcode_coverage_inc_amplicon

for folder in $(ls folder location); do awk '{print $1, $3}' $file > output; done

##### Get the VAF/variant information  #####
# bcftools query -f '%CHROM:%POS %AF\n' analysis/20190727_183445/18F-199_S1/default/variants/63f8c42f.VD.filtered.vcf.gz > S1_VAF
# bcftools query -f '%CHROM:%POS %AF\n' analysis/20190727_183445/18F-208_S2/default/variants/912d7e9c.VD.filtered.vcf.gz > S2_VAF
# bcftools query -f '%CHROM:%POS %AF\n' analysis/20190727_183445/18F-181_S3/default/variants/c5763424.VD.filtered.vcf.gz > S3_VAF
# bcftools query -f '%CHROM:%POS %AF\n' analysis/20190727_183445/18F-381_S4/default/variants/2c132364.VD.filtered.vcf.gz > S4_VAF
