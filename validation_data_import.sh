# Get the coverage metrics
csplit --digits=2 --quiet --prefix=outfile_S1 analysis/20190727_183445/18F-199_S1/default/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S2 analysis/20190727_183445/18F-208_S2/default/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S3 analysis/20190727_183445/18F-181_S3/default/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S4 analysis/20190727_183445/18F-381_S4/default/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S1_10M analysis/20190727_183445/18F-199_S1/10M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S2_10M analysis/20190727_183445/18F-208_S2/10M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S3_10M analysis/20190727_183445/18F-181_S3/10M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S4_10M analysis/20190727_183445/18F-381_S4/10M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S1_5M analysis/20190727_183445/18F-199_S1/5M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S2_5M analysis/20190727_183445/18F-208_S2/5M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S3_5M analysis/20190727_183445/18F-181_S3/5M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S4_5M analysis/20190727_183445/18F-381_S4/5M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S1_2.5M analysis/20190727_183445/18F-199_S1/2.5M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S2_2.5M analysis/20190727_183445/18F-208_S2/2.5M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S3_2.5M analysis/20190727_183445/18F-181_S3/2.5M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S4_1.25M analysis/20190727_183445/18F-381_S4/2.5M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S1_1.25M analysis/20190727_183445/18F-199_S1/1.25M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S2_1.25M analysis/20190727_183445/18F-208_S2/1.25M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S3_1.25M analysis/20190727_183445/18F-181_S3/1.25M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"
csplit --digits=2 --quiet --prefix=outfile_S4_1.25M analysis/20190727_183445/18F-381_S4/1.25M_resample/metrics/*.metrics.targetcoverage "/#/+1" "{*}"

for file in $(ls ~/work/analysis/20190727_183445/); do echo $file; done

# Get the coverage per primer pair metrics
# awk '{print $1, $3}' analysis/20190727_183445/18F-199_S1/default/metrics/63f8c42f.counts.PRIMER > sample_1_primer_barcode_coverage_inc_amplicon
# awk '{print $1, $3}' analysis/20190727_183445/18F-208_S2/default/metrics/912d7e9c.counts.PRIMER > sample_2_primer_barcode_coverage_inc_amplicon
# awk '{print $1, $3}' analysis/20190727_183445/18F-181_S3/default/metrics/c5763424.counts.PRIMER > sample_3_primer_barcode_coverage_inc_amplicon
# awk '{print $1, $3}' analysis/20190727_183445/18F-381_S4/default/metrics/2c132364.counts.PRIMER > sample_4_primer_barcode_coverage_inc_amplicon

for folder in $(ls folder location); do awk '{print $1, $3}' $file > output; done

# Get the VAF/variant information
# bcftools query -f '%CHROM:%POS %AF\n' analysis/20190727_183445/18F-199_S1/default/variants/63f8c42f.VD.filtered.vcf.gz > S1_VAF
# bcftools query -f '%CHROM:%POS %AF\n' analysis/20190727_183445/18F-208_S2/default/variants/912d7e9c.VD.filtered.vcf.gz > S2_VAF
# bcftools query -f '%CHROM:%POS %AF\n' analysis/20190727_183445/18F-181_S3/default/variants/c5763424.VD.filtered.vcf.gz > S3_VAF
# bcftools query -f '%CHROM:%POS %AF\n' analysis/20190727_183445/18F-381_S4/default/variants/2c132364.VD.filtered.vcf.gz > S4_VAF
