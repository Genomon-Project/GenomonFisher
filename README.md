# GenomonFisher
Genomon fisher exact test mutation caller

## Dependency
Python (>= 2.7, >= 3.7), pysam, scipy, builtins,
samtools

## Install

```
git clone https://github.com/Genomon-Project/GenomonFisher.git
cd GenomonFisher
python setup.py build
python setup.py install
```

## Run
Disease sample vs. Control sample Comparison
```
fisher comparison [-h] -1 Disease.bam -2 Control.bam -o Output_text
                  -r Reference_genome_in_fasta_format -s Samtools_path
                  [-q Mapping_quality_threshold (30)] [-Q Base_quality_threshold (15)]
                  [-m Minimum_amount_of_disease_allele_frequency (0.07)]
                  [-M Maximum_amount_of_control_allele_frequency (0.1)]
                  [-f Fisher_thres_hold (0.05)] [-d Minimum_depth (10)] 
                  [-v Minimum_variant_read (4)] 
                  [-e (Print_header)] [-g Log_file] [-l Log_level (DEBUG)]

```
Single sample mutation calling
```
fisher single [-h] -1 target.bam -o Output.text 
              -r Reference_genome_in_fasta_format -s Samtools_path
              [-q Mapping_quality_threshold (30)] [-Q Base_quality_threshold (15)]
              [-m Minimum_amount_of_allele_frequency (0.07)]
              [-p 10_percent_posterior_quantile_threshold (0.02)]
              [-d Minimum_depth (10)]
              [-v Minimum_variant_read (4)] 
              [-e (Print_header)] [-g Log_file] [-l Log_level (DEBUG)]

```

