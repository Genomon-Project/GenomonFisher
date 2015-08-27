# GenomonFisher
Genomon fisher exact test mutation caller

## Dependency
Python (>= 2.7), pysam

## Install

```
git clone https://github.com/Genomon-Project/GenomonFisher.git
cd GenomonFisher
python setup.py build
python setup.py install
```

**target somatic mutation candidats**: the somatic mutation candidates (should be .tsv or .vcf format).  
**target tumor bam**: the indexed bam file of the target tumor sample.  
**target normal bam**: the indexed bam file of the target normal sample.  

## Run

```
fisher comparison [-h] -1 Disease.bam -2 Control.bam -o Output_text
                  -r Reference_genome_in_fasta_format -s Samtools_path
                  [-q Mapping_quality_threshold] [-Q Base_quality_threshold]
                  [-m Minimum_amount_of_disease_allele_frequency]
                  [-M Maximum_amount_of_control_allele_frequency]
                  [-f Fisher_thres_hold] [-d Minimum_depth] 
                  [-v Minimum_variant_read] 
                  [-e (Print_header)] [-g Log_file] [-l Log_level]

```
```
fisher single [-h] -1 target.bam -o Output.text 
              -r Reference_genome_in_fasta_format -s Samtools_path
              [-q Mapping_quality_threshold] [-Q Base_quality_threshold]
              [-m Minimum_amount_of_allele_frequency]
              [-p 10_percent_posterior_quantile_threshold]
              [-d Minimum_depth]
              [-e (Print_header)] [-g Log_file] [-l Log_level]

```

