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
usage: fisher comparison [-h] -1 BAM1 -2 BAM2 [-a SAMPLE1] [-b SAMPLE2] -o
                         OUTPUT -r REF_FA -s SAMTOOLS_PATH
                         [-S SAMTOOLS_PARAMS] [-Q BASE_QUALITY]
                         [-m MIN_ALLELE_FREQ] [-M MAX_ALLELE_FREQ]
                         [-f FISHER_VALUE] [-d MIN_DEPTH]
                         [-v MIN_VARIANT_READ] [-R REGION] [-L REGIONS]
                         [-O {vcf,anno}] [-e] [-g LOG_FILE] [-l LOG_LEVEL]

```
Single sample mutation calling
```
usage: fisher single [-h] -1 BAM1 [-a SAMPLE1] -o OUTPUT -r REF_FA -s
                     SAMTOOLS_PATH [-S SAMTOOLS_PARAMS] [-Q BASE_QUALITY]
                     [-m MIN_ALLELE_FREQ] [-p POST_10_Q] [-d MIN_DEPTH]
                     [-v MIN_VARIANT_READ] [-R REGION] [-L REGIONS]
                     [-O {vcf,anno}] [-e] [-g LOG_FILE] [-l LOG_LEVEL]
```

