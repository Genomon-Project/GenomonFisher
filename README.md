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
fisher [-h] [--version] [-1 control.bam] [-2 disease.bam] [-e]
            [-q mapping quality threshold] [-Q base quality threshold]
            [-m minimum amount of disease allele frequency]
            [-n maximum amount of control allele frequency]
            [-p 10% posterior quantile threshold] [-f fisher threshold]
            [-d minimum depth]
            [-v minimum amount of variant reads [disease]]
            [-g log file name] [-l logging level] -r reference genome in
            fasta format -s samtools_path
            output.txt

```

