# Rscripts

This repository store different Rscripts to plot different kind of biological data. The example files can be found on example_files directory.

## frequency.R

This script use the ggplot package to create a stacked barplot of genomic positions with different allele frequencies, e.g. data of intrahost variants. The tabular file should have six columns:

- genome: The name of genome, will be used as title of barplot
- position: The genomic positions
- major_depth: The depth of major allele
- minor_depth: The depth of minor allele
- major_nuc: The symbol present as major allele
- minor_nuc: The symbol present as minor allele

### usage

> Rscript frequency.R frequency_example.tsv

### result

![](example_figures/frequency_example.png)