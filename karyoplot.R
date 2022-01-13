#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd(".")

bam_file <- paste(args[1])
cytoband_file <- paste(args[2])
genome_name <- paste(args[3])
cytoband_file_table <- read.table(cytoband_file)
reference_name <- cytoband_file_table[2,1]

output_pdf <- paste(genome_name,".pdf")
output_pdf <- gsub(" ","", output_pdf, fixed = TRUE)


library(karyoploteR)
library(Rsamtools)


pdf(file = output_pdf,
    width = 8,
    height = 3)

custom.genome <- toGRanges(data.frame(chr=reference_name, start=1, end=29903))
custom.cytobands <- toGRanges(cytoband_file)
zoom.region <- toGRanges(data.frame(reference_name,1,29903))
kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands,plot.type = 1, labels.plotter = NULL, zoom = zoom.region, clipping=TRUE)
kp <- kpPlotBAMCoverage(kp, data=bam_file, ymax=500)
kpAddChromosomeNames(kp, chr.names=genome_name, cex=0.8)
kpAddBaseNumbers(tick.dist=2000,kp,tick.len=5, add.units=TRUE, digits=2, minor.ticks=TRUE, minor.tick.dist=200, minor.tick.len=2, cex=0.5, tick.col=NULL, minor.tick.col='black', clipping=TRUE)
kpAddCytobandLabels(kp, force.all=TRUE,col="black", cex=0.35, srt=90)
kpAxis(kp, ymax=500, cex=0.8)

dev.off()