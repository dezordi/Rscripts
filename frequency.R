#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd(".")

library(ggplot2)

table_file <- paste(args[1])
output_pdf <- paste(table_file,".stacked.pdf")
output_pdf <- gsub(" ","", output_pdf, fixed = TRUE)

table_data <- read.csv(file = table_file,sep='\t')
df <- data.frame(table_data)
df <- df[order(df$position),]

df$min_freq <- (df$minor_depth/(df$minor_depth + df$major_depth))
df$maj_freq <- (df$major_depth/(df$minor_depth + df$major_depth))

pos <- c(rbind(df$position, df$position))
depth <- c(rbind(df$maj_freq, df$min_freq ))
nuc  <- c(rbind(df$major_nuc, df$minor_nuc))
order <- c(1:(2*(nrow(df))))
tags <- c(rbind(rep('major',each=nrow(df)),rep('minor',each=nrow(df))))
genome <- unique(df$genome)

df2 <- data.frame(order,pos,depth,nuc,tags)
df2$pos <- as.character(as.numeric(df2$pos))

pdf(file = output_pdf,
    width = 8, 
    height = 3)

ggplot(df2, aes(fill=tags, y=depth, x=reorder(pos,order), label = nuc)) + 
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), angle = 90) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Genomic position", y = "Allele frequency") +
  ggtitle(genome)
dev.off()