#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd(".")

library(lubridate)
library(dplyr)


date_file_path <- paste(args[1])
output_name <- paste(args[2])

date_file <- read.csv(date_file_path, sep="\t", header = FALSE)
colnames(date_file)[1] = "name"
colnames(date_file)[2] = "date"

date_file <- date_file %>%
  mutate(date = decimal_date(ymd(date)))

write.csv(date_file, output_name, row.names = FALSE, quote = FALSE)
