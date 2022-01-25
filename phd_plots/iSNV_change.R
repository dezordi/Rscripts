library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(stringr)

setwd('.')
brute_intrahost <- read.csv('all_intrahosts.csv', sep = ',')
sample_list <- read.csv('amostra.txt', header = FALSE)
coinfec <- c("IAM1812","IAM3269","IAM3285","IAM3301","IAM3337","IAM3349","IAM3691",
             "IAM3693","IAM3695","IAM3706","IAM3712","IAM3714","IAM3828")

brute_intrahost$target <- brute_intrahost$GENOME %in% unlist(sample_list)

target_intrahost <- brute_intrahost %>%
  filter(target == TRUE)
target_intrahost$coinfec <- target_intrahost$GENOME %in% unlist(coinfec)

write.csv(target_intrahost, "./target_intrahosts.csv")

target_intrahost_mod <- target_intrahost %>%
  filter(coinfec == FALSE, MAJOR != '-', MINOR != '-', str_length(MAJOR) == 1, str_length(MINOR) == 1) %>%
  mutate(MAJOR = replace(MAJOR, MAJOR=="T","U"),
         MINOR = replace(MINOR, MINOR=="T","U")) %>%
  unite(change,c(MAJOR,MINOR),sep='>') %>%
  select(GENOME, POS, change)

change_order <- c("A>C","C>A","A>G","G>A","A>U","U>A","C>U","U>C","C>G","G>C","G>U","U>G")

change_count_mvs <-target_intrahost_mod %>%
  mutate(change_count = length(change)) %>%
  group_by(change) %>%
  mutate(change_freq_codetec = (n()/change_count)*100)

change_count_mvs <- unique(select(change_count_mvs,change, change_freq_codetec))
change_count_mvs <- change_count_mvs %>%
  arrange(factor(change, levels = change_order))

ggplot(data=change_count_mvs, mapping = aes(x=factor(change, level = change_order), change_freq_codetec)) +
  geom_bar(stat = 'identity', alpha = 0.6) + 
  labs(x = 'Mutação',
       y="Frequência das Mutações (%)")