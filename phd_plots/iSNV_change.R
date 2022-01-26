library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(stringr)

setwd('/media/dezordi/34ECABD5066A98F6/01_Trabalho/projetos/PhD/plots/intrahost/')
brute_intrahost <- read.csv('all_intrahosts.csv', sep = ',')
sample_list <- read.csv('amostra.txt', header = FALSE)
coinfec <- c("IAM1812","IAM3269","IAM3285","IAM3301","IAM3337","IAM3349","IAM3691",
             "IAM3693","IAM3695","IAM3706","IAM3712","IAM3714","IAM3828")

target_lineages <- list("Gamma (P.1)","Gamma-sub (P.1.x)",
                        "Delta (B.1.617.2)","Delta-sub (AY.x)",
                        "P.2","B.1.1.28","B.1.1.33",
                        "Alpha (B.1.1.7)","B.1.1")

cols_lineages <- c("B.1.1"='#FF66B2',"Others"='#C0C0C0',
                   "B.1.1.28"='#66FFFF',"B.1.1.33"='#66B2FF',
                   "Alpha (B.1.1.7)"='#6666FF',"P.2"='#FFFF66',
                   "Gamma (P.1)"="#FF6666","Gamma-sub (P.1.x)"='#FFB266',
                   "Delta (B.1.617.2)"='#B2FF66',"Delta-sub (AY.x)"='#66FF66',
                   "B.1.1.28"='#66FFFF',"B.1.1.33"='#66B2FF')

intrahost_table <- read.csv('./phd_tables-metricas.tsv', sep = '\t', dec = ',')
intrahost_table <- intrahost_table %>%
  filter(cob_iSNV >= 95, lineage_up_major != 'None', coinfection == FALSE) %>%
  mutate(lineage_up_major = str_replace(lineage_up_major,"P.1..*","Gamma-sub"),
         lineage_up_major = str_replace(lineage_up_major,"P.1","Gamma (P.1)"),
         lineage_up_major = str_replace(lineage_up_major,"Gamma-sub","Gamma-sub (P.1.x)"),
         lineage_up_major = str_replace(lineage_up_major,"AY.*","Delta-sub (AY.x)"),
         lineage_up_major = str_replace(lineage_up_major,"B.1.617.2","Delta (B.1.617.2)"),
         lineage_up_major = str_replace(lineage_up_major,"B.1.1.7","Alpha (B.1.1.7)"),
         lineage_up_major = case_when(lineage_up_major %in% unlist(target_lineages) ~ lineage_up_major),
         lineage_up_major = replace_na(lineage_up_major,"Others")) %>%
  select(amostra, Genero, lineage_up_major, collection_date, idade, iSNVs, cob_iSNV, depth_iSNV)

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

target_intrahost_mod <- target_intrahost_mod %>%
  left_join(intrahost_table, c("GENOME" = "amostra"))

target_intrahost_mod <- target_intrahost_mod %>% 
  filter(is.na(lineage_up_major) == FALSE)


# plot geral
change_count_mvs <-target_intrahost_mod %>%
  filter(lineage_up_major != c("Delta-sub (AY.x)","Delta (B.1.617.2)")) %>%
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

##plots por linhagem

target_lineages <- list("Gamma (P.1)","Gamma-sub (P.1.x)",
                        "Delta (B.1.617.2)","Delta-sub (AY.x)",
                        "P.2","B.1.1.28","B.1.1.33",
                        "Alpha (B.1.1.7)","B.1.1")

cols_lineages <- c("B.1.1"='#FF66B2',"Others"='#C0C0C0',
                   "B.1.1.28"='#66FFFF',"B.1.1.33"='#66B2FF',
                   "Alpha (B.1.1.7)"='#6666FF',"P.2"='#FFFF66',
                   "Gamma (P.1)"="#FF6666","Gamma-sub (P.1.x)"='#FFB266',
                   "Delta (B.1.617.2)"='#B2FF66',"Delta-sub (AY.x)"='#66FF66')

# plot gamma
change_count_mvs <-target_intrahost_mod %>%
  filter(lineage_up_major == c("Gamma (P.1)","Gamma-sub (P.1.x)")) %>%
  mutate(change_count = length(change)) %>%
  group_by(change) %>%
  mutate(change_freq_codetec = (n()/change_count)*100)

change_count_mvs <- unique(select(change_count_mvs,change, change_freq_codetec))
change_count_mvs <- change_count_mvs %>%
  arrange(factor(change, levels = change_order))

plot_gamma <- ggplot(data=change_count_mvs, 
                     mapping = aes(x=factor(change, level = change_order),change_freq_codetec)) +
  geom_bar(stat = 'identity', alpha = 0.6,  fill = "#FF6666") + 
  labs(x = 'Mutação',
       y="") +
  scale_x_discrete(limits = change_order)

# plot delta
change_count_mvs <-target_intrahost_mod %>%
  filter(lineage_up_major == c("Delta-sub (AY.x)","Delta (B.1.617.2)")) %>%
  mutate(change_count = length(change)) %>%
  group_by(change) %>%
  mutate(change_freq_codetec = (n()/change_count)*100)

change_count_mvs <- unique(select(change_count_mvs,change, change_freq_codetec))
change_count_mvs <- change_count_mvs %>%
  arrange(factor(change, levels = change_order))

plot_delta <- ggplot(data=change_count_mvs, 
                         mapping = aes(x=factor(change, level = change_order),change_freq_codetec)) +
  geom_bar(stat = 'identity', alpha = 0.6,  fill = "#B2FF66") + 
  labs(x = 'Mutação',
       y="")+
  scale_x_discrete(limits = change_order)

# plot B11
change_count_mvs <-target_intrahost_mod %>%
  filter(lineage_up_major == 'B.1.1') %>%
  mutate(change_count = length(change)) %>%
  group_by(change) %>%
  mutate(change_freq_codetec = (n()/change_count)*100)

change_count_mvs <- unique(select(change_count_mvs,change, change_freq_codetec))
change_count_mvs <- change_count_mvs %>%
  arrange(factor(change, levels = change_order))

plot_b11 <- ggplot(data=change_count_mvs, 
                     mapping = aes(x=factor(change, level = change_order),change_freq_codetec)) +
  geom_bar(stat = 'identity', alpha = 0.6,  fill = "#FF66B2") + 
  labs(x = '',
       y="")+
  scale_x_discrete(limits = change_order)

# plot outras
change_count_mvs <-target_intrahost_mod %>%
  filter(lineage_up_major == "Others") %>%
  mutate(change_count = length(change)) %>%
  group_by(change) %>%
  mutate(change_freq_codetec = (n()/change_count)*100)

change_count_mvs <- unique(select(change_count_mvs,change, change_freq_codetec))
change_count_mvs <- change_count_mvs %>%
  arrange(factor(change, levels = change_order))

plot_outras <- ggplot(data=change_count_mvs, 
                     mapping = aes(x=factor(change, level = change_order),change_freq_codetec)) +
  geom_bar(stat = 'identity', alpha = 0.6,  fill = "#C0C0C0") + 
  labs(x = '',
       y="")+
  scale_x_discrete(limits = change_order)

# plot b1128
change_count_mvs <-target_intrahost_mod %>%
  filter(lineage_up_major == "B.1.1.28") %>%
  mutate(change_count = length(change)) %>%
  group_by(change) %>%
  mutate(change_freq_codetec = (n()/change_count)*100)

change_count_mvs <- unique(select(change_count_mvs,change, change_freq_codetec))
change_count_mvs <- change_count_mvs %>%
  arrange(factor(change, levels = change_order))

plot_b1128 <- ggplot(data=change_count_mvs, 
                     mapping = aes(x=factor(change, level = change_order),change_freq_codetec)) +
  geom_bar(stat = 'identity', alpha = 0.6,  fill = "#66FFFF") + 
  labs(x = "",
       y="")+
  scale_x_discrete(limits = change_order)

# plot b1133
change_count_mvs <-target_intrahost_mod %>%
  filter(lineage_up_major == "B.1.1.33") %>%
  mutate(change_count = length(change)) %>%
  group_by(change) %>%
  mutate(change_freq_codetec = (n()/change_count)*100)

change_count_mvs <- unique(select(change_count_mvs,change, change_freq_codetec))
change_count_mvs <- change_count_mvs %>%
  arrange(factor(change, levels = change_order))

plot_b1133 <- ggplot(data=change_count_mvs, 
                     mapping = aes(x=factor(change, level = change_order),change_freq_codetec)) +
  geom_bar(stat = 'identity', alpha = 0.6,  fill = "#66B2FF") + 
  labs(x = '',
       y="")+
  scale_x_discrete(limits = change_order)


ggarrange(plot_outras, plot_b11,
          plot_b1128, plot_b1133,
          plot_gamma, plot_delta, nrow = 3, ncol = 2, align = 'v', labels = c('A','B','C','D','E','F'))

