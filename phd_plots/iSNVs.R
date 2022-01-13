library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd('/media/dezordi/34ECABD5066A98F6/01_Trabalho/projetos/PhD/plots/intrahost/')

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
  filter(cob_iSNV >= 95, lineage_up != 'None') %>%
  mutate(lineage_up = str_replace(lineage_up,"P.1..*","Gamma-sub"),
         lineage_up = str_replace(lineage_up,"P.1","Gamma (P.1)"),
         lineage_up = str_replace(lineage_up,"Gamma-sub","Gamma-sub (P.1.x)"),
         lineage_up = str_replace(lineage_up,"AY.*","Delta-sub (AY.x)"),
         lineage_up = str_replace(lineage_up,"B.1.617.2","Delta (B.1.617.2)"),
         lineage_up = str_replace(lineage_up,"B.1.1.7","Alpha (B.1.1.7)"),
         lineage_up = case_when(lineage_up %in% unlist(target_lineages) ~ lineage_up),
         lineage_up = replace_na(lineage_up,"Others")) %>%
  select(amostra, Genero, lineage_up, idade, iSNVs, cob_iSNV, depth_iSNV)


##boxplot de iSVNs por genero
intrahost_table_genero <- intrahost_table %>%
  filter(Genero != "Desconhecido")

intrahost_table_genero_valores <- intrahost_table_genero %>%
  group_by(Genero) %>%
  summarize(Mean = mean(iSNVs, na.rm=TRUE),
            Stdev = sd(iSNVs, na.rm=TRUE)) %>%
  mutate_at(2:3, round, 2)

masc_values <- intrahost_table_genero_valores %>%
  filter(Genero == "Masculino")
  
fem_values <- intrahost_table_genero_valores %>%
  filter(Genero == "Feminino")

ggplot(intrahost_table_genero, aes(x = Genero, y = iSNVs)) +
  geom_point(aes(colour = lineage_up, size = iSNVs),
             position = position_jitterdodge(dodge.width = 0.0, jitter.height = 0.5, jitter.width = 0.5),
             alpha = 0.7)+ 
  stat_summary(fun.data='mean_sdl', fun.args = list(mult=1), 
               geom="pointrange") +
  theme_bw() + 
  labs(x='Gênero',
       y='iSNVs')  +
  annotate("text", x='Feminino', y=43, label= paste("Média: ",fem_values$Mean,sep=' ')) +
  annotate("text", x='Masculino', y=43, label= paste("Média: ",masc_values$Mean,sep=' ')) +
  annotate("text", x='Feminino', y=39, label= paste("Desvio: ",fem_values$Stdev,sep=' ')) +
  annotate("text", x='Masculino', y=39, label= paste("Desvio: ",masc_values$Stdev,sep=' ')) +
  scale_color_manual(values = cols_lineages)


##boxplot de iSNVs por linhagem

intrahost_table_valores <- intrahost_table %>%
  group_by(lineage_up) %>%
  summarize(Mean = mean(iSNVs, na.rm=TRUE),
            Stdev = sd(iSNVs, na.rm=TRUE)) %>%
  mutate_at(2:3, round, 2)

alpha_values <- intrahost_table_valores %>%
  filter(lineage_up == "Alpha (B.1.1.7)")

b11_values <- intrahost_table_valores %>%
  filter(lineage_up == "B.1.1")

b1128_values <- intrahost_table_valores %>%
  filter(lineage_up == "B.1.1.28")

b1133_values <- intrahost_table_valores %>%
  filter(lineage_up == "B.1.1.33")

gamma_values <- intrahost_table_valores %>%
  filter(lineage_up == "Gamma (P.1)")

gammasub_values <- intrahost_table_valores %>%
  filter(lineage_up == "Gamma-sub (P.1.x)")

p2_values <- intrahost_table_valores %>%
  filter(lineage_up == "P.2")

delta_values <- intrahost_table_valores %>%
  filter(lineage_up == "Delta (B.1.617.2)")

deltasub_values <- intrahost_table_valores %>%
  filter(lineage_up == "Delta-sub (AY.x)")

outras_values <- intrahost_table_valores %>%
  filter(lineage_up == "Others")

ggplot(intrahost_table, aes(x = lineage_up, y = iSNVs)) +
  geom_point(aes(colour = lineage_up, size = iSNVs),
             position = position_jitterdodge(dodge.width = 0.0, jitter.height = 0.5, jitter.width = 0.5),
             alpha = 0.7)+ 
  stat_summary(fun.data='mean_sdl', fun.args = list(mult=1), 
               geom="pointrange") +
  theme_bw() + 
  labs(x='Linhagem',
       y='iSNVs')  +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95)) +
  scale_color_manual(values = cols_lineages) +
  annotate("text", x='Alpha (B.1.1.7)', y=43, label= paste("Média: ",alpha_values$Mean,sep=' '))+
  annotate("text", x='Alpha (B.1.1.7)', y=41, label= paste("Desvio: ",alpha_values$Stdev,sep=' ')) +
  annotate("text", x='B.1.1', y=43, label= paste("Média: ",b11_values$Mean,sep=' '))+
  annotate("text", x='B.1.1', y=41, label= paste("Desvio: ",b11_values$Stdev,sep=' ')) +
  annotate("text", x='B.1.1.28', y=43, label= paste("Média: ",b1128_values$Mean,sep=' '))+
  annotate("text", x='B.1.1.28', y=41, label= paste("Desvio: ",b1128_values$Stdev,sep=' ')) +
  annotate("text", x='B.1.1.33', y=43, label= paste("Média: ",b1133_values$Mean,sep=' '))+
  annotate("text", x='B.1.1.33', y=41, label= paste("Desvio: ",b1133_values$Stdev,sep=' ')) +
  annotate("text", x='Delta (B.1.617.2)', y=43, label= paste("Média: ",delta_values$Mean,sep=' '))+
  annotate("text", x='Delta (B.1.617.2)', y=41, label= paste("Desvio: ",delta_values$Stdev,sep=' ')) +
  annotate("text", x='Delta-sub (AY.x)', y=43, label= paste("Média: ",deltasub_values$Mean,sep=' '))+
  annotate("text", x='Delta-sub (AY.x)', y=41, label= paste("Desvio: ",deltasub_values$Stdev,sep=' ')) +
  annotate("text", x='Gamma (P.1)', y=43, label= paste("Média: ",gamma_values$Mean,sep=' '))+
  annotate("text", x='Gamma (P.1)', y=41, label= paste("Desvio: ",gamma_values$Stdev,sep=' ')) +
  annotate("text", x='Gamma-sub (P.1.x)', y=43, label= paste("Média: ",gammasub_values$Mean,sep=' '))+
  annotate("text", x='Gamma-sub (P.1.x)', y=41, label= paste("Desvio: ",gammasub_values$Stdev,sep=' ')) +
  annotate("text", x='Others', y=43, label= paste("Média: ",outras_values$Mean,sep=' '))+
  annotate("text", x='Others', y=41, label= paste("Desvio: ",outras_values$Stdev,sep=' ')) +
  annotate("text", x='P.2', y=43, label= paste("Média: ",p2_values$Mean,sep=' '))+
  annotate("text", x='P.2', y=41, label= paste("Desvio: ",p2_values$Stdev,sep=' ')) 

##boxplot de iSNVs por faixa etária

intrahost_table_idade <- intrahost_table %>%
  filter(idade!="Desconhecido",
         idade !="",
         idade !="#VALUE!") %>%
  mutate(grupo = case_when(as.numeric(idade) <= 15 ~ '0-15',
                           as.numeric(idade) <= 65 ~ '16-65',
                           as.numeric(idade) > 65 ~ '66+'))

intrahost_idade_valores <- intrahost_table_idade %>%
  group_by(grupo) %>%
  summarize(Mean = mean(iSNVs, na.rm=TRUE),
            Stdev = sd(iSNVs, na.rm=TRUE)) %>%
  mutate_at(2:3, round, 2)

values_015 <- intrahost_idade_valores %>%
  filter(grupo == "0-15")

values_1665 <- intrahost_idade_valores %>%
  filter(grupo == "16-65")

values_66mais <- intrahost_idade_valores %>%
  filter(grupo == "66+")

ggplot(intrahost_table_idade, aes(x = grupo, y = iSNVs)) +
  geom_point(aes(colour = lineage_up, size = iSNVs),
             position = position_jitterdodge(dodge.width = 0.0, jitter.height = 0.5, jitter.width = 0.5),
             alpha = 0.7)+ 
  stat_summary(fun.data='mean_sdl', fun.args = list(mult=1), 
               geom="pointrange") +
  theme_bw() + 
  labs(x='Faixa Etária',
       y='iSNVs')  +
  scale_color_manual(values = cols_lineages) +
  annotate("text", x='0-15', y=43, label= paste("Média: ",values_015$Mean,sep=' '))+
  annotate("text", x='0-15', y=41, label= paste("Desvio: ",values_015$Stdev,sep=' '))+
  annotate("text", x='16-65', y=43, label= paste("Média: ",values_1665$Mean,sep=' '))+
  annotate("text", x='16-65', y=41, label= paste("Desvio: ",values_1665$Stdev,sep=' '))+
  annotate("text", x='66+', y=43, label= paste("Média: ",values_66mais$Mean,sep=' '))+
  annotate("text", x='66+', y=41, label= paste("Desvio: ",values_66mais$Stdev,sep=' '))

##correlacao iSNV profundidade

ggscatter(intrahost_table, x = "depth_iSNV", y = "iSNVs", 
          add = "reg.line", conf.int = TRUE, 
          size = "iSNVs", alpha = 0.5,
          xlab = "Profundidade Media", ylab = "Número de iSNVs") +
  stat_cor(method = "pearson", label.x = 3500, label.y = 45) +
  theme(legend.position = 'none')
  

