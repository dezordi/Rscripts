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

#checando normalidade por gênero
intrahost_table_genero_normalidade <- intrahost_table_genero %>%
  group_by(Genero) %>%
  mutate(shapiro_pvalue = shapiro.test(iSNVs)$p.value) %>%
  select(Genero, shapiro_pvalue) %>%
  unique() %>%
  mutate(normalidade = case_when(shapiro_pvalue > 0.05 ~ 'TRUE',
                                 shapiro_pvalue < 0.05 ~ 'FALSE'))

ggplot(intrahost_table_genero, aes(x = Genero, y = iSNVs)) +
  geom_point(aes(colour = lineage_up_major, size = iSNVs),
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
  scale_color_manual(values = cols_lineages) +
  stat_compare_means(method = 'wilcox.test', label = "p.signif", label.y = 30, label.x=1.5)


##boxplot de iSNVs por linhagem

intrahost_table_valores <- intrahost_table %>%
  group_by(lineage_up_major) %>%
  summarize(Mean = mean(iSNVs, na.rm=TRUE),
            Stdev = sd(iSNVs, na.rm=TRUE)) %>%
  mutate_at(2:3, round, 2)

alpha_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "Alpha (B.1.1.7)")

b11_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "B.1.1")

b1128_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "B.1.1.28")

b1133_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "B.1.1.33")

gamma_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "Gamma (P.1)")

gammasub_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "Gamma-sub (P.1.x)")

p2_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "P.2")

delta_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "Delta (B.1.617.2)")

deltasub_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "Delta-sub (AY.x)")

outras_values <- intrahost_table_valores %>%
  filter(lineage_up_major == "Others")

#checando normalidade por linhagem
intrahost_table_lineage_normalidade <- intrahost_table %>%
  group_by(lineage_up_major) %>%
  mutate(shapiro_pvalue = shapiro.test(iSNVs)$p.value) %>%
  select(lineage_up_major, shapiro_pvalue) %>%
  unique() %>%
  mutate(normalidade = case_when(shapiro_pvalue > 0.05 ~ 'TRUE',
                                 shapiro_pvalue < 0.05 ~ 'FALSE'))
shapiro.test(intrahost_table$iSNVs)

ggplot(intrahost_table, aes(x = lineage_up_major, y = iSNVs)) +
  geom_point(aes(colour = lineage_up_major, size = iSNVs),
             position = position_jitterdodge(dodge.width = 0.0, jitter.height = 0.5, jitter.width = 0.5),
             alpha = 0.7)+ 
  stat_summary(fun.data='mean_sdl', fun.args = list(mult=1), 
               geom="pointrange") +
  theme_bw() + 
  labs(x='Linhagem',
       y='iSNVs')  +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95)) +
  scale_color_manual(values = cols_lineages) +
  annotate("text", x='Alpha (B.1.1.7)', y=33, label= paste("Média: ",alpha_values$Mean,sep=' '))+
  annotate("text", x='Alpha (B.1.1.7)', y=31, label= paste("Desvio: ",alpha_values$Stdev,sep=' ')) +
  annotate("text", x='B.1.1', y=33, label= paste("Média: ",b11_values$Mean,sep=' '))+
  annotate("text", x='B.1.1', y=31, label= paste("Desvio: ",b11_values$Stdev,sep=' ')) +
  annotate("text", x='B.1.1.28', y=33, label= paste("Média: ",b1128_values$Mean,sep=' '))+
  annotate("text", x='B.1.1.28', y=31, label= paste("Desvio: ",b1128_values$Stdev,sep=' ')) +
  annotate("text", x='B.1.1.33', y=33, label= paste("Média: ",b1133_values$Mean,sep=' '))+
  annotate("text", x='B.1.1.33', y=31, label= paste("Desvio: ",b1133_values$Stdev,sep=' ')) +
  annotate("text", x='Delta (B.1.617.2)', y=33, label= paste("Média: ",delta_values$Mean,sep=' '))+
  annotate("text", x='Delta (B.1.617.2)', y=31, label= paste("Desvio: ",delta_values$Stdev,sep=' ')) +
  annotate("text", x='Delta-sub (AY.x)', y=33, label= paste("Média: ",deltasub_values$Mean,sep=' '))+
  annotate("text", x='Delta-sub (AY.x)', y=31, label= paste("Desvio: ",deltasub_values$Stdev,sep=' ')) +
  annotate("text", x='Gamma (P.1)', y=33, label= paste("Média: ",gamma_values$Mean,sep=' '))+
  annotate("text", x='Gamma (P.1)', y=31, label= paste("Desvio: ",gamma_values$Stdev,sep=' ')) +
  annotate("text", x='Gamma-sub (P.1.x)', y=33, label= paste("Média: ",gammasub_values$Mean,sep=' '))+
  annotate("text", x='Gamma-sub (P.1.x)', y=31, label= paste("Desvio: ",gammasub_values$Stdev,sep=' ')) +
  annotate("text", x='Others', y=33, label= paste("Média: ",outras_values$Mean,sep=' '))+
  annotate("text", x='Others', y=31, label= paste("Desvio: ",outras_values$Stdev,sep=' ')) +
  annotate("text", x='P.2', y=33, label= paste("Média: ",p2_values$Mean,sep=' '))+
  annotate("text", x='P.2', y=31, label= paste("Desvio: ",p2_values$Stdev,sep=' ')) +
  stat_compare_means(method = "kruskal.test", label.y = 45) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 40)
  
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

#checando normalidade por idade
intrahost_table_idade_normalidade <- intrahost_table_idade %>%
  group_by(grupo) %>%
  mutate(shapiro_pvalue = shapiro.test(iSNVs)$p.value) %>%
  select(grupo, shapiro_pvalue) %>%
  unique() %>%
  mutate(normalidade = case_when(shapiro_pvalue > 0.05 ~ 'TRUE',
                                 shapiro_pvalue < 0.05 ~ 'FALSE'))
shapiro.test(intrahost_table_idade$iSNVs)

ggplot(intrahost_table_idade, aes(x = grupo, y = iSNVs)) +
  geom_point(aes(colour = lineage_up_major, size = iSNVs),
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
  annotate("text", x='66+', y=41, label= paste("Desvio: ",values_66mais$Stdev,sep=' ')) +
  stat_compare_means(method = "kruskal.test", label.y = 45) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 30)

##correlacao iSNV profundidade

ggscatter(intrahost_table, x = "depth_iSNV", y = "iSNVs", 
          add = "reg.line", conf.int = TRUE, 
          size = "iSNVs", alpha = 0.5,
          xlab = "Profundidade Media", ylab = "Número de iSNVs") +
  stat_cor(method = "pearson", label.x = 3500, label.y = 45) +
  theme(legend.position = 'none')
  
##By period

x_order <- c("2020_01", "2020_02", "2020_03", "2020_04", "2020_05", "2020_06", 
             "2020_07", "2020_08", "2020_09", "2020_10", "2020_11", "2020_12", 
             "2021_01", "2021_02", "2021_03", "2021_04", "2021_05", "20210_6",
             "2021_07","2021_08", "2021_09", "2021_10")

lineages_period <- intrahost_table %>%
  mutate(month = format(as.Date(collection_date), "%m"), year = format(as.Date(collection_date), "%Y")) %>%
  unite(period,c(year,month))

data_lineages <- lineages_period %>%
  filter(period %in% x_order) %>%
  arrange(factor(period, levels = x_order))

#checando normalidade por periodo
intrahost_table_data_normalidade <- data_lineages %>%
  group_by(period) %>%
  mutate(shapiro_pvalue = shapiro.test(iSNVs)$p.value) %>%
  select(period, shapiro_pvalue) %>%
  unique() %>%
  mutate(normalidade = case_when(shapiro_pvalue > 0.05 ~ 'TRUE',
                                 shapiro_pvalue < 0.05 ~ 'FALSE'))
shapiro.test(data_lineages$iSNVs)

intrahost_data_valores <- data_lineages %>%
  group_by(period) %>%
  summarize(Mean = mean(iSNVs, na.rm=TRUE),
            Stdev = sd(iSNVs, na.rm=TRUE)) %>%
  mutate_at(2:3, round, 2)

values_2020_04 <- intrahost_data_valores %>%
  filter(period == "2020_04")

values_2020_05 <- intrahost_data_valores %>%
  filter(period == "2020_05")

values_2020_06 <- intrahost_data_valores %>%
  filter(period == "2020_06")

values_2020_07 <- intrahost_data_valores %>%
  filter(period == "2020_07")

values_2020_08 <- intrahost_data_valores %>%
  filter(period == "2020_08")

values_2020_09 <- intrahost_data_valores %>%
  filter(period == "2020_09")

values_2020_10 <- intrahost_data_valores %>%
  filter(period == "2020_10")

values_2020_11 <- intrahost_data_valores %>%
  filter(period == "2020_11")

values_2020_12 <- intrahost_data_valores %>%
  filter(period == "2020_12")

values_2021_01 <- intrahost_data_valores %>%
  filter(period == "2021_01")

values_2021_02 <- intrahost_data_valores %>%
  filter(period == "2021_02")

values_2021_03 <- intrahost_data_valores %>%
  filter(period == "2021_03")

values_2021_04 <- intrahost_data_valores %>%
  filter(period == "2021_04")

values_2021_05 <- intrahost_data_valores %>%
  filter(period == "2021_05")

values_2021_06 <- intrahost_data_valores %>%
  filter(period == "2021_05")

values_2021_07 <- intrahost_data_valores %>%
  filter(period == "2021_07")

values_2021_08 <- intrahost_data_valores %>%
  filter(period == "2021_08")

values_2021_09 <- intrahost_data_valores %>%
  filter(period == "2021_09")

values_2021_10 <- intrahost_data_valores %>%
  filter(period == "2021_10")

ggplot(data_lineages, aes(x=factor(period, levels = x_order),y=iSNVs)) +
  geom_point(aes(colour = lineage_up_major, size = iSNVs),
             position = position_jitterdodge(dodge.width = 0.0, jitter.height = 0.5, jitter.width = 0.5),
             alpha = 0.7)+ 
  stat_summary(fun.data='mean_sdl', fun.args = list(mult=1), 
               geom="pointrange") +
  theme_bw() + 
  labs(x='Período',
       y='iSNVs')  +
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = cols_lineages)  +
  stat_compare_means(method = "kruskal.test", label.y = 45) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 30) +
  annotate("text", x='2020_04', y=37, label= paste("Média: ",values_2020_04$Mean,sep=' '), size = 1.8)+
  annotate("text", x='2020_04', y=36, label= paste("Desvio: ",values_2020_04$Stdev,sep=' '), size = 1.8) +
  annotate("text", x='2020_05', y=37, label= paste("Média: ",values_2020_05$Mean,sep=' '), size = 1.8)+
  annotate("text", x='2020_05', y=36, label= paste("Desvio: ",values_2020_05$Stdev,sep=' '), size = 1.8) +
  annotate("text", x='2020_09', y=37, label= paste("Média: ",values_2020_09$Mean,sep=' '), size = 1.8)+
  annotate("text", x='2020_09', y=36, label= paste("Desvio: ",values_2020_09$Stdev,sep=' '), size = 1.8) +
  annotate("text", x='2020_10', y=37, label= paste("Média: ",values_2020_10$Mean,sep=' '), size = 1.8)+
  annotate("text", x='2020_10', y=36, label= paste("Desvio: ",values_2020_10$Stdev,sep=' '), size = 1.8) +
  annotate("text", x='2020_11', y=37, label= paste("Média: ",values_2020_11$Mean,sep=' '), size = 1.8)+
  annotate("text", x='2020_11', y=36, label= paste("Desvio: ",values_2020_11$Stdev,sep=' '), size = 1.8) +
  annotate("text", x='2021_07', y=37, label= paste("Média: ",values_2021_07$Mean,sep=' '), size = 1.8)+
  annotate("text", x='2021_07', y=36, label= paste("Desvio: ",values_2021_07$Stdev,sep=' '), size = 1.8) +
  annotate("text", x='2021_08', y=37, label= paste("Média: ",values_2021_08$Mean,sep=' '), size = 1.8)+
  annotate("text", x='2021_08', y=36, label= paste("Desvio: ",values_2021_08$Stdev,sep=' '), size = 1.8)
