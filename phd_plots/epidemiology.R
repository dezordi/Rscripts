### vacinas_fora do ar = https://opendatasus.saude.gov.br/dataset/covid-19-vacinacao/resource/ef3bd0b8-b605-474b-9ae5-c97390c197a8

library(dplyr)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(geobr)
library(ggplot2)
library(scatterpie)

##### AMBIENTE E VARIAVEIS #####

setwd('/media/dezordi/34ECABD5066A98F6/01_Trabalho/projetos/PhD/plots/epidemiology/')
target_lineages <- list("Gamma (P.1)","Gamma-sub (P.1.x)",
                        "Delta (B.1.617.2)","Delta-sub (AY.x)",
                        "P.2","B.1.1.28","B.1.1.33",
                        "Alpha (B.1.1.7)","B.1.1")

x_order <- c("2020_19", "2020_20", "2020_21", "2020_22", "2020_23", "2020_24", 
             "2020_25", "2020_26", "2020_27", "2020_28", "2020_29", "2020_30", 
             "2020_31", "2020_32", "2020_33", "2020_34", "2020_35", "2020_36",
             "2020_37","2020_38", "2020_39", "2020_40", "2020_41", "2020_42",
             "2020_43", "2020_44", "2020_45", "2020_46", "2020_47", "2020_48",
             "2020_49", "2020_50", "2020_51", "2020_52", "2020_53", "2021_1", 
             "2021_2", "2021_3", "2021_4", "2021_5", "2021_6", "2021_7","2021_8",
             "2021_9", "2021_10", "2021_11", "2021_12", "2021_13", "2021_14",
             "2021_15", "2021_16", "2021_17", "2021_18","2021_19", "2021_20",
             "2021_21", "2021_22", "2021_23", "2021_24", "2021_25", "2021_26",
             "2021_27", "2021_28", "2021_29", "2021_30", "2021_31", "2021_32",
             "2021_33", "2021_34", "2021_35", "2021_36", "2021_37","2021_38",
             "2021_39", "2021_40", "2021_41", "2021_42", "2021_43", "2021_44",
             "2021_45", "2021_46", "2021_47", "2021_48")

cols_vac <- c('1st shot'='red', '2nd shot or Unique'='blue')
cols_lineages <- c("B.1.1"='#FF66B2',"Others"='#C0C0C0',
                   "B.1.1.28"='#66FFFF',"B.1.1.33"='#66B2FF',
                   "Alpha (B.1.1.7)"='#6666FF',"P.2"='#FFFF66',
                   "Gamma (P.1)"="#FF6666","Gamma-sub (P.1.x)"='#FFB266',
                   "Delta (B.1.617.2)"='#B2FF66',"Delta-sub (AY.x)"='#66FF66',
                   "B.1.1.28"='#66FFFF',"B.1.1.33"='#66B2FF')

###### PERNAMBUCO ######

PE_cases_death_table <- read.csv('basegeral_PE_2021.csv', sep = ';') ### https://dados.seplag.pe.gov.br/apps/corona.html update: 21 dez. 2021
PE_vaccine_table <- read.csv('basegeral_vacina_PE_2021.csv', sep=';') ### https://www.pecontracoronavirus.pe.gov.br/ update: 21 dez. 2021
PE_lineage_table <- read.csv('gisaid_PE_2021_2.tsv', sep='\t') ### gisaid metadada download | grep "Pernambuco" update: 21 dez. 2021

### MAPA ###

lineages_to_map <- PE_lineage_table %>%
  filter(as.Date(collection_date) < '2021-12-01',
         as.Date(submission_date) < '2021-12-01', 
         lineage != 'None',
         location != 'Fernando de noronha',
         location != 'Fernando De Noronha',
         location != 'Fernando de Noronha') %>%
  mutate(lineage = str_replace(lineage,"P.1..*","Gamma-sub"),
         lineage = str_replace(lineage,"P.1","Gamma (P.1)"),
         lineage = str_replace(lineage,"Gamma-sub","Gamma-sub (P.1.x)"),
         lineage = str_replace(lineage,"AY.*","Delta-sub (AY.x)"),
         lineage = str_replace(lineage,"B.1.617.2","Delta (B.1.617.2)"),
         lineage = str_replace(lineage,"B.1.1.7","Alpha (B.1.1.7)"),
         lineage = case_when(lineage %in% unlist(target_lineages) ~ lineage),
         lineage = replace_na(lineage,"Others"))

  lineages_to_map_sum <- lineages_to_map %>%
  group_by(lat, long) %>%
  summarise(B.1.1 = sum(lineage == "B.1.1"),
            Others = sum(lineage == "Others"),
            B.1.1.28 = sum(lineage == "B.1.1.28"),
            B.1.1.33 = sum(lineage == "Gamma-sub (P.1.x)"),
            Alpha = sum(lineage == "Alpha (B.1.1.7)"),
            P.2 = sum(lineage == "P.2"),
            Gamma = sum(lineage == "Gamma (P.1)"),
            Gamma.sub = sum(lineage == "Gamma-sub (P.1.x)"),
            Delta = sum(lineage == "Delta (B.1.617.2)"),
            Delta.sub = sum(lineage == "Delta-sub (AY.x)")) 

lineages_to_map_sum <- lineages_to_map_sum %>%
  group_by(lat,long) %>%
  mutate(count =sum(B.1.1, Others, B.1.1.28, B.1.1.33, Alpha, P.2, Gamma, Gamma.sub, Delta, Delta.sub)) %>%
  mutate(count_norm = case_when(count <= 5 ~ 0.8,
                                count <= 30 ~ 1.2,
                                count <= 100 ~ 1.6,
                                count > 100 ~ 2))

PE_map <- read_municipality(code_muni= "PE", year=2020)
PE_map <- PE_map %>%
  filter(name_muni !="Fernando De Noronha")

plot_map_PE <- ggplot() +
  geom_sf(data=PE_map, size=.2, show.legend = FALSE, alpha=0.7) +
  geom_scatterpie(data=lineages_to_map_sum, size=.3,
                  aes(x=long,y=lat, r= 0.05 * count_norm, group=lat),
                  cols = colnames(lineages_to_map_sum[,c(3:12)]),
                  alpha=1) +
  theme(legend.position = 'top') +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Lineage") + 
  guides(fill=guide_legend(ncol=10))
plot_map_PE
### FREQUENCIAS ###

lineages_period <- PE_lineage_table%>%
  filter(as.Date(collection_date) < '2021-12-01',
         as.Date(submission_date) < '2021-12-01',
         lineage != 'None') %>%
  mutate(week = epiweek(as.Date(collection_date)), year = format(as.Date(collection_date), "%Y")) %>%
  mutate(lineage = str_replace(lineage,"P.1..*","Gamma-sub"),
         lineage = str_replace(lineage,"P.1","Gamma (P.1)"),
         lineage = str_replace(lineage,"Gamma-sub","Gamma-sub (P.1.x)"),
         lineage = str_replace(lineage,"AY.*","Delta-sub (AY.x)"),
         lineage = str_replace(lineage,"B.1.617.2","Delta (B.1.617.2)"),
         lineage = str_replace(lineage,"B.1.1.7","Alpha (B.1.1.7)"),
         lineage = case_when(lineage %in% unlist(target_lineages) ~ lineage),
         lineage = replace_na(lineage,"Others"))%>%
  unite(period,c(year,week)) %>%
  group_by(period) %>%
  mutate(genome_count = n()) %>%
  group_by(period,lineage) %>%
  mutate(lineage_freq = n()/genome_count, lineage_count = n())
lineages_period <- unique(select(lineages_period, period, lineage, lineage_freq, lineage_count, genome_count))

data_lineages <- lineages_period %>%
  filter(period %in% x_order) %>%
  arrange(factor(period, levels = x_order))

plot_lineage_PE <- ggplot(data_lineages, aes(x=factor(period, levels = x_order),y=lineage_freq, fill=lineage)) +
  geom_bar(position = 'fill', stat = 'identity', colour='black',alpha=0.7) +
  geom_line(aes(y=genome_count/200, group=1), colour = 'black', size = 2, alpha=0.6)+
  scale_y_continuous(
    name = 'Frequência',
    sec.axis = sec_axis(~.*200, name ='Número de Genomas (n)')
  ) +
  scale_x_discrete(limits = x_order) +
  labs(x = 'Período (Ano/Semana Epidemiológica)') +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 90)) +
  guides(fill=guide_legend(ncol=12)) +
  scale_fill_manual(values = cols_lineages)

plot_lineage_PE
### GENOMAS/CASOS ###

cases_table_period <- PE_cases_death_table %>%
  filter(is.na(dt_notificacao) != TRUE & as.Date(dt_notificacao) < '2021-12-01') %>%
  mutate(week = epiweek(as.Date(dt_notificacao)), year = format(as.Date(dt_notificacao), "%Y")) %>%
  unite(period,c(year,week)) %>%
  count(period) %>%
  rename(case_number = n)

genome_by_period <- lineages_period %>%
  group_by(period) %>%
  summarise(genome_count)
genome_by_period <- unique(select(genome_by_period,period,genome_count))

cases_genome_by_period <- reduce(list(cases_table_period, genome_by_period), full_join, by='period')
cases_genome_by_period[is.na(cases_genome_by_period)] <- 0
cases_genome_by_period <- cases_genome_by_period %>%
  mutate(genomesBycases = (genome_count / case_number)*100)

plot_cases_PE <- ggplot(cases_genome_by_period, aes(x=factor(period, levels = x_order),y=case_number)) +
  geom_bar(stat = 'identity', colour='black',alpha=0.3) +
  geom_line(aes(y=genomesBycases*25000, group=1), colour = 'red', size = 2, alpha=0.6)+
  scale_y_continuous(
    name = 'Número de Casos (n)',
    sec.axis = sec_axis(~./25000, name ='Genomas Sequenciados / Número de Casos (%)')
  ) +
  scale_x_discrete(limits = x_order) +
  labs(x = '') +
  theme(axis.text.x = element_blank())

plot_map_PE
plot_cases_PE
plot_lineage_PE


###### RIO GRANDE DO SUL ######

RS_cases_death_table <- read.csv('basegeral_RS_2021.csv', sep = ';') ### https://dados.seplag.pe.gov.br/apps/corona.html update: 21 dez. 2021
RS_cases_death_table <- RS_cases_death_table %>%
  mutate(DATA_CONFIRMACAO = dmy(DATA_CONFIRMACAO))
RS_vaccine_table <- read.csv('basegeral_vacina_PE_2021.csv', sep=';') ### https://www.pecontracoronavirus.pe.gov.br/ update: 21 dez. 2021
RS_lineage_table <- read.csv('gisaid_RS_2021_2.tsv', sep='\t') ### gisaid metadada download | grep "Pernambuco" update: 21 dez. 2021

##Lineages
IsDate <- function(mydate, date.format = "%Y-%m-%d") {
  tryCatch(!is.na(as.Date(mydate, date.format)),  
           error = function(err) {FALSE})  
}

RS_lineage_table$complete_collection_DATE <- IsDate(RS_lineage_table$collection_date)
RS_lineage_table$complete_submission_DATE <- IsDate(RS_lineage_table$submission_date)

### MAPA ###

lineages_to_map <- RS_lineage_table %>%
  filter(complete_collection_DATE == TRUE,
         complete_submission_DATE == TRUE) %>%
  filter(as.Date(collection_date) < '2021-12-01',
         as.Date(submission_date) < '2021-12-01', 
         lineage != 'None') %>%
  mutate(lineage = str_replace(lineage,"P.1..*","Gamma-sub"),
         lineage = str_replace(lineage,"P.1","Gamma (P.1)"),
         lineage = str_replace(lineage,"Gamma-sub","Gamma-sub (P.1.x)"),
         lineage = str_replace(lineage,"AY.*","Delta-sub (AY.x)"),
         lineage = str_replace(lineage,"B.1.617.2","Delta (B.1.617.2)"),
         lineage = str_replace(lineage,"B.1.1.7","Alpha (B.1.1.7)"),
         lineage = case_when(lineage %in% unlist(target_lineages) ~ lineage),
         lineage = replace_na(lineage,"Others"))

lineages_to_map_sum <- lineages_to_map %>%
  group_by(lat, long) %>%
  summarise(B.1.1 = sum(lineage == "B.1.1"),
            Others = sum(lineage == "Others"),
            B.1.1.28 = sum(lineage == "B.1.1.28"),
            B.1.1.33 = sum(lineage == "Gamma-sub (P.1.x)"),
            Alpha = sum(lineage == "Alpha (B.1.1.7)"),
            P.2 = sum(lineage == "P.2"),
            Gamma = sum(lineage == "Gamma (P.1)"),
            Gamma.sub = sum(lineage == "Gamma-sub (P.1.x)"),
            Delta = sum(lineage == "Delta (B.1.617.2)"),
            Delta.sub = sum(lineage == "Delta-sub (AY.x)")) 

lineages_to_map_sum <- lineages_to_map_sum %>%
  group_by(lat,long) %>%
  mutate(count =sum(B.1.1, Others, B.1.1.28, B.1.1.33, Alpha, P.2, Gamma, Gamma.sub, Delta, Delta.sub)) %>%
  mutate(count_norm = case_when(count <= 5 ~ 0.8,
                                count <= 30 ~ 1.2,
                                count <= 100 ~ 1.6,
                                count > 100 ~ 2))

RS_map <- read_municipality(code_muni= "RS", year=2020)

plot_map_RS <- ggplot() +
  geom_sf(data=RS_map, size=.2, fill='white', show.legend = FALSE, alpha=0.7) +
  geom_scatterpie(data=lineages_to_map_sum, size=.3,
                  aes(x=long,y=lat, r= 0.09 * count_norm, group=lat),
                  cols = colnames(lineages_to_map_sum[,c(3:12)]),
                  alpha=1) +
  theme(legend.position = 'top') +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Lineage") + 
  guides(fill=guide_legend(ncol=10))

plot_map_RS

### FREQUENCIAS ###

lineages_period <-RS_lineage_table%>%
  filter(complete_collection_DATE == TRUE,
         complete_submission_DATE == TRUE) %>%
  filter(as.Date(collection_date) < '2021-12-01',
         as.Date(submission_date) < '2021-12-01',
         lineage != 'None') %>%
  mutate(week = epiweek(as.Date(collection_date)), year = format(as.Date(collection_date), "%Y")) %>%
  mutate(lineage = str_replace(lineage,"P.1..*","Gamma-sub"),
         lineage = str_replace(lineage,"P.1","Gamma (P.1)"),
         lineage = str_replace(lineage,"Gamma-sub","Gamma-sub (P.1.x)"),
         lineage = str_replace(lineage,"AY.*","Delta-sub (AY.x)"),
         lineage = str_replace(lineage,"B.1.617.2","Delta (B.1.617.2)"),
         lineage = str_replace(lineage,"B.1.1.7","Alpha (B.1.1.7)"),
         lineage = case_when(lineage %in% unlist(target_lineages) ~ lineage),
         lineage = replace_na(lineage,"Others"))%>%
  unite(period,c(year,week)) %>%
  group_by(period) %>%
  mutate(genome_count = n()) %>%
  group_by(period,lineage) %>%
  mutate(lineage_freq = n()/genome_count, lineage_count = n())
lineages_period <- unique(select(lineages_period, period, lineage, lineage_freq, lineage_count, genome_count))

data_lineages <- lineages_period %>%
  filter(period %in% x_order) %>%
  arrange(factor(period, levels = x_order))

plot_lineage_RS <- ggplot(data_lineages, aes(x=factor(period, levels = x_order),y=lineage_freq, fill=lineage)) +
  geom_bar(position = 'fill', stat = 'identity', colour='black',alpha=0.7) +
  geom_line(aes(y=genome_count/120, group=1), colour = 'black', size = 2, alpha=0.6)+
  scale_y_continuous(
    name = 'Frequência',
    sec.axis = sec_axis(~.*120, name ='Número de Genomas (n)')
  ) +
  scale_x_discrete(limits = x_order) +
  labs(x = 'Período (Ano/Semana Epidemiológica)') +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 90)) +
  guides(fill=guide_legend(ncol=12)) +
  scale_fill_manual(values = cols_lineages)

plot_lineage_RS

### GENOMAS/CASOS ###

cases_table_period <- RS_cases_death_table %>%
  filter(is.na(DATA_CONFIRMACAO) != TRUE & as.Date(DATA_CONFIRMACAO) < '2021-12-01') %>%
  mutate(week = epiweek(as.Date(DATA_CONFIRMACAO)), year = format(as.Date(DATA_CONFIRMACAO), "%Y")) %>%
  unite(period,c(year,week)) %>%
  count(period) %>%
  rename(case_number = n)

genome_by_period <- lineages_period %>%
  group_by(period) %>%
  summarise(genome_count)
genome_by_period <- unique(select(genome_by_period,period,genome_count))

cases_genome_by_period <- reduce(list(cases_table_period, genome_by_period), full_join, by='period')
cases_genome_by_period[is.na(cases_genome_by_period)] <- 0
cases_genome_by_period <- cases_genome_by_period %>%
  mutate(genomesBycases = (genome_count / case_number)*100)

x_order <- c("2020_5","2020_6","2020_7","2020_8","2020_9","2020_10","2020_11","2020_12",
             "2020_13","2020_14","2020_15","2020_16","2020_17",
             "2020_18","2020_19", "2020_20", "2020_21", "2020_22", "2020_23", "2020_24", 
             "2020_25", "2020_26", "2020_27", "2020_28", "2020_29", "2020_30", 
             "2020_31", "2020_32", "2020_33", "2020_34", "2020_35", "2020_36",
             "2020_37","2020_38", "2020_39", "2020_40", "2020_41", "2020_42",
             "2020_43", "2020_44", "2020_45", "2020_46", "2020_47", "2020_48",
             "2020_49", "2020_50", "2020_51", "2020_52", "2020_53", "2021_1", 
             "2021_2", "2021_3", "2021_4", "2021_5", "2021_6", "2021_7","2021_8",
             "2021_9", "2021_10", "2021_11", "2021_12", "2021_13", "2021_14",
             "2021_15", "2021_16", "2021_17", "2021_18","2021_19", "2021_20",
             "2021_21", "2021_22", "2021_23", "2021_24", "2021_25", "2021_26",
             "2021_27", "2021_28", "2021_29", "2021_30", "2021_31", "2021_32",
             "2021_33", "2021_34", "2021_35", "2021_36", "2021_37","2021_38",
             "2021_39", "2021_40", "2021_41", "2021_42", "2021_43", "2021_44",
             "2021_45", "2021_46", "2021_47", "2021_48")

plot_cases_RS <- ggplot(cases_genome_by_period, aes(x=factor(period, levels = x_order),y=case_number)) +
  geom_bar(stat = 'identity', colour='black',alpha=0.3) +
  geom_line(aes(y=genomesBycases*5000, group=1), colour = 'red', size = 2, alpha=0.6)+
  scale_y_continuous(
    name = 'Número de Casos (n)',
    sec.axis = sec_axis(~./5000, name ='Genomas Sequenciados / Número de Casos (%)')
  ) +
  scale_x_discrete(limits = x_order) +
  labs(x = '') +
  theme(axis.text.x = element_blank())

plot_map_RS
plot_cases_RS
plot_lineage_RS

#####################################OLD PLOTS##############################

##CASES & DEPTH
death_table_period <- PE_cases_death_table %>%
  filter(evolucao == 'OBITO' & is.na(dt_notificacao) != TRUE & as.Date(dt_notificacao) < '2021-12-01') %>%
  mutate(month = format(as.Date(dt_notificacao), "%m"), year = format(as.Date(dt_notificacao), "%Y")) %>%
  unite(period,c(year,month)) %>%
  count(period,sort = FALSE)%>%
  rename(death_number = n) %>%
  rename(date = period)



rm(PE_cases_death_table)

##VACCINE
vacine_1st_table_period_PE <- PE_vaccine_table %>%
  filter(is.na(vacina_dataaplicacao) != TRUE & as.Date(vacina_dataaplicacao) < '2021-12-01') %>% 
  filter(grepl("1",vacina_descricao_dose_nova)) %>%
  mutate(month = format(as.Date(vacina_dataaplicacao), "%m"), year = format(as.Date(vacina_dataaplicacao), "%Y")) %>%
  unite(period,c(year,month)) %>% 
  count(period) %>%
  rename(one_shot_number = n) %>%
  rename(date = period)

vacine_2nd_table_period_PE <- PE_vaccine_table %>%
  filter(is.na(vacina_dataaplicacao) != TRUE & as.Date(vacina_dataaplicacao) < '2021-12-01') %>% 
  filter(grepl("2",vacina_descricao_dose_nova) | grepl("unica",vacina_descricao_dose_nova)) %>%
  mutate(month = format(as.Date(vacina_dataaplicacao), "%m"), year = format(as.Date(vacina_dataaplicacao), "%Y")) %>%
  unite(period,c(year,month)) %>% 
  count(period) %>%
  rename(two_shot_number = n) %>%
  rename(date = period)

rm(PE_vaccine_table)

##Lineages
##Inserir genomas IAM

##Merging
data_list <- list(cases_table_period,death_table_period,vacine_1st_table_period_PE,vacine_2nd_table_period_PE)
data_merged <- reduce(data_list, full_join, by = 'date')
data_merged <- data_merged %>%
  mutate(one_shot_number = replace_na(one_shot_number,0), two_shot_number = replace_na(two_shot_number,0))

data_merged <- data_merged %>%
  filter(date %in% x_order) %>%
  arrange(factor(date, levels = x_order))

data_lineages <- lineages_period %>%
  filter(period %in% x_order) %>%
  arrange(factor(period, levels = x_order))

#pop PE elegivel para vacinacao 2021 7.692.429
first_plot <- ggplot(data_merged, aes(x=factor(date, level = x_order))) +
  geom_bar(aes(y=case_number), stat="identity",fill='#FF6666',alpha=0.7,colour='black') +
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits = x_order) +
  geom_line(aes(y=cumsum(((one_shot_number/7692429)*100)*2100), colour='1st shot',group=1), size = 2, alpha=0.5) +
  geom_line(aes(y=cumsum(((two_shot_number/7692429)*100)*2100), colour='2nd shot or Unique',group=1), size = 2, alpha=0.5) +
  scale_y_continuous(
    name = "Number of Cases (n)",
    sec.axis = sec_axis(~./210000, name = "Percentage of vaccinated population (%)")
  ) +
  labs(x='Period',
       colour='Vaccine doses')+
  scale_colour_manual(values = cols) +
  theme(legend.position = 'top')
first_plot

second_plot <- ggplot(data_lineages, aes(x=factor(period, levels = x_order),y=lineage_freq, fill=lineage)) +
  geom_bar(position = 'fill', stat = 'identity', colour='black',alpha=0.7) +
  geom_line(aes(y=genome_count/300, group=1), colour = 'black', size = 2, alpha=0.6)+
  scale_y_continuous(
    name = 'Lineage Frequency',
    sec.axis = sec_axis(~.*300, name ='Number of Genomes (n)')
  ) +
  theme(axis.text.x = element_blank())+
  scale_x_discrete(limits = x_order) +
  labs(x = '') +
  theme(legend.position = 'top') +
  guides(fill=guide_legend(ncol=12)) +
  scale_fill_manual(values = cols_lineages)

ggarrange(first_plot,second_plot, heights = c(5,6), nrow = 2, align = 'v', labels = c('A','B'))

rm(first_plot,second_plot,third_plot,cases_table_period,data_lineages,data_merged,death_table_period,
   uti_occ_period,vacine_1st_table_period_PE,vacine_2nd_table_period_PE, lineages_period)

### RIO GRANDE DO SUL ###

RS_cases_death_table <- read.csv('basegeral_RS_2021.csv', sep = ';') ### https://ti.saude.rs.gov.br/covid19/api update: 21 dez. 2021
RS_cases_death_table <- RS_cases_death_table %>%
  mutate(DATA_CONFIRMACAO = dmy(DATA_CONFIRMACAO))
RS_vaccine_table <- read.csv('basegeral_vacina_RS_2021.csv', sep=';') ### https://vacina.saude.rs.gov.br/ update: 21 dez. 2021
RS_vaccine_table <- RS_vaccine_table %>%
  mutate(Data = dmy(Data))
RS_lineage_table <- read.csv('gisaid_RS_2021.tsv', sep='\t') ### gisaid metadada download | grep "Rio Grande do Sul" update: 21 dez. 2021
#####################################MONTH PLOTS##############################

##CASES & DEPTH
death_table_period <- RS_cases_death_table %>%
  filter(EVOLUCAO == 'OBITO' & is.na(DATA_CONFIRMACAO) != TRUE & as.Date(DATA_CONFIRMACAO) < '2021-12-01') %>%
  mutate(month = format(as.Date(DATA_CONFIRMACAO), "%m"), year = format(as.Date(DATA_CONFIRMACAO), "%Y")) %>%
  unite(period,c(year,month)) %>%
  count(period,sort = FALSE)%>%
  rename(death_number = n) %>%
  rename(date = period)

cases_table_period <- RS_cases_death_table %>%
  filter(is.na(DATA_CONFIRMACAO) != TRUE & as.Date(DATA_CONFIRMACAO) < '2021-12-01') %>%
  mutate(month = format(as.Date(DATA_CONFIRMACAO), "%m"), year = format(as.Date(DATA_CONFIRMACAO), "%Y")) %>%
  unite(period,c(year,month)) %>%
  count(period) %>%
  rename(case_number = n) %>%
  rename(date = period)

rm(RS_cases_death_table)

##VACCINE

vacine_1st_table_period_RS <- RS_vaccine_table %>%
  filter(is.na(Data) != TRUE & as.Date(Data) < '2021-12-01') %>%
  mutate(month = format(as.Date(Data), "%m"), year = format(as.Date(Data), "%Y")) %>%
  unite(period,c(year,month)) %>%
  select(period,X1.dose) %>%
  group_by(period) %>%
  summarise(n = sum(X1.dose)) %>%
  rename(one_shot_number = n) %>%
  rename(date = period)

vacine_2nd_table_period_RS <- RS_vaccine_table %>%
  filter(is.na(Data) != TRUE & as.Date(Data) < '2021-12-01') %>%
  mutate(month = format(as.Date(Data), "%m"), year = format(as.Date(Data), "%Y")) %>%
  unite(period,c(year,month)) %>%
  select(period,X2.dose,Dose.unica) %>%
  group_by(period) %>%
  summarise(n = sum(X2.dose,Dose.unica)) %>%
  rename(two_shot_number = n) %>%
  rename(date = period)

rm(RS_vaccine_table)



##Inserir genomas IAM
lineages_period <- RS_lineage_table %>%
  filter(complete_DATE == TRUE) %>%
  filter(as.Date(collection_date) < '2021-12-01', lineage != 'None') %>%
  mutate(month = format(as.Date(collection_date), "%m"), year = format(as.Date(collection_date), "%Y")) %>%
  mutate(lineage = str_replace(lineage,"P.1..*","Gamma-sub"),
         lineage = str_replace(lineage,"P.1","Gamma (P.1)"),
         lineage = str_replace(lineage,"Gamma-sub","Gamma-sub (P.1.x)"),
         lineage = str_replace(lineage,"AY.*","Delta-sub (AY.x)"),
         lineage = str_replace(lineage,"B.1.617.2","Delta (B.1.617.2)"),
         lineage = str_replace(lineage,"B.1.1.7","Alpha (B.1.1.7)"),
         lineage = case_when(lineage %in% unlist(target_lineages) ~ lineage),
         lineage = replace_na(lineage,"Others"))%>%
  unite(period,c(year,month)) %>%
  group_by(period) %>%
  mutate(genome_count = n()) %>%
  group_by(period,lineage) %>%
  mutate(lineage_freq = n()/genome_count, lineage_count = n())
lineages_period <- unique(select(lineages_period, period, lineage, lineage_freq, lineage_count, genome_count))
rm(RS_lineage_table)  

##Merging
data_list <- list(cases_table_period,death_table_period,vacine_1st_table_period_RS,vacine_2nd_table_period_RS)
data_merged <- reduce(data_list, full_join, by = 'date')
data_merged <- data_merged %>%
  mutate(one_shot_number = replace_na(one_shot_number,0), two_shot_number = replace_na(two_shot_number,0))

data_merged <- data_merged %>%
  filter(date %in% x_order) %>%
  arrange(factor(date, levels = x_order))

data_lineages <- lineages_period %>%
  filter(period %in% x_order) %>%
  arrange(factor(period, levels = x_order))

#pop RS elegivel vacinacao 2021 9.750.642
first_plot <- ggplot(data_merged, aes(x=factor(date, level = x_order))) +
  geom_bar(aes(y=case_number), stat="identity",fill='#66B2FF',alpha=0.7,colour='black') +
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits = x_order) +
  geom_line(aes(y=cumsum(((one_shot_number/9750642)*100)*2300), colour='1st shot',group=1), size = 2, alpha=0.5) +
  geom_line(aes(y=cumsum(((two_shot_number/9750642)*100)*2300), colour='2nd shot or Unique',group=1), size = 2, alpha=0.5) +
  scale_y_continuous(
    name = "Number of Cases (n)",
    sec.axis = sec_axis(~./230000, name = "Percentage of vaccinated population (%)")
  ) +
  labs(x='Period',
       colour='Vaccine doses')+
  scale_colour_manual(values = cols) +
  theme(legend.position = 'top')
first_plot

second_plot <- ggplot(data_lineages, aes(x=factor(period, levels = x_order),y=lineage_freq, fill=lineage)) +
  geom_bar(position = 'fill', stat = 'identity', colour='black',alpha=0.7) +
  geom_line(aes(y=genome_count/300, group=1), colour = 'black', size = 2, alpha=0.6)+
  scale_y_continuous(
    name = 'Lineage Frequency',
    sec.axis = sec_axis(~.*300, name ='Number of Genomes (n)')
  ) +
  theme(axis.text.x = element_blank())+
  scale_x_discrete(limits = x_order) +
  labs(x = '') +
  theme(legend.position = 'top') +
  guides(fill=guide_legend(ncol=12)) +
  scale_fill_manual(values = cols_lineages)
second_plot
ggarrange(first_plot,second_plot, heights = c(4,6), nrow = 2, align = 'v', labels = c('A','B'))