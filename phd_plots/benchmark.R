library(ggplot2)

setwd('/media/dezordi/34ECABD5066A98F6/01_Trabalho/projetos/PhD/plots/benchmark/')

time_table <- read.csv('tools.tsv', sep = '\t')

ggplot(time_table, aes(group,time,fill=group)) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(size=6, alpha=0.7) +
  scale_color_manual(values=c("#FF6666",
                              "#66FF66",
                              "#66B2FF")) +
  labs(x = "Versão Workflow",
       y = "Tempo de Análise (minutos)") +
  ylim(50,200) +
  scale_x_discrete(labels=c("v.0.0.2","v.0.0.2.5","v.0.0.3")) +
  theme(legend.position = "none")