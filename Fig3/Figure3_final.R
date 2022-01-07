library(ggplot2)
library(ggExtra)
library(gridExtra)
library(scales)
library(deSolve)
library(dbplyr)
library(dplyr) 
library(plotly)
library(viridis)
library(ggpubr)

LCMV_sub_file = 'LCMV_figure3_3_8_21.csv'
LASV_sub_file = 'LASV_figure3_3_8_21.csv'

LCMV_sub <- read.csv(LCMV_sub_file, header=TRUE)
LASV_sub <- read.csv(LASV_sub_file, header=TRUE)

NumPlot1 <- ggplot(data=LASV_sub, aes(x=time)) +
  geom_ribbon(aes(ymin=minP/500, ymax=maxP/500), fill = "grey53",) +
  geom_line(aes(y=meanP/500), color="orange", linetype = 2) +
  geom_vline(xintercept = LASV_sub$ReducMin[1], linetype = 3, color='grey60') +
  geom_vline(xintercept = LASV_sub$ReducMax[1], linetype = 3, color='grey60') +
  geom_vline(xintercept = LASV_sub$ReducMean[1], linetype = 3, color='grey60') +
  xlim(0, 1750) +
  ylim(0, 0.10) +
  labs(x="Time (days)", y="Pathogen prevalence",
       title="Lassa virus", color = "MCMV transmission") +
  theme_bw() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title.align=0.5) +
  guides(color=guide_legend(title.position = "top")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
NumPlot1

NumPlot2 <- ggplot(data=LCMV_sub, aes(x=time)) +
  geom_ribbon(aes(ymin=minP/500, ymax=maxP/500), fill = "grey53",) +
  geom_line(aes(y=meanP/500), color="orange", linetype = 2) +
  geom_vline(xintercept = LCMV_sub$ReducMin[1], linetype = 3, color='grey60') +
  geom_vline(xintercept = LCMV_sub$ReducMax[1], linetype = 3, color='grey60') +
  geom_vline(xintercept = LCMV_sub$ReducMean[1], linetype = 3, color='grey60') +
  xlim(0, 1750) +
  ylim(0, 0.10) +
  labs(x="Time (days)", y="Pathogen prevalence",
       title="Lymphocytic choriomeningitis virus") +
  theme_bw() +
  guides(color=guide_legend(title.position = "top")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
NumPlot2

FinalFig <- ggarrange(NumPlot1, NumPlot2, 
                      labels = c("a", "b"),
                      nrow = 1, ncol = 2,
                      common.legend = TRUE)
FinalFig

ggsave('',
       plot = FinalFig,
       height = 100,
       width = 200,
       units = "mm",
       dpi = 300,
       device="pdf")
