library(reshape2)
library(plotly)
library(ggplot2)
library(viridis)
library(ggrepel)

final.df <- read.csv('/Users/tannervarrelman/Documents/MCMV_project_5_20/DurationPathData/Figure2Sim_2_10_21.csv', header=TRUE)
path.file <- '/Users/tannervarrelman/Documents/MCMV_project_5_20/PaperFigures/Figure1_pathogen_ests.csv'
path.df <- read.csv(path.file, header=TRUE)

# Find the reduction time for specific parameters
sub1 <- subset(final.df, Rop %in% c(1.50))
sub2 <- subset(sub1, (1/delt) %in% c(10))

# Find the reduction time for specific parameters
sub3 <- subset(final.df, Rop %in% unique(final.df$Rop)[150])
sub4 <- subset(sub3, (1/delt) %in% unique(1/final.df$delt)[356])

# Find the reduction time for specific parameters
LASVsub1 <- subset(final.df, Rop %in% c(1.50))
LASVsub2 <- subset(LASVsub1, (1/delt) %in% c(1/unique(LASVsub1$delt)[13]))

# Find the reduction time for specific parameters
LCMVsub1 <- subset(final.df, Rop %in% c(1.10))
LCMVsub2 <- subset(LCMVsub1, (1/delt) %in% unique(1/final.df$delt)[356])

reducFig2 <- ggplot(data=subset(final.df, betav %in% c(0.033))) +
  geom_raster(aes(x=1/delt, y=Rop, fill=time)) +
  geom_contour(aes(x=1/delt, y=Rop, z=time), color="grey", linetype="dashed") +
  geom_point(data=path.df, aes(x=Duration, R0), color='black') +
  geom_text_repel(data=path.df, aes(x=Duration , y=R0, label=Virus), 
                  show.legend=FALSE, color='black', nudge_y=0.15, angle=90,
                  min.segment.length = unit(0, 'lines')) +
  labs(x="Average infectious period (days)", y=expression(Pathogen~R[0]), fill="95% Pathogen reduction (days)") +
  theme(legend.position="top", legend.key.width = unit(1.5, "cm"), legend.title.align=0.5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(direction=-1, breaks=c(0,250,500,750, 1000, 1250),
                     labels=c(0, 250, 500, 750, 1000, 1250)) +
  guides(fill=guide_colorbar(title.position="top"))
reducFig2

ggsave('',
       plot = reducFig2,
       height = 130,
       width = 120,
       units = "mm",
       dpi = 300,
       device="pdf")

