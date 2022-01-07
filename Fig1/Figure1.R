library(ggplot2)
library(ggExtra)

file.path <- '../data/'
file <- 'TruePosterior_1_29_21.csv'
my.file <- paste(file.path, file, sep="")
data <- read.csv(my.file, header=TRUE)
names(data) <- c('Beta','Sigma_1','Sigma_2')
# Convert transmission to freq. dependent
data$betaFreq <- data$Beta*22

p <- ggplot(data, aes(x=betaFreq, y=Sigma_1)) +
  geom_hex(bins = 80, aes(fill=..count..)) +
  geom_point(alpha = 0, show.legend=FALSE) +
  scale_fill_continuous(type = "viridis") +
  xlab(expression(paste("Transmission rate ", (beta[v])))) +  
  ylab(expression(paste("Rate of becoming infectious ", (sigma)))) +
  labs(fill='Count')+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.9, 0.70),
        legend.background = element_rect(fill = "gainsboro",size=0.5, linetype="solid", 
                                         colour ="black"),
        panel.background = element_rect(size=1, 
                                        linetype="solid", color="black"))

p2 <- ggMarginal(p,data = data, x=Beta, y=Sigma_1, type =  "histogram", color='black',margins = "both", size = 5)

ggsave('',
       plot = p2,
       height = 100,
       width = 135,
       units = "mm",
       dpi = 1200,
       device="pdf")

