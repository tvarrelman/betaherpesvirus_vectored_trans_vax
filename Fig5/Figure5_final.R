library(ggplot2)
library(ggExtra)
library(gridExtra)
library(scales)
library(deSolve)
library(dbplyr)
library(plotly)
library(parallel)
library(viridis)

setwd('')

CMV_mod <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    dS <- b - (d*betap*S*P)/b - (d*betav*S*V)/b - d*S
    dE <- (d*betav*S*V)/b - sig*E - d*E
    dV <- sig*E - d*V
    dP <-  (d*betap*S*P)/b - delt*P - d*P
    list(c(dS, dE, dV, dP))
  })
}

file <- '../data/GormanBetavEsts.csv'
beta.df <- read.csv(file, header=TRUE)

# d is chosen based on a lifespan of 365 days
d <- 0.002739726
# b is chosen to reflect a constant population size of 500 (found in lit)
b <- 500*d
# This estimate comes from our MCMV ABC parameter estimates 
sig <- 0.099

times <- seq(0, 30000, by = 1)
df_final = data.frame(matrix(ncol = 7, nrow = 0))

duration.list = c(75)
Rop.list <- c(1.5, 2)
df_final = data.frame()
for (j in 1:length(Rop.list)){
  Rop <- Rop.list[j]
  for (duration in duration.list){
    delt <- 1/duration
    betap <- Rop*(d+delt)
    for (i in 1:nrow(beta.df)){
      # Pathogen steady state (Vaccine absent)
      SSS <- (b*(d+delt))/(d*betap)
      #SSS <- 320
      SSP <- b*((1/(d+delt))-(1/betap))
      # Number of individuals initially vaccinated
      E0 <- SSS*.10
      state1 <- c(S = SSS, E = E0, V = 0, P = SSP)
      betav <- beta.df$betavEst[i]
      location <- beta.df$location[i]
      strain <- beta.df$strain[i]
      parameters <- c(b=b, betap=betap, betav=betav, d=d, sig=sig, delt=delt)
      out <- ode(y = state1, times = times, func = CMV_mod, parms = parameters)
      df <- as.data.frame(out)
      filtered.sub <- dplyr::filter(df, P <= SSP*0.05)
      timetoreduc <- filtered.sub$time[1]
      #print(timetoreduc)
      Vaccine.beta <- paste("\U03b2:",sprintf("%0.4f", round(betav, digits = 3)))
      Vaccine.strain <- paste("MCMV strain:",strain)
      reduc.df <- data.frame("betap"=betap, "betav"=betav, "time"=timetoreduc, 
                             "d"=d, "delt"=delt, "Rop"=Rop, "VaccStrain"=Vaccine.strain,
                             "Location"=location)
      df_final = rbind(df_final, reduc.df)
    }
  } 
}

df_final$Rop[df_final$Rop == 1.5] = "Pathogen~R[0]:~1.5"
df_final$Rop[df_final$Rop == 2] = "Pathogen~R[0]:~2"
CannotReduce <- c(NA, NA, NA, NA, NA, NA, NA, NA, 
                  NA, NA, NA, NA, NA, NA, NA, "R[0~p]>R[0~v]")
Coord <- c(NA, NA, NA, NA, NA, NA, NA, NA, 
           NA, NA, NA, NA, NA, NA, NA, 2000)
df_final$CannotReduce <- CannotReduce
df_final$Coord <- Coord

reducFig <- ggplot(data=df_final, aes(x=Location, y=time)) +
  #geom_point(aes(shape=VaccStrain), colour = "black", size = 4, stroke=2) +
  geom_point(aes(shape=VaccStrain, fill=Location),size=3, stroke = 1, position = position_dodge(0.5)) +
  scale_shape_manual(values=c(21, 24)) + 
  guides(fill = FALSE) +
  scale_color_manual(values = c("Boullanger Island"="#F8766D",
                                "Canberra"="#7CAE00",
                                "Macquarie Island"="#00BFC4",
                                "Walpeup"="#C77CFF")) +
  scale_x_discrete(limits=c("Boullanger Island",
                            "Walpeup",
                            "Canberra",
                            "Macquarie Island")) +
  geom_text(aes(x=Location, y=Coord, label=CannotReduce, color=Location), 
            angle=90, parse=TRUE, fontface="bold", color="grey30") +
  #scale_color_viridis(option="cividis", discrete=TRUE) +
  labs(x="Location", y="Time to 95% pathogen reduction (days)", shape="", fill="") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, hjust=1),
        legend.position = "top", legend.direction = "horizontal",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
new.lab <- as_labeller(c(a="A", b="B", c="italic(C)"), label_parsed)
Roplabs <- c(expression(Pathogen~R[0]: 1.5), expression(Pathogen~R[0]: 2))
Fig2 <- reducFig + facet_wrap(~Rop, labeller=label_parsed)
Fig2

ggsave('', 
       plot = Fig2,
       height = 5.5,
       width = 8.7,
       units = "cm",
       dpi = 1000,
       device="pdf",
       scale = 2)
