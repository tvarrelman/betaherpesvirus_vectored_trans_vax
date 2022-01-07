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

CMV_mod <- function(t, state, parameters){
  with(as.list(c(state,  parameters)),{
    dS <- b - (d*betap*S*P)/b - (d*betav*S*V)/b - d*S
    dE <- (d*betav*S*V)/b - sig*E - d*E
    dV <- sig*E - d*V
    dP <-  (d*betap*S*P)/b - delt*P - d*P
    list(c(dS, dE, dV, dP))
  })
}

delayed_inf <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dS <- b - (d*betap*S*P)/b  - ((d*betap*S*Ep)/b) - (d*betav*S*V)/b - (d*betav*S*Vp)/b - (d*betav*S*Vr)/b - d*S
    dE <- (d*betav*S*V)/b + (d*betav*S*Vp)/b + (d*betav*S*Vr)/b - (d*betap*E*P)/b - ((d*betap*E*Ep)/b) - sig*E - d*E
    dEp <- (d*betap*E*P)/b + ((d*betap*E*Ep)/b) - sig*Ep - delt*Ep - d*Ep
    dEr <- delt*Ep - sig*Er - d*Er
    dV <- sig*E - d*V
    dVp <- sig*Ep - delt*Vp - d*Vp
    dVr <- delt*Vp + sig*Er - d*Vr
    dP <- (d*betap*S*P)/b + ((d*betap*S*Ep)/b) - delt*P - d*P
    dR <- delt*P - d*R
    list(c(dS, dE, dEp, dEr, dV, dVp, dVr, dP, dR))
  })
}

times <- seq(0, 10000, by = 1)
# d is chosen based on a lifespan of 365 days
d <- 0.002739726
# b is chosen to reflect a constant population size of 500 (found in lit)
b <- 500*d
betap.list <- c(0.072659585, 0.006052883)
delt.list <- c(0.0457, 0.002739726)
path.label <- c("Lassa virus", "Lymphocytic choriomeningitis virus")
betav <- 0.033
sig <- 0.099

SSP <- b*((1/(d+delt))-(1/betap))

df1_total = data.frame()
df2_total = data.frame()
for (i in 1:length(betap.list)){
  betap <- betap.list[i]
  delt <- delt.list[i]
  label <- path.label[i]
  SSS <- (b*(d+delt))/(d*betap)
  SSP <- b*((1/(d+delt))-(1/betap))
  SSR <- (b*delt*((1/(d+delt))-(1/betap)))/(d)
  E0 <- SSS*0.1
  SSS2 <- SSS-E0
  df_total = data.frame()
  state1 <- c(S = SSS2, E=E0, Ep=0, Er=0, V=0, Vp=0, Vr=0, P=SSP, R=SSR)
  state2 = c(S = SSS2, E = E0, V = 0, P= SSP)
  parameters <- c(b=b, betap=betap, betav=betav, d=d, sig=sig, delt=delt)
  out <- ode(y = state1, times = times, func = delayed_inf, parms = parameters, method="ode45")
  out2 <- ode(y = state2, times = times, func = CMV_mod, parms = parameters, method="ode45")
  df <- as.data.frame(out)
  df2 <- as.data.frame(out2)
  df$Path.label <- label
  df$betap <- betap
  df$delt <- delt
  
  mod2 = rep("Orignal Model", nrow(df2))
  mod1 = rep("Delayed Immunity Model", nrow(df))
  df$model <- mod1
  df2$model <- mod2
  
  df2$Path.label <- label
  df2$betap <- betap
  df2$delt <- delt
  
  df1_total <- rbind(df1_total,df)
  df2_total <- rbind(df2_total, df2)
}

lasv_sub <- subset(df1_total, Path.label %in% c('Lassa virus'))
lcmv_sub <- subset(df1_total, Path.label %in% c('Lymphocytic choriomeningitis virus'))

lasv_sub2 <- subset(df2_total, Path.label %in% c('Lassa virus'))
lcmv_sub2 <- subset(df2_total, Path.label %in% c('Lymphocytic choriomeningitis virus'))

lasv_ssp <- b*((1/(d+0.0457))-(1/0.072659585))
lcmv_ssp <- b*((1/(d+0.002739726))-(1/0.006052883))

filtered.sub1 <- dplyr::filter(lasv_sub, (P+Vp+Ep) <= lcmv_ssp*0.05)
filtered.sub2 <- dplyr::filter(lasv_sub2, P <= lcmv_ssp*0.05)
timetoreduc <- filtered.sub1$time[1]
timetoreduc2 <- filtered.sub2$time[1]


lasvFig <- ggplot() +
  geom_line(data=lasv_sub, aes(x=time, y=(P+Vp+Ep)/500, color=model)) +
  geom_line(data=lasv_sub2, aes(x=time, y=(P)/500, color=model)) +
  #geom_line(data=lasv_sub2, aes(x=time, y=(P+Vp+Ep)/500, linetype=rho), color='black', show.legend=FALSE) +
  xlab("Time (days)") +
  ylab("Pathogen prevalence") +
  labs(color = "") +
  theme_bw() +
  theme(legend.position = "left", legend.direction = "horizontal",
        legend.title.align=0.5) +
  guides(color=guide_legend(title.position = "top", nrow = 1)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlim(0, 2000) +
  ylim(0, .10) +
  ggtitle("Lassa virus")
lasvFig

lcmvFig <- ggplot() +
  geom_line(data = lcmv_sub, aes(x=time, y=(P+Vp+Ep)/500, color=model), show.legend=TRUE) +
  geom_line(data = lcmv_sub2, aes(x=time, y=(P)/500, color=model), show.legend=TRUE) +
  #geom_line(data=lcmv_sub2, aes(x=time, y=(P+Vp+Ep)/500, linetype=rho), color='black', show.legend=FALSE) +
  xlab("Time (days)") +
  ylab("Pathogen prevalence") +
  labs(color="") +
  theme_bw() +
  #theme(legend.position = "top", legend.direction = "horizontal",
  #      legend.title.align=0.5) +
  guides(color=guide_legend(title.position = "left")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlim(0, 2000) +
  ylim(0, .10) +
  ggtitle("Lymphocytic choriomeningitis virus")
lcmvFig

FinalEffFig <- ggarrange(lasvFig, lcmvFig, labels = c("a", "b"),
                         nrow = 1, ncol = 2,
                         common.legend = TRUE)
FinalEffFig

ggsave('',
       plot = FinalEffFig,
       height = 100,
       width = 200,
       units = "mm",
       dpi = 300,
       device="pdf")


