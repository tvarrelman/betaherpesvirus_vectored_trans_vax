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

# original model
CMV_mod <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    dS <- b - (d*betap*S*P)/b - (d*betav*S*V)/b - d*S
    dE <- (d*betav*S*V)/b - sig*E - d*E
    dV <- sig*E - d*V
    dP <-  (d*betap*S*P)/b - delt*P - d*P
    dR <- delt*P - d*R
    list(c(dS, dE, dV, dP, dR))
  })
}

# partial transmission blocking
CMV_mod2 <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dS <- b - (d*betap*S*P)/b - (1-rho)*((d*betap*S*Vp)/b) - (1-rho)*((d*betap*S*Ep)/b) - (d*betav*S*V)/b - (d*betav*S*Vp)/b - (d*betav*S*Vr)/b - d*S
    dE <- (d*betav*S*V)/b + (d*betav*S*Vp)/b + (d*betav*S*Vr)/b - (d*betap*E*P)/b -(1-rho)*((d*betap*E*Vp)/b)-(1-rho)*((d*betap*E*Ep)/b) - sig*E - d*E
    dEp <- (d*betap*E*P)/b + (1-rho)*((d*betap*E*Vp)/b) +(1-rho)*((d*betap*E*Ep)/b) - sig*Ep - delt*Ep - d*Ep
    dEr <- delt*Ep - sig*Er - d*Er
    dV <- sig*E - (d*betap*V*P)/b - (1-rho)*((d*betap*V*Vp)/b)-(1-rho)*((d*betap*V*Ep)/b) - d*V
    dVp <- (d*betap*V*P)/b + (1-rho)*((d*betap*V*Vp)/b)+ (1-rho)*((d*betap*V*Ep)/b)+ sig*Ep - delt*Vp - d*Vp
    dVr <- delt*Vp + sig*Er - d*Vr
    dP <- (d*betap*S*P)/b + (1-rho)*((d*betap*S*Vp)/b) + (1-rho)*((d*betap*S*Ep)/b) - delt*P - d*P
    dR <- delt*P - d*R
    list(c(dS, dE, dEp, dEr, dV, dVp, dVr, dP, dR))
  })
}
# partial infection blocking
CMV_mod3 <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dS <- b - (d*betap*S*(P+Ep+Vp))/b - (d*betav*S*(V+Vp+Vr))/b - d*S
    dE <- (d*betav*S*(V+Vp+Vr))/b - ((1-rho)*d*betap*E*(P+Ep+Vp))/b - (sig+d)*E
    dEp <- ((1-rho)*d*betap*E*(P+Ep+Vp))/b - (sig+delt+d)*Ep
    dEr <- delt*Ep - (sig+d)*Er
    dV <- sig*E - ((1-rho)*d*betap*V*(P+Ep+Vp))/b - d*V
    dVp <- sig*Ep + ((1-rho)*d*betap*V*(P+Ep+Vp))/b - (d+delt)*Vp
    dVr <- delt*Vp + sig*Er - d*Vr
    dP <- (d*betap*S*(P+Ep+Vp))/b - (delt+d)*P
    dR <- delt*P - d*R
    list(c(dS, dE, dEp, dEr, dV, dVp, dVr, dP, dR))
  })
}

#delayed immunity model
CMV_mod4 <- function(t, state, parameters){
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

# to show the difference between main Fig 4 and SI Fig 3, we need to overlay CMV_mod2 and CMV_mod3

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

rho.list <- c(0.5, 0.75, 1)
df_final = data.frame()
df_thresh = data.frame()
for (i in 1:length(betap.list)){
  betap <- betap.list[i]
  delt <- delt.list[i]
  label <- path.label[i]
  SSS <- (b*(d+delt))/(d*betap)
  SSP <- b*((1/(d+delt))-(1/betap))
  SSR <- (b*delt*((1/(d+delt))-(1/betap)))/(d)
  E0 <- SSS*0.01
  #SSS2 <- SSS-E0
  df_total = data.frame()
  rho_thresh = (sig*(d+delt-betap)*betav)/(betap*(d*(d+sig)-sig*betav))
  print(rho_thresh)
  for (j in 1:length(rho.list)){
    rho <- rho.list[j]
    state1 <- c(S = SSS, E=E0, Ep=0, Er=0, V=0, Vp=0, Vr=0, P=SSP, R=SSR)
    parameters <- c(b=b, betap=betap, betav=betav, d=d, sig=sig, delt=delt, rho=rho)
    # transmission blocking
    out <- ode(y = state1, times = times, func = CMV_mod2, parms = parameters, method="ode45")
    df <- as.data.frame(out)
    df$Path.label <- label
    df$betap <- betap
    df$delt <- delt
    df$rho <- as.numeric(rho)
    df$model <- rep('Transmission blocking', nrow(df))
    # In the model of partial efficacy (transmission blocking), only 1-rho individuals transmit the pathogen
    # We need to account for this when evaluating pathogen reduction. Therefore, adjust Vp and Ep to only include those that transmit pathogen
    ep_ad <- list()
    vp_ad <- list()
    for (k in 1:nrow(df)){
      vp_ad <- append(vp_ad, df$Vp[k]*(1-as.numeric(rho)))
      ep_ad <- append(ep_ad, df$Ep[k]*(1-as.numeric(rho)))
    }
    #print(vp_ad)
    df$Vp <- unlist(vp_ad)
    df$Ep <- unlist(ep_ad)
    # infection blocking
    out2 <- ode(y = state1, times = times, func = CMV_mod3, parms = parameters, method="ode45")
    df2 <- as.data.frame(out2)
    df2$Path.label <- label
    df2$betap <- betap
    df2$delt <- delt
    df2$rho <- as.numeric(rho)
    df2$model <- rep('Infection blocking', nrow(df2))

    df_total <- rbind(df_total, df, df2)
  }
  df_final <- rbind(df_final, df_total)
}

lasv_sub <- subset(df_final, Path.label %in% c('Lassa virus'))
lcmv_sub <- subset(df_final, Path.label %in% c('Lymphocytic choriomeningitis virus'))

lasv_ssp <- b*((1/(d+0.0457))-(1/0.072659585))
lcmv_ssp <- b*((1/(d+0.002739726))-(1/0.006052883))

lasv <- lasv_sub[(lasv_sub$rho==0.5) & (lasv_sub$model=='Infection blocking'),]
lasv2 <- lasv_sub[(lasv_sub$rho==0.5) & (lasv_sub$model=='Transmission blocking'),]

lcmv <- lcmv_sub[(lcmv_sub$rho==0.5) & (lcmv_sub$model=='Infection blocking'),]
lcmv2 <- lcmv_sub[(lcmv_sub$rho==0.5) & (lcmv_sub$model=='Transmission blocking'),]

lasv_sub_infection <- dplyr::filter(lasv, (P+Vp+Ep) <= lasv_ssp*0.05)$time[1]
lasv_sub_trans <- dplyr::filter(lasv2, (P+Vp+Ep) <= lasv_ssp*0.05)$time[1]

lcmv_sub_infection <- dplyr::filter(lcmv, (P+Vp+Ep) <= lcmv_ssp*0.05)$time[1]
lcmv_sub_trans <- dplyr::filter(lcmv2, (P+Vp+Ep) <= lcmv_ssp*0.05)$time[1]

lasvFig <- ggplot(data=lasv_sub) +
  geom_line(aes(x=time, y=(P+Vp+Ep)/500, color=as.factor(rho), linetype=model)) +
  xlab("Time (days)") +
  ylab("Pathogen prevalence") +
  labs(color = "Vaccine efficacy", linetype='MCMV model') +
  theme_bw() +
  theme(legend.position = "left", legend.direction = "horizontal",
        legend.title.align=0.5) +
  guides(color=guide_legend(title.position = "top", nrow = 1),
         linetype=guide_legend(title.position = "top", nrow = 1)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlim(0, 2000) +
  ylim(0, .10) +
  ggtitle("Lassa virus")
lasvFig

lcmvFig <- ggplot(lcmv_sub) +
  geom_line(data=lcmv_sub, aes(x=time, y=(P+Vp+Ep)/500, color=as.factor(rho), linetype=model), show.legend=FALSE) +
  xlab("Time (days)") +
  ylab("Pathogen prevalence") +
  theme_bw() +
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
       dpi = 1200,
       device="pdf")




