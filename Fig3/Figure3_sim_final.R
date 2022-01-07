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

setwd('')
posterior_file = '../data/TruePosterior_1_29_21.csv'
posterior_df = read.csv(posterior_file, header=TRUE)
names(posterior_df) <- c('beta', 'sigma1', 'simga2')
set.seed(1223)
post_pars <- sample_n(posterior_df, 100, replace = TRUE)
# d is chosen based on a lifespan of 365 days
d <- 0.002739726
# b is chosen to reflect a constant population size of 500 (found in lit)
b <- 500*d
rho <- 0.7

deltap.list <- c(0.0457, 0.002739726)
# Pathogen transmission estimates from the literature
betap.list <- c( 0.072659585, 0.006052883)
betap.label <- c('Lassa virus', 'Lymphocytic choriomeningitis virus')

file <- 'GormanBetavEsts_1_5_21.csv'
beta.df <- read.csv(file, header=TRUE)

# MCMV frequency depedent model
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

times <- seq(0, 3000, by = 1)
df_final = data.frame(matrix(ncol = 7, nrow = 0))
df_final2 = data.frame(matrix(ncol = 7, nrow = 0))
for (j in 1:length(betap.list)){
  betap <- betap.list[j]
  lab <- betap.label[j]
  delt <- deltap.list[j]
  print(betap/(d+delt))
  df_total = data.frame(matrix(ncol = 7, nrow = 0))
  df_total2 = data.frame(matrix(ncol = 7, nrow = 0))
  for (i in 1:nrow(post_pars)){
    betav <- post_pars$beta[i]*22
    sig <- post_pars$sigma1[i]
    SSS <- (b*(d+delt))/(d*betap)
    SSP <- b*((1/(d+delt))-(1/betap))
    SSR <- (b*delt*((1/(d+delt))-(1/betap)))/(d)
    # Number of individuals initially vaccinated
    E0 <- SSS*0.1
    state1 <- c(S = SSS, E = E0, V = 0, P = SSP, R=SSR)
    state2 <- c(S = SSS, E=E0, Ep=0, Er=0, V=0, Vp=0, Vr=0, P=SSP, R=SSR)
    parameters <- c(b=b, betap=betap, betav=betav, d=d, sig=sig, delt=delt)
    parameters2 <- c(b=b, betap=betap, betav=betav, d=d, sig=sig, delt=delt, rho=rho)
    out <- ode(y = state1, times = times, func = CMV_mod, parms = parameters)
    out2 <- ode(y = state2, times = times, func = CMV_mod2, parms = parameters2)
    df <- as.data.frame(out)
    df2 <- as.data.frame(out2)
    Vaccine.beta <- rep(paste("\U03b2:",sprintf("%0.4f", round(betav, digits = 3))), times=nrow(df))
    df$pathReduc <- rep(SSP*0.05, nrow(df))
    df2$pathReduc <- rep(SSP*0.05, nrow(df))
    df$Path.label <- lab 
    df2$Path.label <- lab 
    df$betav <- betav
    df2$betav <- betav
    df$simga1 <- sig
    df2$simga1 <- sig
    df$SimNum <- i
    df2$SimNum <- i
    df1 <- cbind(df, Vaccine.beta)
    df3 <- cbind(df2, Vaccine.beta)
    df_total <- rbind(df_total,df1)
    df_total2 <- rbind(df_total2,df3) 
  } 
  Pathogen.beta <- rep(paste("\U03b2:",betap), times=nrow(df_total))
  df_total$Pathogen.beta <- Pathogen.beta
  df_total2$Pathogen.beta <- Pathogen.beta
  df_final <- rbind(df_final, df_total)
  df_final2 <- rbind(df_final2, df_total2)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

print(getmode(df_final$betav))

write.csv(df_final, 'Figure3Sim_3_8_21.csv')
# The second dataframe isn't used in the paper, but looks at varrying efficacy w/ sampling from the posterior
#write.csv(df_final2, 'Figure3Sim_reduc_3_8_21.csv')

