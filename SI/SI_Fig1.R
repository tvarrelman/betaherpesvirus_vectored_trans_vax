library(ggplot2)
library(ggExtra)
library(gridExtra)
library(scales)
library(deSolve)
library(dbplyr)
library(plotly)
library(parallel) 
library(viridis)
library(ggrepel)
library(pbmcapply)
library(gtools)
library(ggpubr)

setwd('')
#Define the model
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

times <- seq(0, 100000, by = 1)
max_local <- max(times) - 1
# d is chosen based on a lifespan of 365 days
d <- 0.002739726
# b is chosen to reflect a constant population size of 500 (found in lit)
b <- 500*d

delt <- 0.002739726
# Beta_p's chosen based on R0
betap.list <- c(0.0109589, 0.01643836, 0.02191781, 0.02739726)

rho_list <- seq(0.5, 1, 0.05)
betav <- 0.033
sig <- 0.099
param.grid <- expand.grid(betap.list,rho_list)

path_reduc_fun <- function(iter){
  #betap <- betap.list[iter]
  betap <- param.grid$Var1[iter]
  rho <- param.grid$Var2[iter]
  Rov <- (betav*sig)/(d*(d+sig))
  Rop <- betap/(d+delt)
  SSS <- (b*(d+delt))/(d*betap)
  SSP <- b*((1/(d+delt))-(1/betap))
  # Number of individuals initially vaccinated
  E0 <- SSS*0.1
  SSR <- (b*delt*((1/(d+delt))-(1/betap)))/(d)
  SSS2 <- SSS-E0
  #crit_rho = (Rov/Rop)*((1-Rop)/(1-Rov)) 
  crit_rho <- (sig*(d+delt-betap)*betav)/(betap*(d*(d+sig)-sig*betav))
  state1 <- c(S = SSS2, E=E0, Ep=0, Er=0, V=0, Vp=0, Vr=0, P=SSP, R=SSR)
  parameters1 <- c(b=b, betap=betap, betav=betav, d=d, sig=sig, delt=delt, rho=rho)
  out <- ode(y = state1, times = times, func = CMV_mod2, parms = parameters1)
  df <- as.data.frame(out)
  p <- df$P[max_local]
  ep <- df$Ep[max_local]*(1-rho)
  vp <- df$Vp[max_local]*(1-rho)
  prev <- (p+ep+vp)/500
  Rop_label <- Rop
  reduc.df <- data.frame("betap"=betap, "betav"=betav, 
                         "rho"=rho, "d"=d, "delt"=delt, 
                         "Rop"=Rop_label, "Rov"=Rov, "crit_rho"=crit_rho,
                         "Prevalence"=prev)
  rbind(reduc.df)
}

iter <- seq(1, nrow(param.grid), by=1)
system.time({results <- pbmclapply(iter, path_reduc_fun, mc.cores=2)})
final.df <- do.call("rbind", results)

erad_plot <- ggplot(data = final.df) +
  geom_point(aes(x = rho, y = Prevalence), show.legend = TRUE) +
  geom_vline(aes(xintercept=crit_rho), linetype="dashed", color="red") +
  xlab('Vaccine efficacy') +
  ylab('Pathogen prevalence') +
  labs(color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
erad_plot

new.labs <- c("0.0109589"="R[0~p]==2", "0.01643836"="R[0~p]==3", 
              "0.02191781"="R[0~p]==4", "0.02739726"="R[0~p]==5")
final.fig <- erad_plot + facet_wrap(~betap, ncol = 2, 
                       labeller = labeller(betap = as_labeller(new.labs, label_parsed)))

final.fig
ggsave('',
       plot = final.fig,
       height = 4.5,
       width = 4.5,
       units = "in",
       dpi = 300,
       device="pdf")