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
CMV_mod <- function(t, state, parameters){
  with(as.list(c(state,  parameters)),{
    dS <- b - (d*betap*S*P)/b - (d*betav*S*V)/b - d*S
    dE <- (d*betav*S*V)/b - sig*E - d*E
    dV <- sig*E - d*V
    dP <-  (d*betap*S*P)/b - delt*P - d*P
    list(c(dS, dE, dV, dP))
  })
}

# List of pathogen duration values
deltap.list <- seq(10, 370, by=1)
Rop.list <- seq(1.01, 2.5, by= 0.01)

# All possible combinations of the two lists
param.grid <- expand.grid(deltap.list, Rop.list)


path_reduc_fun <- function(iter){
  # this takes a while
  times <- seq(0, 10000, by = 2)
  sig <- 0.099
  # MCMV transmission rate
  betav <- 0.033
  #betav <- param.grid$Var3[iter]
  # Path. transmission
  delt <- 1/param.grid$Var1[iter]
  #delt <- 1/deltap.list[iter]
  duration <- param.grid$Var1[iter]
  Rop <- param.grid$Var2[iter]
  #print(Rop)
  #print(delt)
  d <- 0.002739726
  betap <- Rop*(d+delt)
  b <- d*500
  #SSP <- 130
  # Susceptible steady state (Vaccine absent)
  SSS <- (b*(d+delt))/(d*betap)
  #SSS <- 320
  SSP <- b*((1/(d+delt))-(1/betap))
  # Number of individuals initially vaccinated
  # Initial fraction vaccinated
  E0 <- SSS*0.10
  state1 <- c(S = SSS, E = E0, V = 0, P = SSP)
  # This estimate comes from our MCMV ABC parameter estimates 
  parameters <- c(b=b, betap=betap, betav=betav, d=d, sig=sig, delt=delt)
  out <- ode(y = state1, times = times, func = CMV_mod, parms = parameters, method="lsoda")
  df <- as.data.frame(out)
  filtered.sub <- dplyr::filter(df, P <= SSP*0.05)
  timetoreduc <- filtered.sub$time[1]
  reduc.df <- data.frame("betap"=betap, "betav"=betav, "time"=timetoreduc, "d"=d, "delt"=delt, "Rop"=Rop)
}

# Set up the number of simluations that we want to run
# This is set to the length of the parameter grid that we will be indexing
iter <- seq(1,nrow(param.grid), by=1)
system.time({results <- pbmclapply(iter, path_reduc_fun, mc.cores=3)})
final.df <- do.call("rbind", results)

write.csv(final.df, 'Figure2Sim.csv')

