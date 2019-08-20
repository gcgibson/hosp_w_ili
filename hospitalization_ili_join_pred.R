library(rjags)
library(R2jags)
# read ILI + hosp data
hosp_data <- read.csv("/Users/gcgibson/final_delay/data-raw//hosp_data.csv")
nat_data <- read.csv("/Users/gcgibson/final_delay/data-raw/flu_data.csv")
nat_data <- nat_data[nat_data$region == "National" & nat_data$year > 2009  | (nat_data$region == "National" & nat_data$year == 2009 &nat_data$week >= 35),]

# get common data
common_dates <- intersect(hosp_data$time,nat_data$time)
nat_data <- nat_data[nat_data$time %in% common_dates,]
hosp_data <- hosp_data[hosp_data$time %in% common_dates,]


# subset to training data

nat_data_train <- nat_data$weighted_ili[1:200]

# remap to [0,1]
nat_data_train <- nat_data_train/(max(nat_data_train)+.00001)

hosp_data_train <- hosp_data[hosp_data$age == 2,]$weeklyrate[1:200]
# map to poisson
hosp_data_train <- hosp_data_train*10


jags.data = list("ILI"=nat_data_train,"N"=length(nat_data_train),
                 "HOSP" =hosp_data_train)
## joint estimation of hosp and wILI

model.loc=("ss_model.txt")
jagsscript = cat("
                 model {  
                 # priors on parameters
                 mu ~ dnorm(0, 0.0001); 
                 tau.pro ~ dgamma(0.001,0.001); 

                 v <-1
                 phi_m1 ~ dnorm(0,1);
                 
                 TRUE_INFLUENZA_M1[1] <- mu;
                 beta_mean[1] <- exp(TRUE_INFLUENZA_M1[1])/(1+exp(TRUE_INFLUENZA_M1[1]))
                 ILI[1] ~ dbeta( v*beta_mean[1],(1-beta_mean[1])*v );
                 HOSP[1] ~ dpois(exp(TRUE_INFLUENZA_M1[1]));
                 
                for(i in 2:N) {
                   TRUE_INFLUENZA_M1[i] ~ dnorm(phi_m1*TRUE_INFLUENZA_M1[i-1],tau.pro); 
                    beta_mean[i] <- exp(TRUE_INFLUENZA_M1[i])/(1+exp(TRUE_INFLUENZA_M1[i]))

                    ILI[i] ~ dbeta( v*beta_mean[i],(1-beta_mean[i])*v);
                    HOSP[i] ~ dpois(exp(TRUE_INFLUENZA_M1[i]) +.00001);
                   }
                 }  
                 ",file=model.loc)


jags.params=c("TRUE_INFLUENZA_M1","phi_m1")
mod_ss = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, n.chains = 3, 
n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE) 
phi_m1 <- mean(mod_ss$BUGSoutput$sims.array[,,"phi_m1"])
#plot(colMeans(TRUE_INFLUENZA))

#plot(exp(colMeans(TRUE_INFLUENZA))/(1+exp(colMeans(TRUE_INFLUENZA))))

#plot(exp(colMeans(TRUE_INFLUENZA)))



### INDEPENDENT ESTIMATION

#### wILI


model.loc=("ss_model.txt")
jagsscript = cat("
                 model {  
                 # priors on parameters
                 mu ~ dnorm(0, 0.0001); 
                 tau.pro ~ dgamma(0.001,0.001); 
                 
                 v <-1
                 phi_m2 ~ dnorm(0,1);
                 TRUE_INFLUENZA_M2[1] <- mu;
                 beta_mean[1] <- exp(TRUE_INFLUENZA_M2[1])/(1+exp(TRUE_INFLUENZA_M2[1]))
                 ILI[1] ~ dbeta( v*beta_mean[1],(1-beta_mean[1])*v );

                 for(i in 2:N) {
                 TRUE_INFLUENZA_M2[i] ~ dnorm(phi_m2*TRUE_INFLUENZA_M2[i-1],tau.pro); 
                 beta_mean[i] <- exp(TRUE_INFLUENZA_M2[i])/(1+exp(TRUE_INFLUENZA_M2[i]))
                 
                 ILI[i] ~ dbeta( v*beta_mean[i],(1-beta_mean[i])*v);
                 }
                 }  
                 ",file=model.loc)

jags.params=c("TRUE_INFLUENZA_M2","phi_m2")
mod_ss_independent_wILI = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, n.chains = 3, 
              n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE) 
phi_m2 <- mean(mod_ss_independent_wILI$BUGSoutput$sims.array[,,"phi_m2"])


### HOSP Independent


model.loc=("ss_model.txt")
jagsscript = cat("
                 model {  
                 # priors on parameters
                 mu ~ dnorm(0, 0.0001); 
                 tau.pro ~ dgamma(0.001,0.001); 
                 
                 v <-1
                 phi_m3 ~ dnorm(0,1);
                 TRUE_INFLUENZA_M3[1] <- mu
                 HOSP[1] ~ dpois(exp(TRUE_INFLUENZA_M3[1]));
                 
                for(i in 2:N) {
                   TRUE_INFLUENZA_M3[i] ~ dnorm(phi_m3*TRUE_INFLUENZA_M3[i-1],tau.pro); 

                    HOSP[i] ~ dpois(exp(TRUE_INFLUENZA_M3[i]) +.00001);
                   }
                 }  
                 ",file=model.loc)

jags.params=c("TRUE_INFLUENZA_M3","phi_m3")
mod_ss_independent_hosp = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, n.chains = 3, 
              n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE) 
phi_m3 <- mean(mod_ss_independent_hosp$BUGSoutput$sims.array[,,"phi_m3"])


### AR(1) Evaluation

nat_data_test <- nat_data$weighted_ili[201:276]
nat_data_test <- nat_data_test/(max(nat_data_test)+.00001)
hosp_data_test <- hosp_data[hosp_data$age == 2,]$weeklyrate[201:276]
hosp_data_test <- hosp_data_test*10

mse_m1_ili <- c()
mse_m1_hosp <- c()
mse_m2 <- c()
mse_m3 <- c()

for (i in 1:(length(nat_data_test)-1)){
  ili_hat_under_joint <- phi_m1*nat_data_test[i]
  ili_hat_under_indp <- phi_m2*nat_data_test[i]
  
  hosp_hat_under_joint <- phi_m1*hosp_data_test[i]
  hosp_hat_under_ind   <- phi_m3*hosp_data_test[i]
  
  truth_ili <- nat_data_test[i+1]
  truth_hosp <- hosp_data_test[i+1]
  
  mse_m1_ili <- c(mse_m1_ili,(ili_hat_under_joint-truth_ili)^2)
  mse_m1_hosp <- c(mse_m1_hosp,(hosp_hat_under_joint-truth_hosp)^2)
  mse_m2 <- c(mse_m2,(ili_hat_under_indp-truth_ili)^2)
  mse_m3 <- c(mse_m3,(hosp_hat_under_ind-truth_hosp)^2)
}

mean(mse_m1_ili)
mean(mse_m2)

mean(mse_m1_hosp)
mean(mse_m3)

