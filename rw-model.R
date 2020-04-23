setwd("/")

library(R.matlab)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)

# Import Data -------------------------------------------------------------

Ss  <- list(
  S1  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_3_2.mat")$rec.task,
  S2  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_4_1.mat")$rec.task,
  S3  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_5_1.mat")$rec.task,
  S4  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_6_2.mat")$rec.task,
  S5  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_7_2.mat")$rec.task,
  S6  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_8_1.mat")$rec.task,
  S7  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_9_2.mat")$rec.task,
  S8  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_10_2.mat")$rec.task,
  S9  <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_11_1.mat")$rec.task,
  S10 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_12_1.mat")$rec.task,
  S11 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_13_2.mat")$rec.task,
  S12 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_14_1.mat")$rec.task,
  S13 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_15_1.mat")$rec.task,
  S14 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_17_2.mat")$rec.task,
  S15 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_18_1.mat")$rec.task,
  S16 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_19_1.mat")$rec.task,
  S17 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_20_2.mat")$rec.task,
  S18 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_21_2.mat")$rec.task,
  S19 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_22_1.mat")$rec.task,
  S20 <- readMat("Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_23_1.mat")$rec.task
)




f_alpha <- function(mu, sigma){
  alpha <- (((1-mu)/(sigma^2)) - (1/mu))*(mu^2)
  return(alpha)
}

f_beta <- function(alpha, mu){
  beta <- alpha*((1/mu) - 1)
  return(beta)
}




# Define Model (Rescorla Wagner)-------------------------------------------------------------------

rw <- function(params){
  for(s in 1:length(freqs)){
    for (t in 1:freqs[s]){
      if(t == 1){
        L[t, s] <- .5
        V[t, s] <- .5 + params[1]*(f[t, s] - .5)
      } else{
        A       <- 1 / (1 + exp((2*V[t-1, s] - 1) * -params[2]))
        L[t, s] <- ifelse(r[t, s] == 1, A, 1 - A)
        V[t, s] <- V[t-1, s] + params[1]*(f[t, s] - V[t-1, s])
      }
    }
  }
  log.likelihoods <- log(L)
  deviance <- -2 * sum(log.likelihoods, na.rm = TRUE)
  return(deviance)
}


# Fit Model (Rescorla Wagner)----------------------------------------------

results <- matrix(NA, nrow=length(Ss), ncol=9,
                  dimnames = list(c(), c("Subject",
                                         "Threat_Alpha",
                                         "Threat_Beta",
                                         "Threat_LL",
                                         "Threat_Conv",
                                         "Safety_Alpha",
                                         "Safety_Beta",
                                         "Safety_LL",
                                         "Safety_Conv")))

for(s in 1:length(Ss)){
  S <- as.data.frame(Ss[[s]])
  S <- S[which(S$V1 != 3 & !grepl("NaN", S$V7)), c(1, 4, 6, 7)]
  colnames(S) <- c("group", "stim", "feedback", "choice")
  S$choice <- ifelse(S$choice == 1, 1, 0)
  stimOld <- as.integer(names(table(S$stim)))
  stimNew <- 1:length(table(S$stim))
  for(i in 1:length(stimOld)){
    S$stim[which(S$stim == stimOld[i])] <- stimNew[i]
  }
  for(g in 1:max(S$group)){
    data  <- S[which(S$group == g),]
    stim  <- data$stim
    freqs <- as.vector(table(data$stim))
    V     <- matrix(NA, nrow = max(freqs), ncol = length(freqs))
    L     <- matrix(NA, nrow = max(freqs), ncol = length(freqs))
    f     <- matrix(NA, nrow = max(freqs), ncol = length(freqs))
    r     <- matrix(NA, nrow = max(freqs), ncol = length(freqs))
    for(j in 1:length(freqs)){
      fData <- data$feedback[which(stim==j)]
      rData <- data$choice[which(stim==j)]
      if(length(fData) > 0){
        f[1:length(fData), j] <- fData
        r[1:length(rData), j] <- rData
      }
    }
    
    rw.fits <- optim(par = c(.1, 0), fn = rw, hessian = T, lower = c(0, 0), upper = c(1, 50), method="L-BFGS-B")
    if(g==1){
      results[s, 1] <- s
      results[s, 2] <- rw.fits$par[1]
      results[s, 3] <- rw.fits$par[2]
      results[s, 4] <- rw.fits$value
      results[s, 5] <- rw.fits$convergence
    } else if(g==2){
      results[s, 6] <- rw.fits$par[1]
      results[s, 7] <- rw.fits$par[2]
      results[s, 8] <- rw.fits$value
      results[s, 9] <- rw.fits$convergence
    }
  }
}


# Simulated Prediction Errors ---------------------------------------------

for(s in 1:length(Ss)){
  S <- as.data.frame(Ss[[s]])
  S <- S[, c(1, 4, 6, 7)]
  colnames(S) <- c("group", "stim", "feedback", "choice")
  S$choice <- ifelse(S$choice == 1, 1, 0)
  stimOld <- as.integer(names(table(S$stim)))
  stimNew <- 1:length(table(S$stim))
  for(i in 1:length(stimOld)){
    S$stim[which(S$stim == stimOld[i])] <- stimNew[i]
  }
  for(g in 1:max(S$group)){
    data  <- S[which(S$group == g),]
    stim  <- data$stim
    freqs <- as.vector(table(data$stim))
    V     <- matrix(NA, nrow = max(freqs), ncol = length(freqs))
    L     <- matrix(NA, nrow = max(freqs), ncol = length(freqs))
    f     <- matrix(NA, nrow = max(freqs), ncol = length(freqs))
    r     <- matrix(NA, nrow = max(freqs), ncol = length(freqs))
    for(j in 1:length(freqs)){
      fData <- data$feedback[which(stim==j)]
      rData <- data$choice[which(stim==j)]
      if(length(fData) > 0){
        f[1:length(fData), j] <- fData
        r[1:length(rData), j] <- rData
      }
    }
    
    rw.fits <- optim(par = c(.1, 0), fn = rw, hessian = T, lower = c(0, 0), upper = c(1, 50), method="L-BFGS-B")
    if(g==1){
      results[s, 1] <- s
      results[s, 2] <- rw.fits$par[1]
      results[s, 3] <- rw.fits$par[2]
      results[s, 4] <- rw.fits$value
      results[s, 5] <- rw.fits$convergence
    } else if(g==2){
      results[s, 6] <- rw.fits$par[1]
      results[s, 7] <- rw.fits$par[2]
      results[s, 8] <- rw.fits$value
      results[s, 9] <- rw.fits$convergence
    }
  }
}



# Simulation Model --------------------------------------------------------
rw_PE <- function(params){
  for(s in 1:length(freqs)){
    for (t in 1:freqs[s]){
      if(t == 1){
        V[t, s] <- .5 + params[1]*(f[t, s] - .5)
        Err[t, s] <- .5
      } else{
        V[t, s] <- V[t-1, s] + params[1]*(f[t, s] - V[t-1, s])
        Err[t, s] <- f[t, s] - V[t-1, s]
      }
    }
  }
  return(Err)
}

# Simulation Response -----------------------------------------------------
SimOutput <- data.frame(Subject=integer(), Condition=integer(), Stimulus=integer())
for(s in 1:length(Ss)){
  S <- as.data.frame(Ss[[s]])
  S <- S[which(S$V1 != 3), c(1, 4, 6, 7)]
  names <- c("group", "stim", "feedback", "choice")
  colnames(S) <- names
  S$choice <- ifelse(S$choice == 1, 1, 0)
  for(g in 1:max(S$group)){
    data  <- S[which(S$group == g),]
    stim  <- data$stim
    freqs <- as.vector(table(data$stim))
    V     <- matrix(NA, nrow = max(freqs), ncol = 10)
    Err   <- matrix(NA, nrow = max(freqs), ncol = 10)
    f     <- matrix(NA, nrow = max(freqs), ncol = 10)
    for(j in 1:length(freqs)){
      fData <- data$feedback[which(stim==j)]
      if(length(fData) > 0){
        f[1:length(fData), j] <- fData
      }
    }
    col1  <- ifelse(g == 1, results[s, 2], results[s, 6])
    SimData <- as.data.frame(rw_PE(col1))
    SimData <- as.data.frame(cbind(1:10, rbind(t(SimData[1:10]))))
    for(k in 1:10){
      SimOutput[length(SimOutput$Stimulus) + 1, 1:(ncol(SimData) + 2)] <- cbind(s, g, SimData[k,])
    }
  }
}

S <- as.data.frame(Ss[[s]])
S <- S[, c(1, 4, 6, 7, 9, 11, 12)]
colnames(S) <- c("Condition", "Stimulus", "Feedback", "Choice", "StimOnset", "FeedbackOnset", "TrialEndTime")
S[which(S$Condition == 3), "TrialEndTime"] <- lead(S[which(S$Condition == 3), "TrialEndTime"], 3)
S[which(S$Condition == 3), "FeedbackOnset"] <- S[which(S$Condition == 3), "StimOnset"]
S$Keep <- NA
S[which(S$Condition == 3), "Keep"] <- rep(c(1, 0, 0, 0), 12)
S[which(S$Condition < 3), "Keep"] <- 1
S <- S[which(S$Keep == 1), ]
S$PredictionError <- NA
S[which(S$Condition == 3), "PredictionError"] <- 0
S$Duration <- S$TrialEndTime - S$FeedbackOnset
for(i in 1:max(S$Stimulus, na.rm = T)){
  for(j in 1:max(S$Condition)){
    PE <- SimOutput[which(SimOutput$Stimulus == i & SimOutput$Condition == j), 4:18]
    PE <- as.vector(unlist(PE[which(!is.na(PE))]))
    S[which(S$Stimulus == i & S$Condition == j), "PredictionError"] <- PE
  }
}

# Create 3 Col Text File --------------------------------------------------
S6_L1 <- cbind(S$FeedbackOnset[1:72], S$Duration[1:72], S$PredictionError[1:72])
write.table(S6_L1, file = "S6_L1_3Col.txt", sep = "\t",
            row.names = F, col.names = F)
