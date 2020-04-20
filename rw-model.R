setwd("C:/Users/David/Documents/GitHub/affect-and-multidimensional-rl/")

library(R.matlab)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)

change to document 2 56
another change 47
another line of code
# Import Data -------------------------------------------------------------

Ss  <- list(
  S1  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_3_2.mat")$rec.task,
  S2  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_4_1.mat")$rec.task,
  S3  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_5_1.mat")$rec.task,
  S4  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_6_2.mat")$rec.task,
  S5  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_7_2.mat")$rec.task,
  S6  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_8_1.mat")$rec.task,
  S7  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_9_2.mat")$rec.task,
  S8  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_10_2.mat")$rec.task,
  S9  <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_11_1.mat")$rec.task,
  S10 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_12_1.mat")$rec.task,
  S11 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_13_2.mat")$rec.task,
  S12 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_14_1.mat")$rec.task,
  S13 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_15_1.mat")$rec.task,
  S14 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_17_2.mat")$rec.task,
  S15 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_18_1.mat")$rec.task,
  S16 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_19_1.mat")$rec.task,
  S17 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_20_2.mat")$rec.task,
  S18 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_21_2.mat")$rec.task,
  S19 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_22_1.mat")$rec.task,
  S20 <- readMat("C:/Users/David/Desktop/Data_Analysis_Scripts/Data_Analysis_Scripts/AM_MRI_analysis/AM_MRI_dataonly/AM_MRI_LEARN_23_1.mat")$rec.task
)




f_alpha <- function(mu, sigma){
  alpha <- (((1-mu)/(sigma^2)) - (1/mu))*(mu^2)
  return(alpha)
}

f_beta <- function(alpha, mu){
  beta <- alpha*((1/mu) - 1)
  return(beta)
}




# Define Model (Pearce Hall)-------------------------------------------------------------------

rw <- function(params){
  Priors  <- c(params[4], (1-params[4]))
  alpha   <- (((1-params[1])/(params[2]^2)) - (1/params[1]))*(params[1]^2)
  beta    <- alpha*((1/params[1]) - 1)
  A       <- rep(alpha, 4)
  B       <- rep(beta, 4)
  A.Pairs <- rep(alpha, 6)
  B.Pairs <- rep(beta, 6)
  
  for (t in 1:length(data$group)){
    
    # Define Stimuli
    X <- c(
      as.integer(data[t, "Stim1"] == 1 | data[t, "Stim2"] == 1),
      as.integer(data[t, "Stim1"] == 2 | data[t, "Stim2"] == 2),
      as.integer(data[t, "Stim1"] == 3 | data[t, "Stim2"] == 3),
      as.integer(data[t, "Stim1"] == 4 | data[t, "Stim2"] == 4))
    X.Pairs <- c(
      as.integer(data[t, "Stim1"] == 1 & data[t, "Stim2"] == 2),
      as.integer(data[t, "Stim1"] == 1 & data[t, "Stim2"] == 3),
      as.integer(data[t, "Stim1"] == 1 & data[t, "Stim2"] == 4),
      as.integer(data[t, "Stim1"] == 2 & data[t, "Stim2"] == 3),
      as.integer(data[t, "Stim1"] == 2 & data[t, "Stim2"] == 4),
      as.integer(data[t, "Stim1"] == 3 & data[t, "Stim2"] == 4))
    
    # Compute expectations conditioned on each model
    Cond.Exp <- c(((X%*%A)/sum(X)) / (((X%*%A)/sum(X)) + ((X%*%B)/sum(X))), X.Pairs%*%A.Pairs / (X.Pairs%*%A.Pairs + X.Pairs%*%B.Pairs))
    Cond.Exp[is.na(Cond.Exp)] <- 0
    
    # Compute expectation given trial type (single vs. double stimulus)
    v <- Cond.Exp[1]^as.integer(sum(X) == 1) * (Cond.Exp%*%Priors)^as.integer(sum(X) == 2)
    
    # Compute likelihood of usbject's choice
    P.Choice <- (v^data[t, "choice"])*((1-v)^(1-data[t, "choice"]))
    L[t] <- exp(params[3]*P.Choice) / (exp(params[3]*v) + exp(params[3]*(1-v)))
    
    # Update Alpha and Beta
    A[which(X == 1)] <- A[which(X == 1)] + data[t, "choice"]*as.integer(sum(X) == 1) + (data[t, "choice"]/2)*as.integer(sum(X) == 2)
    A.Pairs[which(X.Pairs == 1)] <- A.Pairs[which(X.Pairs == 1)] + data[t, "choice"]*as.integer(sum(X) == 2)
    B[which(X == 1)] <- B[which(X == 1)] + (1-data[t, "choice"])*as.integer(sum(X) == 1) + ((1-data[t, "choice"])/2)*as.integer(sum(X) == 2)
    B.Pairs[which(X.Pairs == 1)] <- B.Pairs[which(X.Pairs == 1)] + (1 - data[t, "choice"])*as.integer(sum(X) == 2)
    
    # Update Priors
    JointLs <- (Cond.Exp^data[t, "choice"])*((1-Cond.Exp)^(1-data[t, "choice"]))*Priors
    Priors <- Priors^as.integer(sum(X) == 1)*(JointLs / sum(JointLs))^as.integer(sum(X) == 2)
  }
  log.likelihoods <- log(L)
  MAP <- -2 * sum(log.likelihoods, na.rm = TRUE)
  return(MAP)
}


# Fit Model (Pearce Hall)----------------------------------------------

results <- matrix(NA, nrow=length(Ss), ncol=9,
                  dimnames = list(
                    c(), 
                    c("Subject", "Threat_Alpha", "Threat_Beta", "Threat_LL", "Threat_Conv", "Safety_Alpha", "Safety_Beta", "Safety_LL", "Safety_Conv")))

for(s in 1:length(Ss)){
  
  # Prep Data
  S <- Ss[[s]]
  S <- as.data.frame(S)
  S <- S[which(S$V1 != 3 & !grepl("NaN", S$V7)), c(1, 4, 6, 7)]
  names <- c("group", "stim", "feedback", "choice")
  colnames(S) <- names
  S$choice <- ifelse(S$choice == 1, 1, 0)
  S$win <- ifelse(S$choice == S$feedback, 1, 0)
  pairs <- c("1","2","3","4","1,2","1,3","1,4","2,3","2,4","3,4")
  S$pairs <- sapply(S$stim, function(x) pairs[x])
  S <- separate(S, "pairs", into = c("Stim1", "Stim2"), sep = ",")
  S[which(is.na(S$Stim2)), "Stim2"] <- 99
  
  # Set Priors
  for(g in 1:max(S$group)){
    data  <- S[which(S$group == g),]
    L <- NA
    rw.fits <- optim(par = c(.5, .2, 0, .5), fn = rw, hessian = T, lower = c(.01, .01, 0, 0), upper = c(.99, .2886751, 50, 1), method="L-BFGS-B")
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



names(RealOutput)[6:20] <- 1:15
RealOutput$Simulated <- 0
names(SimOutput)[6:20] <- 1:15
SimOutput$Simulated <- 1

Stack <- rbind(RealOutput, SimOutput)

SimLong <- melt(Stack,
                id.vars = c("Simulated", "Condition", "Stimulus"),
                measure.vars = 6:20,
                variable.name = "Trial", value.name = "Sunny")
SimLong$Condition <- as.factor(SimLong$Condition)
levels(SimLong$Condition) <- c("Threat", "Safety")
SimLong$Simulated <- as.factor(SimLong$Simulated)
levels(SimLong$Simulated) <- c("Subject", "Model")
SimLong$Stimulus <- as.factor(SimLong$Stimulus)
levels(SimLong$Stimulus) <- c("Stim 1", "Stim 2", "Stim 3", "Stim 4", "Stim 5",
                              "Stim 6", "Stim 7", "Stim 8", "Stim 9", "Stim 10")


tiff("test2.tiff", units="in", width=8, height=3, res=300)
ggplot(data=SLsim) +
  #  theme_bw() +
  theme(panel.spacing.x=unit(.5, "lines") , panel.spacing.y=unit(.5,"lines"))+
  scale_color_brewer(palette = "Set1") +
  theme(legend.title=element_blank()) +
  geom_line(aes(x=Trial, y=Sunny), group = 1, size = .75, color = "blue") + 
  geom_point(aes(x=Trial, y=Sunny2), size = 1) +
  facet_grid(rows = vars(Condition), cols = vars(Stimulus)) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1))
dev.off()


curve(dbeta(x,81,219),col = "blue", xlab = "Batting Average", ylab = "Probability", xlim=c(0.1,0.4), ylim=c(0,25))
par(new=TRUE)
curve(dbeta(x,(10),(10)),col = "red", xlab = "Batting Average", ylab = "Probability", xlim=c(0.1,0.4),ylim=c(0,25))