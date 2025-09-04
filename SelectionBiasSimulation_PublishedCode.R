#########################################################
##### Bias Simulations                              #####
##### Created by: Jennifer Head                     #####
##### Last Updated: August 15, 2025                   #####
#########################################################

#### SECTION 0. PREP ####

library(tidyverse)
library(cowplot)
set.seed(1234)  # ensure consistency across simulations
nreps <- 1000   # 1,000 simulations
ndraws <- c(1000, 10000) # sample size: 1,000 vs 10,000

#### SECTION 1.0 CONFOUNDING BIAS ####

# Range of values for Beta
BA <- seq(0, 0.5, 0.1)
BY <- seq(0, 0.5, 0.1)
# Note - code makes the assumption that BA = BY = B1. Can add an extra loop to 
# have BA and BY vary independently

## SIMULATION 1.1: SIMPLE CONFOUNDING DAG WITH POSITIVES EDGES
## A (+) <- C -> (+) Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(BA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
      
      # simulate Z with prevalence of 50%
      Z <- rbinom(ndraws[k], 1, 0.5)
      
      # Simulate A and Y as a function of Z
      A <- rbinom(ndraws[k], 1, prob = plogis(0.01 + BA[j]*Z))
      Y <- rbinom(ndraws[k], 1, prob = plogis(-0.5 + BY[j]*Z))
      
      # Compute components of the unbiased risk difference
      EY1 <- mean(Y[A == 1 & Z == 1] )* mean(Z==1) + mean(Y[A == 1 & Z == 0]) * mean(Z==0)
      EY0 <- mean(Y[A == 0 & Z == 1]) * mean(Z==1) + mean(Y[A == 0 & Z == 0]) * mean(Z==0)
      
      # Compute components of the biased risk difference
      S1 <- mean(Y[A == 1])
      S0 <- mean(Y[A == 0])
      
      # Save the values of B1 and the bias parameters
      save_RDj[rep] <- (S1-S0) - (EY1-EY0) # Biased RD - Unbiased RD
      save_RRj[rep] <- (S1/S0) - (EY1/EY0) # Biased RR - Unbiased RR
      B1j <- rep(BA[j], nreps)
      nj <- rep(ndraws[k], nreps)
      
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }
    
}
  
# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding1 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding1$n <- paste0("Sample size: ", bias_confounding1$n)

a <- ggplot(bias_confounding1) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (NOT adjusting for Z)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

b <- ggplot(bias_confounding1) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (NOT adjusting for Z)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(a, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(b, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
             x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)


## SIMULATION 1.2: SIMPLE CONFOUNDING DAG WITH POSITIVE & NEGATIVE EDGES
## A (+) <- C -> (-) Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(BA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
      
    # simulate Z with prevalence of 50%
    Z <- rbinom(ndraws[k], 1, 0.5)
    
    # Simulate A and Y as a function of Z
    A <- rbinom(ndraws[k], 1, prob = plogis(0.01 + BA[j]*Z))
    Y <- rbinom(ndraws[k], 1, prob = plogis(-0.5 - BY[j]*Z))
    
    # Compute components of the unbiased risk difference
    EY1 <- mean(Y[A == 1 & Z == 1] )* mean(Z==1) + mean(Y[A == 1 & Z == 0]) * mean(Z==0)
    EY0 <- mean(Y[A == 0 & Z == 1]) * mean(Z==1) + mean(Y[A == 0 & Z == 0]) * mean(Z==0)
    
    # Compute components of the biased risk difference
    S1 <- mean(Y[A == 1])
    S0 <- mean(Y[A == 0])
    
    # Save the values of B1 and the bias parameters
    save_RDj[rep] <- (S1-S0) - (EY1-EY0) # Biased RD - Unbiased RD
    save_RRj[rep] <- (S1/S0) - (EY1/EY0) # Biased RR - Unbiased RR
    B1j <- rep(BA[j], nreps)
    nj <- rep(ndraws[k], nreps)
    
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }
  
}


# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding2 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding2$n <- paste0("Sample size: ", bias_confounding2$n)

c <- ggplot(bias_confounding2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (NOT adjusting for Z)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

d <- ggplot(bias_confounding2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (NOT adjusting for Z)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(c, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(d, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)


## SIMULATION 1.3: SIMPLE CONFOUNDING DAG WITH BOTH NEGATIVE EDGES
## A (-) <- C -> (-) Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(BA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
    
    # simulate Z with prevalence of 50%
    Z <- rbinom(ndraws[k], 1, 0.5)
    
    # Simulate A and Y as a function of Z
    A <- rbinom(ndraws[k], 1, prob = plogis(0.01 - BA[j]*Z))
    Y <- rbinom(ndraws[k], 1, prob = plogis(-0.5 - BY[j]*Z))
    
    # Compute components of the unbiased risk difference
    EY1 <- mean(Y[A == 1 & Z == 1] )* mean(Z==1) + mean(Y[A == 1 & Z == 0]) * mean(Z==0)
    EY0 <- mean(Y[A == 0 & Z == 1]) * mean(Z==1) + mean(Y[A == 0 & Z == 0]) * mean(Z==0)
    
    # Compute components of the biased risk difference
    S1 <- mean(Y[A == 1])
    S0 <- mean(Y[A == 0])
    
    # Save the values of B1 and the bias parameters
    save_RDj[rep] <- (S1-S0) - (EY1-EY0) # Biased RD - Unbiased RD
    save_RRj[rep] <- (S1/S0) - (EY1/EY0) # Biased RR - Unbiased RR
    B1j <- rep(BA[j], nreps)
    nj <- rep(ndraws[k], nreps)
    
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
  
  }

}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding3 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding3$n <- paste0("Sample size: ", bias_confounding3$n)

e <- ggplot(bias_confounding3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (NOT adjusting for Z)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

f <- ggplot(bias_confounding3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (NOT adjusting for Z)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(e, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(f, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)

#### SECTION 2.0 SELECTION BIAS ####

## SIMULATION 2.1: SIMPLE CONFOUNDING DAG WITH POSITIVES EDGES
## A (+) <- C -> (+) Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(BA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
    
    # simulate A and Y with prevalence of 50% and 40%
    A <- rbinom(ndraws[k], 1, 0.5)
    Y <- rbinom(ndraws[k], 1, 0.4)
    
    # Simulate for C as a function of both UA and UY
    C <- rbinom(ndraws[k], 1, prob = plogis(0.01 + BA[j]*A + BY[j]*Y))
    
    # Compute components of the biased risk difference
    EY1 <- mean(Y[A == 1 & C == 1] )* mean(C==1) + mean(Y[A == 1 & C == 0]) * mean(C==0)
    EY0 <- mean(Y[A == 0 & C == 1]) * mean(C==1) + mean(Y[A == 0 & C == 0]) * mean(C==0)
    
    # Compute components of the unbiased risk difference
    S1 <- mean(Y[A == 1])
    S0 <- mean(Y[A == 0])
    
    # Save the values of B1 and the bias parameters
    save_RDj[rep] <- (EY1-EY0) - (S1-S0) # Biased RD - Unbiased RD
    save_RRj[rep] <- (EY1/EY0) - (S1/S0) # Biased RR - Unbiased RR
    B1j <- rep(BA[j], nreps)
    nj <- rep(ndraws[k], nreps)
    
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }
  
}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding4 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding4$n <- paste0("Sample size: ", bias_confounding4$n)

g <- ggplot(bias_confounding4) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

h <- ggplot(bias_confounding4) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(g, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(h, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)

## SIMULATION 2.2: SIMPLE COLLIDER DAG WITH ONE POSITIVE EDGE AND ONE NEGATIVE
## A (+) -> C <- (-) Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(BA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
    
      # simulate A and Y with prevalence of 50% and 40%
      A <- rbinom(ndraws[k], 1, 0.5)
      Y <- rbinom(ndraws[k], 1, 0.4)
      
      # Simulate for C as a function of both UA and UY
      C <- rbinom(ndraws[k], 1, prob = plogis(0.01 + BA[j]*A - BY[j]*Y))
      
      # Compute components of the biased risk difference
      EY1 <- mean(Y[A == 1 & C == 1] )* mean(C==1) + mean(Y[A == 1 & C == 0]) * mean(C==0)
      EY0 <- mean(Y[A == 0 & C == 1]) * mean(C==1) + mean(Y[A == 0 & C == 0]) * mean(C==0)
      
      # Compute components of the unbiased risk difference
      S1 <- mean(Y[A == 1])
      S0 <- mean(Y[A == 0])
      
      # Save the values of B1 and the bias parameters
      save_RDj[rep] <- (EY1-EY0) - (S1-S0) # Biased RD - Unbiased RD
      save_RRj[rep] <- (EY1/EY0) - (S1/S0) # Biased RR - Unbiased RR
      B1j <- rep(BA[j], nreps)
      nj <- rep(ndraws[k], nreps)
      
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
  
  }

}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding5 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding5$n <- paste0("Sample size: ", bias_confounding5$n)

i <- ggplot(bias_confounding5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

j <- ggplot(bias_confounding5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(i, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(j, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)

## SIMULATION 2.3: SIMPLE COLLIDER DAG WITH BOTH NEGATIVE EDGES
## A (-) -> C <- (-) Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(BA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
    
      # simulate A and Y with prevalence of 50% and 40%
      A <- rbinom(ndraws[k], 1, 0.5)
      Y <- rbinom(ndraws[k], 1, 0.4)
      
      # Simulate for C as a function of both UA and UY
      C <- rbinom(ndraws[k], 1, prob = plogis(0.01 - BA[j]*A - BY[j]*Y))
      
      # Compute components of the biased risk difference
      EY1 <- mean(Y[A == 1 & C == 1] )* mean(C==1) + mean(Y[A == 1 & C == 0]) * mean(C==0)
      EY0 <- mean(Y[A == 0 & C == 1]) * mean(C==1) + mean(Y[A == 0 & C == 0]) * mean(C==0)
      
      # Compute components of the unbiased risk difference
      S1 <- mean(Y[A == 1])
      S0 <- mean(Y[A == 0])
      
      # Save the values of B1 and the bias parameters
      save_RDj[rep] <- (EY1-EY0) - (S1-S0) # Biased RD - Unbiased RD
      save_RRj[rep] <- (EY1/EY0) - (S1/S0) # Biased RR - Unbiased RR
      B1j <- rep(BA[j], nreps)
      nj <- rep(ndraws[k], nreps)
    
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }
  
}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding6 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding6$n <- paste0("Sample size: ", bias_confounding6$n)

k <- ggplot(bias_confounding6) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

l <- ggplot(bias_confounding6) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = BA, labels = BA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(k, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(l, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)


#### SECTION 3.0 MORE COMPLEX (FIVE NODE) SELECTION BIAS ####


# Define new set of values for the beta's (now we have 4)
B_UaA <- seq(0, 1., 0.2)
B_UyY <- seq(0, 1., 0.2)
B_UaC <- seq(0, 1., 0.2)
B_UyC <- seq(0, 1., 0.2)

# Note - code makes the assumption that B_UaA = B_UyY = B_UaC = B_UyC = B1. Can add an extra loop to 
# have beta's vary independently


## SIMULATION 3.1: MULTI-NODE COLLIDER DAG WITH POSITIVES EDGES
## A <- (+) UA (+) -> C <- (+) UY (+) -> Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(B_UaA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
    
      # simulate the unknown factors with prevalence of 50% and 40%
      UA <- rbinom(ndraws[k], 1, 0.5)
      UY <- rbinom(ndraws[k], 1, 0.4)
      
      # Simulate A as a function of UA and Y as a function of UY 
      A <- rbinom(ndraws[k], 1, prob = plogis(0.2 + B_UaA[j]*UA))
      Y <- rbinom(ndraws[k], 1, prob = plogis(-0.2 + B_UyY[j]*UY))
      
      # Simulate  C as a function of both UA and UY
      C <- rbinom(ndraws[k], 1, prob = plogis(0.01 + B_UaC[j]*UA + B_UyC[j]*UY))
      
      # Compute components of the biased risk difference
      EY1 <- mean(Y[A == 1 & C == 1] )* mean(C==1) + mean(Y[A == 1 & C == 0]) * mean(C==0)
      EY0 <- mean(Y[A == 0 & C == 1]) * mean(C==1) + mean(Y[A == 0 & C == 0]) * mean(C==0)
      
      # Compute components of the unbiased risk difference
      S1 <- mean(Y[A == 1])
      S0 <- mean(Y[A == 0])
      
      # Save the values of B1 and the bias parameters
      save_RDj[rep] <- (EY1-EY0) - (S1-S0) # Biased RD - Unbiased RD
      save_RRj[rep] <- (EY1/EY0) - (S1/S0) # Biased RR - Unbiased RR
      B1j <- rep(B_UaA[j], nreps)
      nj <- rep(ndraws[k], nreps)
      
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }
  
}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding7 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding7$n <- paste0("Sample size: ", bias_confounding7$n)

m <- ggplot(bias_confounding7) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = B_UaA, labels = B_UaA) +
  facet_wrap(~n)

n <- ggplot(bias_confounding7) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = B_UaA, labels = B_UaA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(m, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(n, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)

## SIMULATION 3.2: MULTI-NODE COLLIDER DAG WITH ONE NEGATIVE EDGE
## A <- (-) UA (+) -> C <- (+) UY (+) -> Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(B_UaA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
    
      # simulate the unknown factors with prevalence of 50% and 40%
      UA <- rbinom(ndraws[k], 1, 0.5)
      UY <- rbinom(ndraws[k], 1, 0.4)
      
      # Simulate A as a function of UA and Y as a function of UY 
      A <- rbinom(ndraws[k], 1, prob = plogis(0.2 - B_UaA[j]*UA))
      Y <- rbinom(ndraws[k], 1, prob = plogis(-0.2 + B_UyY[j]*UY))
      
      # Simulate  C as a function of both UA and UY
      C <- rbinom(ndraws[k], 1, prob = plogis(0.01 + B_UaC[j]*UA + B_UyC[j]*UY))
      
      # Compute components of the biased risk difference
      EY1 <- mean(Y[A == 1 & C == 1] )* mean(C==1) + mean(Y[A == 1 & C == 0]) * mean(C==0)
      EY0 <- mean(Y[A == 0 & C == 1]) * mean(C==1) + mean(Y[A == 0 & C == 0]) * mean(C==0)
      
      # Compute components of the unbiased risk difference
      S1 <- mean(Y[A == 1])
      S0 <- mean(Y[A == 0])
      
      # Save the values of B1 and the bias parameters
      save_RDj[rep] <- (EY1-EY0) - (S1-S0) # Biased RD - Unbiased RD
      save_RRj[rep] <- (EY1/EY0) - (S1/S0) # Biased RR - Unbiased RR
      B1j <- rep(B_UaA[j], nreps)
      nj <- rep(ndraws[k], nreps)
      
      }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }
  
}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding8 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding8$n <- paste0("Sample size: ", bias_confounding8$n)

o <- ggplot(bias_confounding8) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = B_UaA, labels = B_UaA) +
  facet_wrap(~n)

p <- ggplot(bias_confounding8) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (Selection on C)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = B_UaA, labels = B_UaA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(o, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(p, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)


#### SECTION 4.0 INFORMATION BIAS (Exposure misclassification) ####

# Define new set of values for the beta's (now we have 4)
B_UA <- seq(0, 1.5, 0.25)
B_UY <- seq(0, 1.5, 0.25)
# Note - code makes the assumption that B_UA = B_UY = B1. Can add an extra loop to 
# have beta's vary independently

## SIMULATION 4.1: DIFFERENTIAL EXPOSURE MISCLASSIFICATION, POSITIVE EDGES
## A -> (+) A* <- (+) U (+) -> Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(B_UA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
    
      # simulate A and U with prevalence of 50% and 40%
      A <- rbinom(ndraws[k], 1, 0.5)
      U <- rbinom(ndraws[k], 1, 0.4)
  
      # Simulate Y as a function of U 
      Y <- rbinom(ndraws[k], 1, prob = plogis(-0.2 + B_UY[j]*U))
      
      # Simulate A* as a function of A and U 
      A_star <- rbinom(ndraws[k], 1, prob = plogis(2*(A==1) - 2*(A==0) + B_UA[j]*U))
  
      # Compute components of the observed risk difference
      EY1_star <- mean(Y[A_star == 1])
      EY0_star <- mean(Y[A_star == 0])
      
      # Compute components of the true risk difference
      EY1 <- mean(Y[A == 1])
      EY0 <- mean(Y[A == 0])
      
      # Save the values of B1 and the bias parameters
      save_RDj[rep] <- (EY1_star - EY0_star) - (EY1 - EY0) # Biased RD - Unbiased RD
      save_RRj[rep] <- (EY1_star / EY0_star) - (EY1 / EY0) # Biased RR - Unbiased RR
      B1j <- rep(B_UA[j], nreps)
      nj <- rep(ndraws[k], nreps)
    
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }
  
}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding9 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding9$n <- paste0("Sample size: ", bias_confounding9$n)

q <- ggplot(bias_confounding9) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (Misclassification of A)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = B_UA, labels = B_UA) +
  facet_wrap(~n)

r <- ggplot(bias_confounding9) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (Misclassification of A)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = B_UA, labels = B_UA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(q, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(r, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)


## SIMULATION 4.1: DIFFERENTIAL EXPOSURE MISCLASSIFICATION, POSITIVE & NEGATIVE EDGES
## A -> (+) A* <- (-) U (+) -> Y

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(B_UA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
    
      # simulate A and U with prevalence of 50% and 40%
      A <- rbinom(ndraws[k], 1, 0.5)
      U <- rbinom(ndraws[k], 1, 0.4)
      
      # Simulate Y as a function of U 
      Y <- rbinom(ndraws[k], 1, prob = plogis(-0.2 - B_UY[j]*U))
      
      # Simulate A* as a function of A and U 
      A_star <- rbinom(ndraws[k], 1, prob = plogis(2*(A==1) - 2*(A==0) + B_UA[j]*U))
      
      # Compute components of the observed risk difference
      EY1_star <- mean(Y[A_star == 1])
      EY0_star <- mean(Y[A_star == 0])
      
      # Compute components of the true risk difference
      EY1 <- mean(Y[A == 1])
      EY0 <- mean(Y[A == 0])
      
      # Save the values of B1 and the bias parameters
      save_RDj[rep] <- (EY1_star - EY0_star) - (EY1 - EY0) # Biased RD - Unbiased RD
      save_RRj[rep] <- (EY1_star / EY0_star) - (EY1 / EY0) # Biased RR - Unbiased RR
      B1j <- rep(B_UA[j], nreps)
      nj <- rep(ndraws[k], nreps)
      
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }

}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding10 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding10$n <- paste0("Sample size: ", bias_confounding10$n)

s <- ggplot(bias_confounding10) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = B1)) +
  theme_bw() +
  ylab("Bias (Misclassification of A)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = B_UA, labels = B_UA) +
  facet_wrap(~n)

t <- ggplot(bias_confounding10) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RR, group = B1)) +
  theme_bw() +
  ylab("Bias (Misclassification of A)") +
  xlab(expression(~beta)) +
  scale_x_continuous(breaks = B_UA, labels = B_UA) +
  facet_wrap(~n)

ggdraw() +
  draw_plot(s, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(t, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(label = c("A. Risk Difference", "B. Risk Ratio"),
                  x = c(0,0), y = c(1, 0.5), size = 13, hjust = -0.3)


### Make the plot for the main text now

BA <- seq(0, 1, 0.2)
BY <- seq(0, 1, 0.2)

# Initialize storage vectors
save_RD <- save_RR <- B1 <- n <- c()

# Loop through each value of B1
for (k in 1:length(ndraws)){
  
  for (j in 1:length(BA)){
    
    save_RDj <- save_RRj <- rep(NA, nreps) # initialize storage for single simulation
    
    for (rep in 1:nreps){
      
      # simulate A and Y with prevalence of 50% and 40%
      A <- rbinom(ndraws[k], 1, 0.5)
      Y <- rbinom(ndraws[k], 1, 0.4)
      
      # Simulate for C as a function of both UA and UY
      C <- rbinom(ndraws[k], 1, prob = plogis(0.01 + BA[j]*A + BY[j]*Y))
      
      # Compute components of the biased risk difference
      EY1 <- mean(Y[A == 1 & C == 1] )* mean(C==1) + mean(Y[A == 1 & C == 0]) * mean(C==0)
      EY0 <- mean(Y[A == 0 & C == 1]) * mean(C==1) + mean(Y[A == 0 & C == 0]) * mean(C==0)
      
      # Compute components of the unbiased risk difference
      S1 <- mean(Y[A == 1])
      S0 <- mean(Y[A == 0])
      
      # Save the values of B1 and the bias parameters
      save_RDj[rep] <- (EY1-EY0) - (S1-S0) # Biased RD - Unbiased RD
      save_RRj[rep] <- (EY1/EY0) - (S1/S0) # Biased RR - Unbiased RR
      B1j <- rep(BA[j], nreps)
      nj <- rep(ndraws[k], nreps)
      
    }
    
    # concatenate results from each simulation
    save_RR <- c(save_RR, save_RRj)
    save_RD <- c(save_RD, save_RDj)
    B1 <- c(B1, B1j)
    n <- c(n, nj)
    
  }
  
}

# Save results as a dataframe and create boxplot of bias by Beta1 and n
bias_confounding11 <- data.frame(cbind(save_RD, save_RR, B1, n))
bias_confounding11$n <- paste0("Sample size: ", bias_confounding11$n)

u <- ggplot(bias_confounding11) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = paste0(B1, n), fill = n)) +
  theme_bw() +
  ylab("Bias (Controlling for C)") +
  xlab(expression(~beta)) +
  scale_fill_manual("Sample size",
                    breaks = c("Sample size: 1000",
                               "Sample size: 10000"),
                    labels = c("1,000", "10,000"),
                    values = c("lightblue2", "darkseagreen3")) +
  ylim(c(-0.07, 0.01))

v <- ggplot(bias_confounding7) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = B1, y = save_RD, group = paste0(B1, n), fill = n)) +
  theme_bw() +
  ylab("Bias (Controlling for C)") +
  xlab(expression(~beta)) +
  scale_fill_manual("Sample size",
                    breaks = c("Sample size: 1000",
                               "Sample size: 10000"),
                    labels = c("1,000", "10,000"),
                    values = c("lightblue2", "darkseagreen3")) +
  ylim(c(-0.07, 0.01))

ggpubr::ggarrange(u,v, common.legend = T, legend = "right",
                  labels = c("A)", "B)"))

ggsave("C:\\Users\\jrhead\\University of Michigan Dropbox\\Jennifer Head\\Teaching\\BiasMS\\Figures\\Simulations.jpg",
       dpi = 600, width = 10, height = 7)
