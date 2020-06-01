#script to compare conventional paired t test with permutation paired t test over all subjects (PPTT)
#nsim and nperm is number of simulations/permutations
start_time <- Sys.time() #start time simulation
library(multcomp)
permutation_fun <- function(nsim, nperm, rho, n, distribution, meanDiff = 0) {
  n1 <- n + 1
  n2 <- 2 * n
  tcrit <- qt(0.95, n - 1)
  PPTT1 <- c()
  #-------------------Generate Permutation Matrices------------------------------#
  P1 <- matrix(0, nrow = n2, ncol = nperm)
  for (h in 1:nperm) {
    P1[, h] <- sample(1:n2)
  }
  #-----------------------Data Generation------------------------------#
  t_test_vec <- numeric(nsim) #vector t test
  x <- matrix(0, ncol = nsim, nrow = n2)
  for (h in 1:nsim) {
    x11 <- rnorm(n)
    x22 <- rho * x11 + sqrt(1 - rho^2) * rnorm(n) #formula to make x11 and x22 correlated (pre and post)
    #homogenous distributions
    if (distribution == "Normal") { #generate normal data
      x[, h] <- c(x11, x22)
    }
    if (distribution == "Exp") {
      x[, h] <- c(qexp(pnorm(x11)), qexp(pnorm(x22))) #transform normal data to exponential
    }
    if (distribution == "LNorm") {
      x[, h] <- c(qlnorm(pnorm(x11)), qlnorm(pnorm(x22))) #transform to lognorm
    }
    if (distribution == "Gamma") {
      x[, h] <- c(qgamma(pnorm(x11),shape = 2), qgamma(pnorm(x22),shape = 2)) #transform to gamma, shape = 2 this resembles gamma in psychology but value can be changed
    }
    
#-------- add mean diff, add t test, compute means, variances and calculate t test value-------#    
    
    #adds mean difference for power calculation
    x[1:n,h] <- x[1:n,h] + meanDiff
    
    #adds t-test
    t_test_vec[h] <- t.test(x[1:n,h],x[n1:n2,h], paired = TRUE)$p.value < 0.05
  }
  x1 <- x[1:n, ]
  x2 <- x[n1:n2, ]
  
  diffs <- x1 - x2
  mdiff <- colMeans(diffs)
  vdiff <- (colSums(diffs^2) - n * mdiff^2) / (n - 1)
  Tpar <- sqrt(n) * (mdiff) / sqrt(vdiff) #t test formula
  NaT <- is.na(Tpar)
  Tpar[NaT] <- 5000 #value in case no value is identified
  #--------------Start of Simulation Loop----------------------------------#
  for (s in 1:nsim) {
    xx <- x[, s]
    #--------------------Permutation test------------------------------#
    xP1 <- matrix(xx[P1], ncol = nperm)
    xP11 <- xP1[1:n, ]
    xP12 <- xP1[n1:n2, ]
    DP1 <- xP11 - xP12
    mDP1 <- colMeans(DP1)
    vDP1 <- (colSums(DP1^2) - n * mDP1^2) / (n - 1)
    TP1 <- sqrt(n) * mDP1 / sqrt(vDP1) #permutation t test formula
    NAP1 <- is.na(TP1) 
    TP1[NAP1] <- 5000 #value in case no value is identified
    PPTTx <- 2*min(c(mean(Tpar[s] <= TP1),mean(Tpar[s] >= TP1)))
    PPTT1[s] <- (PPTTx < 0.05)
  }
  result <- data.frame(
    Distribution = distribution,
    rho = rho,
    tTest = mean(t_test_vec),
    PPTT1 = mean(PPTT1))
  return(result) #return results from simulation
}
sample_size <- c(10,20,50,100) #sample sizes used, can be changed
distribution <- c("Normal","Exp","LNorm","Gamma") #data generation distributions used
rho <- c(0.3,0.5,0.9) #within pair correlations, can be changed
meanDiff <- c(0,0.5) #0 for type I error rate, 1 for power 0.5 but any value different than 1 works
Design <- expand.grid(samp=sample_size,dist=distribution,rho=rho,meanDiff=meanDiff) #design parameter vector
n_permutations <- 10000 #10000 for accurate simulation results
n_simulations <- 10000 
type1_Pow_PPTT <- numeric(nrow(Design)) #create PPTT results vector
type1_Pow_conventional_t <- numeric(nrow(Design)) #create t test results vector
for (i in 1:nrow(Design)){
  tmp_results <- permutation_fun(n_simulations, n_permutations, Design$rho[i], Design$samp[i], distribution = Design$dist[i], meanDiff = Design$meanDiff[i])
  type1_Pow_conventional_t[i] <- tmp_results[['tTest']]
  type1_Pow_PPTT[i] <- tmp_results[['PPTT1']]
}
results <- cbind(Design,type1_Pow_PPTT,type1_Pow_conventional_t) #collect results into vector from design parameter vector, PPTT vector and t vector
end_time <- Sys.time() #end time simulation

