library(coin)
library(extraDistr) #Laplace distribution
library(ggplot2) #Plotting barplots
#----------------------------------------------------------------------------------------------
#Code written by Toon
median.test.two_sided <- function(vec1, vec2, N=1000, exact=FALSE) {
  # Calculate observed test statistic
  observed_test_statistic <- abs(median(vec1) - median(vec2))
  
  # Combine both vectors
  full_vec <- c(vec1, vec2)
  
  # Check if the number of combinations is less than 10,000
  if (choose(length(full_vec), length(vec1)) < N | exact) {
    
    # Generate the null distribution by exact combinations
    t_star <- combn(length(full_vec), length(vec1), function(x) {
      group1 <- full_vec[x]
      group2 <- full_vec[-x]  # all elements not in x
      test_stats <- abs(median(group1) - median(group2))
      test_stats
    })
  } else {
    
    #message("Using approximation method with 10,000 random permutations.")
    
    # Approximate null distribution by permutation if combinations exceed 10,000
    #N <- 500
    Groupvec <- c(rep(1, length(vec1)), rep(2, length(vec2)))  # Create group vector
    perm_data <- data.frame(value = full_vec, group = Groupvec)  # Create data frame for permutations
    t_star <- replicate(N, {
      # Shuffle the group labels:
      perm_data$group <- sample(perm_data$group)
      # Calculate the test statistic for shuffled groups:
      test_stats <- abs(median(perm_data$value[perm_data$group == 1]) - median(perm_data$value[perm_data$group == 2]))
      test_stats
    })
  }
  
  # Calculate the p-value
  p_value <- mean(t_star >= observed_test_statistic)
  
  p_value
}


#----------------------------------------------------------------------------------------------

# Sample simulation study of Q2 for t3:


set.seed(100) #For reproductibility
nsim <- 500
#t3 distribution
delta <- sqrt(3)/2 #delta = sqrt(variance)/2
p.t <- p.wmw <- p.med <- numeric(nsim) #empty vector of p values for each test


sample_sizes <- c(5, 10, 20, 50, 100) # sample sizes for which we do the test. 
n_sample_sizes <- length(sample_sizes)


results_t3 <- data.frame(
  sample_size = sample_sizes,
  t_test_power = numeric(n_sample_sizes),
  wmw_test_power = numeric(n_sample_sizes),
  med_test_power = numeric(n_sample_sizes)
)

for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim)
  {
    Y1 <- rt(n, 3)
    Y2 <- rt(n, 3) + delta
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    #Wilcoxon-Mann-Whitney test
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    #median test
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  # Monte-Carlo approximation of the power
  # for alpha = .05
  results_t3$t_test_power[j] <- mean(p.t < 0.05)
  results_t3$wmw_test_power[j] <- mean(p.wmw < 0.05)
  results_t3$med_test_power[j] <- mean(p.med < 0.05)
}
  


#----------------------------------------------------------------------------------------------

# Sample simulation study of Q2 for Exponential distribution

results_exp <- data.frame(
  sample_size = sample_sizes,
  t_test_power = numeric(n_sample_sizes),
  wmw_test_power = numeric(n_sample_sizes),
  med_test_power = numeric(n_sample_sizes)
)

#variance = 1/rate**2
delta <- sqrt(1)/2 #delta = sqrt(variance)/2
for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for (i in 1:nsim) {
    Y1 <- rexp(n, rate = 1)
    Y2 <- rexp(n, rate = 1)+delta
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2, exact=TRUE)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_exp$t_test_power[j] <- mean(p.t < 0.05)
  results_exp$wmw_test_power[j] <- mean(p.wmw < 0.05)
  results_exp$med_test_power[j] <- mean(p.med < 0.05)
}




#----------------------------------------------------------------------------------------------

# Sample simulation study of Q2 for Laplace distribution

results_lap <- data.frame(
  sample_size = sample_sizes,
  t_test_power = numeric(n_sample_sizes),
  wmw_test_power = numeric(n_sample_sizes),
  med_test_power = numeric(n_sample_sizes)
)
# variance = 2*sigma**2
delta <- sqrt(2)/2 #delta = sqrt(variance)/2

for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim){
    Y1 <- rlaplace(n, mu = 0, sigma = 1)
    Y2 <- rlaplace(n, mu = delta, sigma = 1)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_lap$t_test_power[j] <- mean(p.t < 0.05)
  results_lap$wmw_test_power[j] <- mean(p.wmw < 0.05)
  results_lap$med_test_power[j] <- mean(p.med < 0.05)
}

#----------------------------------------------------------------------------------------------

# Sample simulation study of Q2 for t5 distribution

results_t5 <- data.frame(
  sample_size = sample_sizes,
  t_test_power = numeric(n_sample_sizes),
  wmw_test_power = numeric(n_sample_sizes),
  med_test_power = numeric(n_sample_sizes)
)

#variance t-distribution = dof/(dof-2)
#with dof the degrees of freedom
delta <- sqrt(5/3)/2 #delta = sqrt(variance)/2

for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim)
  {
    Y1 <- rt(n, 5)
    Y2 <- rt(n, 5) + delta
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_t5$t_test_power[j] <- mean(p.t < 0.05)
  results_t5$wmw_test_power[j] <- mean(p.wmw < 0.05)
  results_t5$med_test_power[j] <- mean(p.med < 0.05)
}

#----------------------------------------------------------------------------------------------

# Sample simulation study of Q2 for Logistic distribution


results_log <- data.frame(
  sample_size = sample_sizes,
  t_test_power = numeric(n_sample_sizes),
  wmw_test_power = numeric(n_sample_sizes),
  med_test_power = numeric(n_sample_sizes)
)

#variance = scale**2 * pi**2 / 3 
delta <- sqrt(pi**2/3)/2 #delta = sqrt(variance)/2
for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim)
  {
    Y1 <- rlogis(n, location=0, scale=1)
    Y2 <- rlogis(n, location = delta, scale=1)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_log$t_test_power[j] <- mean(p.t < 0.05)
  results_log$wmw_test_power[j] <- mean(p.wmw < 0.05)
  results_log$med_test_power[j] <- mean(p.med < 0.05)
}



#----------------------------------------------------------------------------------------------

# Sample simulation study of Q2 for normal distribution

results_norm <- data.frame(
  sample_size = sample_sizes,
  t_test_power = numeric(n_sample_sizes),
  wmw_test_power = numeric(n_sample_sizes),
  med_test_power = numeric(n_sample_sizes)
)
#variance = sqrt(sigma)
delta <- sqrt(1)/2 #delta = sqrt(variance)/2
for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim)
  {
    Y1 <- rnorm(n, 0, 1)
    Y2 <- rnorm(n, delta, 1)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_norm$t_test_power[j] <- mean(p.t < 0.05)
  results_norm$wmw_test_power[j] <- mean(p.wmw < 0.05)
  results_norm$med_test_power[j] <- mean(p.med < 0.05)
}

#----------------------------------------------------------------------------------------------

# Sample simulation study of Q2 for uniform distribution

results_unif <- data.frame(
  sample_size = sample_sizes,
  t_test_power = numeric(n_sample_sizes),
  wmw_test_power = numeric(n_sample_sizes),
  med_test_power = numeric(n_sample_sizes)
)
#variance = 1/12 * (max-min)**2
delta <- sqrt(1/12)/2 #delta = sqrt(variance)/2
for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
    
  for(i in 1:nsim)
  {
    Y1 <- runif(n, min=0, max=1)
    Y2 <- runif(n, min=0, max=1)+delta
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_unif$t_test_power[j] <- mean(p.t < 0.05)
  results_unif$wmw_test_power[j] <- mean(p.wmw < 0.05)
  results_unif$med_test_power[j] <- mean(p.med < 0.05)
}

# Save files
save(results_unif, file = "results_unif.RData")
save(results_exp, file = "results_exp.RData")
save(results_t3, file = "results_t3.RData")
save(results_t5, file = "results_t5.RData")
save(results_lap, file = "results_lap.RData")
save(results_log, file = "results_log.RData")
save(results_norm, file = "results_norm.RData")
