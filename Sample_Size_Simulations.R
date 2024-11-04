library(coin)
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

# Sample imulation study of Q2 for t3:


set.seed(100) #For reproductibility
nsim <- 500
#t3 distribution
delta <- sqrt(3)/2 #delta = sqrt(variance)/2
p.t <- p.wmw <- p.med <- numeric(nsim) #empty vector of p values for each test


sample_sizes <- c(5, 10, 20, 50, 100) # sample sizes for which we do the test. 
n_sample_sizes <- length(sample_sizes)


results_t3 <- data.frame(
  sample_size = sample_sizes,
  t_test_power = numeric(n_levels),
  wmw_test_power = numeric(n_levels),
  med_test_power = numeric(n_levels)
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
  

