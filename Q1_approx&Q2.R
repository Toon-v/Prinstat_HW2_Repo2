##############################################################################
### Include approximation to Question 1

median.test.two_sided <- function(vec1, vec2) {
  # Calculate observed test statistic
  observed_test_statistic <- abs(median(vec1) - median(vec2))
  
  # Combine both vectors
  full_vec <- c(vec1, vec2)
  
  # Check if the number of combinations is less than 10,000
  if (choose(length(full_vec), length(vec1)) < 500) {
    
    #message("Using exact method with all combinations.") # deleted this because in the simulation, it would generate 1000 times
    
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
    N <- 500
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
  
  return(p_value)
}

######################################################################
# Question 2

library(coin)
p.t <- p.wmw <- p.med <- c()
n <- 20
delta <- sqrt(3)/2
for(i in 1:1000)
{
  Y1 <- rt(n, 3)
  Y2 <- rt(n, 3) + delta
  X <- factor(c(rep("A",n),rep("B",n)))
  Y <- c(Y1, Y2)
  #permutation t-test
  p.t[i] <- pvalue(oneway_test(Y ~ X,
                               distribution=approximate(B=9999)))
  p.wmw[i] <- wilcox.test(Y1,Y2,
                          exact = TRUE)$p.value
  p.med[i] <- median.test.two_sided(Y1, Y2)
}
# Monte-Carlo approximation of the power
# for alpha = .05
mean(p.t < .05)
mean(p.wmw < .05)
mean(p.med < .05)
