

#_____________________________________________________________________________
# Code supplementary to HW2 Vancraeynest Vanhauwe

#_____________________________________________________________________________

# read library
library(coin)
library(extraDistr) #Laplace distribution
library(ggplot2) #Plotting barplots
library(tidyverse)

#_____________________________________________________________________________

# Question 1

median.test.two_sided <- function(vec1, vec2, N=1000, exact=FALSE) {
  
  observed_test_statistic <- abs(median(vec1) - median(vec2))
  
  
  full_vec <- c(vec1, vec2)
  
  # Check if the number of combinations is less than N - default N = 1000
  if (choose(length(full_vec), length(vec1)) < N | exact) {
    
    # Generate the null distribution by exact combinations
    t_star <- combn(length(full_vec), length(vec1), function(x) {
      group1 <- full_vec[x]
      group2 <- full_vec[-x]  # all elements not in x
      test_stats <- abs(median(group1) - median(group2))
      test_stats
    })
  } else {
    
    
    # Approximate null distribution by permutation if combinations exceed N
    Groupvec <- c(rep(1, length(vec1)), 
                  rep(2, length(vec2)))  # Create group vector
    
    perm_data <- data.frame(value = full_vec,
                            group = Groupvec)  # Create data frame for permutations
    t_star <- replicate(N, {
      # Shuffle the group labels:
      perm_data$group <- sample(perm_data$group)
      # Calculate the test statistic for shuffled groups:
      test_stats <- abs(median(perm_data$value[perm_data$group == 1])
                        - median(perm_data$value[perm_data$group == 2]))
      test_stats
    })
  }
  
  # Calculate the p-value
  p_value <- mean(t_star >= observed_test_statistic)
  
  p_value
}

#_____________________________________________________________________________

# Question 2

#----------------------------------------------------------------------------------

# set-up for Power simulation 

nsim <- 500
p.t <- p.wmw <- p.med <- numeric(nsim) #empty vector of p values for each test

sample_sizes <- c(5, 10, 20, 50, 100) # sample sizes for which we do the test. 
n_sample_sizes <- length(sample_sizes)

#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for t3:

set.seed(1234) #For reproductibility
delta <- sqrt(3)/2 #delta = sqrt(variance)/2

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



#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for Exponential distribution

set.seed(10234) #For reproductibility

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




#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for Laplace distribution

set.seed(100234) #For reproductibility

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

#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for t5 distribution

set.seed(1000234) #For reproductibility

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

#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for Logistic distribution

set.seed(10000234) #For reproductibility

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



#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for normal distribution

set.seed(100000234) #For reproductibility

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

#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for uniform distribution

set.seed(1000000234) #For reproductibility

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

# Add distribution column
results_unif$distribution <- "uniform"
results_exp$distribution <- "exponential"
results_t3$distribution <- "t with dof 3"
results_t5$distribution <- "t with dof 5"
results_lap$distribution <- "laplace"
results_log$distribution <- "logistic"
results_norm$distribution <- "normal"

# Combine all datasets into one
results_Power <- rbind(
  results_unif,
  results_exp,
  results_t3,
  results_t5,
  results_lap,
  results_log,
  results_norm
)

#_____________________________________________________________________________

# Question 3

Power_n20 <- results_Power %>% filter(sample_size == 20)

long_Power_n20 <- Power_n20 %>%
  pivot_longer(cols = c(t_test_power, wmw_test_power, med_test_power),
               names_to = "test_type",
               values_to = "power")

ggplot(long_Power_n20, aes(x = distribution, y = power, fill = test_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Power Comparison by Test Type (Sample Size = 20)",
       x = "Distribution",
       y = "Power",
       fill = "Test Type") +
  coord_flip() + 
  theme_minimal()

#_____________________________________________________________________________

# Question 4

# set-up for Type I error rate simulation

nsim <- 500
p.t <- p.wmw <- p.med <- numeric(nsim) #empty vector of p values for each test

sample_sizes <- c(5, 10, 20, 50, 100) # sample sizes for which we do the test. 
n_sample_sizes <- length(sample_sizes)

#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for t3:

set.seed(1234) #For reproductibility

results_t3 <- data.frame(
  sample_size = sample_sizes,
  t_test_TypeI = numeric(n_sample_sizes),
  wmw_test_TypeI = numeric(n_sample_sizes),
  med_test_TypeI = numeric(n_sample_sizes)
)

for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim)
  {
    Y1 <- rt(n, 3)
    Y2 <- rt(n, 3)
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
  results_t3$t_test_TypeI[j] <- mean(p.t < 0.05)
  results_t3$wmw_test_TypeI[j] <- mean(p.wmw < 0.05)
  results_t3$med_test_TypeI[j] <- mean(p.med < 0.05)
}



#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for Exponential distribution

set.seed(10234) #For reproductibility

results_exp <- data.frame(
  sample_size = sample_sizes,
  t_test_TypeI = numeric(n_sample_sizes),
  wmw_test_TypeI = numeric(n_sample_sizes),
  med_test_TypeI = numeric(n_sample_sizes)
)

for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for (i in 1:nsim) {
    Y1 <- rexp(n, rate = 1)
    Y2 <- rexp(n, rate = 1)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2, exact=TRUE)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_exp$t_test_TypeI[j] <- mean(p.t < 0.05)
  results_exp$wmw_test_TypeI[j] <- mean(p.wmw < 0.05)
  results_exp$med_test_TypeI[j] <- mean(p.med < 0.05)
}




#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for Laplace distribution

set.seed(100234) #For reproductibility

results_lap <- data.frame(
  sample_size = sample_sizes,
  t_test_TypeI = numeric(n_sample_sizes),
  wmw_test_TypeI = numeric(n_sample_sizes),
  med_test_TypeI = numeric(n_sample_sizes)
)


for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim){
    Y1 <- rlaplace(n, mu = 0, sigma = 1)
    Y2 <- rlaplace(n, mu = 0, sigma = 1)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_lap$t_test_TypeI[j] <- mean(p.t < 0.05)
  results_lap$wmw_test_TypeI[j] <- mean(p.wmw < 0.05)
  results_lap$med_test_TypeI[j] <- mean(p.med < 0.05)
}

#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for t5 distribution

set.seed(1000234) #For reproductibility

results_t5 <- data.frame(
  sample_size = sample_sizes,
  t_test_TypeI = numeric(n_sample_sizes),
  wmw_test_TypeI = numeric(n_sample_sizes),
  med_test_TypeI = numeric(n_sample_sizes)
)


for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim)
  {
    Y1 <- rt(n, 5)
    Y2 <- rt(n, 5)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_t5$t_test_TypeI[j] <- mean(p.t < 0.05)
  results_t5$wmw_test_TypeI[j] <- mean(p.wmw < 0.05)
  results_t5$med_test_TypeI[j] <- mean(p.med < 0.05)
}

#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for Logistic distribution

set.seed(10000234) #For reproductibility

results_log <- data.frame(
  sample_size = sample_sizes,
  t_test_TypeI = numeric(n_sample_sizes),
  wmw_test_TypeI = numeric(n_sample_sizes),
  med_test_TypeI = numeric(n_sample_sizes)
)

for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim)
  {
    Y1 <- rlogis(n, location=0, scale=1)
    Y2 <- rlogis(n, location=0, scale=1)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_log$t_test_TypeI[j] <- mean(p.t < 0.05)
  results_log$wmw_test_TypeI[j] <- mean(p.wmw < 0.05)
  results_log$med_test_TypeI[j] <- mean(p.med < 0.05)
}



#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for normal distribution

set.seed(100000234) #For reproductibility

results_norm <- data.frame(
  sample_size = sample_sizes,
  t_test_TypeI = numeric(n_sample_sizes),
  wmw_test_TypeI = numeric(n_sample_sizes),
  med_test_TypeI = numeric(n_sample_sizes)
)

for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  for(i in 1:nsim)
  {
    Y1 <- rnorm(n, 0, 1)
    Y2 <- rnorm(n, 0, 1)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_norm$t_test_TypeI[j] <- mean(p.t < 0.05)
  results_norm$wmw_test_TypeI[j] <- mean(p.wmw < 0.05)
  results_norm$med_test_TypeI[j] <- mean(p.med < 0.05)
}

#----------------------------------------------------------------------------------

# Sample simulation study of Q2 for uniform distribution

set.seed(1000000234) #For reproductibility

results_unif <- data.frame(
  sample_size = sample_sizes,
  t_test_TypeI = numeric(n_sample_sizes),
  wmw_test_TypeI = numeric(n_sample_sizes),
  med_test_TypeI = numeric(n_sample_sizes)
)

for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  
  for(i in 1:nsim)
  {
    Y1 <- runif(n, min=0, max=1)
    Y2 <- runif(n, min=0, max=1)
    X <- factor(c(rep("A",n),rep("B",n)))
    Y <- c(Y1, Y2)
    #permutation t-test
    p.t[i] <- pvalue(oneway_test(Y ~ X,
                                 distribution=approximate(nresample=1000)))
    p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
    p.med[i] <- median.test.two_sided(Y1, Y2)
  }
  results_unif$t_test_TypeI[j] <- mean(p.t < 0.05)
  results_unif$wmw_test_TypeI[j] <- mean(p.wmw < 0.05) 
  results_unif$med_test_TypeI[j] <- mean(p.med < 0.05)
}



# Add distribution variable
results_unif$distribution <- "uniform"
results_exp$distribution <- "exponential"
results_t3$distribution <- "t with dof 3"
results_t5$distribution <- "t with dof 5"
results_lap$distribution <- "laplace"
results_log$distribution <- "logistic"
results_norm$distribution <- "normal"

results_TypeI <- rbind(
  results_unif,
  results_exp,
  results_t3,
  results_t5,
  results_lap,
  results_log,
  results_norm
)

TypeI_n20 <- results_TypeI %>% filter(sample_size == 20)

long_TypeI_n20 <- TypeI_n20 %>%
  pivot_longer(cols = c(t_test_TypeI, wmw_test_TypeI, med_test_TypeI),
               names_to = "test_type",
               values_to = "Type_I")

ggplot(long_TypeI_n20, aes(x = distribution, y = Type_I, fill = test_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Type I error rate Comparison by Test Type (Sample Size = 20)",
       x = "Distribution",
       y = "Type I error rate",
       fill = "Test Type") +
  coord_flip() + 
  geom_hline(yintercept = 0.05, 
             color = "red", 
             linetype = "dashed", 
             linewidth = 1) +
  theme_minimal()



#______________________________________________________________________________

# Appendix 1 : QQ Plots


#QQPlot for t3 distribution

n <- 20
set.seed(1234) #For reproductibility
delta <- sqrt(3)/2 #delta = sqrt(variance)/2


Y1 <- rt(n, 3)
Y2 <- rt(n, 3) + delta

# Create Q-Q plot
qqplot(Y1, Y2, main = "Q-Q Plot of for t3 distribution",
       xlab = "Quantiles of group1", ylab = "Quantiles of group2",
       pch = 19, col = "blue")

# Add a diagonal line
abline(0, 1, col = "red", lty = 2, lwd = 2)


#----------------------------------------------------------------------------------

#QQPlot for Exponential distribution

set.seed(10234) #For reproductibility

#variance = 1/rate**2

delta <- sqrt(1)/2 #delta = sqrt(variance)/2

Y1 <- rexp(n, rate = 1)
Y2 <- rexp(n, rate = 1)+delta


# Create Q-Q plot
qqplot(Y1, Y2, main = "Q-Q Plot of for exponential distribution",
       xlab = "Quantiles of group1", ylab = "Quantiles of group2",
       pch = 19, col = "blue")

# Add a diagonal line
abline(0, 1, col = "red", lty = 2, lwd = 2)



#----------------------------------------------------------------------------------


#QQPlot for Laplace distribution

set.seed(100234) #For reproductibility

delta <- sqrt(2)/2 #delta = sqrt(variance)/2

Y1 <- rlaplace(n, mu = 0, sigma = 1)
Y2 <- rlaplace(n, mu = delta, sigma = 1)

# Create Q-Q plot
qqplot(Y1, Y2, main = "Q-Q Plot of for Laplace distribution",
       xlab = "Quantiles of group1", ylab = "Quantiles of group2",
       pch = 19, col = "blue")

# Add a diagonal line
abline(0, 1, col = "red", lty = 2, lwd = 2)

#----------------------------------------------------------------------------------

# QQPlot for t5 distribution

set.seed(1000234) #For reproductibility

delta <- sqrt(5/3)/2 #delta = sqrt(variance)/2


Y1 <- rt(n, 5)
Y2 <- rt(n, 5) + delta

# Create Q-Q plot
qqplot(Y1, Y2, main = "Q-Q Plot of for t5 distribution",
       xlab = "Quantiles of group1", ylab = "Quantiles of group2",
       pch = 19, col = "blue")

# Add a diagonal line
abline(0, 1, col = "red", lty = 2, lwd = 2)
#----------------------------------------------------------------------------------

# QQPlot for Logistic distribution

set.seed(10000234) #For reproductibility

#variance = scale**2 * pi**2 / 3 
delta <- sqrt(pi**2/3)/2 #delta = sqrt(variance)/2

Y1 <- rlogis(n, location=0, scale=1)
Y2 <- rlogis(n, location = delta, scale=1)

# Create Q-Q plot
qqplot(Y1, Y2, main = "Q-Q Plot of for logistic distribution",
       xlab = "Quantiles of group1", ylab = "Quantiles of group2",
       pch = 19, col = "blue")

# Add a diagonal line
abline(0, 1, col = "red", lty = 2, lwd = 2)

#----------------------------------------------------------------------------------

# QQPlot for normal distribution

set.seed(100000234) #For reproductibility


#variance = sqrt(sigma)
delta <- sqrt(1)/2 #delta = sqrt(variance)/2

Y1 <- rnorm(n, 0, 1)
Y2 <- rnorm(n, delta, 1)

# Create Q-Q plot
qqplot(Y1, Y2, main = "Q-Q Plot of for normal distribution",
       xlab = "Quantiles of group1", ylab = "Quantiles of group2",
       pch = 19, col = "blue")

# Add a diagonal line
abline(0, 1, col = "red", lty = 2, lwd = 2)

#----------------------------------------------------------------------------------

# QQPlot for uniform distribution

set.seed(1000000234) #For reproductibility

#variance = 1/12 * (max-min)**2
delta <- sqrt(1/12)/2 #delta = sqrt(variance)/2

Y1 <- runif(n, min=0, max=1)
Y2 <- runif(n, min=0, max=1)+delta

# Create Q-Q plot
qqplot(Y1, Y2, main = "Q-Q Plot of for uniform distribution",
       xlab = "Quantiles of group1", ylab = "Quantiles of group2",
       pch = 19, col = "blue")

# Add a diagonal line
abline(0, 1, col = "red", lty = 2, lwd = 2)

