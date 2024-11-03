##############################################################################
library(coin) #package needed for question two 
library(extraDistr)
library(ggplot2)
### Include approximation to Question 1

median.test.two_sided <- function(vec1, vec2, N=1000, exact=FALSE) {
  # Calculate observed test statistic
  observed_test_statistic <- abs(median(vec1) - median(vec2))
  
  # Combine both vectors
  full_vec <- c(vec1, vec2)
  
  # Check if the number of combinations is less than 10,000
  if (choose(length(full_vec), length(vec1)) < N | exact) {
    
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

######################################################################
# Question 2
nsim <- 500
output <- list()

#t3 distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
n <- 20
delta <- sqrt(3)/2
for(i in 1:nsim)
{
  Y1 <- rt(n, 3)
  Y2 <- rt(n, 3) + delta
  X <- factor(c(rep("A",n),rep("B",n)))
  Y <- c(Y1, Y2)
  #permutation t-test
  p.t[i] <- pvalue(oneway_test(Y ~ X,
                               distribution=approximate(nresample=1000)))
  p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
  p.med[i] <- median.test.two_sided(Y1, Y2)
}
# Monte-Carlo approximation of the power
# for alpha = .05
output <- append(output, list(t3_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

######################################################################
# Question 3
#Exponential distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
delta <- sqrt(1)/2
for (i in 1:nsim) {    # Or: replicate(N, { })
  Y1 <- rexp(n, rate = 1)
  Y2 <- rexp(n, rate = 1)+delta
  X <- factor(c(rep("A",n),rep("B",n)))
  Y <- c(Y1, Y2)
  
  p.t[i] <- pvalue(oneway_test(Y ~ X,
                               distribution=approximate(nresample=1000)))
  p.wmw[i] <- wilcox.test(Y1,Y2, exact=TRUE)$p.value
  p.med[i] <- median.test.two_sided(Y1, Y2)
}
mean(p.wmw<0.05)
output <- append(output, list(exp_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))


#Laplace distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
delta <- sqrt(2)/2
for(i in 1:nsim)
{
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
output <- append(output, list(laplace_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

#t5 distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
delta <- sqrt(5/3)/2
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
output <- append(output, list(t5_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

#Logistic distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
delta <- sqrt(pi**2/3)/2
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
output <- append(output, list(logistic_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

#normal distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
delta <- sqrt(1)/2
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
output <- append(output, list(normal_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

#uniform distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
delta <- sqrt(1/12)/2
for(i in 1:nsim)
{
  Y1 <- runif(n, 0, 1)
  Y2 <- runif(n, 0, 1)+delta
  X <- factor(c(rep("A",n),rep("B",n)))
  Y <- c(Y1, Y2)
  #permutation t-test
  p.t[i] <- pvalue(oneway_test(Y ~ X,
                               distribution=approximate(nresample=1000)))
  p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
  p.med[i] <- median.test.two_sided(Y1, Y2)
}
output <- append(output, list(uniform_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

output


# create a dataset
distributions <- rep(names(output), each=3)
tests <- rep(c("permutation t" , "wilcox" , "median") , 7)
power <- unlist(output)
data <- data.frame(distributions,tests,power)

# Grouped
ggplot(data, aes(fill=tests, y=power, x=distributions)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip()

######################################################################
# Question 4
#Type I error rate
nsim <- 500
output2 <- list()

#t3 distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
n <- 20
delta <- 0
for(i in 1:nsim)
{
  Y1 <- rt(n, 3)
  Y2 <- rt(n, 3) + delta
  X <- factor(c(rep("A",n),rep("B",n)))
  Y <- c(Y1, Y2)
  #permutation t-test
  p.t[i] <- pvalue(oneway_test(Y ~ X,
                               distribution=approximate(nresample=1000)))
  p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
  p.med[i] <- median.test.two_sided(Y1, Y2)
}
# for alpha = .05
output2 <- append(output2, list(t3_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))


#Exponential distribution
p.t <- p.wmw <- p.med <- numeric(nsim)
for (i in 1:nsim) {    # Or: replicate(N, { })
  Y1 <- rexp(n, rate = 1)
  Y2 <- rexp(n, rate = 1)+delta
  X <- factor(c(rep("A",n),rep("B",n)))
  Y <- c(Y1, Y2)
  
  p.t[i] <- pvalue(oneway_test(Y ~ X,
                               distribution=approximate(nresample=1000)))
  p.wmw[i] <- wilcox.test(Y1,Y2, exact=TRUE)$p.value
  p.med[i] <- median.test.two_sided(Y1, Y2)
}
mean(p.wmw<0.05)
output2 <- append(output2, list(exp_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))


#Laplace distribution
p.t <- p.wmw <- p.med <- numeric(nsim)

for(i in 1:nsim)
{
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
output2 <- append(output2, list(laplace_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

#t5 distribution
p.t <- p.wmw <- p.med <- numeric(nsim)

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
output2 <- append(output2, list(t5_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

#Logistic distribution
p.t <- p.wmw <- p.med <- numeric(nsim)

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
output2 <- append(output2, list(logistic_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

#normal distribution
p.t <- p.wmw <- p.med <- numeric(nsim)

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
output2 <- append(output2, list(normal_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

#uniform distribution
p.t <- p.wmw <- p.med <- numeric(nsim)

for(i in 1:nsim)
{
  Y1 <- runif(n, 0, 1)
  Y2 <- runif(n, 0, 1)+delta
  X <- factor(c(rep("A",n),rep("B",n)))
  Y <- c(Y1, Y2)
  #permutation t-test
  p.t[i] <- pvalue(oneway_test(Y ~ X,
                               distribution=approximate(nresample=1000)))
  p.wmw[i] <- wilcox.test(Y1,Y2)$p.value
  p.med[i] <- median.test.two_sided(Y1, Y2)
}
output2 <- append(output2, list(uniform_distribution=c(mean(p.t<0.05), mean(p.wmw<0.05), mean(p.med<0.05))))

output2
library(MKinfer)
perm.t.test(Y~X, R=1000)$p.value

# create a dataset
distributions <- rep(names(output2), each=3)
tests <- rep(c("permutation t" , "wilcox" , "median") , 7)
type1 <- unlist(output2)
data <- data.frame(distributions,tests,type1)

# Grouped
ggplot(data, aes(fill=tests, y=type1, x=distributions)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip()