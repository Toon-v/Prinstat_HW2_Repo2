library(extraDistr) #Laplace distribution

#----------------------------------------------------------------------------------

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
    
