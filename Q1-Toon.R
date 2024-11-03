##############################################################################
# One-sided mean test
rm(list=ls())
mean.test.one_sided <- function(vec1, vec2) {
  # Calculate observed test statistic
  observed_test_statistic <- mean(vec1) - mean(vec2)
  
  # Combine both vectors
  full_vec <- c(vec1, vec2)
  
  # Generate all combinations of length(vec1) elements as a matrix
  comb1 <- combn(full_vec, length(vec1), simplify = TRUE)
  
  # Create a matrix where each column is the complement of each comb1 subset
  comb2 <- apply(comb1, 2, function(x) {
    setdiff(full_vec, x)
  })
  
  # Convert `comb2` to a matrix with length(vec2) rows for consistent structure
  comb2 <- matrix(unlist(comb2), nrow = length(vec2), byrow = FALSE)
  
  # Calculate the test statistics for each combination
  distribution = apply(comb1, 2 , mean) - apply(comb2, 2, mean)
  
  # Calculate the p-value
  p_value <- sum(distribution >= observed_test_statistic)/length(distribution)
}

# Example from the course 
vec1 <- c(37, 49, 55, 57)  
vec2 <- c(23, 31, 46)
outcome <- mean.test.one_sided(vec1, vec2) # Same as coursebook
outcome # Same as coursebook

##############################################################################
# Change to two-sided mean

#   to change to two-sided, we take the absolute value of the distribution,
#   by which we assume the null distribution is symmetric (and H0 : F1 = F2
#   can be shown to imply that the null distribtuion is symmetric..)


rm(list=ls())
mean.test.two_sided <- function(vec1, vec2) {
  # Calculate observed test statistic
  observed_test_statistic <- mean(vec1) - mean(vec2)
  
  # Combine both vectors
  full_vec <- c(vec1, vec2)
  
  # Generate all combinations of length(vec1) elements as a matrix
  comb1 <- combn(full_vec, length(vec1), simplify = TRUE)
  
  # Create a matrix where each column is the complement of each comb1 subset
  comb2 <- apply(comb1, 2, function(x) {
    setdiff(full_vec, x)
  })
  
  # Convert `comb2` to a matrix with length(vec2) rows for consistent structure
  comb2 <- matrix(unlist(comb2), nrow = length(vec2), byrow = FALSE)
  
  # Calculate the test statistics for each combination
  distribution = abs(apply(comb1, 2 , mean) - apply(comb2, 2, mean))
  
  # Calculate the p-value
  p_value <- sum(distribution >= observed_test_statistic)/length(distribution)
}

# Example from the course 
vec1 <- c(37, 49, 55, 57)  
vec2 <- c(23, 31, 46)
outcome <- mean.test.two_sided(vec1, vec2) 
outcome


# Change to two sided median

rm(list=ls())
median.test.two_sided <- function(vec1, vec2) {
  # Calculate observed test statistic
  observed_test_statistic <- median(vec1) - median(vec2)
  
  # Combine both vectors
  full_vec <- c(vec1, vec2)
  
  # Generate all combinations of length(vec1) elements as a matrix
  comb1 <- combn(full_vec, length(vec1), simplify = TRUE)
  
  # Create a matrix where each column is the complement of each comb1 subset
  comb2 <- apply(comb1, 2, function(x) {
    setdiff2(full_vec, x)
  })
  
  # Convert `comb2` to a matrix with length(vec2) rows for consistent structure
  comb2 <- matrix(unlist(comb2), nrow = length(vec2), byrow = FALSE)
  
  # Calculate the test statistics for each combination
  distribution = abs(apply(comb1, 2 , median) - apply(comb2, 2, median))
  
  # Calculate the p-value
  p_value <- sum(distribution >= observed_test_statistic)/length(distribution)
}

# Example from the course 
vec1 <- c(37, 46, 55, 57)  
vec2 <- c(23, 31, 46)

outcome <- median.test.two_sided(vec1, vec2) 
outcome


#Test code for the right setdiff function
g <- c(vec1, vec2)

#Right setdiff function when x contains y
setdiff2 <- function(x, y){
  tabx <- table(x)
  taby <- table(y)
  
  out <- numeric(length(tabx))
  names(out) <- names(tabx)
  
  for (i in names(tabx)){
    if (i %in% names(taby)){
      out[i] <- tabx[i]-taby[i]
    }
    else {
      out[i] <- tabx[i]
    }
  }
  rep(as.numeric(names(out)), as.vector(out))
}
setdiff(g, vec2)

