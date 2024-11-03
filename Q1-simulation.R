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

median.test_sim <- function(x, y, nsim = 1000){
  distribution <- numeric(nsim)
  observed_t <- median(x)-median(y)
  full_vec <- c(x, y)
  for (i in 1:nsim){
    comb1 <- sample(full_vec, length(x))
    comb2 <- setdiff2(full_vec, comb1)
    distribution[i] <- median(comb1)-median(comb2)
  }
  hist(distribution)
  p_value <- sum(distribution >= observed_t)/length(distribution)
  p_value
}

set.seed(100)
x <- rt(20, 3)
y <- rt(20, 3)+3
p <- median.test_sim(x, y)
