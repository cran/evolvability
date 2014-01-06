randomBeta = function(n = 1, k = 2){
  X = matrix(rnorm(n*k), ncol = n)
  X = apply(X, 2, function(x) x/sqrt(sum(x^2)))
  rownames(X) = paste("dim", 1:k, sep="")
  X
}

