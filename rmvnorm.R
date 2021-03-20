#q4

# args load

args = commandArgs(trailingOnly = TRUE)
n = as.integer(args[1])
p = as.numeric(args[2])


#current load
#n =10
#p = 0.5
num_of_sam = 5e5

times = n*log(p)

### making Sigma/cov matrix
S = matrix(NA,ncol = n,nrow = n)
for (h in 1:n){
  for (j in 1:n){
    l1= (h-1)/n
    l2 = (j-1)/n
    S[h,j] = exp(times*(abs(l1-l2)^1.99)-abs(cos(l1))-abs((cos(l2))))
    
  }
}

#Sample from a multivariate normal distribution based on cholesky decomposition
rmvnorm_chol <- function(n,mu=rep(0,nrow(S)),Sigma= S){
  p = length(mu)
  Z <- matrix(rnorm(n*p),nrow=p,ncol=n)
  U <- chol(Sigma)
  X <- mu + crossprod(U, Z)
  return(t(X))
}

### est
x_chol <- rmvnorm_chol(1e5)
###

firstE = mean(apply(x_chol, 1, max))

###

secondE = mean(sqrt(apply(x_chol^2, 1, sum)))

###
dim1 = x_chol[,1]
dim2 = x_chol[,2]
dim12 = dim1*dim2
pr = length(dim12[dim12>0.5*p])/length(dim12)


# print output
output = c(firstE,secondE,pr)
cat(paste(formatC(output,digits = 1,flag="-"),"\n",collapse=""))