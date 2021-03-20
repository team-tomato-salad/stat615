#HW1 Q2
#read input values from the argument
args = commandArgs(trailingOnly = TRUE)


degree = as.integer(args[1])
ys = as.numeric(args[2:length(args)])
n = length(ys)

# test2 = read.table("Q2_testcase.txt.args",header = FALSE,fill = TRUE,sep = " ")
# degree = as.integer(test2[1])
# ys = as.numeric(test2[2:length(test2)])
# n = length(ys)


#make x matrix
originalx = 0

for (i in 1:n){
  originalx[i] =i/n 
  
}

processedx = 0
for (i in 1:n){
  eachProcessedvec= 0
  for (j in 1:degree){
    eachProcessedvec[j] = (1/factorial(j))*(originalx[i])^j
    
  }
  
  processedx[i] = sum(eachProcessedvec)

}

#Xmatrix = cbind(rep(1, n),processedx)

#run linear reg on x matrix for beta0 and beta 1


fastSimpleLinearRegression <- function(x, y) {
  y <- y - mean(y)
  x <- x - mean(x)
  n <- length(y)
  stopifnot(length(x) == n)        # for error handling
  s2y <- sum( y * y ) / ( n - 1 )  # \sigma_y^2
  s2x <- sum( x * x ) / ( n - 1 )  # \sigma_x^2
  sxy <- sum( x * y ) / ( n - 1 )  # \sigma_xy
  rxy <- sxy / sqrt( s2y * s2x )   # \rho_xy
  b <- rxy * sqrt( s2y / s2x )

  return(beta1 = b)
}


b1 = fastSimpleLinearRegression(processedx,ys)
b0 = mean(ys)-(b1*mean(processedx))
# finding following betas
followb = 0

for (i in 1:degree-1){
  followb[i] = b1/(factorial(i+1))
  
}

bs = c(b0,b1,followb)
#cat(bs,sep=" ",digits = 8)
cat(paste(formatC(bs,digits = 8,flag="-")," ",collapse=""))



