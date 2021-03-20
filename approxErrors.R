#hw3 q2

#load args

args = commandArgs(trailingOnly = TRUE)
load(args[1])
#load current
#load("~/Desktop/615/approxErrors.test2.RData")

len =length(X)
### find the true curve
coefs = coef(lm(Y ~ poly(X, 2,raw=T)))
quad = function(each){
  coefs[1]+coefs[2]*each+coefs[3]*(each^2)}
trueYs = sapply(Z,quad)

###

#Piecewise Linear Interpolation
piecewise.linear.interpolation  = function(x,y){
  n <- length (x)
  od_x = order(x)
  y <- y[od_x]
  x <- x[od_x]
  #coefs <- sapply(1:(n-1L),function(i) linear.interpolation(x[i], y[i], x[i + 1], y[i + 1]))
  m = diff(y)/diff(x)
  b = y[-1] - m*x[-1]
  return(cbind(x=x[-1],b,m))
}

eval.poly <- function(x,coefs){ y <- rep (0, length (x))
for(i in length(coefs):1L) y <- coefs [i] + x * y
return (y) }

eval.piecewise.poly <- function(x,coefs){
  n = nrow(coefs)
  x_bound = c(-Inf,coefs[-n,"x"],Inf)
  y = rep(NA,length=length(x))
  for(i in 1:n){
    idx = which((x <= x_bound[i+1]) & (x > x_bound[i]))
    y[idx] = eval.poly(x[idx],coefs[i,c("b","m")])
  }
  return(y)
}

Z_hat <- eval.piecewise.poly(Z,piecewise.linear.interpolation(X,Y))

#compute error
s = 0

for (loc in 1:length(Z)){
  
  diff = (Z_hat[loc]-trueYs[loc])^2
  s = s+diff
  
}

err = s/length(Z)

cat(formatC(err, digits=5,flag =" -"))
