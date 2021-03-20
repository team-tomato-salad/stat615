#q1

### load from command

#args = commandArgs(trailingOnly = TRUE)
#test <- read.csv(args[1])
#X = test$X
#Y = test$Y
#Z = test$Z

### current load
start = Sys.time()
test1 <- read.csv("Desktop/615/hw4/logit_test1.csv")
X = test1$X
Y = test1$Y
Z = test1$Z
len = length(X)
#alpha = -1.7175

###
#function summing up deritives


base_func = function(alpha,x,z){
  return(alpha^2*z^2-alpha*x)
}

times_func = function(alpha,x,z){
  return(2*alpha*z^2-x)
}


long_func = function(alpha,x,y,z){
  v = base_func(alpha,x,z)
  b = times_func(alpha,x,z)
  expv = exp(v)
  expnv = exp(-v)
  return(-y*b*expv/(1+expv)+(y-1)*-b*expnv/(1+expnv))
}

sum_deritive= function(alpha,X,Y,Z){
  sums = 0
  for (i in 1:length(X)){
    xi = X[i]
    yi = Y[i]
    zi = Z[i]

    sums = sums + long_func(alpha,xi,yi,zi)

  }
  # each_part = sapply(1:length(X),function(i){
  #   
  #   long_func(alpha,X[i],Y[i],Z[i])
  # })
  
  return(sums)
}


##### chunk -5,5 into smaller pieces

original_maximum = function(base_func,alpha,X,Y,Z){
  sums = 0
  for (i in 1:length(X)){
    xi = X[i]
    yi = Y[i]
    zi = Z[i]
    bf = base_func(alpha,xi,zi)
    loglikeli = -yi*log(1+exp(bf))+(yi-1)*log(1+exp(-bf))
    sums = sums + loglikeli
    
  }
  return(sums)
}

#interval_score = list()
possible_pos = seq(-5.25, 5.25, by=0.5)

# for (pos in 1:length(possible_pos)){
#   alph = possible_pos[pos]
#   
#   interval_score = rbind(interval_score,original_maximum(base_func,alph,X,Y,Z))
# }


interval_score = sapply(possible_pos,function(alph){
  
  original_maximum(base_func,alph,X,Y,Z)
})


where = which.max(interval_score)
pos0 = possible_pos[where]-0.5
pos1 = possible_pos[where]+0.5

#0.25 increment will time out
######plug in 

secant <- function(f,X,Y,Z,x0,x1,tol=1e-10,max_iter=1000){
  convergence = 1
  f0 = f(x0,X,Y,Z); f1 = f(x1,X,Y,Z)
  if (abs(f0 -f1)<tol ){
    warning ( " Expect a huge jump ! " )
    break
  }
  x12 <- -f1/(f1-f0)*(x1-x0); x2 <- x1 + x12 
  for(iter in 1:max_iter){
    if (abs(x12)<tol ){ 
      convergence = 0
      break
    }
    f0 <- f1; x1 <- x2; f1 <- f(x2,X,Y,Z); f01 <- f1 - f0 
    if (abs(f01)<tol ){
      warning ( " Expect a huge jump ! " ) 
      break
    }
    x12 <- -f1/f01*x12
    x2 <- x1 + x12 }
  return(list(root=x2,f_root = f(x2,X,Y,Z),iter=iter ,convergence=convergence)) }



### output
a = secant(sum_deritive,X,Y,Z,pos0,pos1)

output = a$root

cat(formatC(output, digits=5,flag =" -"))


end = Sys.time()

end - start









