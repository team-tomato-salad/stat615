#q3 

# args load

args = commandArgs(trailingOnly = TRUE)
alpha = as.numeric(args[1])
beta = as.numeric(args[2])


#current load
#alpha =2
#beta = 2
num_of_sam = 10e6


#sample from uniform
U = runif(num_of_sam)

#sample from beta
betas = rbeta(num_of_sam,alpha,beta)

#calc M (exp^(-x2))
M = 1

# compare and reject

dev = exp(-1*(betas^2))/U
acc = dev[dev >= M]

acept_rate =  length(acc)/length(dev)
l= quantile(betas[which(dev >= M)],c(0.01,0.25,0.5,0.75,0.99))


#prints

cat(paste(formatC(acept_rate,digits = 1,flag="-"),"\n",collapse=""))
cat(formatC(l, digits=1,flag =" -"))




