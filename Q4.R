# HW1 Q4

args = commandArgs(trailingOnly = TRUE)
xbinary = args[1]
ybinary = args[2]
lam = as.numeric(args[3])

#read data from test case 9



########

# read X matrix

xfile = file(xbinary, "rb")  # need to change to arg1
x_row = readBin(xfile, integer())
x_col = readBin(xfile, integer())

X = matrix(nrow = x_row,ncol = x_col)

for (col in 1:x_col) {
 
    X[,col] =  readBin(xfile,numeric(),x_row)
  
}
close(xfile)

#read y matrix
yfile = file(ybinary, "rb")  # need to change to arg1
y_row = readBin(yfile, integer())
y_col = readBin(yfile, integer())

Y = matrix(nrow = y_row,ncol = y_col)
for (row in 1:x_row){
  Y[row,] =  readBin(yfile,numeric())
}
close(yfile)


#calculate closed form

I = diag(x_col)
amI = lam*I  #lamI matrix

xt = t(X) #transpose X

beforeinv = (xt %*% X)+ amI

# invsvd = svd(beforeinv)
# U = invsvd$u
# diag = invsvd$d
# Vt = invsvd$v
# invdiag = 1/diag
# invD = invdiag*I
# 
# inved = t(Vt) %*% invD %*% t(U)
# 
# betas = inved %*% xt %*% Y
# 
# 
# inved = chol2inv(beforeinv,size = NCOL(beforeinv))
# betas = inved %*% xt %*% Y
ch = chol(beforeinv)
b =  xt %*% Y
z = forwardsolve(t(ch),b,upper.tri=FALSE)
betas = round(backsolve(ch,z,upper.tri=TRUE))
output = NULL 
index = NULL
for (i in 1:length(betas)){
  if (betas[i]!= 0 ){
    index= c(index,i)
    output =c(output,betas[i]) 
  }
}

#inved = chol2inv(ch)
#inved = solve(beforeinv)
#betas = inved %*% xt %*% Y
#b = round(betas)

#output = NULL 
#index = NULL
#for (i in 1:length(b)){
#  if (b[i]!= 0 ){
#    index= c(index,i)
#    output =c(output,b[i]) 
#  }
#}

cat(paste(index,output,"\n",collapse=""))
