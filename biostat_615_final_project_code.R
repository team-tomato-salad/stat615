library(scrime)
library(doParallel)
library(gglasso)
library(far)
registerDoParallel(cores=10)
library("glmnet")
set.seed(112620)

#This file contains our permutation method for tuning the penalty parameter in the group lasso, our modified method for tuning the penalty parameter in the regular lasso, and a simulation to assess the performance of our methods on data with a group structure. The simulation takes about 2 minutes to run and produces three visuals. 

################################################################################
#Extension to group lasso - novel work
################################################################################

#fit_group_lasso_FDR_control fits a group lasso model using a lambda value chosen to control the type-I error rate under the null (see details of helper functions)

#Input:
#X: SNP Matrix
#Y: Response
#Group: Group membership of SNPs
#Repeats: Number of permutation used to simulate draws from the null
#Alpha: Desired level of test

#Output:
#Group Lasso coefficients

fit_group_lasso_FDR_control <-function(X,Y,groups,repeats=10,alpha=.01){
  perm_group_lasso <- gglasso(X,Y,groups,lambda=select_lambda_glasso(X,Y,groups,alpha,repeats))
  beta <- perm_group_lasso$beta
}

#select_lambda_glasso uses a permutation test to identify the value of lambda which guarantees no more than alpha*p false discoveries under the null
#Y is permuted a given number of times
#For each permutation, bisection_search_group is called to identify the lambda value that yields alpha*p false discoveries
#Parallel processing is used to simultaneously process each permutation
#The median of the permutation distribution is returned

#Input:
#X: SNP matrix
#Y: Response
#Groups: Grouping of SNPs
#Alpha: Desired Level of Test
#Repeats: Number of Permutations

#Output:
#FDR-Controlling Lambda

select_lambda_glasso <- function(X,Y,groups,alpha=0.01,repeats=50){
  n <- length(Y)
  p <- dim(X)[2]
  nonzero <- round(alpha*p)
  max_group_size <- max(table(groups))
  ind=sample(1:n,n*repeats,replace=T)
  y.perm=matrix(Y[ind],n,repeats)
  perm_distribution <- foreach(i=1:repeats, .combine=c) %dopar% {
    lambda <- bisection_search_group(X,y.perm[,i],groups,nonzero,max_group_size)
    lambda
  }
  return(median(perm_distribution))
}

#Bisect_search_group employs a bisection search to efficiently find the lambda value that yields a desired number of nonzero group lasso coefficients 
#If it is impossible to have exactly the desired number of nonzero coefficients, then bisect_search_group will return the smallest lambda that yields no more than the desired number of nonzero coefficients under the null

#Input:
#X: SNP matrix
#Y: Response
#Nonzero: Desired number of nonzero coefficients
#Upper_bound: Used to set upper bound for lambda path (if zero, then 'nonzero' is the upper bound)
#Tolerance: Used to check for convergence in the event that it is impossible to have exactly the desired number of nonzero coefficients
#Maxiter: Maximum number of bisections

#Output:
#Lambda value

bisection_search_group=function(X,Y,groups,nonzero=1,upper_bound=1,tol=.00001,maxiter=10){
  #Initial lambda path
  fit=gglasso(X,Y,groups,dfmax=nonzero+2*upper_bound)
  lambda_path=cbind(fit$df,fit$lambda)
  lambda_min = lambda_path[,2][min(which(lambda_path[,1]>=nonzero))]
  #Check for a solution on the initial lambda path
  if(nonzero%in%lambda_path[,1]){
    return(lambda_min)
  }
  lambda_max =lambda_path[,2][max(which(lambda_path[,1]>=0& lambda_path[,1]<nonzero))]
  
  #Set bounds of interval
  a= lambda_min
  b= lambda_max
  c=NA
  
  fitL=gglasso(X,Y,groups, lambda=a)
  fitU=gglasso(X,Y,groups, lambda=b)
  
  for(i in 1:maxiter){
  
  #check for convergence 
  if(fitL$df==fitU$df){
    return(fitU$lambda)
  }
  if(abs(a-b)<tol){
    return(fitU$lambda)
  }
  
  #Otherwise, bisect interval
  if(fitU$df<nonzero&nonzero<fitL$df){
    c= 0.5*(a+b)
    fitM=gglasso(X,Y,groups,lambda=c)
    if(fitM$df==nonzero){
      return(fitM$lambda)
    } else if(fitM$df<nonzero){
      b=c
    } else if(fitM$df>nonzero){
      a=c
    }
  }
    
  fitL=gglasso(X,Y,groups, lambda=a)
  fitU=gglasso(X,Y,groups, lambda=b)
  
  #check for an exact solution 
  if(fitL$df==nonzero){
    return(fitL$lambda)
  }
  if(fitU$df==nonzero){
    return(fitU$lambda)
  }
  
  }
  #Return best value after maxiter
  lambda=fitU$lambda
  return(lambda)
}

################################################################################
#Modification of original permutation method
################################################################################
#Introduces sparse matrices
#Implements parallel processing for permutation tests
#Offers a new approximate method for finding lambda that avoids the bisection search while still guaranteeing a type-I error rate no greater than alpha (the desired level)

lassoPL_sparse_parallel = function(Xmat, Y,lambda=NULL,repeats=10,alpha=0.05,approx=T){
  start=Sys.time()
  Xmat <- as(Xmat,"sparseMatrix")
  Y=as.vector(Y)
  permLams=NULL
  
  if(!is.null(lambda)){
    model = glmnet(Xmat, Y,  alpha=1,lambda = lambda)
    lambda.final=lambda
  }
  
  if(is.null(lambda)){
    temp=permlam_par(Xmat,Y,repeats=repeats,alpha=alpha,approx)
    lambda.final=median(temp[,1],na.rm=T)
    permLams=temp
    model = glmnet(Xmat, Y,alpha=1,lambda = lambda.final)
  }
  
  beta = coef(model)
  beta=as.vector(beta)
  names(beta)=c("Intercept",colnames(Xmat))
  ending=Sys.time()
  time= ending-start
  #print(time)
  
  significant_predictors= beta[-1][which(beta[-1]!=0)]
  
  stats = list(betahat = beta,significant_predictors=significant_predictors ,DF=length(which(model$beta!=0)),
               lambda.final=lambda.final,permLams=permLams,time=time)
  
  return(stats)
}

permlam_par=function(Xmat,Y,alpha=1/ncol(Xmat),repeats=3,approx=T){
  registerDoParallel(cores=5)
  ##### alpha >= 1/p
  if(alpha>=1/ncol(Xmat)){
    nonzero= round(alpha*ncol(Xmat))
    lams=matrix(ncol=2)
    #count=0
    #start=Sys.time()
    while(sum(!lams[,2]%in% c(nonzero))>0){
      lams <- foreach(i=1:repeats, .combine=rbind) %dopar% {
        ind=sample(1:length(Y),length(Y))
        y.perm=Y[ind]
        if(approx==T){
          a=bisect_efficient(Xmat,y.perm,nonzero=nonzero) 
        }else{
          a=bisect(Xmat,y.perm,nonzero=nonzero)
        }
        a
      }
    }
    colnames(lams)=c("lambda","perm DF")
    return(lams)
  }
  
  #### alpha < 1/p
  if(alpha< 1/ncol(Xmat)){
    if(as.integer(round(1/alpha,0))%%ncol(Xmat)!=0) {
      stop("If alpha < 1/p, then this program can only handle alpha values of the following form: alpha= 1/(p*k) for some integer k, and p=ncol(X)")
    }
    nonzero=1
    reps2= 1/(alpha*ncol(Xmat))
    maxlams=vector()
    for(w in 1:repeats){
      lams=vector()
      count=0
      start=Sys.time()
      while(count<reps2){
        ind=sample(1:length(Y),length(Y))
        y.perm=Y[ind]
        a=bisect(Xmat,y.perm,nonzero=nonzero)
        if(a[2]==nonzero){
          count=count+1
          lams=rbind(lams,a)
        }
      }
      colnames(lams)=c("lambda","perm DF")
      end=Sys.time()
      maxlams=rbind(maxlams,c(max(lams[,1]),round(end-start,3),all(lams[,2]==1)))
    }
    colnames(maxlams)=c("lambda","time","DF Check")
    return(maxlams)
  }
}

bisect_efficient=function(Xmat,Y,nonzero=1){
  fit=glmnet(Xmat,Y,alpha=1,dfmax=nonzero*1.5)
  LDF=cbind(fit$df,fit$lambda)
  return(c(LDF[,2][max(which(LDF[,1]>=0& LDF[,1]<=nonzero))],nonzero))
  
}

is_TE = function(x){  #see if "object" has error or not
  inherits(x, "try-error")
}

bisect=function(Xmat,Y,nonzero=1){
  
  fit=glmnet(Xmat,Y,alpha=1,dfmax=nonzero+100)
  LDF=cbind(fit$df,fit$lambda)
  tmax=LDF[,2][max(which(LDF[,1]>=0& LDF[,1]<nonzero))]
  tmin= try(LDF[,2][min(which(LDF[,1]>=nonzero))])
  if(nonzero%in%LDF[,1]){
    return(c(tmin,nonzero))
  }
  #if(is_TE(tmin)){
  #print("error: unable to calculate tmin, try increasing dfmax or nlam?")
  #}
  if(tmin<=0 |is.na(tmin)){
    tmin=0.00000000000001
  }
  if(tmax>100){
    count3=0
    dfmax2=nonzero+100
    while(count3<1){
      #print(paste("warning: bad starting lambda values (too high), will try again"))
      dfmax2=dfmax2+5
      fit=glmnet(Xmat,Y,alpha=1,dfmax=dfmax2)
      LDF=cbind(fit$df,fit$lambda)
      #print(LDF)
      tmax=LDF[,2][max(which(LDF[,1]>=0& LDF[,1]<nonzero))]
      tmin= try(LDF[,2][min(which(LDF[,1]>=nonzero))])
      #if(is_TE(tmin)){
      #print("error: unable to calculate tmin, try increasing dfmax")
      #}
      if(tmin<=0 |is.na(tmin)){
        tmin=0.00000000000001
      }
      if(tmax<10){
        count3=1
      }
    }
  }
  
  tL=tmin
  tU=tmax
  tm=NA
  tm2=NA
  finalDF=NA
  lam=NA
  mdf=NA
  
  count=0
  count2=0
  while(count<1){
    if(!is.na(tm)){tm2=tm}
    fitL=glmnet(Xmat,Y,lambda=tL)
    fitU=glmnet(Xmat,Y,lambda=tU)
    if(nonzero>fitL$df){#print("problem, all possible lam values result in < desired nonzero coeffs")
      tL=0.00000000000001
    }
    
    if(nonzero<fitU$df){#print("problem, all possible lam values result in > desired nonzero coeffs")
      tU= tU+quantile(abs(Y),0.1)
    }
    if(fitL$df==nonzero){lambda=fitL$lambda
    count=1
    finalDF=nonzero
    }
    if(fitU$df==nonzero){lambda=fitU$lambda
    count=1
    finalDF=nonzero
    }
    if(fitU$df<nonzero&nonzero<fitL$df){tm= 0.5*(tL+tU)
    fitM=glmnet(Xmat,Y,lambda=tm)
    lam=fitM$lambda
    mdf=fitM$df
    if(fitM$df==nonzero){lambda=fitM$lambda
    count=count+1
    finalDF=fitM$df
    } else if(fitM$df<nonzero){
      tU=tm
    } else if(fitM$df>nonzero){
      tL=tm
    }
    }
    
    if(!is.na(tm)& !is.na(tm2) & tm==tm2){
      #count2= count2+1
      count2=5
    }
    else if(!is.na(tm)& !is.na(tm2) & tm!=tm2){
      count2=0
    }
    if(count2==5){
      count=1
      lambda=fitM$lambda
      finalDF=fitM$df
      print("tm has not changed for 2 iterations, will try again")
    }
  }
  return(c(lambda,finalDF))
}


################################################################################
#Simulation
#150 groups of 3, 150 groups of 2, 150 groups of a single SNP
#1 group of 3 and 1 group of 2 are casual
#We test the precision and recall of our methods as well as 2 methods for selecting lambda through cross-validation
#We also compare the type-I error rate. In this simulation, we suppose that we want to make no more than 2 false rejections
#The code produces three visuals
#The first compares the average type-I error rate, the second compares the average precision, and the last compares the average recall across 10 simulations with the same underlying truth
################################################################################ 

perm_group_lasso_precision <- c()
perm_group_lasso_recall <- c()
perm_lasso_precision <- c()
perm_lasso_recall <- c()
cv_precision <- c()
cv_recall <- c()
cvlse_precision <- c()
cvlse_recall <- c()
perm_group_lasso_type_1_err <- c()
perm_lasso_type_1_err <- c()
cv_type_1_err <- c()
cvlse_type_1_err <- c()

#Repeat simulation 10 times
for(i in 1:10){
#Generate simulation SNP data
lisprecision <- list(c(2, 0, 1), c(0, 2))
list2 <- list(c(1, 0, 1), c(1, 0))
sim <- simulateSNPs(600, 903, c(3, 2), list.ia.val = lisprecision,
                       list.equal = list2, vec.ia.num = c(50, 75), maf = 0.25)

#Training and validation sets for CV methods generated from same underlying truths
train <- simulateSNPs(400, 903, c(3, 2), list.ia.val = lisprecision,
                              list.equal = list2, vec.ia.num = c(round(.67*50), round(.67*75)), maf = 0.25)
val <- simulateSNPs(200, 903, c(3, 2), list.ia.val = lisprecision,
                      list.equal = list2, vec.ia.num = c(round(.33*50), round(.33*75)), maf = 0.25)

#Set group membership
#Causal groups are 1 and 2
groups <- c(rep(1,3),rep(2,2),rep(seq(3,151),3),rep(seq(152,301),2),seq(301,451))

X <- sim$data
Y <- sim$cl

X_train <- train$data
Y_train <- train$cl
X_val <- val$data
Y_val <- val$cl

#Test performance of our lambda-selection method for group lasso
#Permit no more than 2 false rejection
beta <- fit_group_lasso_FDR_control(X,Y,groups,repeats=10,alpha=2/903)
perm_group_lasso_precision <- c(perm_group_lasso_precision,sum(beta[1:5]!=0)/max(sum(beta!=0),1))
perm_group_lasso_recall <- c(perm_group_lasso_recall,sum(beta[1:5]!=0)/5)
perm_group_lasso_type_1_err <- c(perm_group_lasso_type_1_err,sum(beta[6:903]!=0))

#Test the modified permutation method for selecting lambda without accounting for groups
perm_lasso <- lassoPL_sparse_parallel(Xmat=X,Y=Y,repeats=10,alpha=2/903)
perm_lasso_precision <- c(perm_lasso_precision,sum(names(perm_lasso$significant_predictors)%in%c("SNP1","SNP2","SNP3","SNP4","SNP5"))/max(length(names(perm_lasso$significant_predictors)),1))
perm_lasso_type_1_err <- c(perm_lasso_type_1_err,sum(!names(perm_lasso$significant_predictors)%in%c("SNP1","SNP2","SNP3","SNP4","SNP5")))
perm_lasso_recall <- c(perm_lasso_recall,sum(names(perm_lasso$significant_predictors)%in%c("SNP1","SNP2","SNP3","SNP4","SNP5"))/5)

#Compare to cross-validation methods for choosing lambda

#Choose lambda for group lasso to minimize CV MSE
cv_lasso <- gglasso(X_train,Y_train,groups,lambda=cv.gglasso(X_val,Y_val,groups)$lambda.min)$beta
cv_precision <- c(cv_precision,sum(cv_lasso[1:5]!=0)/max(sum(cv_lasso!=0),1))
cv_type_1_err <- c(cv_type_1_err,sum(cv_lasso[6:903]!=0))
cv_recall <- c(cv_recall,sum(cv_lasso[1:5]!=0)/5)

#Choose lambda for group lasso 1se away from cv.min
cv_lasso_1se <- gglasso(X_train,Y_train,groups,lambda=cv.gglasso(X_val,Y_val,groups)$lambda.1se)$beta
cvlse_precision <- c(cvlse_precision,sum(cv_lasso_1se[1:5]!=0)/max(sum(cv_lasso_1se!=0),1))
cvlse_type_1_err <- c(cvlse_type_1_err,sum(cv_lasso_1se[6:903]!=0))
cvlse_recall <- c(cvlse_recall,sum(cv_lasso_1se[1:5]!=0)/5)

}

#Data visualizations 

#Type-I Error Rate Control
df <- data.frame(method=c("Permutation G-Lasso","Permutation Lasso","cv.min gglasso","cv.1se gglasso"),type1err=c(mean(perm_group_lasso_type_1_err[1:10]),mean(perm_lasso_type_1_err[1:10]),mean(cv_type_1_err[1:10]),mean(cvlse_type_1_err[1:10])),type1err_sd=c(sd(perm_group_lasso_type_1_err[1:10]),sd(perm_lasso_type_1_err[1:10]),sd(cv_type_1_err[1:10]),sd(cvlse_type_1_err[1:10])),type2err=c(mean(perm_group_lasso_recall[1:10]),mean(perm_lasso_recall[1:10]),mean(cv_recall[1:10]),mean(cvlse_recall[1:10])),type2err_sd=c(sd(perm_group_lasso_recall[1:10]),sd(perm_lasso_recall[1:10]),sd(cv_recall[1:10]),sd(cvlse_recall[1:10])))
df$method <- factor(df$method,levels = c("Permutation G-Lasso", "Permutation Lasso", "cv.min gglasso", "cv.1se gglasso"))

t1<-ggplot(df,aes(x=method,y=type1err,fill=method))+ geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=type1err, ymax=type1err+type1err_sd), width=.2,position=position_dodge(.9))+ geom_hline(yintercept=2, linetype="dashed", color = "red", size=1.5)+ theme(legend.position='none',plot.title = element_text(hjust = 0.5))+ylab("Number of Type-I Errors")+xlab("Method")+ggtitle("Comparison of FDR Control")

t1

#Precision and Recall
df2 <- data.frame(method=c("Permutation G-Lasso","Permutation Lasso","cv.min gglasso","cv.1se gglasso"),type1err=c(mean(perm_group_lasso_precision[1:10]),mean(perm_lasso_precision[1:10]),mean(cv_precision[1:10]),mean(cvlse_precision[1:10])),type1err_sd=c(sd(perm_group_lasso_precision[1:10]),sd(perm_lasso_precision[1:10]),sd(cv_precision[1:10]),sd(cvlse_precision[1:10])),type2err=c(mean(perm_group_lasso_recall[1:10]),mean(perm_lasso_recall[1:10]),mean(cv_recall[1:10]),mean(cvlse_recall[1:10])),type2err_sd=c(sd(perm_group_lasso_recall[1:10]),sd(perm_lasso_recall[1:10]),sd(cv_recall[1:10]),sd(cvlse_recall[1:10])))
#,type2err=c((perm_group_lasso_recall[1:10]),(perm_lasso_recall[1:10]),(cv_recall[1:10]),(cvlse_recall[1:10]) )
df2$method <- factor(df$method,levels = c("Permutation G-Lasso", "Permutation Lasso", "cv.min gglasso", "cv.1se gglasso"))

precision <-ggplot(df2,aes(x=method,y=type1err,fill=method))+ geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=type1err, ymax=type1err+type1err_sd), width=.2,position=position_dodge(.9))+ theme(legend.position='none',plot.title = element_text(hjust = 0.5))+ylab("")+xlab("Method")+ggtitle("Average Precision")

precision

recall<-ggplot(df2,aes(x=method,y=type2err,fill=method))+ geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=type2err, ymax=type2err+type2err_sd), width=.2,position=position_dodge(.9))+ theme(legend.position='none',plot.title = element_text(hjust = 0.5))+ylab("")+xlab("Method")+ggtitle("Average Recall")

recall
