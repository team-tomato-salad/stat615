### Readme for the final project

Sparse SNP group Lasso regression is an R based algorithm for SNP data analysis. 
We reimplement group Lasso method by taking advantage of sparse matrices and parallel processing to improve the efficiency, particularly for datasets with fewer than 100,000 SNPs.  This package successfully identify grouped SNPs for specific genetic related disease.

Installation and dependencies
	This package is implemented with R studio. R (>= 2.14.0)
	> llibrary(scrime)
	>library(doParallel)
	>library(gglasso)
	>library(far)
		registerDoParallel(cores=10)
	>library("glmnet")

Usage
	Sample code:
		set.seed(112620)
		perm_lasso <-lassoPL_sparse_parallel(Xmat=X,Y=Y,repeats=10,alpha=2/903). # X is the SNP matrix. # Y is permuted response 
		beta <- fit_group_lasso_FDR_control(X,Y,groups,repeats=10,alpha=2/903)
	Main functions:	
		fit_group_lasso_FDR_control
		select_lambda_glasso
		lassoPL_sparse_parallel
Processed data and metadata
	Raw data matrix will be sample by SNP sized matrix with 3 possible entries, 0, 1, and 2. During the process, the matrix will be converted to a sparse matrix for efficiency. 
	Sample data can be simulated with R package
Release Notes
	Version 1.0.0 - First available 
Authors:
	Wenjia Cao, Yizhou Qiu, Natasha Stewart 
		

