#################################################################
# Job: Fit Exponentiated Weibull distribution for bathtub hazard
#
# Name: Reubenadat
# Date: June 03 2021
# 
##############################################################################
rm(list = ls())
libs = c("utils", "MASS", "copula","numDeriv","scatterplot3d", "reliaR","extRemes","Rmpfr")
lapply(libs, require, character.only = TRUE) #"rmultil" not available for 4.02

### Hazard function and associated plot ####
hexpweibull <- function(t,lambda, kappa, alpha, log=FALSE){
  if ((!is.numeric(t)) || (!is.numeric(lambda)) || (!is.numeric(kappa)) || (!is.numeric(alpha))){ 
    stop("non-numeric argument to mathematical function")}
  if( (length(lambda)!= 1) || (length(kappa)!= 1) || (length(alpha)!= 1) ){
    stop("Non-t parameters must be atomic")}
  if (  (lambda <= 0) || (kappa <= 0) || (alpha <= 0) ){ 
    stop("Invalid arguments.  qp,lp,kp ap must be > 0 ")}
  log.pdf <-  log(alpha) + (alpha-1)*pweibull(t,scale = lambda,shape=kappa,log.p=TRUE) + 
    dweibull(t,scale = lambda,shape=kappa,log=TRUE)
  cdf <- exp(alpha*pweibull(t,scale = lambda,shape = kappa,log.p=TRUE) )
  log.h <- log.pdf - log(1-cdf)
  ifelse(log, return(log.h), return(exp(log.h)))
}   

h <- Vectorize( function(t) hexpweibull(t,3,3,.1) )
curve(h,0,6,ylim=c(0,4),col= "black",lwd=2,n = 1000,
      xlab="t",ylab="Hazard function",cex.axis=1.5, cex.lab = 1.5, 
      main = expression(paste("Exponentiated Weibull  Hazard for ",
                              lambda, " = 3, ", kappa, " = 4, ", alpha, " = 0.1")) )
##############################################################################
##############################################################################
# Cumulative distribution function
pexpweibull<- function(q,lambda,kappa,alpha,log.p=FALSE){
  if ((!is.numeric(q)) || (!is.numeric(lambda))|| (!is.numeric(kappa)) || (!is.numeric(alpha))){ 
    stop("non-numeric argument to mathematical function")}
  if((length(lambda)!= 1) || (length(kappa)!=1) || (length(alpha)!=1)){
    stop("Non-q parameters must be atomic")}
  #if ((min(q) <= 0) || (lambda <= 0) || (kappa <= 0) || (alpha <= 0) ){ 
  # stop("Invalid arguments.  qp,lp,kp ap must be > 0 ")}
  log.cdf <- alpha*pweibull(q,scale = lambda,shape = kappa,log.p=TRUE)
  mpfr(ifelse(log.p, return(log.cdf), return(exp(log.cdf))),128L)
}  
################################################## 
# Probability density function
dexpweibull<- function(t,lambda,kappa,alpha,log=FALSE){
  if ((!is.numeric(t)) || (!is.numeric(lambda))|| (!is.numeric(kappa)) || (!is.numeric(alpha))){ 
    stop("non-numeric argument to mathematical function")}
  if((length(lambda)!= 1) || (length(kappa)!= 1) || (length(alpha)!= 1)){
    stop("Non-x parameters must be atomic")}
  #if ((min(t) <= 0) || (lambda <= 0) || (kappa <= 0) || (alpha <= 0)  ){ 
  # stop("Invalid arguments. td,ld,kd,ad must be > 0 ")}
  log.pdf <-  log(alpha) + (alpha - 1)*pweibull(t,scale = lambda, shape = kappa,log.p  = TRUE)  + 
    dweibull(t,scale = lambda,shape = kappa, log = TRUE)
  mpfr(return(exp(log.pdf)), 128L)
  #ifelse(log, return(log.pdf), return(exp(log.pdf)))
}
##################################################   
# Quantile function
qexpweibull<- function(p,lambda,kappa,alpha){
  quant <-  qweibull(p^(1/alpha),scale = lambda,shape = kappa)
  return(quant)
} 
##################################################   
# Random number generation
rexpweibull<- function(n,lambda,kappa,alpha,log=FALSE){
  if ((!is.numeric(n)) || (!is.numeric(lambda)) ||(!is.numeric(kappa)) || (!is.numeric(alpha))) 
    stop("non-numeric argument to mathematical function")
  if ((n <= 0) || (min(lambda) <= 0) || (min(kappa) <= 0) || (min(alpha) <= 0)) 
    stop("Invalid arguments")
  u = runif(n)
  sim <-  qweibull(u^(1/alpha),scale = lambda,shape=kappa)
  return(sim)
} 
############################################################
smart_df = function(...){
  data.frame(...,stringsAsFactors = FALSE)
}
smart_WT = function(...){
  write.table(...,row.names = FALSE,quote = FALSE)
}
smart_hist = function(xx,...){
  hist(xx,freq = FALSE,...); lines(density(xx),lwd = 2,lty = 2,col = "red")
}
####################### Event time #################################
# true parameters
theta = 4/2;alpha1 = 0.1;alpha2 = 2;lambda1 = 3;lambda2 = 2.1;kappa = 4
true_mean = c(lambda1, kappa,alpha1)

set.seed(223); n.obs = 5000; num_rep = 500
datagum = list()
for (jj in seq(num_rep)) {
  gumbel_cop = gumbelCopula(param = theta, dim = 2,use.indepC =c("message","TRUE","FALSE"))
  bvd = mvdc(copula = gumbel_cop, margins = c("expweibull","weibull"),paramMargins = list(list( lambda = lambda1, kappa = kappa, alpha = alpha1), 
                                                                                          list( shape = alpha2, scale = lambda2)),check = TRUE  )
  all_data = rMvdc(n.obs,bvd)
  datagum[[jj]] = all_data
}
########################################################################
## likelihood ##
ns_llike = function(data_curr,init_params) {
  
  params = init_params; lambda = params[1]; kappa = params[2];  alpha = params[3]; theta = theta
  
  tmppdf1 = dexpweibull(data_curr$X, lambda, kappa, alpha, log = FALSE)
  tmppdf2 = dweibull(data_curr$X, alpha2, lambda2, log = FALSE)
  
  tmpcdf1 = pexpweibull(data_curr$X, lambda, kappa, alpha, log.p = FALSE)
  tmpcdf2 = pweibull(data_curr$X, alpha2, lambda2, lower.tail = TRUE, log.p = FALSE)
  
  tmpu1 = -log( pmax(1e-6,tmpcdf1) )
  tmpu2 = -log( pmax(1e-6,tmpcdf2) )
  
  tmpFxy = exp(-( tmpu1^theta + tmpu2^theta )^(1/theta))
  tmpfxy = as.double( tmpFxy*( (tmpu1^theta + tmpu2^theta)^( (1/theta) - 1) )*
                        ( ( tmpu1^(theta - 1) )*( tmppdf1/tmpcdf1 ) + ( tmpu2^(theta - 1) )*(tmppdf2/tmpcdf2) ) )
  
  tmppdf = pmax(1e-6,tmppdf1 + tmppdf2 - tmpfxy)
  tmpSxy = pmax(1e-6,1 -  tmpcdf1 - tmpcdf2 +  tmpFxy)
  
  ll = sum( event*log(tmppdf) , (1 - event)*log(tmpSxy)  ) 
  ll
}

ns_ll = function(init_params){
  ns_llike(data_curr,init_params)
}
func_hess = function(X){+          tryCatch({
  hess1 = diag(ginv(X))  }, error=function(e){ warning(e)
    error.flag = TRUE })
  hess1
}
###########################optimization##########################
my_optim = function(data_curr) {
  out_opt = optim(init_params, ns_ll, control = list(fnscale = -1, maxit = 1e6), hessian = TRUE)
  hess = -out_opt$hessian; curr_var = func_hess(hess);   curr_params = as.numeric(out_opt$par)
  unc_params = rep(NA,7); unc_params =  c(curr_params,curr_var,propC, init_params )
  all_estimates = rbind(all_estimates, smart_df(S = what_rep, N = sample_size, lambda = unc_params[1],kappa = unc_params[2],alpha = unc_params[3], 
                                                lambda.var = unc_params[4], kappa.var = unc_params[5],alpha.var = unc_params[6], propC = unc_params[7]) )
  all_estimates
}
######################Execution###########################
all_estimates = c()
sample_sizes = 1000#c(1000, 2500,5000) 
for (sample_size in sample_sizes ) {
  cat(paste0("N = ",sample_size,"\n"))
  for (what_rep in seq(num_rep)) {
    if (what_rep %% 5 == 0) cat(paste0(what_rep," "))
    curr_data = datagum[[what_rep]][seq(sample_size),]
    data_cur  = cbind(T = pmin(curr_data[,1], curr_data[,2]), C = runif(sample_size, 0, 6) )
    
    data_curr = transform(data_cur, X = pmin(T, C) )
    event = smart_df(E = ifelse(data_curr$T < data_curr$C, 1, 0)  )#1==yes, 0==no
    propC = table(event)[1]/nrow(event)    
    
    init_params = exp(c(runif(2,-5,2), runif(1,-5,0)) ) #runif(2,1,1.2)#c(runif(1,0,1), runif(1,2,7)) #, runif(1,1/5,1/3) )#rep(1,2)#
    
    all_estimates = my_optim(data_curr)
  }
  cat("\n") #print each line
}
#####################################################################

all_estimates$CP_lambda = all_estimates$lambda - qnorm(0.975)* sqrt(all_estimates$lambda.var) < true_mean[1] & 
  true_mean[1] < all_estimates$lambda + qnorm(0.975) * sqrt(all_estimates$lambda.var)
all_estimates$CP_kappa = all_estimates$kappa - qnorm(0.975)* sqrt(all_estimates$kappa.var) < true_mean[2] & 
  true_mean[2] < all_estimates$kappa + qnorm(0.975) * sqrt(all_estimates$kappa.var)
all_estimates$CP_alpha = all_estimates$alpha - qnorm(0.975)* sqrt(all_estimates$alpha.var) < true_mean[3] & 
  true_mean[3] < all_estimates$alpha + qnorm(0.975) * sqrt(all_estimates$alpha.var)

table_func = function(xx){
  as.numeric(sum(xx)/nrow(all_estimates[which(all_estimates$N== N  ),] ) )
}

summary = c()
for (N in sample_sizes ) {
  tmp_mean = apply(all_estimates[which(all_estimates$N == N), c("lambda","kappa","alpha")],2,mean)
  tmp_df = smart_df(N,params = c("lambda","kappa","alpha"), mean.est = as.numeric(tmp_mean))
  tmp_df$bias = (tmp_df$mean.est) - true_mean
  tmp_df$Mod_var = apply(all_estimates[which(all_estimates$N == N), c("lambda.var","kappa.var","alpha.var")],2, mean)
  tmp_df$Emp_var = apply(all_estimates[which(all_estimates$N == N), c("lambda","kappa","alpha")],2,var)
  tmp_df$CP = apply(all_estimates[which(all_estimates$N == N), c("CP_lambda","CP_kappa","CP_alpha")],2,table_func)
  summary = rbind(summary,tmp_df)
}
summary[,c(3:7)] = lapply(summary[c(3:7)],round,3) ;summary 
summary(all_estimates["propC"])




