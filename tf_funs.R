library(glmnet)
library(glmgen)
library(MASS)
library(remotes)
library(devtools)
library(glmgen)
library(cplm)
library(genlasso)
#library(trendfiltering)
library(pracma)
library(parallel)

ts_sd <- function(x) { sqrt(mean((x[-1] - x[-n])^2)/2)}
# From LASSO solution, identify where knots lie by evaluating discrete derivatives and finding which are non-zero. 
find_knots <- function(fit_y,x,k=1,tol=0.01){
  first_deriv = diff(fit_y)/diff(x)
  second_deriv = diff(first_deriv)/diff(x)[2:(length(x)-1)]
  third_deriv = diff(second_deriv)/diff(x)[3:(length(x)-1)]
  if(k==1) {
    knots = which(abs(second_deriv)>tol) +2
  }
  if(k==2){
    knots = which(abs(third_deriv)>tol) +3
  }  
  return(knots - 1)
}

BM_helper <- function(t,sigmasq, start=0 ) {
  ## first, simulate a set of random deviates
  x <- rnorm(n = length(t) - 1, sd = sqrt(sigmasq))
  ## now compute their cumulative sum
  x <- c(0, cumsum(x))
  x=x+start
  return(x)  
}


BM <- function(max_iter, step_size, sigmasq,niter=1000,start_vec=NULL){
  sig2 <- sigmasq*step_size
  t <- 0:(max_iter/step_size) *step_size
  if(is.null(start_vec)) {
    df <- data.frame(sapply(1:niter,function(x) BM_helper(t,sig2) ))
  }
  else if(length(start_vec) == 1){
    df <- data.frame(sapply(1:niter,function(x) BM_helper(t,sig2,start=start_vec) ))  
  }
  else{
    niter = length(start_vec)
    df <- data.frame(sapply(1:niter,function(x) BM_helper(t,sig2,start=start_vec[x]) ))  
  }
  
  
  rownames(df) <- t  
  return(df)
}

# Construct time series trend to evaluate trend filtering on, as a function of number of knots, slope, and noise. 
time_seq = function(n, slope, p, sigma){
  x = rep(0, n); v = rep(0, n)
  x[1] = 0; v[1] = 0; knots = c()
  #runif(1, min = -0.5, max = 0.5)
  for(i in 2:n){
    ind = rbinom(1, 1, 1 - p)
    v[i] = v[i-1] + ind*(runif(1, min = -slope, max = slope) - v[i-1])
    x[i] = x[i-1] + v[i-1]
    if (ind) {knots = c(knots, i)}
  }
  y = x + rnorm(n, sd = sigma)
  return(list(EY = x, Y = y, knots = knots,t = 1:n))
}


generate_linear <- function(n,p,K,signal,sd,rho){
  Sigma = toeplitz(rho^(0:(p-1)))
  X = mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
  beta <- rep(0,p)
  beta[0] <- signal
  beta[1] <- -signal
  EY = X %*% beta
  Y = rnorm(n,0,sd) +  EY
  return(list(X=X,Y=Y,EY=EY,beta=beta,signal=signal,sd=sd,rho=rho))
}

tf_fit <- function(sel_Y,inf_Y,t,EY,onesd_rule = TRUE,type="SURE",deg=1 ){
  if(type=="SURE"){
    tf = sure_trendfilter(1:n,sel_Y, weights =rep(1,n),k=deg)
    lambda_min = tf$lambda_min
    lambda_1se= tf$lambda_1se
    
    
    if(onesd_rule) {
      fit_y = tf$fitted_values[,which(tf$lambda == lambda_1se)]
    }else{
      fit_y = tf$fitted_values[,which(tf$lambda == lambda_min)]
    }
    
  } else{
    
    select_model = genlasso::trendfilter(sel_Y)
    cv = cv.trendfilter(select_model)
    
    lambda_min = as.character(round(cv$lambda.min,3))
    lambda_1se = as.character(round(cv$lambda.1se,3))
    
    
    if(onesd_rule) {
      fit_y = select_model$fit[,lambda_min]
    }else{
      fit_y = select_model$fit[,lambda_1se]
    }
    select_model$fit[,lambda_min]
    select_model$fit[,lambda_1se]
  }
  
  knots = find_knots(fit_y,t,k=deg)
  k= length(knots)+deg
  basis = tp(t, knots = knots, degree=deg,k=k)
  
  bs_X = cbind(rep(1, n), basis$X, basis$Z)
  n_x = ncol(bs_X) 
  
  reg_inf <- lm(inf_Y~ 0 + bs_X)
  reg_sel <- lm(sel_Y~0 + bs_X)
  reg_truth <- lm(EY~0 + bs_X)
  return(list(reg_inf=reg_inf,reg_sel=reg_sel,true_curve=predict(reg_truth),basis=bs_X,knots =knots))
}





fission_conditional_tf <- function(f_y,Y,t,BM_t,EY,gamma,sigmasq,onesd_rule=TRUE,type="SURE",deg=1,alpha=0.1){
  
  tfs <- tf_fit(f_y,Y,t,EY,onesd_rule = onesd_rule,type=type,deg=deg)
  
  reg_inf = tfs$reg_inf
  reg_sel = tfs$reg_sel
  X = tfs$basis
  
  var_y <- ts_sd(Y)**2 
  var_f <- ts_sd(f_y)**2 
  
  sigmasq_y <- var_y*t/(var_y - BM_t)
  sigmasq_f <- var_f - BM_t
  gamma_est <- sigmasq_f / var_f
  
  
  inv_gram <- pinv(t(X) %*% X)
  oracle_reg <- (coef(reg_inf) - gamma*coef(reg_sel))/(1-gamma)
  oracle_cov <- inv_gram *sigmasq/(1-gamma)
  reg_high <- (coef(reg_inf) - gamma_est*coef(reg_sel))/(1-gamma_est)
  reg_low <- coef(reg_inf) 
  reg_covar <- inv_gram*sigmasq_f/(1-gamma_est)
  
  oracle_pred <- get_conf_ols(oracle_reg,X,oracle_cov,alpha=alpha)
  high_pred <- get_conf_ols(reg_high,X,reg_covar,alpha=alpha)
  low_pred <- get_conf_ols(reg_low,X,reg_covar,alpha=alpha)
  
  emp_low <- apply(cbind(low_pred$pred_low,high_pred$pred_low),1,min)
  emp_high <- apply(cbind(low_pred$pred_high,high_pred$pred_high),1,max)
  
  
  
  
  
  CIs <- cbind(tfs$true_curve,oracle_pred$pred_low,oracle_pred$pred_high, emp_low, emp_high)
  colnames(CIs) <- c("Truth","Lower_Oracle","Upper_Oracle","Lower_Empirical","Upper_Empirical")
  avg_CI_oracle <- mean(CIs[,"Upper_Oracle"] - CIs[,"Lower_Oracle"])
  avg_CI_empirical <- mean(CIs[,"Upper_Empirical"] - CIs[,"Lower_Empirical"])
  oracle_coverage <- mean((CIs[,"Truth"] > CIs[,"Lower_Oracle"])&(CIs[,"Truth"] < CIs[,"Upper_Oracle"]))
  emp_coverage <- mean((CIs[,"Truth"] > CIs[,"Lower_Empirical"])&(CIs[,"Truth"] < CIs[,"Upper_Empirical"]))
  metrics <- c(oracle_coverage,emp_coverage,avg_CI_oracle,avg_CI_empirical,length(tfs$knots),sigmasq_f,gamma_est)
  return(list(CIs=CIs,metrics=metrics,knots=tfs$knots))
}




get_conf_ols <- function(beta,X,covar,alpha=0.1){
  pred <- X %*% beta
  sd <- sqrt(diag(X %*% covar %*% t(X)))
  pred_low =  pred+ qnorm(alpha/2)*sd
  pred_high =  pred- qnorm(alpha/2)*sd
  return(list(pred=pred,sd=sd,pred_high=pred_high,pred_low=pred_low))
}


get_metrics <- function(tfs) {
  CI <- cbind(tfs$true_curve,predict(tfs$reg_inf,level=1-alpha,interval="confidence"))
  coverage = mean((CI[,1] > CI[,3]) & (CI[,1] < CI[,4]))
  length = mean(CI[,4] - CI[,3])
  num_knots = length(tfs$knots)
  metrics = c(coverage,length,num_knots)
  colnames(CI) <- c("truth","fit","lwr","upper")
  return(list(CI=CI,metrics=metrics))
}


singlerun_tf <-function(iter,n,slope,prob_knot,sd,type="SURE",onesd_rule=TRUE,deg=1,alpha = 0.1,step_size = 0.05,max_iter=2){
  print(iter)
  sigmasq = sd**2
  
  dat = time_seq(n,slope,prob_knot,sd)
  tfs <- tf_fit(dat$Y,dat$Y,dat$t,dat$EY,onesd_rule = onesd_rule,type=type,deg=deg)
  naive_run <- get_metrics(tfs)
  Z <- BM(max_iter,step_size,sigmasq,niter=n,start_vec = 0)
  
  CIs = list()
  for(wh in 2:nrow(Z)) {
    BM_t = as.numeric(rownames(Z)[wh])
    z <- as.vector(data.matrix(Z)[wh,])
    
    gamma <- sigmasq /(sigmasq+ BM_t)
    sigma_hat <- ts_sd(dat$Y)
    f_y <- dat$Y + z
    g_y <- dat$Y - (sigmasq/BM_t)*z
    g_y_est <- dat$Y - (sigma_hat**2/BM_t)*z
    
    
    tf_cond <- fission_conditional_tf(f_y,dat$Y,dat$t,BM_t,dat$EY,gamma,sigmasq,type=type)
    tfs_fission_est <- tf_fit(f_y,g_y_est,dat$t,dat$EY,onesd_rule = onesd_rule,type=type,deg=deg)
    tfs_fission <- tf_fit(f_y,g_y,dat$t,dat$EY,onesd_rule = onesd_rule,type=type,deg=deg)
    
    est_run <- get_metrics(tfs_fission_est)
    fission_run <- get_metrics(tfs_fission)
    metrics = c(sigmasq,BM_t,gamma,length(dat$knots),sigma_hat,fission_run$metrics,est_run$metrics,tf_cond$metrics)
    CIs[BM_t] = list(CI_conditional = tf_cond$CIs, CI_fission = fission_run$CI, CI_fission_est = est_run$CI)
    if(wh == 2){
      results = metrics
    }
    else {
      results = rbind(results,metrics)    
    }
  }
  
  results <- data.frame(results)
  
  colnames(results) <- c("sigmasq","t","gamma","num_knots","hat_sigma","cover_oracle","length_fission",
                         "nk_fission","cover_est","length_est","nk_fissionest","cover_cond","cover_condest","length_cond","length_condest","nk_cond","sigma_hat","gamma_est")
  return(results)
  
}


