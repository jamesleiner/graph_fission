library(MASS)
library(matrixcalc)
library(pracma)
library(dplyr)
library(tidyr)
library(plotly)
library(genlasso)


##########################################################################################################
# Simple helper functions
##########################################################################################################


proj <- function(X,Y){
  X %*% pinv(t(X)%*%X) %*% t(X) %*% Y
}


proj_matrix <- function(X){
  X %*% pinv(t(X)%*% X ) %*% t(X)
}


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


get_active <- function(graph, row,x_sep=4,y_sep=5) {
  flag <- (graph[row[1],1] == x_sep)*(graph[row[2],1] == x_sep + 1)
  flag2 <- (graph[row[1],2] == y_sep)*(graph[row[2],2] == y_sep + 1)
  return(flag | flag2)
}



##########################################################################################################
# connected_component: Takes in a graph as a list of edges and identified all connected components
#
# Inputs: G - Graph represented as an edge list (vector)
#         start - Starting node for the connected component search (integer)
#
# Outputs: List of connected component
############################################################################################################

connected_component <- function(G, start) {
  last_component = c()
  component = start
  
  while(!identical(component,last_component)){
    last_component = component 
    component = unique(c(component,G[G[,1] %in% component,2]))
    component = unique(c(component,G[G[,2] %in% component,1]))
    component = component[order(component)]
  }
  return(component)
}

##########################################################################################################
# run_tf: Runs graph trend filtering given ege matrix and order
#
# Inputs: response - Response variable (numeric vector)
#         edge_mat - Graph Laplacian matrix (matrix)
#         k - Order of the penalty term (integer)
#
# Outputs: Fitted trend using genlasso
##########################################################################################################

run_tf <-function(response,edge_mat, k= 0) {
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0 ) {
    penalty = edge_mat %*% matrix.power(graph_laplacian,k/2)
  }
  if(k%%2 == 1){
    penalty = matrix.power(graph_laplacian,(k+1)/2)
  }
  h = genlasso(y=response,D=penalty)
  return(h)
  
}
##########################################################################################################
# elastic_loss: Calculates square loss with elastic net loss for fitted trend
#
# Inputs: beta - Fitted trend (numeric)
#         y - Response variable (numeric)
#         edge_mat - Graph Laplacian matrix (matrix)
#         k - Order of the penalty term (integer)
#         lambda_1 - L1 penalty parameter (numeric)
#         lambda_2 - L2 penalty parameter (numeric)
#
# Outputs: Loss 
##########################################################################################################

elastic_loss <- function(beta,y,edge_mat,k=0,lambda_1=1,lambda_2=1){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0 ) {
    penalty = edge_mat %*% matrix.power(graph_laplacian,k/2)
  }
  if(k%%2 == 1){
    penalty = matrix.power(graph_laplacian,(k+1)/2)
  }
  sum((beta - y)**2) + lambda_1 * sum(abs(penalty%*% beta)) + lambda_2 * sum((penalty%*% beta)**2 )
}




##########################################################################################################
# multi_split: Splits the input vector into K components using decomposition rules
#
# Inputs: x - Assumed Gaussian  (numeric)
#         K - Number of splitsinteger)
#         sd_level - Standard deviation of error term, assumed known (integer)
#
# Outputs: Matrix containing K vectors, one fore each split
##########################################################################################################


multi_split <- function(x,K,sd_level) {
  cov <- matrix(data = -1/(K**2),nrow=K,ncol=K)
  diag(cov) <- rep((1/K)*(1-1/K),K)
  cov <- cov*sd_level
  y <- mvrnorm(n=1, mu = rep(x/K,K),Sigma=cov)
  return(y)
}

##########################################################################################################
# generate_data: Generates synthetic data based on the given graph structure, edges, and parameters
#
# Inputs: graph - Graph structure (matrix)
#         edges - Edge in the form of a list with (node1, node2) entries (matrix)
#         edge_mat - Graph Laplacian matrix (matrix)
#         active_set - Active set of edges (vector)
#         k - Order of basis (integer)
#         sd_level - Standard deviation level for error term (numeric)
#         tau - SProportion of information to be generated for selection (numeric)
#         scale - Scaling of coefficients for generated data (numeric)
#         type - Type of distribution for generating response variable (character)
#         err_type - Type of distribution for generating noise (either normal, skew normal, t, or Laplacian)
#         est_var - Flag to estimate noise variance from the data (logical)
#
# Outputs: List containing generated data and related information
##########################################################################################################


generate_data <- function(graph,edges,edge_mat, active_set,k=0,sd_level=1,tau=1,scale=10,type="Gaussian",err_type="normal",est_var=FALSE){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  consider_set = unique(c(edges[,1],edges[,2]))
  G = edges[-active_set,]
  num=nrow(graph)
  
  
  if(k%%2 ==0 ) {
    diff = edge_mat %*% matrix.power(graph_laplacian,k/2)
    
    c = 1
    while(length(consider_set) > 0 ){
      component = connected_component(G,consider_set[1])
      span_vec <- rep(0,nrow(graph))
      span_vec[component] = 1
      consider_set = setdiff(consider_set, component)
      if(c==1){
        basis = span_vec
      }
      else{
        basis = cbind(basis,span_vec)
      }
      c = c + 1
      
    }
    
    if(sum(basis)==length(basis)){
      print("yo")
      basis = as.matrix(basis)
    }
    else {
      basis = cbind(rep(1,nrow(basis)),matrix.power(pinv(graph_laplacian),k/2) %*% basis)
    }
    
  } 
  if(k%%2 ==1) {
    diff = matrix.power(graph_laplacian,(k+1)/2)
    basis = matrix.power(pinv(graph_laplacian),(k+1)/2)
    basis = cbind(rep(1,nrow(basis)), basis[,active_set])
  }
  
  beta <- scale*rnorm(ncol(basis))
  truth = basis %*% beta
  get_basis(truth,edge_mat,k)
  
  err = rnorm(num,0,sd_level)
  
  
  if(err_type == "t"){
    df = 5
    error = rt(num,df) 
    err = sd_level*error/sqrt(df/(df-2))
  }
  if(err_type == "laplace"){
    err = rlaplace(num,0,sd_level/sqrt(2))
  }
  if(err_type == "sn"){
    omega = 1
    alpha = 5
    error = rsn(num,0,omega,alpha)
    mu = omega*(alpha/(sqrt(1+alpha**2)))*sqrt(2/pi)
    sd = sqrt(omega**2 *(1-2*(alpha**2)/(1+alpha**2)/pi))
    err = sd_level*(error-mu)/sd
  }
  
  
  response = truth + err
  
  
  
  if(est_var == FALSE){
    Z = rnorm(num,0,sd_level)
    sd_est = sd_level
  }
  else{
    
    h = run_tf(response,edge_mat,k=k)
    lambda = sqrt(log(num)/num)
    ind = max(which(h$lambda > 1))
    sd_est = sqrt(sum((h$beta[,ind] - response)**2)/h$df[ind])
    Z = rnorm(num,0,sd_est)
  }
  f_Y = response + tau*Z
  g_Y = response - (1/tau)*Z
  
  return(list(beta=beta,basis=basis,mean=truth,Y=response,f_Y=f_Y,g_Y=g_Y,sd_est=sd_est))
}




get_basis <- function(points,edge_mat, k,precision = 8){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0) {
    mat = round(matrix.power(graph_laplacian,k/2) %*% points,precision)
    c=1
    for(vals in unique(mat)){
      print(vals)
      vec = as.numeric(mat == vals)
      if(c==1){
        basis = vec
      }
      else{
        basis = cbind(basis,vec)
      }
      c= c+1
    }
    
    if(sum(basis)==length(basis)){
      print("yo")
      basis = as.matrix(basis)
    }
    else {
      basis = cbind(rep(1,nrow(basis)),matrix.power(pinv(graph_laplacian),k/2) %*% basis)
    }
  }
  else{
    mat = matrix.power(graph_laplacian,(k+1)/2) %*% points
    active_set = which(round(mat,precision) != round(getmode(mat),precision))
    basis = matrix.power(pinv(graph_laplacian),(k+1)/2)
    basis = cbind(rep(1,nrow(basis)), basis[,active_set])
  }
  return(basis)
}




CI_robust <-function(basis,comp_df,select_df,true_mean,sd_high,sd_low,sd_0,cv=0.9){
  alpha = 1- cv
  gamma_low = sd_low **2/ (sd_low**2 + sd_0**2)
  gamma_high = sd_high**2/ (sd_high**2 + sd_0**2)
  eta_Y = proj(basis,comp_df)
  eta_fY =proj(basis,select_df)
  est_high = (eta_Y - gamma_high *eta_fY)/(1-gamma_high)
  est_low = (eta_Y - gamma_low *eta_fY)/(1-gamma_low)
  true_mean = proj(basis,true_mean)
  width_high = (sd_high**2/(1-gamma_high))
  covar_high =width_high* proj_matrix(basis)
  A1 <- apply(cbind(est_low,est_high),1,min)
  A2 <- apply(cbind(est_low,est_high),1,max)
  CI_low = A1 + qnorm(alpha/2)*sqrt(diag(covar_high))
  CI_high = A2 - qnorm(alpha/2)*sqrt(diag(covar_high))
  return(list(pred=cbind(CI_low,CI_high),gamma_low=gamma_low,gamma_high=gamma_high))
}




get_fit <-function(basis,comp_df,select_df,infer_df,true_mean,sd_low,sd_high,sd_0,cv=0.9) {
  CI_rob = CI_robust(basis,comp_df,select_df,true_mean,sd_high,sd_low,sd_0)
  projected_mean <- lm(true_mean~0+basis)$fitted.values
  reg <- lm(infer_df~0+basis)
  pred <- predict(reg,interval="confidence",level = cv)
  fit <- cbind(true_mean,projected_mean,pred,CI_rob$pred,CI_rob$gamma_high,CI_rob$gamma_low)
  
  df <- ncol(basis)
  colnames(fit) <- c("true_mean","proj_mean","fit","lwr","upr","lwr_rob","upr_rob","gamma_low","gamma_high")
  errtrue = mean((fit[,"fit"] - fit[,"true_mean"])**2)
  errproj = mean((fit[,"fit"] - fit[,"proj_mean"])**2)
  length = mean(fit[,"upr"] - fit[,"lwr"])
  length_robust = mean(fit[,"upr_rob"]-fit[,"lwr_rob"])
  cover <- mean((fit[,"proj_mean"] >fit[,"lwr"]) & (fit[,"proj_mean"] <fit[,"upr"]))
  cover_robust <- mean((fit[,"proj_mean"] >fit[,"lwr_rob"]) & (fit[,"proj_mean"] <fit[,"upr_rob"]))
  metrics <- c(df,errtrue,errproj,length,length_robust,cover,cover_robust)
  
  
  return(list(fit=fit,metrics=metrics))
}

get_CI2 <- function(select_df,infer_df,comp_df,edge_mat,true_mean,sd_0,cv=0.9,k=0,K=5){
  alpha = 1-cv
  #Use graph-cross val to select lambda
  splits <- t(sapply(select_df, function(x) multi_split(x,K,sd_0)))
  for(fold in 1:K){
    train <- (rowSums(splits) - splits[,fold])/(1-1/K)
    test <- splits[,fold]*K
    h = run_tf(train,edge_mat,k=k)
    res <- sapply(1:ncol(h$fit),function(x) h$fit[,x] - test)
    mse <- apply(res**2,2,mean)
    if(fold ==1 ){
      cv_mse = mse
    }
    else{
      cv_mse = rbind(cv_mse,mse)
    }
  }
  mse <- apply(cv_mse,2,mean)
  sd <- apply(cv_mse,2,sd)
  
  ind <- which(mse == min(mse))
  ind_1se <- min(which(mse < (mse[ind] + sd[ind])))
  lambda_1se = h$lambda[ind_1se]
  lambda_min = h$lambda[ind]
  
  #SURE to select lambda on f(Y)
  select_tf = run_tf(select_df,edge_mat,k=k)
  resid = sapply(1:dim(select_tf$beta)[2], function(x) select_tf$beta[,x] - select_df)
  mse = colMeans(resid**2)
  sure = mse+2*select_tf$df*sd_0/nrow(select_tf$beta)
  lambda_sure = select_tf$lambda[which(sure == min(sure))]    
  
  fit_1se = coef(select_tf,lambda=lambda_1se)$beta
  basis_1se <- get_basis(fit_1se,edge_mat, k)
  
  fit_min = coef(select_tf,lambda=lambda_min)$beta
  basis_min <- get_basis(fit_min,edge_mat, k)
  
  fit_sure = coef(select_tf,lambda=lambda_sure)$beta
  basis_sure <- get_basis(fit_sure,edge_mat, k)
  
  full_tf = run_tf(comp_df,edge_mat,k=k)
  fit_full <-coef(full_tf,lambda=lambda_min)$beta
  basis_full <- get_basis(fit_full,edge_mat, k)
  
  
  reg_high <- lm(infer_df~0+basis_1se)
  reg_low <- lm(comp_df~0+basis_full)
  sd_low <- summary(reg_low)$sigma
  sd_high <- summary(reg_high)$sigma
  
  if(is.nan(sd_low)){
    sd_low = 0 
  }
  if(is.nan(sd_high)){
    sd_high = sd(comp_df)
  }
  
  one_se <- get_fit(basis_1se,comp_df,select_df,infer_df,true_mean,sd_low,sd_high,sd_0,cv=cv)
  min <- get_fit(basis_min,comp_df,select_df,infer_df,true_mean,sd_low,sd_high,sd_0,cv=cv)
  sure <- get_fit(basis_sure,comp_df,select_df,infer_df,true_mean,sd_low,sd_high,sd_0,cv=cv)
  metrics <- cbind(rbind(one_se$metrics,min$metrics,sure$metrics),sd_low,sd_high,sd_0)
  colnames(metrics) <- c("df","errtrue","errproj","length","length_robust","coverage","coverage_robust","sd_low","sd_high","sd_0")
  rownames(metrics) <- c("1 SE Rule","Min CV","SURE")
  return(list(fit_1se=one_se$fit,fit_min=min$fit,fit_sure=sure$fit,metrics=metrics))
}




elastic_estimates <-function(train,test,edge_mat,k=2,lambda_list = c(0.01,0.1,1,2,3,5,10,100)){
  mse_vec = c()
  for(lambda in lambda_list){
    o <- optim(rep(0,100),function(x) elastic_loss(x,train,edge_mat,k=k,lambda_1=0,lambda_2 =lambda),method="L-BFGS")
    mse = mean((o$par - test)**2)
    mse_vec = c(mse,mse_vec)
    if(lambda == lambda_list[1]){
      pred = o$par
    }
    else{
      pred =cbind(pred,o$par)
    }
  }
  return(list(pred=pred,mse=mse_vec,lambda=lambda_list)) 
}

run_experiment2 <- function(graph,edges,edge_mat,k=0,cv=0.9,sd_level=1,tau=1,scale=10,err_type="normal",est_var=FALSE,K_list =c(2,3,4,5,6,7,8,9,10,15,20,25,30) ){
  if(k %*% 2 == 0){
    active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
  }
  else{
    active_set = c(5,20,37,43)
  }
  dat = generate_data(graph,edges,edge_mat, active_set,scale=scale,k=k,sd_level=sd_level,tau=1,err_type=err_type,est_var=est_var)
  sd_est = dat$sd_est
  for(K in K_list){
    print(K)
    CI <- get_CI(dat$f_Y,dat$g_Y,edge_mat,dat$mean,k=k,sd_level=sd_est,K=K,cv=cv)
    comp <- setNames(data.frame(cbind(CI$num_components,CI$lambda,CI$errproj,CI$errtrue,CI$length,CI$cover)),c("df","lambda","errproj","errtrue","length","cover"))
    comp$K = K
    comp$method = c("1se","min","SURE")
    comp$split = "fission"
    
    CI2 <- get_CI(dat$Y,dat$Y,edge_mat,dat$mean,k=k,sd_level=sd_est,K=K,cv=cv)
    comp2 <- setNames(data.frame(cbind(CI2$num_components,CI2$lambda,CI2$errproj,CI2$errtrue,CI2$length,CI2$cover)),c("df","lambda","errproj","errtrue","length","cover"))
    comp2$K = K
    comp2$method = c("1se","min","SURE")
    comp2$split = "full"
    if(K==K_list[1]){
      results = rbind(comp,comp2)
    }
    else{
      results = rbind(results,comp,comp2)
    }
    
  }
  results= cbind(results,sd_est,sd_level)
  return(results)
}

ridge_soln <-function(Y,edge_mat,k=0,lambda=1){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0 ) {
    penalty = edge_mat %*% matrix.power(graph_laplacian,k/2)
  }
  if(k%%2 == 1){
    penalty = matrix.power(graph_laplacian,(k+1)/2)
  }
  n=dim(penalty)[2]
  gram = diag(n) + lambda*n*t(penalty) %*% penalty
  beta_hat = pinv(gram) %*% Y
  df = sum(diag(pinv(gram)))
  return(list(beta = beta_hat,df=df))
}


run_experiment_ridge <- function(graph,edges,edge_mat,k=0,cv=0.9,sd_level=1,tau=1,scale=1,err_type="normal",est_var=FALSE,K_list =c(2,3,4,5,6,7,8,9,10,15,20,25,30),lambda_list = c(1:100/10000,2:100/1000,2:10/100)){
  if(k %*% 2 == 0){
    active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
  }
  else{
    active_set = c(5,20,37,43)
  }
  dat = generate_data(graph,edges,edge_mat, active_set,k=k,sd_level=sd_level,scale=scale,tau=tau,err_type=err_type,est_var=est_var)
  sd_est = dat$sd_est
  for(K in K_list){
    splits <- t(sapply(dat$Y, function(x) multi_split(x,K,sd_est)))
    for(fold in 1:K){
      train <- (rowSums(splits) - splits[,fold])/(1-1/K)
      test <- splits[,fold]*K
      a=sapply(lambda_list, function(x) ridge_soln(train,edge_mat,k=k,lambda=x))  
      beta= as.matrix(data.frame(a[1,])) 
      df = as.vector(data.frame(a[2,]))
      mse = colMeans((beta- test)**2)
      if(fold ==1 ){
        cv_mse = mse
        cv_df = df
      }
      else{
        cv_mse = rbind(cv_mse,mse)
        cv_df = rbind(cv_df,df)
      }
    }
    
    mse <- apply(cv_mse,2,mean)
    sd <- apply(cv_mse,2,sd)
    
    ind <- which(mse == min(mse))
    ind_1se <- max(which(mse < (mse[ind] + sd[ind])))
    lambda_1se = lambda_list[ind_1se]
    lambda_min = lambda_list[ind]
    
    a = sapply(lambda_list, function(x) ridge_soln(dat$Y,edge_mat,k=k,lambda=x))
    beta= as.matrix(data.frame(a[1,])) 
    df = unlist(a[2,])
    mse = colMeans((beta- as.vector(dat$Y))**2)
    err_true = colMeans((beta- as.vector(dat$mean))**2)
    sure = mse+2*df*sd_est/nrow(beta)
    ind_sure = which(sure == min(sure))
    lambda_sure = lambda_list[ind_sure]    
    res = data.frame(type=c("1se","min","SURE"),
                     df=c(as.matrix(cbind(df[ind_1se],df[ind],df[ind_sure]))),
                     lambda = c(as.matrix(cbind(lambda_1se,lambda_min,lambda_sure))),
                     mse =c(as.matrix(cbind(mse[ind_1se],mse[ind],mse[ind_sure]))),
                     errtrue =c(as.matrix(cbind(err_true[ind_1se],err_true[ind],err_true[ind_sure]))))
    res$k = k 
    res$K = K
    res$scale = scale
    res$sd_level = sd_level
    res$sd_est = dat$sd_est
    res$err_type = err_type
    res$est_var = est_var
    if(K==K_list[1]){
      results = res
    }
    else{
      results = rbind(results,res)
    }
  }
  return(results)
}

run_experiment3 <- function(graph,edges,edge_mat,k=0,cv=0.9,tau=1,scale=10,err_type="normal",est_var=FALSE,K=K,sd_list=c(0.1,0.5,1,2,3,5)){
  if(k %*% 2 == 0){
    active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
  }
  else{
    active_set = c(5,20,37,43)
  }
  for(sd in sd_list){
    print(sd)
    dat = generate_data(graph,edges,edge_mat, active_set,scale=scale,k=k,sd_level=sd,tau=tau,err_type=err_type,est_var=est_var)
    sd_est = dat$sd_est
    res <- get_CI2(dat$f_Y,dat$g_Y,dat$Y,edge_mat,dat$mean,dat$sd_est,cv=cv,k=k,K=K)
    comp <- cbind(res$metrics,k,cv,sd,tau,scale,est_var,K)
    if(sd==sd_list[1]){
      results = comp
    }
    else{
      results = rbind(results,comp)
    }
    method <- rownames(results)
    colnames(results) <- c(colnames(res$metrics),"k","cv","sd","tau","scale","est_var","K")
  }
  return(results)
}



pois_loss <- function(beta,y,edge_mat,k=0,lambda=1){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0 ) {
    penalty = edge_mat %*% matrix.power(graph_laplacian,k/2)
  }
  if(k%%2 == 1){
    penalty = matrix.power(graph_laplacian,(k+1)/2)
  }
  loss <-mean(-y*beta+exp(beta)) + lambda*sum(abs(penalty %*% beta))
  return(loss)
}


run_tf_pois <- function(y,edge_mat,k=0,lambda=1) {
  o <- optim(y,function(x) pois_loss(x,y,edge_mat,k=k,lambda=lambda))
  return(o)
}


split_poisson <- function(y,K){
  
  folds <-sapply(1:length(y),function(i)rmultinom(1,y[i],rep(1,K)/K))
  return(t(folds))
}


get_fit_pois <-function(Y,truth,basis,cv=0.9){
  alpha = 1-cv
  reg <- glm(Y~0+basis,family="poisson")
  proj_trend <- glm(truth~0+basis,family="poisson")$fitted.values
  cov_norm<-vcov(reg)
  cov_robust<-vcovCR(reg,type='CR2',cluster=factor(1:length(Y)))
  cov_pred_norm <- basis %*% cov_norm %*% t(basis)
  cov_pred_robust <- basis %*% cov_robust %*% t(basis)
  pred = exp(basis %*% reg$coefficients)
  lwr = exp(basis%*% reg$coefficients + sqrt(diag(cov_pred_norm))*qnorm(alpha/2) )
  upr = exp(basis%*% reg$coefficients - sqrt(diag(cov_pred_norm))*qnorm(alpha/2) )
  lwr_robust = exp(basis%*% reg$coefficients + sqrt(diag(cov_pred_robust))*qnorm(alpha/2) )
  upr_robust = exp(basis%*% reg$coefficients - sqrt(diag(cov_pred_robust))*qnorm(alpha/2) )
  fit <- cbind(truth, proj_trend,pred,lwr,upr,lwr_robust,upr_robust)
  df <- ncol(basis)
  
  cov_normal <- mean((fit[,4] < fit[,2]) & (fit[,2] < fit[,5]))
  cov_robust  <- mean((fit[,6] < fit[,2]) & (fit[,2] < fit[,7]))
  length_normal = mean(abs(fit[,5]-fit[,4]))
  length_robust = mean(abs(fit[,7]-fit[,6]))
  err_proj = mean(abs(proj_trend - pred))
  err_true =  mean(abs(truth-pred))
  metrics <- c(df,cov_normal,cov_robust,length_normal,length_robust,err_proj,err_true)
  return(list(fit=fit,metrics=metrics))
}


generate_data_pois <- function(graph,edges,edge_mat, active_set,k=0,scale=10){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  consider_set = unique(c(edges[,1],edges[,2]))
  G = edges[-active_set,]
  num=nrow(graph)
  
  
  if(k%%2 ==0 ) {
    diff = edge_mat %*% matrix.power(graph_laplacian,k/2)
    
    c = 1
    while(length(consider_set) > 0 ){
      component = connected_component(G,consider_set[1])
      span_vec <- rep(0,nrow(graph))
      span_vec[component] = 1
      consider_set = setdiff(consider_set, component)
      if(c==1){
        basis = span_vec
      }
      else{
        basis = cbind(basis,span_vec)
      }
      c = c + 1
      
    }
    
    if(sum(basis)==length(basis)){
      basis = as.matrix(basis)
    }
    else {
      basis = cbind(rep(1,nrow(basis)),matrix.power(pinv(graph_laplacian),k/2) %*% basis)
    }
    
  } 
  if(k%%2 ==1) {
    diff = matrix.power(graph_laplacian,(k+1)/2)
    basis = matrix.power(pinv(graph_laplacian),(k+1)/2)
    basis = cbind(rep(1,nrow(basis)), basis[,active_set])
  }
  
  beta <- scale*runif(ncol(basis),0,1)
  truth = basis %*% beta
  get_basis(truth,edge_mat,k)
  
  truth = exp(truth)
  response= rpois(length(truth),truth)
  return(list(beta=beta,basis=basis,mean=truth,Y=response))
}




run_experiment_poisson <- function(graph,edges,ege_mat,scale=1,k=0,cv=0.9,K=5,lambda_list=c(1:10/100,1:10,20,25,30)){
  alpha = 1-cv
  if(k %*% 2 == 0){
    active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
  }
  else{
    active_set = c(5,20,37,43)
  }
  dat<-generate_data_pois(graph,edges,edge_mat, active_set,scale=scale,k=k)
  sel_inf <- split_poisson(dat$Y,2)
  f_Y = sel_inf[,1]
  g_Y = sel_inf[,2]
  splits<-split_poisson(f_Y,K)
  for(fold in 1:K){
    train <- (rowSums(splits) - splits[,fold])
    test <- splits[,fold]
    o = run_tf_pois(train,edge_mat,k=0,lambda=1)
    beta = sapply(lambda_list, function(x) run_tf_pois(train,edge_mat,k=0,lambda=x)$par/(K-1))
    loss = apply(beta,2,function(x) pois_loss(x,test,k=0,edge_mat,lambda=0))
    if(fold ==1 ){
      cv_loss = loss
    }
    else{
      cv_loss = rbind(cv_loss,loss)
    }
  }
  
  
  mean <- apply(cv_loss,2,mean)
  sd <- apply(cv_loss,2,sd)
  
  ind <- min(which(loss == min(loss)))
  ind_1se <- min(which(loss< (mean[ind] + sd[ind])))
  lambda_1se = lambda_list[ind_1se]
  lambda_min = lambda_list[ind]
  beta = sapply(c(lambda_1se,lambda_min), function(x) run_tf_pois(f_Y,edge_mat,k=0,lambda=x)$par)
  
  
  basis_1se = get_basis(beta[,1],edge_mat,k)[,-1]
  basis_min = get_basis(beta[,2],edge_mat,k)[,-1]
  
  
  fit_1se <- get_fit_pois(g_Y,dat$mean/2,basis_1se,cv=cv)
  fit_min <- get_fit_pois(g_Y,dat$mean/2,basis_1se,cv=cv)
  
  return(list(fit_1se = fit_1se$fit,fit_min = fit_min$fit,metrics=rbind(fit_1se$metrics,fit_min$metrics)))
}

