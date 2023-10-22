rm(list=ls())
library(MASS)
library(matrixcalc)
library(pracma)
library(dplyr)
library(tidyr)
library(plotly)
library(genlasso)
source("funs.R")
graph = expand.grid(x=seq(0,9,1),y=seq(0,9,1))


for(node_ind in 1:nrow(graph)){
  print(node_ind)
  node = graph[node_ind,]
  
  dist = apply(graph,1,function(x) sum(abs(x - node)))
  adj = as.numeric(rownames(graph[dist==1,]))
  edge_mat = matrix(data=0,nrow=length(adj),ncol=nrow(graph))
  
  for(i in 1:length(adj)) {
    node_1 = min(node_ind,adj[i])
    node_2 = max(node_ind,adj[i])
    
    edge_mat[i,node_1] = -1 
    edge_mat[i,node_2] = 1
  }
  if(node_ind == 1){
    edge_mat_total = edge_mat
  }
  else{
    edge_mat_total = rbind(edge_mat_total,edge_mat)
  }
}

edge_mat <- as.matrix(distinct(data.frame(edge_mat_total)))
edges <- t(apply(edge_mat,1,function(x) which(x != 0)))
run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=3,scale=1,tau=1,K_list = c(2,3,4,5,6,10),lambda_list = c(0.00001,0.00004,0.00007,0.001,0.002,0.003,0.004,0.1,0.2,0.3,0.4))


run_experiment_ridge <- function(graph,edges,edge_mat,k=0,cv=0.9,sd_level=1,tau=1,scale=10,err_type="normal",est_var=FALSE,K_list =c(2,3,4,5,6,7,8,9,10,15,20,25,30),lambda_list = c(1:100/10000,2:100/1000,2:10/100)){
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


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=3,scale=1,tau=1,K_list = c(2,3,4,5,6,10),lambda_list = c(0.00001,0.00004,0.00007,0.001,0.002,0.003,0.004,0.1,0.2,0.3,0.4))


res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_normal.Rdata")











pois_loss <- function(beta,y,edge_mat,k=0,lambda=1){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0 ) {
    penalty = edge_mat %*% matrix.power(graph_laplacian,k/2)
  }
  if(k%%2 == 1){
    penalty = matrix.power(graph_laplacian,(k+1)/2)
  }
  loss <-mean(-beta*y + exp(beta)) + lambda*sum(abs(penalty %*% beta))
  return(loss)
}


run_tf_pois <- function(y,edge_mat,k=0,lambda=1) {
  o <- optim(y,function(x) pois_loss(x,y,edge_mat,k=k,lambda=lambda))
  return(o)
}



generate_data_pois <- function(graph,edges,edge_mat, active_set,k=0,sd_level=1,tau=1,scale=10,K=5){
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


split_poisson <- function(y,K){
  
  folds <-sapply(1:length(y),function(i)rmultinom(1,y[i],rep(1,K)/K))
  return(t(folds))
}

active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
dat<-generate_data_pois(graph,edges,edge_mat, active_set,scale=2,k=0,sd_level=1,tau=1)
K=6
k=0
Y<-dat$Y
sel_inf <- split_poisson(Y,2)
splits<-split_poisson(Y,K)

sapply(c(0.1,1,2,3), function(x) log(run_tf_pois(train,edge_mat,k=0,lambda=x)$par)/(K-1))

run_tf_pois(dat$Y,edge_mat,k=0,lambda=0.001)

for(fold in 1:K){
  train <- (rowSums(splits) - splits[,fold])
  test <- splits[,fold]
  o = run_tf_pois(train,edge_mat,k=0,lambda=1)
  beta = log(o$par)/(K-1)
  
  
  basis = get_basis(beta,edge_mat,k=k)
  reg_train = glm(train~0+basis,family="poisson")
  loss = pois_loss(beta,test,k=0,edge_mat,lambda=0)/length(theta_hat)
  df = ncol(basis)
  if(fold ==1 ){
    cv_loss = loss
    cv_df = df
  }
  else{
    cv_loss = rbind(cv_loss,loss)
    cv_df = rbind(cv_df,df)
  }
}


mse <- apply(cv_mse,2,mean)
sd <- apply(cv_mse,2,sd)

ind <- which(mse == min(mse))
ind_1se <- max(which(mse < (mse[ind] + sd[ind])))
lambda_1se = lambda_list[ind_1se]
lambda_min = lambda_list[ind]




run_tf_pois(dat$Y,0+edge_mat,k=0,lambda=1)


o = run_tf_pois(train,edge_mat,k=0,lambda=1)
beta = log(o$par)/(K-1)
basis = get_basis(beta,edge_mat,k=k)
reg_train = glm(train~0+basis,family="poisson")
loss_test = pois_loss(beta,test,k=0,edge_mat,lambda=0)/length(theta_hat)



pois_loss(beta,test,k=0,edge_mat,lambda=0)



beta = o$par*K/(K-1)


h = run_tf(train,edge_mat,k=k)
res <- sapply(1:ncol(h$fit),function(x) h$fit[,x] - test)
mse <- apply(res**2,2,mean)
if(fold ==1 ){
  cv_mse = mse
}
else{
  cv_mse = rbind(cv_mse,mse)
}

graph = expand.grid(x=seq(0,9,1),y=seq(0,9,1))


for(node_ind in 1:nrow(graph)){
  print(node_ind)
  node = graph[node_ind,]
  
  dist = apply(graph,1,function(x) sum(abs(x - node)))
  adj = as.numeric(rownames(graph[dist==1,]))
  edge_mat = matrix(data=0,nrow=length(adj),ncol=nrow(graph))
  
  for(i in 1:length(adj)) {
    node_1 = min(node_ind,adj[i])
    node_2 = max(node_ind,adj[i])
    
    edge_mat[i,node_1] = -1 
    edge_mat[i,node_2] = 1
  }
  if(node_ind == 1){
    edge_mat_total = edge_mat
  }
  else{
    edge_mat_total = rbind(edge_mat_total,edge_mat)
  }
}

edge_mat <- as.matrix(distinct(data.frame(edge_mat_total)))
edges <- t(apply(edge_mat,1,function(x) which(x != 0)))



pois_loss <- function(beta,y,edge_mat,k=0,lambda=1){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0 ) {
    penalty = edge_mat %*% matrix.power(graph_laplacian,k/2)
  }
  if(k%%2 == 1){
    penalty = matrix.power(graph_laplacian,(k+1)/2)
  }
  loss <-sum(-beta*y + exp(beta)) + lambda*sum(penalty %*% beta)
  return(loss)
}


run_tf_pois <- function(y,edge_mat,k=0,lambda=1) {
  o <- optim(y,function(x) pois_loss(x,y,edge_mat,k=k,lambda=lambda))
  return(o)
}



truth = exp(dat$mean)
y= rpois(length(truth),truth)
o <- optim(y,function(x) pois_loss(x,y,edge_mat,k=2,lambda=2))
graph_laplacian <- t(edge_mat) %*% edge_mat
k=2
    mat = round(matrix.power(graph_laplacian,k/2) %*% o$par,7)
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

get_basis(o$par,edge_mat,0)

a = run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,lambda_list = c(1:100/10000,2:100/1000,2:10/100))



a= run_experiment3(graph,edges,edge_mat,k=1,cv=0.9,tau=1,scale=0.5,err_type="sn",est_var=TRUE,K=4,sd_list=c(1,2,3,4))

sd= 2


k=0
scale = 0.5
tau= 1
err_type="sn"
est_var=TRUE
cv=0.9
if(k %*% 2 == 0){
  active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
}else{
  active_set = c(5,20,37,43)
}
K=4

dat = generate_data(graph,edges,edge_mat, active_set,scale=0.25,k=2,sd_level=sd,tau=tau,err_type=err_type,est_var=est_var)
sd_est = dat$sd_est
res <- get_CI2(dat$f_Y,dat$g_Y,dat$Y,edge_mat,dat$mean,dat$sd_est,cv=cv,k=k,K=K)


k=0
sd_level =3
scale = 0.5
tau= 1
err_type="sn"
est_var=TRUE
cv=0.9
if(k %*% 2 == 0){
  active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
}else{
  active_set = c(5,20,37,43)
}


dat = generate_data(graph,edges,edge_mat, active_set,k=k,sd_level=sd_level,scale=scale,tau=tau,err_type=err_type,est_var=est_var)

select_df=dat$f_Y
infer_df = dat$g_Y
comp_df=dat$Y
edge_mat = edge_mat
true_mean=dat$Y
sd_0 = dat$sd_est
K=5

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
gamma = sd_level**2/(sd_level**2 + sd_0**2)

one_se <- get_fit(basis_1se,comp_df,select_df,infer_df,true_mean,sd_low,sd_high,sd_0)
min <- get_fit(basis_min,comp_df,select_df,infer_df,true_mean,sd_low,sd_high,sd_0)
sure <- get_fit(basis_sure,comp_df,select_df,infer_df,true_mean,sd_low,sd_high,sd_0)
metrics <- cbind(rbind(one_se$metrics,min$metrics,sure$metrics),sd_low,sd_high,sd_0)
colnames(metrics) <- c("df","errtrue","errproj","length","length_robust","coverage","coverage_robust","sd_low","sd_high","sd_0")
rownames(metrics) <- c("1 SE Rule","Min CV","SURE")



CI_robust <-function(basis,comp_df,select_df,sd_high,sd_low,sd_0){
  gamma = sd_level**2/(sd_level**2 + sd_0**2)
  gamma_low = sd_low **2/ (sd_low**2 + sd_0**2)
  gamma_high = sd_high**2/ (sd_high**2 + sd_0**2)
  eta_Y = proj(basis,comp_df)
  eta_fY =proj(basis,select_df)
  est_high = (eta_Y - gamma_high *eta_fY)/(1-gamma_high)
  est_low = (eta_Y - gamma_low *eta_fY)/(1-gamma_low)
  true_mean = proj(basis,dat$mean)
  width_true = (sd_level**2/(1-gamma))
  width_high = (sd_high**2/(1-gamma_high))
  covar =width_true* proj_matrix(basis)
  covar_high =width_high* proj_matrix(basis)
  A1 <- apply(cbind(est_low,est_high),1,min)
  A2 <- apply(cbind(est_low,est_high),1,max)
  CI_low = A1 + qnorm(alpha/2)*sqrt(diag(covar_high))
  CI_high = A2 - qnorm(alpha/2)*sqrt(diag(covar_high))
  return(list(pred=cbind(CI_low,CI_high),gamma_low=gamma_low,gamma_high=gamma_high))
}


get_fit <-function(basis,comp_df,select_df,infer_df,true_mean,sd_low,sd_high,sd_0) {
  browser()
  CI_rob = CI_robust(basis,comp_df,select_df,sd_high,sd_low,sd_0)
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




res <- get_CI2(dat$f_Y,dat$g_Y,dat$Y,edge_mat,dat$mean,dat$sd_est,cv=cv,k=k,K=K)




res=mclapply(1:3, function(x) run_experiment3(graph,edges,edge_mat,k=0,cv=0.9,sd_level=1,tau=1,scale=3,err_type="normal",est_var=TRUE,K=2,sd_list=c(1,2)),mc.cores=4)

a= run_experiment3(graph,edges,edge_mat,k=0,cv=0.9,sd_level=1,tau=1,scale=3,err_type="normal",est_var=TRUE,K=2,sd_list=c(1,2))


results <- data.frame(results)
results$method = method
results$err_type = err_type


fit = fit_sure

errtrue = mean((fit[,"fit"] - fit[,"true_mean"])**2)
errproj = mean((fit[,"fit"] - fit[,"proj_mean"])**2)
length = mean(fit[,"upr"] - fit[,"lwr"])
length_robust = mean(fit[,"upr_rob"]-fit[,"lwr_rob"])
cover <- mean((fit[,"proj_mean"] >fit[,"lwr"]) & (fit[,"proj_mean"] <fit[,"upr"]))
cover_rob <- mean((fit[,"proj_mean"] >fit[,"lwr_rob"]) & (fit[,"proj_mean"] <fit[,"upr_rob"]))


basis = basis_1se


reg <- lm(infer_df~0+basis)
pred <- predict(reg,interval="confidence",level = cv)
projected_mean <- lm(true_mean~0+basis)$fitted.values
fit = setNames(cbind(true_mean,projected_mean,pred),c("proj_mean","fit","lwr","upr"))
errtrue = mean((pred[,1] - true_mean)**2)
errproj = mean((pred[,1] - projected_mean)**2)



full_tf = run_tf(comp_df,edge_mat,k=k)
fit_full <-coef(full_tf,lambda=lambda_min)$beta
basis_full <- get_basis(fit_full,edge_mat, k)



reg_high <- lm(comp_df~0+basis_1se)
reg_low <- lm(comp_df~0+basis_full)
sd_low <- summary(reg_low)$sigma
sd_high <- summary(reg_high)$sigma



reg_high <- lm(comp_df~0+basis)
reg_low <- lm(comp_df~0+basis_full)
sd_low <- summary(reg_low)$sigma
sd_high <- summary(reg_high)$sigma
if(is.nan(sd_low)){
  sd_low = 0 
}
if(is.nan(sd_high)){
  sd_high = sd(comp_df)
}
gamma = sd_level**2/(sd_level**2 + sd_0**2)
CI_rob = CI_robust(basis,comp_df,select_df,sd_high,sd_low)
CI_naive = predict(reg_high,level=1-alpha,interval="confidence")






CI_metric_trad <-function(infer_df,basis,true_mean,cv=cv){
  reg <- lm(infer_df~0+basis)
  pred <- predict(reg,interval="confidence",level = cv)
  projected_mean <- lm(true_mean~0+basis)$fitted.values
  fit = setNames(cbind(projected_mean,pred),c("proj_mean","fit","lwr","upr"))
  errtrue = mean((pred[,1] - true_mean)**2)
  errproj = mean((pred[,1] - projected_mean)**2)
  length = mean(fit[,5] - fit_sure[,4])
  cover <- mean((fit[,1] >fit_sure[,4]) & (fit_sure[,1] <fit_sure[,5]))
  return(list(pred=pred,errtrue=errtrue,errproj=errproj,length=length,cover=cover))
}





get_CI <- function(select_df,infer_df,edge_mat,true_mean,cv=0.9,k=0,sd_level=1,K=5){
  
  h = run_tf(select_df,edge_mat,k=k)
  
  
  resid = sapply(1:dim(h$beta)[2], function(x) h$beta[,x] - select_df)
  mse = colMeans(resid**2)
  sure = mse+2*h$df*sd_level/nrow(h$beta)
  lambda_sure = h$lambda[which(sure == min(sure))]    
  
  
  
  splits <- t(sapply(select_df, function(x) multi_split(x,K,sd_level)))
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
  
  fit_1se = coef(h,lambda=lambda_1se)$beta
  basis_1se <- get_basis(fit,edge_mat, k)
  
  fit_min = coef(h,lambda=lambda_min)$beta
  basis_min <- get_basis(fit,edge_mat, k)
  
  fit_sure = coef(h,lambda=lambda_sure)$beta
  basis_sure <- get_basis(fit,edge_mat, k)
  
  
  
  fit = coef(h,lambda=lambda_min)$beta
  basis <- get_basis(fit,edge_mat, k)
  reg <- lm(infer_df~0+basis)
  pred <- predict(reg,interval="confidence",level = cv)
  projected_mean <- lm(true_mean~0+basis)$fitted.values
  fit_min = setNames(cbind(projected_mean,fit,pred),c("proj_mean","fit_1","fit_2","lwr","upr"))
  errtrue_min = mean((pred[,1] - true_mean)**2)
  errproj_min = mean((pred[,1] - projected_mean)**2)
  df_min <- ncol(basis)
  length_min = mean(fit_min[,5] - fit_min[,4])
  cover_min <- mean((fit_min[,1] >fit_min[,4]) & (fit_min[,1] <fit_min[,5]))
  

  df_sure <- ncol(basis)
  
  return(list(fit_1se=fit_1se,fit_min=fit_min,fit_sure = fit_sure,num_components = c(df_1se,df_min,df_sure),lambda=c(lambda_1se,lambda_min,lambda_sure),
              errproj=c(errproj_1se,errproj_min,errproj_sure),errtrue=c(errtrue_1se,errtrue_min,errtrue_sure),
              cover=c(cover_1se,cover_min,cover_sure),length=c(length_1se,length_min,length_sure)))
}







CI_metric_trad <-function(infer_df,basis,true_mean){
  reg <- lm(infer_df~0+basis)
  pred <- predict(reg,interval="confidence",level = cv)
  projected_mean <- lm(true_mean~0+basis)$fitted.values
  fit = setNames(cbind(projected_mean,pred),c("proj_mean","fit_1","fit_2","lwr","upr"))
  errtrue = mean((pred[,1] - true_mean)**2)
  errproj = mean((pred[,1] - projected_mean)**2)
  length = mean(fit_sure[,5] - fit_sure[,4])
  cover <- mean((fit_sure[,1] >fit_sure[,4]) & (fit_sure[,1] <fit_sure[,5]))
  return(list())
}



reg <- lm(infer_df~0+basis)
pred <- predict(reg,interval="confidence",level = cv)
projected_mean <- lm(true_mean~0+basis)$fitted.values
fit = setNames(cbind(projected_mean,fit,pred),c("proj_mean","fit_1","fit_2","lwr","upr"))
errtrue = mean((pred[,1] - true_mean)**2)
errproj = mean((pred[,1] - projected_mean)**2)
length = mean(fit_sure[,5] - fit_sure[,4])
cover <- mean((fit_sure[,1] >fit_sure[,4]) & (fit_sure[,1] <fit_sure[,5]))


alpha = 0.1
sd_level = 2
dat = generate_data(graph,edges,edge_mat, active_set,k=0,sd_level=sd_level,scale=3,tau=1,err_type="normal",est_var=TRUE)
select_df = dat$f_Y
infer_df = dat$g_Y
comp_df = dat$Y
sd_0 = dat$sd_est

h = run_tf(select_df,edge_mat,k=k)
resid = sapply(1:dim(h$beta)[2], function(x) h$beta[,x] - select_df)
mse = colMeans(resid**2)
sure = mse+2*h$df*sd_level/nrow(h$beta)
lambda_sure = h$lambda[which(sure == min(sure))]    



splits <- t(sapply(select_df, function(x) multi_split(x,K,sd_level)))
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


fit = coef(h,lambda=lambda_1se)$beta
basis_1se <- get_basis(fit,edge_mat, k)

fit = coef(h,lambda=lambda_min)$beta
basis <- get_basis(fit,edge_mat, k)


full_tf = run_tf(comp_df,edge_mat,k=k)
resid = sapply(1:dim(h$beta)[2], function(x) h$beta[,x] - select_df)
mse = colMeans(resid**2)
sure = mse+2*h$df*sd_level/nrow(h$beta)
lambda_sure_full = h$lambda[which(sure == min(sure))]   
fit_full <-coef(full_tf,lambda=lambda_sure_full)$beta
basis_full <- get_basis(fit_full,edge_mat, k)


fit = coef(h,lambda=lambda_1se)$beta
basis_1se <- get_basis(fit,edge_mat, k)


reg_high <- lm(comp_df~0+basis)
reg_low <- lm(comp_df~0+basis_full)
sd_low <- summary(reg_low)$sigma
sd_high <- summary(reg_high)$sigma
if(is.nan(sd_low)){
  sd_low = 0 
}
if(is.nan(sd_high)){
  sd_high = sd(comp_df)
}
gamma = sd_level**2/(sd_level**2 + sd_0**2)
CI_rob = CI_robust(basis,comp_df,select_df,sd_high,sd_low)
CI_naive = predict(reg_high,level=1-alpha,interval="confidence")




reg <- lm(infer_df~0+basis)
pred <- predict(reg,interval="confidence",level = cv)
projected_mean <- lm(true_mean~0+basis)$fitted.values
fit_1se = setNames(cbind(projected_mean,fit,pred),c("proj_mean","fit_1","fit_2","lwr","upr"))
errtrue_1se = mean((pred[,1] - true_mean)**2)
errproj_1se = mean((pred[,1] - projected_mean)**2)
length_1se = mean(fit_1se[,5] - fit_1se[,4])
cover_1se <- mean((fit_1se[,1] >fit_1se[,4]) & (fit_1se[,1] <fit_1se[,5]))
df_1se <- ncol(basis)



qnorm(alpha/2)*sqrt(diag(covar))
qnorm(alpha/2)*sqrt(diag(covar_high))


A1 <- apply(cbind(est_low,est_high),1,min)
A2 <- apply(cbind(est_low,est_high),1,max)
CI_low = A1 + qnorm(alpha/2)*sqrt(diag(covar_high))
CI_high = A2 - qnorm(alpha/2)*sqrt(diag(covar_high))
mean((true_mean > CI_high) | (true_mean < CI_low))
