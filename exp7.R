library(MASS)
library(matrixcalc)
library(pracma)
library(dplyr)
library(tidyr)
library(plotly)
library(genlasso)

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

elastic_loss <- function(beta,y,edge_mat,k=0,lambda_1=1,lambda_2=1){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0 ) {
    penalty = edge_mat %*% matrix.power(graph_laplacian,k/2)
  }
  if(k%%2 == 1){
    penalty = matrix.power(graph_laplacian,(k+1)/2)
  }
  sum((beta - y)**2) + lambda_1 * sum(penalty%*% beta) + lambda_2 * sum((penalty%*% beta)**2 )
}






multi_split <- function(x,K,sd_level) {
  cov <- matrix(data = -1/(K**2),nrow=K,ncol=K)
  diag(cov) <- rep((1/K)*(1-1/K),K)
  cov <- cov*sd_level
  y <- mvrnorm(n=1, mu = rep(x/K,K),Sigma=cov)
  return(y)
}



generate_data <- function(graph,edges,edge_mat, active_set,k=0,sd_level=1,tau=1,scale=10,type="Gaussian",err_type="normal",est_var=FALSE){
  
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





getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
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



get_CI <- function(select_df,infer_df,edge_mat,true_mean,cv=0.9,k=0,sd_level=1,SURE=TRUE,K=5){
  
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
  basis <- get_basis(fit,edge_mat, k)
  reg <- lm(infer_df~0+basis)
  pred <- predict(reg,interval="confidence",level = cv)
  projected_mean <- lm(true_mean~0+basis)$fitted.values
  fit_1se = setNames(cbind(projected_mean,fit,pred),c("proj_mean","fit_1","fit_2","lwr","upr"))
  errtrue_1se = mean((pred[,1] - true_mean)**2)
  errproj_1se = mean((pred[,1] - projected_mean)**2)
  length_1se = mean(fit_1se[,5] - fit_1se[,4])
  cover_1se <- mean((fit_1se[,1] >fit_1se[,4]) & (fit_1se[,1] <fit_1se[,5]))
  df_1se <- ncol(basis)
  
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
  
  fit = coef(h,lambda=lambda_sure)$beta
  basis <- get_basis(fit,edge_mat, k)
  reg <- lm(infer_df~0+basis)
  pred <- predict(reg,interval="confidence",level = cv)
  projected_mean <- lm(true_mean~0+basis)$fitted.values
  fit_sure = setNames(cbind(projected_mean,fit,pred),c("proj_mean","fit_1","fit_2","lwr","upr"))
  errtrue_sure = mean((pred[,1] - true_mean)**2)
  errproj_sure = mean((pred[,1] - projected_mean)**2)
  length_sure = mean(fit_sure[,5] - fit_sure[,4])
  cover_sure <- mean((fit_sure[,1] >fit_sure[,4]) & (fit_sure[,1] <fit_sure[,5]))
  
  df_sure <- ncol(basis)
  
  return(list(fit_1se=fit_1se,fit_min=fit_min,fit_sure = fit_sure,num_components = c(df_1se,df_min,df_sure),lambda=c(lambda_1se,lambda_min,lambda_sure),
              errproj=c(errproj_1se,errproj_min,errproj_sure),errtrue=c(errtrue_1se,errtrue_min,errtrue_sure),
              cover=c(cover_1se,cover_min,cover_sure),length=c(length_1se,length_min,length_sure)))
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

get_active <- function(graph, row,x_sep=4,y_sep=5) {
  flag <- (graph[row[1],1] == x_sep)*(graph[row[2],1] == x_sep + 1)
  flag2 <- (graph[row[1],2] == y_sep)*(graph[row[2],2] == y_sep + 1)
  return(flag | flag2)
}


elastic_loss <- function(beta,y,edge_mat,k=0,lambda_1=1,lambda_2=1){
  graph_laplacian <- t(edge_mat) %*% edge_mat
  if(k%%2 ==0 ) {
    penalty = edge_mat %*% matrix.power(graph_laplacian,k/2)
  }
  if(k%%2 == 1){
    penalty = matrix.power(graph_laplacian,(k+1)/2)
  }
  mean((beta - y)**2) + lambda_1 * sum(penalty%*% beta) + lambda_2 * sum((penalty%*% beta)**2 )
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


run_experiment_ridge <- function(graph,edges,edge_mat, active_set,k=2,sd_level=1,scale=1,tau=1,err_type="normal",est_var=FALSE,lambda_list = c(0.01,0.1,1,2,3,5,10,100),K_list =c(2,3,4,5,6,7,8,9,10,15)){
  dat = generate_data(graph,edges,edge_mat, active_set,k=k,sd_level=sd_level,scale=scale,tau=tau,err_type=err_type,est_var=est_var)
  sd_level = dat$sd_est
  
  for(K in K_list){
    splits <- t(sapply(dat$Y, function(x) multi_split(x,K,sd_level)))
    for(fold in 1:K){
      train <- (rowSums(splits) - splits[,fold])/(1-1/K)
      test <- splits[,fold]*K
      h = ridge_estimates(train,test,edge_mat,k=k,lambda_list=lambda_list)
      if(fold ==1 ){
        cv_mse = h$mse
      }
      else{
        cv_mse = rbind(cv_mse,h$mse)
      }
    }
    mse <- apply(cv_mse,2,mean)
    sd <- apply(cv_mse,2,sd)
    
    ind <- which(mse == min(mse))
    ind_1se <- min(which(mse < (mse[ind] + sd[ind])))
    lambda_1se = lambda_list[ind_1se]
    lambda_min = lambda_list[ind]
    vec_1se = rep(0,length(lambda_list))
    vec_min = rep(0,length(lambda_list))
    vec_min[ind] = 1
    vec_1se[ind_1se] = 1
    
    est = ridge_estimates(dat$Y,dat$Y,edge_mat,k=k,lambda_list=lambda_list)
    err_true = sapply(1:ncol(est$pred), function(x) mean((est$pred[,x] - dat$mean)**2))
    df = cbind(K,lambda_list,mse,sd,err_true,vec_min,vec_1se)
    if(K==K_list[1]){
      results = df
    }
    else{
      results = rbind(results,df)
    }
  }
  return(results)
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
  df = tr(pinv(gram) )
  return(list(beta = beta_hat,df=df))
}
run_experiment_ridge <- function(graph,edges,edge_mat,k=0,cv=0.9,sd_level=1,tau=1,scale=10,err_type="normal",est_var=FALSE,K_list =c(2,3,4,5,6,7,8,9,10,15,20,25,30),lambda_list = c(1:100/1000,1:1000/10)){
  if(k %*% 2 == 0){
    active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
  }
  else{
    active_set = c(5,20,37,43)
  }
  dat = generate_data(graph,edges,edge_mat, active_set,k=k,sd_level=sd_level,scale=scale,tau=tau,err_type=err_type,est_var=est_var)
  sd_level = dat$sd_est
  for(K in K_list){
    splits <- t(sapply(dat$Y, function(x) multi_split(x,K,sd_level)))
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
    ind_1se <- min(which(mse < (mse[ind] + sd[ind])))
    lambda_1se = lambda_list[ind_1se]
    lambda_min = lambda_list[ind]
    
    a = sapply(lambda_list, function(x) ridge_soln(dat$Y,edge_mat,k=k,lambda=x))
    beta= as.matrix(data.frame(a[1,])) 
    df = as.vector(data.frame(a[2,]))
    mse = colMeans((beta- as.vector(dat$Y))**2)
    err_true = colMeans((beta- as.vector(dat$mean))**2)
    sure = mse+2*df*sd_level/nrow(beta)
    ind_sure = which(sure == min(sure))
    lambda_sure = lambda_list[ind_sure]    
    res = data.frame(type=c("1se","min","SURE"),
                     df=c(as.matrix(cbind(df[ind_1se],df[ind],df[ind_sure]))),
                     lambda = c(as.matrix(cbind(lambda_1se,lambda_min,lambda_sure))),
                     mse =c(as.matrix(cbind(mse[ind_1se],mse[ind],mse[ind_sure]))),
                     errtrue =c(as.matrix(cbind(err_true[ind_1se],err_true[ind],err_true[ind_sure]))))
    res$K = K
    if(K==K_list[1]){
      results = res
    }
    else{
      results = rbind(results,res)
    }
  }
  return(results)
}



res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=1,tau=1),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_0_normal.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=1,tau=1,est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_0_normal_estvar.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="sn"),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_0_sn.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="sn",est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_0_sn_estvar.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t"),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_0_t.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_0_t_estvar.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="laplace"),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_0_laplace.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="laplace",est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_0_laplace_estvar.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,sd_level=1,tau=1,scale=2),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_normal.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,sd_level=1,tau=1,est_var=TRUE,scale=2),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_normal_estvar.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,sd_level=1,tau=1,err_type="sn",scale=2),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_sn.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,sd_level=1,tau=1,err_type="sn",est_var=TRUE,scale=2),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_sn_estvar.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,sd_level=1,tau=1,err_type="t",scale=2),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_t.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,sd_level=1,tau=1,err_type="t",est_var=TRUE,scale=2),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_t_estvar.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,scale=2,sd_level=1,tau=1,err_type="laplace"),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_laplace.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,scale=2,sd_level=1,tau=1,err_type="laplace",est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_laplace_estvar.Rdata")



res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=2,scale=1,sd_level=1,tau=1),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_2_normal.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=2,scale=1,sd_level=1,tau=1,est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_2_normal_estvar.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=2,scale=1,sd_level=1,tau=1,err_type="sn"),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_2_sn.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=2,scale=1,sd_level=1,tau=1,err_type="sn",est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_2_sn_estvar.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=2,scale=1,sd_level=1,tau=1,err_type="t"),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_2_t.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=2,scale=1,sd_level=1,tau=1,err_type="t",est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_2_t_estvar.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=2,scale=0.5,sd_level=1,tau=1,err_type="laplace"),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_2_laplace.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=2,scale=1,sd_level=1,tau=1,err_type="laplace",est_var=TRUE),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_2_laplace_estvar.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,lambda_list = c(1:1000/100000)))
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_ridge_normal.Rdata")

active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,est_var=TRUE,lambda_list = c(1:1000/100000)))
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_ridge_normal_estvar.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="sn",lambda_list = c(1:1000/100000)))
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_ridge_sn.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="sn",est_var=TRUE,lambda_list = c(1:1000/100000)))
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_ridge_sn_estvar.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",lambda_list = c(1:1000/100000)))
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_ridge_t.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",est_var=TRUE,lambda_list = c(1:1000/100000)))
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_ridge_t_estvar.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",lambda_list = c(1:1000/100000)))
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_ridge_laplace.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",est_var=TRUE,lambda_list = c(1:1000/100000)))
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_ridge_laplace_estvar.Rdata")
