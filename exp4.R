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


multi_split <- function(x,K,sd_level) {
  cov <- matrix(data = -1/(K**2),nrow=K,ncol=K)
  diag(cov) <- rep((1/K)*(1-1/K),K)
  cov <- cov*sd_level
  y <- mvrnorm(n=1, mu = rep(x/K,K),Sigma=cov)
  return(y)
}


generate_data <- function(graph,edges,edge_mat, active_set,k=0,sd_level=1,tau=1){
  
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
  
  #even
  beta <- 10*rnorm(ncol(basis))
  truth = basis %*% beta
  get_basis(truth,edge_mat,k)
  response = truth + rnorm(num,0,sd_level)
  Z = rnorm(num,0,sd_level)
  f_Y = response + tau*Z
  g_Y = response - (1/tau)*Z
  return(list(beta=beta,basis=basis,mean=truth,Y=response,f_Y=f_Y,g_Y=g_Y))
}

get_basis <- function(points,edge_mat, k,precision = 8){
  graph_laplacian <- t(edge_mat) %*% edge_mat
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
  return(basis)
}


get_CI <- function(select_df,infer_df,edge_mat,true_mean,cv=0.9,k=0,sd_level=1,SURE=TRUE,K=5){
  
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
  sd <- apply(cv_mse,2,mean)
  
  ind <- which(mse == min(mse))
  ind_1se <- min(which(mse < (mse[ind] + sd[ind])))
  lambda_1se = h$lambda[ind_1se]
  lambda_min = h$lambda[ind]
  
  h = run_tf(select_df,edge_mat,k=k)
  
  
  resid = sapply(1:dim(h$beta)[2], function(x) h$beta[,x] - select_df)
  mse = colMeans(resid**2)
  sure = mse+2*h$df*sd_level/nrow(h$beta)
  lambda_sure = h$lambda[which(sure == min(sure))]    
  
  
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



run_experiment2 <- function(graph,edges,edge_mat,k=0,cv=0.9,sd_level=1,tau=1,K_list =c(2,3,4,5,6,7,8,9,10,15,20,25,30) ){
  active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
  dat = generate_data(graph,edges,edge_mat, active_set,k=k,sd_level=sd_level,tau=1)
  for(K in K_list){
    print(K)
    CI <- get_CI(dat$f_Y,dat$g_Y,edge_mat,dat$mean,k=k,sd_level=1,K=K,cv=cv)
    comp <- setNames(data.frame(cbind(CI$num_components,CI$lambda,CI$errproj,CI$errtrue,CI$length,CI$cover)),c("df","lambda","errproj","errtrue","length","cover"))
    comp$K = K
    comp$method = c("1se","min","SURE")
    comp$split = "fission"
    
    CI2 <- get_CI(dat$Y,dat$Y,edge_mat,dat$mean,k=k,sd_level=1,K=K,cv=cv)
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
  
  return(results)
}

ntrials = 100


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=0,sd_level=2,tau=1),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval.Rdata")
load("results_graphfission_crossval.Rdata")
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)

df = aggregate(as.numeric(res$df) ~ res$K + res$method + res$split , FUN = mean)
colnames(df) <- c('K',"method","split",'df')
df = df[df$split == "fission",]
df[df$method == "min",2] = "Minimum CV error"
df[df$method == "1se",2] = "1se rule"

cv_df <- df%>%
  ggplot( aes(x=as.numeric(K), y=df,color=as.factor(method))) +
  geom_line(aes(color=as.factor(method)), size = 1.5) +
  geom_point(aes(shape =  as.factor(method), color = as.factor(method)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Number of folds (K)") +
  ylab("Degrees of freedom") +
  geom_hline(yintercept=4, linetype="dashed", 
               color = "black", size=1.5)


ggsave(paste("figures/", "cv_df.pdf",sep=""),cv_df)

df = aggregate(as.numeric(res$errtrue) ~ res$K + res$method + res$split , FUN = mean)
colnames(df) <- c('K',"method","split",'err')
df = df[df$split == "fission",]
df[df$method == "min",2] = "Minimum CV error"
df[df$method == "1se",2] = "1se rule"

cv_err <- df%>%
  ggplot( aes(x=as.numeric(K), y=err,color=as.factor(method))) +
  geom_line(aes(color=as.factor(method)), size = 1.5) +
  geom_point(aes(shape =  as.factor(method), color = as.factor(method)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "bottom", legend.text = element_text(size = 15)) +
  xlab("Number of folds (K)") +
  ylab("Error") 

ggsave(paste("figures/", "cv_error.pdf",sep=""),cv_err)


df = aggregate(as.numeric(res$length) ~ res$K + res$method + res$split , FUN = mean)
colnames(df) <- c('K',"method","split",'length')
df = df[df$split == "fission",]
df[df$method == "min",2] = "Minimum CV error"
df[df$method == "1se",2] = "1se rule"

cv_length <- df%>%
  ggplot( aes(x=as.numeric(K), y=length,color=as.factor(method))) +
  geom_line(aes(color=as.factor(method)), size = 1.5) +
  geom_point(aes(shape =  as.factor(method), color = as.factor(method)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Number of folds (K)") +
  ylab("CI Length") 

ggsave(paste("figures/", "cv_length.pdf",sep=""),cv_length)


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
dat = generate_data(graph,edges,edge_mat, active_set,k=0,sd_level=2,tau=1)
CI <- get_CI(dat$f_Y,dat$g_Y,edge_mat,dat$mean,k=0,sd_level=1,K=10,cv=0.9)
CI$fit_1se



df2 <- cbind(graph,dat$f_Y) 
colnames(df2) <- c("x","y","z")
df3 <- cbind(graph,CI$fit_1se[,5])
colnames(df3) <- c("x","y","z")
df3  <- df3 %>% pivot_wider(names_from =x, values_from = z)
df4 <- cbind(graph,CI$fit_1se[,4])
colnames(df4) <- c("x","y","z")
df4 <- df4 %>% pivot_wider(names_from =x, values_from = z)

df5 <- cbind(graph,CI$fit_1se[,1])
colnames(df5) <- c("x","y","z")
df5 <- df5 %>% pivot_wider(names_from =x, values_from = z)

plot_ly(showscale=FALSE,showlegend=TRUE) %>%  
  add_trace(data = df2, x = ~x, y = ~y, z = ~z, name="Datapoint", mode = "markers", type = "scatter3d", marker = list(size = 5, color = "black", symbol = 104)) %>% 
  add_surface(z=as.matrix(df3)[,-1],opacity=0.2,name="Upper CI", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df4)[,-1],opacity=0.2,name = "Lower CI", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df5)[,-1],opacity=0.8,name = "(Projected) Trend",colorscale = list(c(0, 1), c("blue", "blue")))


df2 <- cbind(graph,dat$f_Y) 
colnames(df2) <- c("x","y","z")
df3 <- cbind(graph,CI$fit_min[,5])
colnames(df3) <- c("x","y","z")
df3  <- df3 %>% pivot_wider(names_from =x, values_from = z)
df4 <- cbind(graph,CI$fit_min[,4])
colnames(df4) <- c("x","y","z")
df4 <- df4 %>% pivot_wider(names_from =x, values_from = z)

df5 <- cbind(graph,CI$fit_min[,1])
colnames(df5) <- c("x","y","z")
df5 <- df5 %>% pivot_wider(names_from =x, values_from = z)

plot_ly(showscale=FALSE,showlegend=TRUE) %>%  
  add_trace(data = df2, x = ~x, y = ~y, z = ~z, name="Datapoint", mode = "markers", type = "scatter3d", marker = list(size = 5, color = "black", symbol = 104)) %>% 
  add_surface(z=as.matrix(df3)[,-1],opacity=0.2,name="Upper CI", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df4)[,-1],opacity=0.2,name = "Lower CI", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df5)[,-1],opacity=0.8,name = "(Projected) Trend",colorscale = list(c(0, 1), c("blue", "blue")))


df2 <- cbind(graph,dat$f_Y) 
colnames(df2) <- c("x","y","z")
df3 <- cbind(graph,CI$fit_sure[,5])
colnames(df3) <- c("x","y","z")
df3  <- df3 %>% pivot_wider(names_from =x, values_from = z)
df4 <- cbind(graph,CI$fit_sure[,4])
colnames(df4) <- c("x","y","z")
df4 <- df4 %>% pivot_wider(names_from =x, values_from = z)

df5 <- cbind(graph,CI$fit_1se[,1])
colnames(df5) <- c("x","y","z")
df5 <- df5 %>% pivot_wider(names_from =x, values_from = z)

plot_ly(showscale=FALSE,showlegend=TRUE) %>%  
  add_trace(data = df2, x = ~x, y = ~y, z = ~z, name="Datapoint", mode = "markers", type = "scatter3d", marker = list(size = 5, color = "black", symbol = 104)) %>% 
  add_surface(z=as.matrix(df3)[,-1],opacity=0.2,name="Upper CI", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df4)[,-1],opacity=0.2,name = "Lower CI", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df5)[,-1],opacity=0.8,name = "(Projected) Trend",colorscale = list(c(0, 1), c("blue", "blue")))




df = aggregate(as.numeric(res$errproj) ~ res$K + res$method + res$split , FUN = mean)



ntrials = 100
k_list =c(0,2,4)
sd_list= c(0.1,0.5,1,2,3)
ntrials = 100
c=1
for(k in k_list) {
  for(sd in sd_list){
      results = cbind(k,sd,ntrials,t(sapply(1:ntrials, function(x) run_experiment(graph,edges,edge_mat,k=k,cv=0.9,sd_level=sd,tau=1))))
      if(c==1){
        res = results
      }
      else{
        res = rbind(res,results)
      }
      c=c+1
  }
}
save(res,file="results_graphfission.Rdata")

colnames(res) <- c('k','sd','ntrials',"cover_fission","cover_full","length_fission","length_full")
res <- data.frame(res)
df = aggregate(res$cover_fission ~ res$k + res$sd , FUN = mean)
df2 = aggregate(res$cover_full ~ res$k + res$sd , FUN = mean)
df <- cbind(df,"Fission")
df2 <- cbind(df2,"Full")
colnames(df) <- c("k","sd","coverage","type")
colnames(df2) <- c("k","sd","coverage","type")
df <- rbind(df,df2)




cover_sd_plot <- df%>%
  ggplot( aes(x=sd, y=coverage,color=as.factor(k))) +
  geom_line(aes(linetype = type,color=as.factor(k)), size = 1.5) +
  geom_point(aes(shape = type, color = as.factor(k)), size = 1) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  geom_hline(yintercept=0.9) + 
  xlab("Amount of Noise (SD)") +
  ylab("Coverage") 


df = aggregate(res$length_fission ~ res$k + res$sd , FUN = mean)
df2 = aggregate(res$length_full ~ res$k + res$sd , FUN = mean)
df <- cbind(df,"Fission")
df2 <- cbind(df2,"Full")
colnames(df) <- c("k","sd","length","type")
colnames(df2) <- c("k","sd","length","type")
df <- rbind(df,df2)


length_sd_plot <- df %>%
  ggplot( aes(x=sd, y=length,color=as.factor(k))) +
  geom_line(aes(linetype = type,color=as.factor(k)), size = 1.5) +
  geom_point(aes(shape = type, color = as.factor(k)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  #theme(legend.key.size = unit(1.5, "cm"))+
  xlab("Amount of Noise (SD)") +
  ylab("CI Length") 


ggsave(paste("figures/", "cover_sd_plot .pdf",sep=""),cover_sd_plot )
ggsave(paste("figures/", "length_sd_plot .pdf",sep=""),length_sd_plot )


