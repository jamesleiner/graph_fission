
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
    print(basis)
    
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


get_CI <- function(select_df,infer_df,edge_mat,true_mean,cv=0.9,k=0){
  h = run_tf(select_df,edge_mat,k=k)
  fit = coef(h,lambda=20)$beta
  basis <- get_basis(fit,edge_mat, k)
  reg <- lm(infer_df~0+basis)
  pred <- predict(reg,interval="confidence",level = cv)
  projected_mean <- lm(true_mean~0+basis)$fitted.values
  return(setNames(cbind(projected_mean,fit,pred),c("proj_mean","fit_1","fit_2","lwr","upr")))
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
graph_laplacian <- t(edge_mat) %*% edge_mat
edges <- t(apply(edge_mat,1,function(x) which(x != 0)))



run_experiment <- function(graph,edges,edge_mat,k=0,pct_active=0.3,cv=0.9,sd_level=1,tau=1){
  active_set = sample(1:nrow(edges))[1:floor(pct_active*nrow(edges))]
  dat = generate_data(graph,edges,edge_mat, active_set,k=k,sd_level=sd_level,tau=tau)
  full_CI <- get_CI(dat$Y,dat$Y,edge_mat,dat$mean,cv=cv,k=k)
  fission_CI <- get_CI(dat$f_Y,dat$g_Y,edge_mat,dat$mean,cv=cv,k=k)
  cover_fission <- mean((fission_CI[,1] >fission_CI[,4]) & (fission_CI[,1] <fission_CI[,5]))
  cover_full <- mean((full_CI[,1] >full_CI[,4]) & (full_CI[,1] <full_CI[,5]))
  fission_length = mean(fission_CI[,5] - fission_CI[,4])
  full_length = mean(full_CI[,5] - full_CI[,4])
  return(c(cover_fission,cover_full,fission_length,full_length))
}

ntrials = 100
k_list =c(0,2,4)
sd_list= c(0.1,0.5,1,2,3)
active_list = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)
ntrials = 100
c=1
for(k in k_list) {
  for(sd in sd_list){
    for(pct_active in active_list){
      results = cbind(k,sd,pct_active,ntrials,t(sapply(1:20, function(x) run_experiment(graph,edges,edge_mat,k=k,pct_active=0.5,cv=0.9,sd_level=1,tau=1))))
      if(c==1){
        res = results
      }
      else{
        res = rbind(res,results)
      }
      c=c+1
    }
  }
}
save(res,file="results_graphfission.Rdata")

colnames(res) <- c('k','sd','pct_active','ntrials',"cover_fission","cover_full","length_fission","length_full")
res <- data.frame(res)
df = aggregate(res$cover_fission ~ res$k + res$sd + res$pct_active, FUN = mean)
df2 = aggregate(res$cover_full ~ res$k + res$sd + res$pct_active, FUN = mean)
df <- cbind(df,"Fission")
df2 <- cbind(df2,"Full")
colnames(df) <- c("k","sd","pct_active","coverage","type")
colnames(df2) <- c("k","sd","pct_active","coverage","type")
df <- rbind(df,df2)

cover_active_plot <- df[df$sd==1,] %>%
  ggplot( aes(x=pct_active, y=coverage,color=as.factor(k))) +
  geom_line(aes(linetype = type,color=as.factor(k)), size = 1.5) +
  geom_point(aes(shape = type, color = as.factor(k)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "bottom", legend.text = element_text(size = 15)) +
  xlab("Percentage of Nodes Active") +
  ylab("Coverage") 



cover_sd_plot <- df[df$pct_active==0.5,] %>%
  ggplot( aes(x=sd, y=coverage,color=as.factor(k))) +
  geom_line(aes(linetype = type,color=as.factor(k)), size = 1.5) +
  geom_point(aes(shape = type, color = as.factor(k)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Amount of Noise (SD)") +
  ylab("Coverage") 



df = aggregate(res$length_fission ~ res$k + res$sd + res$pct_active, FUN = mean)
df2 = aggregate(res$length_full ~ res$k + res$sd + res$pct_active, FUN = mean)
df <- cbind(df,"Fission")
df2 <- cbind(df2,"Full")
colnames(df) <- c("k","sd","pct_active","CI","type")
colnames(df2) <- c("k","sd","pct_active","CI","type")
df <- rbind(df,df2)


length_active_plot <- df[df$sd==1,] %>%
  ggplot( aes(x=pct_active, y=CI,color=as.factor(k))) +
  geom_line(aes(linetype = type,color=as.factor(k)), size = 1.5) +
  geom_point(aes(shape = type, color = as.factor(k)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Percentage of Nodes Active") +
  ylab("CI Length") 

length_sd_plot <- df[df$pct_active==0.5,] %>%
  ggplot( aes(x=sd, y=CI,color=as.factor(k))) +
  geom_line(aes(linetype = type,color=as.factor(k)), size = 1.5) +
  geom_point(aes(shape = type, color = as.factor(k)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "bottom", legend.text = element_text(size = 15)) +
  xlab("Amount of Noise (SD)") +
  ylab("CI Length") 

ggsave(paste("figures/", "cover_active_plot.pdf",sep=""),cover_active_plot )
ggsave(paste("figures/", "cover_sd_plot .pdf",sep=""),cover_sd_plot )
ggsave(paste("figures/", "length_active_plot.pdf",sep=""),length_active_plot )
ggsave(paste("figures/", "length_sd_plot .pdf",sep=""),length_sd_plot )


