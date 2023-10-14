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

##############NEW BATCH############################
ntrials=500

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=0,cv=0.9,tau=1,scale=3,err_type="normal",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_0_normal.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=0,cv=0.9,tau=1,scale=3,err_type="sn",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_0_sn.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=0,cv=0.9,tau=1,scale=3,err_type="t",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_0_t.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=0,cv=0.9,tau=1,scale=3,err_type="laplace",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_0_laplace.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=1,cv=0.9,tau=1,scale=1,err_type="normal",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_1_normal.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=1,cv=0.9,tau=1,scale=1,err_type="sn",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_1_sn.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=1,cv=0.9,tau=1,scale=1,err_type="t",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_1_t.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=1,cv=0.9,tau=1,scale=1,err_type="laplace",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_1_laplace.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=2,cv=0.9,tau=1,scale=0.2,err_type="normal",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_2_normal.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=2,cv=0.9,tau=1,scale=0.2,err_type="sn",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_2_sn.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=2,cv=0.9,tau=1,scale=0.2,err_type="t",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_2_t.Rdata")

res=mclapply(1:ntrials, function(x) run_experiment3(graph,edges,edge_mat,k=2,cv=0.9,tau=1,scale=0.2,err_type="laplace",est_var=TRUE,K=10),mc.cores=128)
save(res,file="results_graphfission_inf_2_laplace.Rdata")


##############INITIAL BATCH############################
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
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_normal.Rdata")

active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,est_var=TRUE,lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_normal_estvar.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="sn",lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_sn.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="sn",est_var=TRUE,lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_sn_estvar.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_t.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",est_var=TRUE,lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_t_estvar.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_laplace.Rdata")


active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=1,tau=1,err_type="t",est_var=TRUE,lambda_list = c(1:100/10000,2:100/1000,2:10/100)),mc.cores=128)
save(res,file="results_graphfission_crossval_ridge_laplace_estvar.Rdata")
