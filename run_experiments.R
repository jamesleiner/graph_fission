##########################################################################################################
#
# Run experiments to verify efficacy of cross-validation and confidence interval construction procedures
#
##########################################################################################################


library(MASS)
library(matrixcalc)
library(pracma)
library(dplyr)
library(tidyr)
library(plotly)
library(genlasso)
source("funs.R")


################ GENERATE GRAPH AND EDGE LIST ####################################################
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

ntrials=100



################ GENERATE SIMULATIONS FOR ASSESSING CONFIDENCE INTERVAL CONSTRUCTION####################################################
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



################ GENERATE SIMULATIONS FOR ASSESSING CROSS-VALIDATION STRENGTH####################################################
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

res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,scale=2,sd_level=1,tau=1,err_type="t",scale=2),mc.cores=128)
res = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2, by = "i", all.x = TRUE),res)
save(res,file="results_graphfission_crossval_1_t.Rdata")


res=mclapply(1:ntrials, function(x) run_experiment2(graph,edges,edge_mat,k=1,scale=2,sd_level=1,tau=1,err_type="t",est_var=TRUE,scale=2),mc.cores=128)
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



################ RUN EXPERIMENTS FOR L2 PENALTIES####################################################
for(err_type in err_type_list){
  for(sd in c(0.5,1,1.5,2,2.5,3,3.5,4,4.5)){
    for(est_var in c(TRUE,FALSE)){
      file = "results_graphfission_crossval_ridge"
      file = paste(file,err_type,sd,est_var,".Rdata",sep="_")
      print(file)
      res=mclapply(1:ntrials, function(x) run_experiment_ridge(graph,edges,edge_mat,k=0,sd_level=sd,scale=1,tau=1,lambda_list = c(1:100/10000,2:100/1000,2:10/100),est_var=est_var,err_type=err_type),mc.cores=128)
      save(res,file=file)
    }
  }
}

################ RUN EXPERIMENTS FOR POISSON LOSS##############################################
ntrials = 100
for(k in c(0,1,2)){
  for(scale in c(0.25,0.5,1,2,3)){
    file = "results/results_graphfission_poisson"
    file = paste(file,scale,k,".Rdata",sep="_")
    print(file)
    res=mclapply(1:ntrials, function(x) run_experiment_poisson(graph,edges,edge_mat,scale=scale,k=k,K=5),mc.cores=128)
    save(res,file=file)
  }
}

