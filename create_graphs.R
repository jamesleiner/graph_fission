library(MASS)
library(matrixcalc)
library(pracma)
library(dplyr)
library(tidyr)
library(plotly)
library(genlasso)
source('funs.R')

k=0


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

k=0
tau=1
active_set = which(apply(edges,1, function(x) get_active(graph,x,x_sep=4))==1)
#active_set = c(5,20,37,43)

dat = generate_data(graph,edges,edge_mat, active_set,scale=10,k=k,sd_level=1,tau=tau,err_type="normal",est_var=FALSE)
res <- get_CI2(dat$f_Y,dat$g_Y,dat$Y,edge_mat,dat$mean,dat$sd_est,cv=0.9,k=k,K=5)
#CI_old <- get_CI(dat$f_Y,dat$g_Y,edge_mat,dat$mean,k=2,sd_level=dat$sd_est,K=5,cv=0.9)


df <- cbind(graph,res$fit_1se[,"proj_mean"]) 
colnames(df) <- c("x","y","truth")
df <- df%>% pivot_wider(names_from =x, values_from = truth)
df2 <- cbind(graph,dat$Y) 
colnames(df2) <- c("x","y","z")
df3 <- cbind(graph,res$fit_1se[,"lwr"])
colnames(df3) <- c("x","y","z")
df3  <- df3 %>% pivot_wider(names_from =x, values_from = z)
df4 <- cbind(graph,res$fit_1se[,"upr"])
colnames(df4) <- c("x","y","z")
df4  <- df4 %>% pivot_wider(names_from =x, values_from = z)

plot_ly(showscale=FALSE,showlegend=TRUE) %>%  
  add_trace(data = df2, x = ~x, y = ~y, z = ~z, name="Datapoint", mode = "markers", type = "scatter3d", marker = list(size = 5, color = "black", symbol = 104)) %>% 
  add_surface(z=as.matrix(df)[,-1],opacity=1,name="Trend", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df3)[,-1],opacity=0.2,name = "Lower CI", colorscale = list(c(0, 1), c("red", "red")))%>% 
  add_surface(z=as.matrix(df4)[,-1],opacity=0.2,name = "Upper CI",colorscale = list(c(0, 1), c("red", "red"))) 


df <- cbind(graph,res$fit_1se[,"proj_mean"]) 
colnames(df) <- c("x","y","truth")
df <- df%>% pivot_wider(names_from =x, values_from = truth)
df2 <- cbind(graph,dat$Y) 
colnames(df2) <- c("x","y","z")
df3 <- cbind(graph,res$fit_1se[,"lwr_rob"])
colnames(df3) <- c("x","y","z")
df3  <- df3 %>% pivot_wider(names_from =x, values_from = z)
df4 <- cbind(graph,res$fit_1se[,"upr_rob"])
colnames(df4) <- c("x","y","z")
df4 <- df4 %>% pivot_wider(names_from =x, values_from = z)

plot_ly(showscale=FALSE,showlegend=TRUE) %>%  
  add_trace(data = df2, x = ~x, y = ~y, z = ~z, name="Datapoint", mode = "markers", type = "scatter3d", marker = list(size = 5, color = "black", symbol = 104)) %>% 
  add_surface(z=as.matrix(df)[,-1],opacity=1,name="Upper CI", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df3)[,-1],opacity=0.2,name = "Lower CI", colorscale = list(c(0, 1), c("blue", "blue")))%>% 
  add_surface(z=as.matrix(df4)[,-1],opacity=0.2,name = "Upper CI",colorscale = list(c(0, 1), c("blue", "blue")))



##########################POISSON GRAPH############################################################
scale_list = c(0.25,0.5,1,2,3)
k_list = c(0,1,2)
first=TRUE

for(k in k_list) {
  for(scale in scale_list){
    file = "results/results_graphfission_poisson"
    file = paste(file,scale,k,".Rdata",sep="_")
    print(file)
    if(file.exists(file)){
      load(file)
      for(i in 1:length(res)){
        if(is.list(res[[i]])){
          temp = data.frame(res[[i]]$metrics)
          temp$k = k
          temp$scale = scale
          if(first==TRUE){
            results = temp
            first=FALSE
          }
            else{
              results = rbind(results,temp)
              }
          }
        }
      }
  }
}


colnames(results) <- c('df','cv','cv','length','length_r','err_true','err_proj','k','scale')
length = aggregate(results$length ~results$scale + results$k , FUN = mean)
results = results[!is.infinite(rowSums(results)),]
cv = aggregate(results$cv ~results$scale + results$k , FUN = mean)



colnames(cv) <- c('scale','k','cv')
colnames(length) <- c('scale','k','length')
plt1 <- cv%>%
  ggplot( aes(x=scale, y=cv,color=as.factor(k),linetype=as.factor(k))) +
  geom_line(aes(color=as.factor(k),linetype=as.factor(k)), linewidth = 1) +
  geom_point(aes(color = as.factor(k),shape=as.factor(k)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("Scale") +
  geom_hline(yintercept=0.8, linetype="dashed",color = "black", size=1.5) + 
  ylab("Coverage")

ggsave(paste("figures/poisson_cv.pdf",sep=""),plt1)


plt2 <- length%>%
  ggplot( aes(x=scale, y=length,color=as.factor(k),linetype=as.factor(k))) +
  geom_line(aes(color=as.factor(k),linetype=as.factor(k)), linewidth = 1) +
  geom_point(aes(color = as.factor(k),shape=as.factor(k)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("Scale") +
  ylab("CI Length") 

ggsave(paste("figures/poisson_length.pdf",sep=""),plt2)


##########################RIDGE CV GRAPHS############################################################
err_type_list = c("normal","t","sn","laplace")
sd_list = c(0.5,1,1.5,2.5,3.5,4.5)
first=TRUE

for(sd in sd_list){
  for(est_var in c(TRUE,FALSE)){
    for(err in err_type_list){
      file = "results_graphfission_crossval_ridge"
      file = paste(file,err,sd,est_var,".Rdata",sep="_")
      if(file.exists(file)){
        load(file)
        for(i in 1:length(res)){
          if(is.data.frame(res[[i]])){
            rownames(res[[i]]) <- 1:nrow(res[[i]])
            temp <- data.frame(res[[i]])
            temp$method = temp$type
            temp$errtype =err
            temp$estvar = est_var
            if(first==TRUE){
              results = temp
              first=FALSE
            }
            else{
              results = rbind(results,temp)
            } 
          }
        }
      }
    
    }
  }
}

results[results$method == "1se","method"] = "1 SE rule"
results[results$method == "min","method"] = "Min CV error"
results[results$errtype== "normal","errtype"] = "Normal"
results[results$errtype== "sn","errtype"] = "Skew-Normal"
results[results$errtype== "laplace","errtype"] = "Laplace"
results[results$errtype== "t","errtype"] = "t"

df = aggregate(results$df ~results$K + results$sd_level +results$method +results$errtype+results$estvar , FUN = mean)
err = aggregate(results$errtrue ~ results$K + results$sd_level+results$method + results$errtype+results$estvar , FUN = mean)

colnames(df) <- c('K',"sd","method","errtype","estvar",'df')
colnames(err) <- c('K',"sd","method","errtype","estvar",'err')

plt1<- subset(err,(errtype=="Normal"&sd==2.5&estvar==FALSE))%>%
  ggplot( aes(x=K, y=err,color=as.factor(method),linetype=as.factor(method))) +
  geom_line(aes(color=as.factor(method),linetype=as.factor(method)), linewidth = 1) +
  geom_point(aes(color = as.factor(method),shape=as.factor(method)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("K") +
  ylab("Error")

ggsave(paste("figures/error_ridge_cv.pdf",sep=""),plt1)


plt2<- subset(df,(errtype=="Normal"&sd==2.5&estvar==FALSE))%>%
  ggplot( aes(x=K, y=df,color=as.factor(method),linetype=as.factor(method))) +
  geom_line(aes(color=as.factor(method),linetype=as.factor(method)), linewidth = 1) +
  geom_point(aes(color = as.factor(method),shape=as.factor(method)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("K") +
  ylab("Degrees of freedom")

ggsave(paste("figures/df_ridge_cv.pdf",sep=""),plt2)





filename = paste("results_graphfission_crossval_ridge","_",err,".Rdata",sep="")
load(filename)

##########################CONFIDENCE INTERVAL GRAPHS############################################################
err_type_list = c("normal","t","sn","laplace")
first=TRUE
for(k in c(0,1)){
    for(err in err_type_list){
      print(filename)
      filename = paste("results/results_graphfission_inf_",k,"_",err,".Rdata",sep="")
      load(filename)
      for(i in 1:length(res)){
        if(is.matrix(res[[i]])){
          method <- rownames(res[[i]])
          rownames(res[[i]]) <- 1:nrow(res[[i]])
          temp <- data.frame(res[[i]])
          temp$method = method
          temp$errtype =err
          if(first==TRUE){
            results = temp
            first=FALSE
          }
          else{
            results = rbind(results,temp)
          } 
        }
      }
  }
}


results[results$errtype== "normal","errtype"] = "Normal"
results[results$errtype== "sn","errtype"] = "Skew-Normal"
results[results$errtype== "laplace","errtype"] = "Laplace"
results[results$errtype== "t","errtype"] = "t"

summary1 = aggregate(results$length ~ results$sd + results$K + results$method + results$k +results$errtype , FUN = mean)
summary2 = aggregate(results$length_robust ~ results$sd + results$K + results$method + results$k +results$errtype , FUN = mean)
summary3 = aggregate(results$coverage ~ results$sd +  results$K + results$method + results$k +results$errtype , FUN = mean)
summary4 = aggregate(results$coverage_robust ~ results$sd + results$K + results$method + results$k +results$errtype , FUN = mean)
summary1$CI_type = "Naive"
summary2$CI_type = "Theorem 1"
summary3$CI_type = "Naive"
summary4$CI_type = "Theorem 1"



colnames(summary1) <- c('sd','K',"method","k","err","length","CI_type")
colnames(summary2) <- c('sd','K',"method","k","err","length","CI_type")
colnames(summary3) <- c('sd','K',"method","k","err","cover","CI_type")
colnames(summary4) <- c('sd','K',"method","k","err","cover","CI_type")

plt1<- subset(rbind(summary1,summary2),(method=="1 SE Rule"&k==1))%>%
  ggplot( aes(x=sd, y=length,color=as.factor(err),linetype=as.factor(CI_type))) +
  geom_line(aes(linetype=as.factor(CI_type),color=as.factor(err)), linewidth = 1) +
  geom_point(aes(color = as.factor(err),shape=as.factor(err)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("Error SD") +
  ylab("CI Length")

ggsave(paste("figures/CI_inf.pdf",sep=""),plt1)




plt2<-subset(rbind(summary3,summary4),(method=="1 SE Rule"&k==1))%>%
  ggplot( aes(x=sd, y=cover,color=as.factor(err),linetype=as.factor(CI_type))) +
  geom_line(aes(linetype=as.factor(CI_type),color=as.factor(err)), linewidth = 1) +
  geom_point(aes(color = as.factor(err),shape=as.factor(err)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("Error SD") +
  ylab("Coverage") + 
  geom_hline(yintercept=0.9,color = "GRAY", size=1.5)

ggsave(paste("figures/cov_inf.pdf",sep=""),plt2)


summary1 = aggregate(results$sd_low ~ results$sd + results$K + results$method + results$k +results$errtype , FUN = mean)
summary2 = aggregate(results$sd_high ~ results$sd + results$K + results$method + results$k +results$errtype , FUN = mean)
summary1$type = "Low"
summary2$type = "High"
colnames(summary1) <- c('sd','K',"method","k","err","sd_est","type")
colnames(summary2) <- c('sd','K',"method","k","err","sd_est","type")


plt3<- subset(rbind(summary1,summary2),(method=="1 SE Rule"&k==1))%>%
  ggplot( aes(x=sd, y=sd_est,color=as.factor(err),shape=as.factor(type))) +
  geom_line(aes(color=as.factor(err)), linewidth = 1) +
  geom_point(aes(color = as.factor(err),shape=as.factor(type)), size = 20) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 25),legend.key.size = unit(1, 'cm'),legend.key = element_rect(fill = "white")) +
  scale_shape_manual(values = c("+", "-")) + 
  xlab("Acutal SD") +
  ylab("Estimated SD")+
  geom_abline(intercept = 0, slope = 1, size = 0.5,linetype="dashed")

ggsave(paste("figures/sd_inf.pdf",sep=""),plt3)



###########OLD CODE#########################
err_type_list = c("normal","t","sn","laplace")

first=TRUE
for(k in c(0,1)){
  for(est_var in c(TRUE,FALSE)){
    for(err in err_type_list){
      filename = paste("results/results_graphfission_crossval_",k,"_",err,sep="")
      if(est_var){
        filename= paste(filename,"_estvar.Rdata",sep="")
      }
      if(!est_var){
        filename= paste(filename,".Rdata",sep="")
      }
      load(filename)
      res$err = err
      res$est_var = est_var
      res$k = k
      
      if(first==TRUE){
        res_comb = res
        first = FALSE
      }
      else{
        res_comb = rbind(res_comb,res)
      }
    }
  }
}
res = res_comb
res[res$method == "1se","method"] = "1se rule"
res[res$method == "min","method"] = "Min CV error"
res[res$err== "normal","err"] = "Normal"
res[res$err== "sn","err"] = "Skew-Normal"
res[res$err== "laplace","err"] = "Laplace"
res[res$err== "t","err"] = "t"
res$K = as.numeric(res$K)
res$df = as.numeric(res$df)
res$errtrue = as.numeric(res$errtrue)
res$errproj = as.numeric(res$errproj)

for(k_choice in c(0,1)){
  for(est_var_choice in c(TRUE,FALSE)){
    for(rule in c("1se rule","Min CV error","SURE")){
      
      summary = aggregate(res$df ~ res$K + res$method + res$split + res$err +res$k +res$est_var , FUN = mean)
      colnames(summary) <- c('K',"method","split",'err','k','est_var','df')
      
      if(est_var_choice == TRUE){
        title = paste("k = ",k_choice,", estimated variance, ",rule,sep="")
      }
      else{
        title = paste("k = ",k_choice,", known variance, ",rule,sep="")
      }
      
      df=subset(summary,(k==k_choice)&(est_var==est_var_choice)&(method==rule)&(split=="fission"))
      cv_df <- df%>%
        ggplot( aes(x=K, y=df,linetype=as.factor(err))) +
        geom_line(aes(linetype=as.factor(err),color=as.factor(err)), linewidth = 1) +
        geom_point(aes(color = as.factor(err),shape=as.factor(err)), size = 4) +
        theme(legend.title = element_blank(),
              panel.background = element_rect(fill = "white", colour = "black"),
              panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
              panel.grid.minor = element_line(colour = "grey"),
              text = element_text(size = 25),
              legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
        xlab("Number of folds (K)") +
        ylab("Degrees of freedom") +
        ggtitle(title) +
        geom_hline(yintercept=4, linetype="dashed",
                   color = "black", size=1.5)
      
      ggsave(paste("figures/",k_choice,"_",est_var_choice,"_",rule,"_df.pdf",sep=""),cv_df)
      
      summary = aggregate(res$errtrue ~ res$K + res$method + res$split + res$err +res$k +res$est_var , FUN = mean)
      colnames(summary) <- c('K',"method","split",'err','k','est_var','errtrue')
      
      
      df=subset(summary,(k==k_choice)&(est_var==est_var_choice)&(method==rule)&(split=="fission"))
      err_df <- df%>%
        ggplot( aes(x=K, y=errtrue,linetype=as.factor(err))) +
        geom_line(aes(linetype=as.factor(err),color=as.factor(err)), linewidth = 1) +
        geom_point(aes(color = as.factor(err),shape=as.factor(err)), size = 4) +
        theme(legend.title = element_blank(),
              panel.background = element_rect(fill = "white", colour = "black"),
              panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
              panel.grid.minor = element_line(colour = "grey"),
              text = element_text(size = 25),
              legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
        xlab("Number of folds (K)") +
        ylab("Error") +
        ylim(c(0,max(df$errtrue)*1.15)) + 
        ggtitle(title) 
      
      ggsave(paste("figures/", k_choice,"_",est_var_choice,"_",rule,"_err.pdf",sep=""),err_df)
    }
  }
}


summary = aggregate(res$errtrue ~ res$K + res$method + res$split + res$err +res$k +res$est_var , FUN = mean)
colnames(summary) <- c('K',"method","split",'err','k','est_var','errtrue')

df=subset(summary,(k==0)&(est_var==TRUE)&(method==rule)&(split=="fission"))

df=subset(summary,(k==0)&(est_var==TRUE)&(method==rule)&(split=="fission"))
err_df <- err%>%
  ggplot( aes(x=K, y=errtrue,linetype=as.factor(err))) +
  geom_line(aes(linetype=as.factor(err),color=as.factor(err)), linewidth = 1) +
  geom_point(aes(color = as.factor(err),shape=as.factor(err)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "bottom", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("Number of folds (K)") +
  ylab("Error") +
  ggtitle(rule) +
  geom_hline(yintercept=4, linetype="dashed",
             color = "black", size=1.5)

ggsave(paste("figures/", "err","_",k_choice,"_",est_var_choice,"_",rule,".pdf",sep=""),err_df)




df = aggregate(as.numeric(res$df) ~ res$K + res$method + res$split , FUN = mean)
colnames(df) <- c('K',"method","split",'df')
cv_df <- df[df$split=='fission',]%>%
  ggplot( aes(x=as.numeric(K), y=df,color=as.factor(method))) +
  geom_line(aes(color=as.factor(method),linetype=as.factor(method)), size = 1.5) +
  geom_point(aes(shape =  as.factor(method), color = as.factor(method)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white"))+
  xlab("Number of folds (K)") +
  ylab("Degrees of freedom") +
  geom_hline(yintercept=4, linetype="dashed",
             color = "black", size=1.5)

ggsave(paste("figures/","cv_df.pdf",sep=""),cv_df)


df = aggregate(as.numeric(res$errtrue) ~ res$K + res$method + res$split , FUN = mean)
colnames(df) <- c('K',"method","split",'df')
err_df <- df[df$split=='fission',]%>%
  ggplot( aes(x=as.numeric(K), y=df,color=as.factor(method))) +
  geom_line(aes(color=as.factor(method),linetype=as.factor(method)), size = 1.5) +
  geom_point(aes(shape =  as.factor(method), color = as.factor(method)), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white"))+
  xlab("Number of folds (K)") +
  ylab("Error") 


ggsave(paste("figures/","cv_error.pdf",sep=""),err_df)

create_graphs <- function(k, err_type, label){
  filename = paste("results_graphfission_crossval_",k,"_",err_type,".Rdata",sep="")
  load(filename)
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
    ggtitle(label) +
    geom_hline(yintercept=4, linetype="dashed",
               color = "black", size=1.5)
  
  
  ggsave(paste("figures/", "cv","_",k,"_",err_type,"_df.pdf",sep=""),cv_df)
  
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
          legend.position = "none", legend.text = element_text(size = 15)) +
    xlab("Number of folds (K)") +
    ylab("Error") +
    ggtitle(label)
  ggsave(paste("figures/", "cv","_",k,"_",err_type,"_error.pdf",sep=""),cv_err)
  
  
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
    ylab("CI Length") +
    ggtitle(label)
  ggsave(paste("figures/", "cv","_",k,"_",err_type,"_length.pdf",sep=""),cv_length)
  
  
}

create_graphs(0,"normal","k=0, Gaussian Errors (w/ known variance)")
create_graphs(0,"normal_estvar","k=0, Gaussian Errors (w/ estimated variance)")
create_graphs(0,"laplace","k=0, Laplace Errors")
create_graphs(0,"t","k=0, t-distributed Errors")
create_graphs(0,"sn","k=0, Skewed Normal Errors")

create_graphs(0,"sn_estvar","k=0, Skewed Normal Errors")
create_graphs(1,"normal","k=1, Gaussian Errors (w/ known variance)")
create_graphs(1,"normal_estvar","k=1, Gaussian Errors (w/ estimated variance)")
create_graphs(1,"laplace","k=1, Laplace Errors")
create_graphs(1,"t","k=1, t-distributed Errors")
create_graphs(1,"sn","k=1, Skewed Normal Errors")
create_graphs(2,"normal","k=2, Gaussian Errors (w/ known variance)")
create_graphs(2,"normal_estvar","k=2, Gaussian Errors (w/ estimated variance)")
create_graphs(2,"laplace","k=2, Laplace Errors")
create_graphs(2,"t","k=2, t-distributed Errors")
create_graphs(2,"sn","k=2, Skewed Normal Errors")
library(MASS)
library(matrixcalc)
library(pracma)
library(dplyr)
library(tidyr)
library(plotly)
library(genlasso)


err_type_list = c("normal","t","sn","laplace")
first=TRUE
for(k in c(0,1)){
  for(err in err_type_list){
    filename = paste("results_graphfission_inf_",k,"_",err,".Rdata",sep="")
    load(filename)
    for(i in 1:length(res)){
      if(is.matrix(res[[i]])){
        method <- rownames(res[[i]])
        rownames(res[[i]]) <- 1:nrow(res[[i]])
        temp <- data.frame(res[[i]])
        temp$method = method
        temp$errtype =err
        if(first==TRUE){
          results = temp
          first=FALSE
        }
        else{
          results = rbind(results,temp)
        } 
      }
    }
  }
}


results[results$errtype== "normal","errtype"] = "Normal"
results[results$errtype== "sn","errtype"] = "Skew-Normal"
results[results$errtype== "laplace","errtype"] = "Laplace"
results[results$errtype== "t","errtype"] = "t"

summary1 = aggregate(results$length ~ results$sd + results$K + results$method + results$k +results$errtype , FUN = mean)
summary2 = aggregate(results$length_robust ~ results$sd + results$K + results$method + results$k +results$errtype , FUN = mean)
summary3 = aggregate(results$coverage ~ results$sd +  results$K + results$method + results$k +results$errtype , FUN = mean)
summary4 = aggregate(results$coverage_robust ~ results$sd + results$K + results$method + results$k +results$errtype , FUN = mean)
summary1$CI_type = "Naive"
summary2$CI_type = "Robust"
summary3$CI_type = "Naive"
summary4$CI_type = "Robust"



colnames(summary1) <- c('sd','K',"method","k","err","length","CI_type")
colnames(summary2) <- c('sd','K',"method","k","err","length","CI_type")
colnames(summary3) <- c('sd','K',"method","k","err","cover","CI_type")
colnames(summary4) <- c('sd','K',"method","k","err","cover","CI_type")

plt1<- subset(rbind(summary1,summary2),(method=="1 SE Rule"&k==0))%>%
  ggplot( aes(x=sd, y=length,color=as.factor(err),linetype=as.factor(CI_type))) +
  geom_line(aes(linetype=as.factor(CI_type),color=as.factor(err)), linewidth = 1) +
  geom_point(aes(color = as.factor(err),shape=as.factor(err)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "bottom", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("Error SD") +
  ylab("CI Length")

ggsave(paste("figures/CI_inf.pdf",sep=""),plt1)




plt2<-subset(rbind(summary3,summary4),(method=="1 SE Rule"&k==0))%>%
  ggplot( aes(x=sd, y=cover,color=as.factor(err),linetype=as.factor(CI_type))) +
  geom_line(aes(linetype=as.factor(CI_type),color=as.factor(err)), linewidth = 1) +
  geom_point(aes(color = as.factor(err),shape=as.factor(err)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("Error SD") +
  ylab("Coverage") + 
  geom_hline(yintercept=0.9, linetype="dashed",color = "black", size=1.5)

ggsave(paste("figures/cov_inf.pdf",sep=""),plt2)

df=subset(summary,(k==k_choice)&(est_var==est_var_choice)&(method==rule)&(split=="fission"))
cv_df <- df%>%
  ggplot( aes(x=K, y=df,linetype=as.factor(err))) +
  geom_line(aes(linetype=as.factor(err),color=as.factor(err)), linewidth = 1) +
  geom_point(aes(color = as.factor(err),shape=as.factor(err)), size = 4) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 25),
        legend.position = "none", legend.text = element_text(size = 15),legend.key.size = unit(2, 'cm'),legend.key = element_rect(fill = "white")) +
  xlab("Number of folds (K)") +
  ylab("Degrees of freedom") +
  ggtitle(title) +
  geom_hline(yintercept=4, linetype="dashed",
             color = "black", size=1.5)
