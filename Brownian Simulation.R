
library(Matrix)
library(MASS)
library(glmnet)
library(dplyr)
library(reshape2)
library(ggplot2)


generate_linear = function(n, p, beta, type, rho,sigmasq=1){
  if (type == "independent") {
    X = matrix(rnorm(n = n*p, sd = 1), ncol = p)
    Sigma = diag(rep(1, p))
  } else {
    p_size = p/5
    Sigma = toeplitz(rho^(0:(p_size-1)))
    Sigma = bdiag(lapply(1:5, function(x) Sigma))
    X = mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
  }
  Y = rnorm(n,0,sd=sigmasq) + X%*%beta
  return(list(X = X, Y = Y, Sigma = Sigma))
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




single_run <- function(n,p,beta,alpha,printst="yo",sigmasq=1,max_iter =50, step_size = 1){
  print(printst)
  truth = which(beta!=0)
  
  dat <- generate_linear(n, p, beta, "independent",0.0)
  Z <- BM(max_iter,step_size,1,niter=n,start_vec = 0)
  t<- as.numeric(rownames(Z))
  apply(Z,1,function(x) x + dat$Y)
  f_y = t(apply(Z,1,function(x) x + as.vector(dat$Y)))
  g_y = t(apply(-Z*sigmasq/t,1,function(x) x + as.vector(dat$Y)))
  
  
  lasso_list = list()
  
  for (i in 2:length(t)) {
    lasso <- cv.glmnet(dat$X,as.vector(as.matrix(f_y)[i,]),family="gaussian")
    lambda <- lasso$glmnet.fit$lambda
    coef <- lasso$glmnet.fit$beta
    selection <- which(coef(lasso, s = 'lambda.1se') != 0)[-1] - 1
    
    #selection performance metrics
    power = sum(beta[selection] != 0 )/ sum(beta != 0)
    precision = sum(beta[selection] != 0)/ max(1,length(selection))
    
    if(length(selection) > 0){
      infer_model = lm(g_y[i,] ~0 + dat$X[,selection])
      conf <- confint(infer_model,level= 1-alpha)
      length = conf[,2] - conf[,1]
      cover = (beta[selection] >= conf[,1]) & (beta[selection] <= conf[,2])
      metrics = c(power,precision,mean(conf),1-mean(cover))
    }
    else{
      metrics = c(power,precision,NA,NA)
    }
    ret = list(lambda = lambda, coef = coef, selection = selection,metrics =metrics)
    lasso_list[[as.character(t[i])]] = ret
  }
  return(lasso_list)
}


niter = 100
n=100
p=100
alpha = 0.1
beta <- c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
#beta <- c(1,0,1,1,0,0,0,0,1,0)
truth = which(beta!=0)
max_iter=10
step_size = 0.1
t <- 0:(max_iter/step_size) *step_size# time

trials <- mclapply(1:100, function(i) single_run(n,p,beta,alpha,max_iter=max_iter, step_size = step_size,printst=i),mc.cores=detectCores())

#Si
#trials = list()
#for(i in 1:niter){
#  trials[[i]] = single_run(n,p,beta,0.9,max_iter=max_iter, step_size = step_size,printst=i)
#}

for(j in 1:niter){
  res <- t(sapply(2:length(t), function(i) c(j,t[i],trials[[j]][[as.character(t[i])]]$metrics)))
  if(j == 1){
    results <- res
  }
  else{
    results <- rbind(results,res) 
  }
}
colnames(results) <- c("Trials","Time","Power","Precision","CI Length","FCR")
results <- data.frame(results)


agg <- aggregate(.~Time,data=results[,-1],FUN=mean)
agg <- melt(agg,id="Time") 

agg %>%
  ggplot( aes(x=Time, y=value, group=variable,color=variable)) +
  geom_line(aes(linetype = variable, color = variable), size = 1.5) + 
  geom_hline(yintercept=alpha) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "bottom", legend.text = element_text(size = 15))













get_lasso_path <- function(coef,k) {
  for(i in 2:length(t)){
    s = t[i]
    res = cbind(coef,t[i],trials[[k]][[as.character(s)]]$lambda, trials[[1]][[as.character(s)]]$coef[coef,])
    if(i==2){
      results =res
    }
    else{
      results=rbind(results,res)
    }
  }
  colnames(results) <- c("Coef","Time","lambda","beta")
  return(results)
}

for(k in 1:niter){
  for(j in 1:p){
    if(j*k == 1){
      results = get_lasso_path(j,k)
    }
    else{
      results = rbind(results,get_lasso_path(j))
    }
  }
}
results = data.frame(results)



results[results$beta ==0,]
agg <- aggregate(lambda~Time +Coef,data=results[results$beta ==0,],FUN=min)
agg[,"type"] =agg$Coef %in% truth
agg[agg$type,"Signal_Type"] = "Non-null"
agg[!agg$type,"Signal_Type"] = "Null"

aggregate(lambda~Time+Signal_Type,data=agg,FUN=mean) %>%
  ggplot( aes(x=Time, y=lambda, group=Signal_Type,color=Signal_Type)) +
  geom_line(aes(x=Time, y=lambda, group=Signal_Type,color=Signal_Type), size = 1.5) + 
  xlab("Time") + 
  ylab("lambda_{t}") + 
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "right", legend.text = element_text(size = 15))




#Figure 1 in document
sigmasq=1
dat <- generate_linear(n, p, beta, "independent",0.0)
f_y <- BM(25,50,1,start_vec = as.vector(dat$Y))
t<- as.numeric(rownames(f_y))
g_y <- sapply(1:ncol(f_y), function(j) f_y[1,j]* (sigmasq+t)/t - f_y[,j]*(sigmasq/t))

plot(t, f_y[,1],type="l",xlab="Time",ylab="")
lines(t,g_y[,1],type="l",col="blue")
abline(h=f_y[1,1],type="l",lty=3)
legend(20,1.5,c("f(Y)","g(Y)","True Y"),lty=c(1,1,3),col=c("black","blue","black"))

max_iter = 0.25
step_size = 25
t <- 0:(max_iter*step_size)/step_size