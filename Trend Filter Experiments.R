source("tf_funs.R")
n <- 100
ntrials <-100
slope = 1
prob_knot = 0.05
sd = 1
step_size = 0.1
alpha = 0.1
resdir = "results"

many_trials <- function(n,ntrials,slope,prob_knot,sd,step_size=0.1,alpha=0.1,save=TRUE,max_iter=5){
  trials <- mclapply(1:ntrials,function(x) singlerun_tf(x,n,slope,1-prob_knot,sd,step_size=step_size,max_iter=max_iter,type="CV"),mc.cores=detectCores())
  df <- do.call("rbind", trials)
  df <- apply(df,2,as.numeric)
  df <- data.frame(df)
  
  df$n = n
  df$ntrials = ntrials
  df$slope = slope
  df$prob_knot = prob_knot
  df$sd = sd
  df$step_size = step_size
  df$alpha = alpha
  
  save(df,file=paste(resdir,"/",paste("trials",n,ntrials,slope,prob_knot,sd,step_size,alpha,sep="_"),".Rdata",sep=""))
  return(df)
}


sigma_seq = c(0.05,0.1,0.2,0.5,1,2)
prob_seq = c(0.05,0.1,0.15,0.2)

for (sigma in sigma_seq) {
  print(sigma)
  for (prob in prob_seq) {
    print(prob)
    df <- many_trials(100,500,1,prob,sigma,step_size=0.1*sigma,max_iter=5*sigma)
    
  }
}

