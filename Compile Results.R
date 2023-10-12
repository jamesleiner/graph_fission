files <- list.files(path="results/",  full.names=TRUE, recursive=FALSE)

i=0
for(file in files){
  load(file)
  if(i==0){
    df_agg = df
  }
  else{
    df_agg = rbind(df_agg,df)
  }
  i=i+1
}

sigma_seq = c(0.05,0.1,0.2,0.5,1,2)
prob_seq = c(0.05,0.1,0.15,0.2)

for (sigma in sigma_seq) {
  print(sigma)
  for (prob in prob_seq) {
    sd_choice = sigma
    prob_choice = prob
    title = paste("SD = ",sd_choice,", Prob. of new knots = ",prob_choice,sep="")
    savefile = paste("figures/",sd_choice,prob_choice,"plot",sep="_")
    
    
    df_summary <- aggregate(. ~ t + sd+ prob_knot, df_agg,mean)
    df_summary$nk_ratio <- df_summary$nk_cond/df_summary$num_knots
    df_summary <- melt(df_summary,id.vars=c("sd","t",'sigmasq',"slope","prob_knot"))
    
    df_subset <- subset(df_summary, (sd==sd_choice)&(prob_knot ==prob_choice)&(variable %in% c("length_est","length_condest")))
    
    df_subset %>%
      ggplot( aes(x=t, y=value, group=variable)) +
      geom_line(aes(linetype = variable, color = variable), size = 1.5) +
      theme(legend.title = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
            panel.grid.minor = element_line(colour = "grey"),
            text = element_text(size = 15),
            legend.position = "bottom", legend.text = element_text(size = 15)) +
      scale_linetype_discrete(labels = c(length_est = "Original Fission", length_condest = "Conditional Fission")) +
      scale_color_discrete(labels = c(length_est = "Original Fission", length_condest = "Conditional Fission"))  +
      xlab("Time")+
      ylab("CI Length") +
      ggtitle(title)
    ggsave(paste("figures/",paste(sd_choice,prob_choice,"CILENGTH.pdf",sep="_"),sep=""))
    
    df_subset <- subset(df_summary, (sd==sd_choice)&(prob_knot ==prob_choice)&(variable %in% c("cover_est","cover_condest")))
    
    df_subset %>%
      ggplot( aes(x=t, y=value, group=variable)) +
      geom_line(aes(linetype = variable, color = variable), size = 1.5) +
      theme(legend.title = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
            panel.grid.minor = element_line(colour = "grey"),
            text = element_text(size = 15),
            legend.position = "bottom", legend.text = element_text(size = 15)) +
      geom_hline(yintercept=0.9)+
      scale_linetype_discrete(labels = c(cover_est = "Original Fission)", cover_condest  = "Conditional")) +
      scale_color_discrete(labels = c(cover_est  = "Original Fission)", cover_condest  = "Conditional")) +
      xlab("Time") +
      ylab("Coverage")+
      ggtitle(title)
    ggsave(paste("figures/",paste(sd_choice,prob_choice,"COVERAGE.pdf",sep="_"),sep=""))
    
    
    df_subset <- subset(df_summary, (sd==sd_choice)&(prob_knot ==prob_choice)&(variable %in% c("nk_ratio","gamma","gamma_est")))
    
    df_subset %>%
      ggplot( aes(x=t, y=value, group=variable)) +
      geom_line(aes(linetype = variable, color = variable), size = 1.5) +
      theme(legend.title = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
            panel.grid.minor = element_line(colour = "grey"),
            text = element_text(size = 15),
            legend.position = "bottom", legend.text = element_text(size = 15)) +
      scale_linetype_discrete(labels = c(gamma = "Gamma", gamma_est = "Gamma (Estimate)",nk_ratio = "Num Knots Selected / Truth")) +
      scale_color_discrete(labels = c(gamma = "Gamma", gamma_est = "Gamma (Estimate)",nk_ratio = "Num Knots Selected / Truth")) +
      xlab("Time") +
      ggtitle(title)
    ggsave(paste("figures/",paste(sd_choice,prob_choice,"GAMMA.pdf",sep="_"),sep=""))
  }
}
