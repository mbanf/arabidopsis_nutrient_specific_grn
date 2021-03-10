


# compute weighted gene coordination
m.gc <- m.de.root #abs(m.de.comb)
#p.de <- colSums(m.gc) / nrow(m.de.bin)
for(j in 1:ncol(m.gc)){
  m.gc[,j] <- m.gc[,j] #- p.de[j]
}
#m.gc[m.gc < 0] <- 0
m.gc = tcrossprod(m.gc)

m.gc.numbers <- m.de.root
m.gc.numbers = tcrossprod(m.gc.numbers)


v.th.likelihood <- numeric(n.sim)
#pb <- txtProgressBar(min = 0, max = n.sim, style = 3)

strt<-Sys.time()
cl<-makeCluster(min(n.sim,n.cpus))
registerDoParallel(cl)
l.th.likelihood <- foreach(s = 1:n.sim, .packages=c("fume", "reshape2")) %dopar% { 
  
  #for(s in 1:n.sim){ # parallel?
  
  #setTxtProgressBar(pb, s)
  
  set.seed((1234 + 123 * s))
  
  m.sim <- abs(m.de.root) #abs(m.de.comb)
  m.sim <- t(apply(m.sim, 1, sample))
  
  #p.de <- colSums(m.sim) / nrow(m.de.root)
  
  for(j in 1:ncol(m.sim)){
    m.sim[,j] <- m.sim[,j] #- p.de[j]
  }
  
  m.sim[m.sim < 0] <- 0
  m.sim = tcrossprod(m.sim)
  
  th.mr <- 1 - th.prob
  
  # v.th.likelihood[s] <- quantile(m.sim, th.mr)
  
  th.likelihood <- quantile(m.sim, th.mr)
  rm(m.sim)
  
  th.likelihood
}

#close(pb)
stopCluster(cl)
print(Sys.time()-strt)

v.th.likelihood <- unlist(l.th.likelihood)
th.min_overlap <- mean(v.th.likelihood)
####

strt<-Sys.time()
cl<-makeCluster(min(length(tfs),n.cpus))
registerDoParallel(cl)
l.regulatoryNetwork <- foreach(r = 1:length(tfs), .packages=c("fume", "reshape2")) %dopar% { 
  
  tf <- tfs[r]
  
  setTxtProgressBar(pb, r)
  
  df.dna.r <- subset(df.dna, df.dna[,1] == tf)
  tgs <- intersect(as.character(df.dna.r[,2]), v.gns)
  
  df.regulatoryNetwork <- c()
  l.treatments <- list()
  
  if(tf %in% v.gns){
    
    tgs.links <- names(which(m.gc[tf, tgs] >= th.min_overlap))
    
    if(length(tgs.links) > 0){
      
      df.regulatoryNetwork <- subset(df.dna.r, df.dna.r$Target %in% tgs.links)
      
      df.regulatoryNetwork["eta"] <- 0
      df.regulatoryNetwork["p.val"] <- 1
      df.regulatoryNetwork["mechanism"] <- 0
      df.regulatoryNetwork["pcc"] <- 0
      
      for(j in 1:length(tgs.links)){
        
        expr.tf <- m.fc.root[tf, ]
        expr.tf.inv <- - m.fc.root[tf, ] 
        expr.tg <- m.fc.root[tgs.links[j], ]
    
        res_orig <- eta_squared_inference(expr.tf, expr.tg, v.treatment_buildingblocks.root, v.repeats.root) 
        res_inv <- eta_squared_inference(expr.tf.inv, expr.tg, v.treatment_buildingblocks.root, v.repeats.root) 
        
        # use eta to select
        if (res_orig[2] < res_inv[2]) { # if (res_orig[2] < res_inv[2]) {
          
          df.regulatoryNetwork$eta[j] <- res_orig[1]
          df.regulatoryNetwork$p.val[j] <- res_orig[2]
          df.regulatoryNetwork$mechanism[j] <- "activation"
          df.regulatoryNetwork$pcc[j] <- cor(expr.tf, expr.tg)
          
        } else {
          
          df.regulatoryNetwork$eta[j] <- res_inv[1]
          df.regulatoryNetwork$p.val[j] <- res_inv[2]
          df.regulatoryNetwork$mechanism[j] <- "repression"
          df.regulatoryNetwork$pcc[j] <- cor(expr.tf, expr.tg)
          
        }
        
      }
    }
  }
  
  #df.regulatoryNetwork$p.val <- p.adjust(df.regulatoryNetwork$p.val,"fdr") 
  #df.regulatoryNetwork <- subset(df.regulatoryNetwork, df.regulatoryNetwork$p.val <= 0.05)
  df.regulatoryNetwork
}
stopCluster(cl)
print(Sys.time()-strt)

