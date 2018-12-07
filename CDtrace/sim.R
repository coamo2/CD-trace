rm(list=ls())
path = 'C:\\Users\\quans\\Desktop\\sim-1206'
setwd(path)


###source files in SpiecEasi-master
###This package is available at https://github.com/zdk123/SpiecEasi by Kurtz et al. (2015)
setwd('SpiecEasi-master\\R')
for(file in list.files()){
  source(file)
}

###source code for gcoda which is available from https://github.com/huayingfang/gCoda
setwd(path)
source("gcoda.R")  

###source code for my method dtrace
source("mycd.R")


###function for AUC computation
ROC_PR = function(pred, true_mat){
  pred = as.matrix(pred)
  TP = sum((pred[upper.tri(pred)]!=0)*(true_mat[upper.tri(pred)]!=0))
  FP = sum((pred[upper.tri(pred)]!=0)*(true_mat[upper.tri(pred)]==0))
  FN = sum((pred[upper.tri(pred)]==0)*(true_mat[upper.tri(pred)]!=0))
  TN = sum((pred[upper.tri(pred)]==0)*(true_mat[upper.tri(pred)]==0))
  TPR = TP/(TP+ FN)
  FPR= FP/(FP+TN)
  precision = TP/(TP + FP)
  recall = TP/(TP + FN)
  return(c(TPR, FPR, precision, recall))
}
summary_ROC_PR = function(path, true_mat){
  res = data.frame(t(sapply(path, function(pred)ROC_PR(pred, true_mat))))
  colnames(res) = c('TPR', 'FPR', 'precision', 'recall')
  res[is.na(res)] = 1
  x = c(0, res$FPR, 1)
  y = c(0, res$TPR, 1)
  roc = sum((y[2:length(y)] + y[1:(length(y)-1)])*(x[2:length(x)] - x[1:(length(x)-1)])/2)
  return(list(res, roc))
}


require(MASS)
set.seed(10)       #seed
p <- 50            #dimension
Times <- 100       #repeat times
Nlambda <- 20      #Number of lambda

#n = 500
#methodtxt<-"block"
for(methodtxt in c("scale_free","band","block")){
  set.seed(666)
  e <- 3*p                               #the number of edges
  graph <- make_graph(methodtxt, p, e)   #make graph
  Omega  <- graph2prec(graph)            #precision matrix
  Sigma <- solve(Omega)                  #The covariance matrix
  require(igraph)
  pdf(paste0(methodtxt, '-network.pdf'))
  gplot = graph_from_adjacency_matrix(graph, mode = 'undirected')
  plot(gplot, vertex.size = 3, vertex.label = NA, vertex.color = 'black',
       edge.color = 'gray', edge.arrow.size = 0.1,
       layout = layout_nicely(gplot))
  dev.off()
  
  for(n in c(100,200,500)){
    #repeat each situation
    mb_list = list()
    gl_list = list()
    gc_list = list()
    dt_list = list()
    roc = NULL
    time = 1
    for(time in 1:Times){
      Y <- exp(mvrnorm(n,rep(0,p),Sigma))    #simulated samples from multivariate normal
      Y <- Y/rowSums(Y)                      #normalization to get the compositional data
      
      #data after clr transformation
      data.clr <- log(Y)%*%(diag(p)-matrix(1,p,p)/p)
      se.mb <- sparseiCov(data.clr, method='mb', lambda.min.ratio=0.1, nlambda=Nlambda)
      se.gl <- sparseiCov(data.clr, method='glasso', lambda.min.ratio=0.1, nlambda=Nlambda)
      
      gc <- gcoda(x = data.clr, counts = F, nlambda=Nlambda, 
                  pseudo=1, lambda.min.ratio = 0.001)
      dc <- CompDtrace(clr.data = data.clr, lambda = NULL, rho = 0.8,
                       nlambda = Nlambda, lambda.min.ratio = 0.001)
      
      #compute roc
      mb_res <- summary_ROC_PR(se.mb$path, Omega)
      gl_res <- summary_ROC_PR(se.gl$path, Omega)
      gc_res <- summary_ROC_PR(gc$path, Omega)
      dt_res <- summary_ROC_PR(dc$path, Omega)
      
      #combine results
      roc = rbind(roc, c(mb_res[[2]], gl_res[[2]], gc_res[[2]], dt_res[[2]]))
      mb_list[[time]] = mb_res[[1]]
      gl_list[[time]] = gl_res[[1]]
      gc_list[[time]] = gc_res[[1]]
      dt_list[[time]] = dt_res[[1]]
      
      #plot(mb_res[[1]]$FPR,mb_res[[1]]$TPR, type='l')
      # lines(dt_res[[1]]$FPR,dt_res[[1]]$TPR, type='l')
      # lines(gc_res[[1]]$FPR,gc_res[[1]]$TPR, type='l',col='red')
      # roc
      cat(n,'   ', methodtxt,'   ', time, '\n')
    }
    #save results
    save.image(file = paste0(methodtxt,n,".rdata"))
  }
}




################summary the results: figures and tables#####
rm(list=ls())
path = 'C:\\Users\\asus\\Desktop\\sim'
setwd(path)


#methodtxt = "scale_free"
#n = 100
all_mean_auc = NULL
all_sd_auc = NULL
for(methodtxt in c("band","scale_free","block")){
  for(n in c(100,200,500)){
    library(MASS)
    load(paste0(methodtxt,n,".rdata"))
    
    #mean and sd of AUC
    mean_auc = data.frame(t(apply(roc,2,mean)), graph = methodtxt, n = n)
    sd_auc = data.frame(t(apply(roc,2,sd)), graph = methodtxt, n = n)
    all_mean_auc = rbind(all_mean_auc, mean_auc)
    all_sd_auc = rbind(all_sd_auc, sd_auc)
    
    #figures
    mb_TPR = sapply(mb_list, function(x)x$TPR)
    mb_FPR = sapply(mb_list, function(x)x$FPR)
    mb_TPR = c(0,rowMeans(mb_TPR),1)
    mb_FPR = c(0,rowMeans(mb_FPR),1)
    
    gl_TPR = sapply(gl_list, function(x)x$TPR)
    gl_FPR = sapply(gl_list, function(x)x$FPR)
    gl_TPR = c(0,rowMeans(gl_TPR),1)
    gl_FPR = c(0,rowMeans(gl_FPR),1)
    
    gc_TPR = sapply(gc_list, function(x)x$TPR)
    gc_FPR = sapply(gc_list, function(x)x$FPR)
    gc_TPR = c(0,rowMeans(gc_TPR),1)
    gc_FPR = c(0,rowMeans(gc_FPR),1)
    
    dt_TPR = sapply(dt_list, function(x)x$TPR)
    dt_FPR = sapply(dt_list, function(x)x$FPR)
    dt_TPR = c(0,rowMeans(dt_TPR),1)
    dt_FPR = c(0,rowMeans(dt_FPR),1)
    
    
    #figures for ROC curve
    pdf(paste(n,methodtxt,"ROC.pdf",sep="-"))
    op <- par(mar = c(5,5,1,1)) 
    plot(dt_FPR,dt_TPR,type="b",lty=1,pch = 1,xlab="1-TN",ylab="TP",cex.axis=1.5,cex.lab=1.5,ylim=c(0,1))
    lines(gc_FPR,gc_TPR,type="b",lty=2,pch = 2)
    lines(gl_FPR,gl_TPR,type="b",lty=3,pch = 3)
    lines(mb_FPR,mb_TPR,type="b",lty=4,pch = 4)
    legend("bottomright",legend=c("CD-trace","gCoda","S-E(glasso)","S-E(mb)"),lty=1:4,pch = 1:4, cex=1.5)
    par(op)
    dev.off()
  }
}

all_mean_auc[,1:4] = all_mean_auc[,4:1]
all_sd_auc[,1:4] = all_sd_auc[,4:1]
colnames(all_mean_auc)[1:4] = c("CD-trace","gCoda","S-E(glasso)","S-E(mb)")
colnames(all_sd_auc)[1:4] = c("CD-trace","gCoda","S-E(glasso)","S-E(mb)")
write.csv(all_mean_auc, "all_mean_auc.csv")
write.csv(all_sd_auc, "all_sd_auc.csv")
  
#barplot for AUC under different situations
a = all_mean_auc[all_mean_auc$graph=='band',]
b = all_mean_auc[all_mean_auc$graph=='block',]
c = all_mean_auc[all_mean_auc$graph=='scale_free',]
a = t(a[,1:4])
colnames(a)<-c("n=100","n=200","n=500")
pdf("band-ROCbar.pdf")
barplot(a-0.6,beside=T,ylim=c(0,0.4),main="AUC under ROC curves",yaxt='n',
        legend.text=TRUE, 
        args.legend=list(x=15, horiz=TRUE,bty='n',y.intersp=58,x.intersp=1))
axis(side=2,at=c(0.60,0.70,0.80,0.90,1.00)-0.6,labels=c('0.6','0.7','0.8','0.9','1.0'))
dev.off()

b = t(b[,1:4])
colnames(b)<-c("n=100","n=200","n=500")
pdf("block-ROCbar.pdf")
barplot(b-0.6,beside=T,ylim=c(0,0.4),main="AUC under ROC curves",yaxt='n',
        legend.text=TRUE, 
        args.legend=list(x=15, horiz=TRUE,bty='n',y.intersp=58,x.intersp=1))
axis(side=2,at=c(0.60,0.70,0.80,0.90,1.00)-0.6,labels=c('0.6','0.7','0.8','0.9','1.0'))
dev.off()

c = t(c[,1:4])
colnames(a)<-c("n=100","n=200","n=500")
pdf("scale-free-ROCbar.pdf")
barplot(c-0.6,beside=T,ylim=c(0,0.4),main="AUC under ROC curves",yaxt='n',
        legend.text=TRUE, 
        args.legend=list(x=15, horiz=TRUE,bty='n',y.intersp=58,x.intersp=1))
axis(side=2,at=c(0.60,0.70,0.80,0.90,1.00)-0.6,labels=c('0.6','0.7','0.8','0.9','1.0'))
dev.off()
