library(MASS)

setwd("/home/youri/Desktop/protocol_hg19/rna/fusion/0_fuma-prostate/extract_breakpoints_in_STAR_fusion_results")


plot_gene <- function(frequencies, frequencies_windowed, entropy, intron_mask) {
  data = read.delim(frequencies)
  #data_w = read.table(frequencies_windowed)
  data_e = read.table(entropy)
  data_i = read.table(intron_mask)
  
  weight_s = median(data$s[data$s != 0])
  weight_e = median(data$e[data$e != 0])
  
  n = length(data$e)
  x = 0:(n-1)
  ymax = log2(max(data$e * data_i$V2,data$s * data_i$V2) + 1)
  
  plot(c(0,n-1),c(0,ymax*2.05),type="n")
  abline(h=ymax*1.05)
  ##------------------------------------------------------------
  ds = c(data$s,data$j)
  #ds2 = ds[ds < 10]
  ds = ds[ds > 0]
  
  #hist(ds,breaks=20)
  fit = fitdistr(ds[ds > 0], 'Poisson')
  simg = dpois(1:10000,fit$estimate)
  simg = simg/max(simg) * 120000*3.6
  lines(1:10000,simg,lwd=2,col="blue")
  #transition1 = qpois(0.950, lambda = fit$estimate[1])
  #transition2 = qpois(0.975, lambda = fit$estimate[1])
  transition3 = qpois(0.990, lambda = fit$estimate[1])
  transition4 = qpois(0.999, lambda = fit$estimate[1])
  transition5 = qpois(0.9999,lambda = fit$estimate[1])
  
  #abline(h=log2(transition1+1),lty=2,col="#880088")
  #abline(h=log2(transition2+1),lty=2,col="#660066")
  abline(h=log2(transition3+1),lty=2,col="#880088")
  abline(h=log2(transition4+1),lty=2,col="#660066")
  abline(h=log2(transition5+1),lty=2,col="#440044")
  
  ## Gamma can't deal with 0 values
  ##fit = fitdistr(ds[ds > 0], 'gamma')
  ##simg = dgamma(1:10000, shape = fit$estimate[1], rate = fit$estimate[2])
  ##simg = simg / max(simg) * 1200
  ##lines(1:10000,simg,lwd=2,col="red") 
  
  fit = fitdistr(ds[ds > 0], 'lognormal')
  #simg = dlnorm(1:10000, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
  #simg = simg / max(simg) * 1200
  #lines(1:10000,simg,lwd=2,col="darkgreen")
  
  #transition1 = qlnorm(0.950, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
  #transition2 = qlnorm(0.975, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
  #transition3 = qlnorm(0.990, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
  
  #abline(h=log2(transition1+1),lty=3,col="#008888")
  #abline(h=log2(transition2+1),lty=3,col="#006666")
  #abline(h=log2(transition3+1),lty=3,col="#004444")
  
## -------------------------------------------------------------
  
  
  #Draw entropy dots first
  points(data_e$V1,(data_e$V2 * ymax) + (ymax*1.05),pch=21,cex=0.001,col="#00000002")

  points(x,log2(data$s+1)*data_i$V2,col="red",pch=21,cex=0.20)
  points(x,log2(data$e+1)/weight_e*weight_s*data_i$V2,col="darkgreen",pch=21,cex=0.20)
  
  # Draw intron/exon boundaries
  abline(v=x[data$j != 0],col="#0000FF55")
  # Draw exon mask
  #points(data_i$V1,(data_i$V2*ymax)+(0.015*ymax),pch=24,cex=0.2)
  
  ## Windows start end plot - not ideal
  ##plot(c(0,11*400),c(0,10000),type="n")
  #x = rep(data_w$V1, times = 1, length.out = NA, each = 2)[2:(length(data_w$V1)*2)]
  #y = rep(data_w$V3, times = 1, length.out = NA, each = 2)[1:(length(data_w$V1)*2)-1]
  ##lines(data_w$V1,data_w$V3)
  #lines(x,log(y),col="red")
  ##abline(v=(0:10)*400,col="gray")
}



samples=c('7046-004-041','7046-004-043')
gene_names=c('erg','tmprss2','ar','sash1','samd5','klk3','gapdh')
#samples=c('7046-004-041')
#gene_names=c('tmprss2')

for(sample in samples) {
  print(sample)
  for(gene in gene_names) {
    print(gene)
    data = paste("data/",sample,"-",gene,"-",sep="")
    figure = paste("plots/",sample,"_",gene,".png",sep="")
    
    v1 = paste(data,"vec1_a.tabular.txt",sep="")
    v2 = paste(data,"vec2.tabular.txt",sep="")
    v3 = paste(data,"vec3.tabular.txt",sep="")
    
    png(figure,width=1250,height=600)
    plot_gene(v1,FALSE,v2,v3)
    dev.off()
    
    #png("043-erg.png",width=1250,height=600)
    #plot_gene("043-erg-vec1_a.tabular.txt","043-erg-vec1_b.tabular.txt","043-erg-vec2.tabular.txt","043-erg-vec3.tabular.txt")
  }
}








## windowed start and end position does not seem to matter that much

