# s-plot like histogram
mcmc_hist=function(x, nbins=30, labels=c(1,1,0,0), freq=T,...){
  # nbins: integer or a character (as breaks in hist function)
  # if freq=T (default) uses frequencies, else uses the densities
  xlim=range(x,na.rm=T)
  if(is.numeric(nbins))
    nbins=seq(xlim[1], xlim[2], len=nbins+1)
  h=hist(x, breaks=nbins, plot=F)
  b=h$breaks; c=h$counts
  if(!freq) c=h$density
  c=rep(c, each=2)
  c=c(0,c,0)
  b=rep(b, each=2)
  df=data.frame(breaks=b, counts=c)
  magplot(df, type='s', side=1:4, las=1,grid=F, labels=labels,...)
  abline(h=0, col='gray')
  return(invisible(df))
}

# display 2D binned image + contours + points
mcmc_plot=function(x, y, gridsize=200, showdensmap=T, cont=c(68,95), 
                   labels=c(1,1,0,0), drawpoints=F,...){
  # gridsize: squared bins of the density image
  # cont: probability contours to be displayed
  # if showdensmap=T display the density image
  # if drawpoints=T plot the points outside the maximum contour
  require(ks)
  require(magicaxis)
  kd=kde(cbind(x,y), gridsize = gridsize, verbose = F)
  im=list(x=kd$eval.points[[1]], y=kd$eval.points[[2]],z=kd$estimate)
  magplot(range(x), range(y), pch='', grid=F, side=1:4, las=1,
          labels=labels,...)
  if(showdensmap){
    magimage(im, col=rev(grey((0:1000)/1000)), add=T, magmap=F)
    magaxis(1:4, labels=c(0,0,0,0),...)
  }
  if(drawpoints){
    require(sp)
    cl=contourLines(im, level=kd$cont[max(cont)])
    pip=point.in.polygon(x,y,cl[[1]]$x, cl[[1]]$y)
    points(x[pip==0], y[pip==0], pch='.', col='gray50')
  }
  for(i in 1:length(cont))
    plot(kd, cont=cont[i], add=T, lwd=2, lty=i, drawlabels=F, col=1)
  #points(mean(x),mean(y),pch=3,col='white')
}

# create corner plot of MCMC chains
cornerplot=function(chains, samples=NULL, labels=NULL, histbins=30, lim=NULL, 
                    gridsize=30, cont=c(.68,.95), ci=c(.16,.84), minorn=2, 
                    drawpoints=F,...){
  # chains: matrix or data.frame where columns corresponds to parameters val.
  # samples: integer giving the number of items to choose from the chains
  # labels: optional names for the parameters (expressions or characters),
  #         if NULL (default), uses the column names
  # histbins: number of bins for the histograms
  # lim: the level of sigma clipping used to cut the chains (from magclip)
  #      if NULL (default), no action is done
  # gridsize: squared bins of the density image
  # cont: probability contours to be displayed
  # ci: confidence interval to be displayed on the histograms, 
  #     also used to display as the intervals of the median values over plots
  # minorn: number of minor-axis divisions (to improve aesthetics)
  # if drawpoints=T plot the points outside the maximum contour
  parnames=colnames(chains)
  if(is.null(labels))
    labels=parnames
  Npar=ncol(chains)
  N=nrow(chains)
  if(Npar <= 1) 
    stop("Need 2+ parameters!")
  if(!is.null(samples))
    chains=chains[sample.int(N,samples),]
  mdvec=apply(chains, 2, median)
  qvec=apply(chains, 2, quantile, ci)
  bqs=vector('list', len=Npar)
  for(i in 1:Npar){
    lbi=labels[i][[1]]
    med=sprintf('%.3f', round(mdvec[i],3))
    lwr=sprintf('%.3f', round(qvec[1,i],3))
    upr=sprintf('%.3f', round(qvec[2,i],3))
    bqs[[i]]=bquote(.(lbi)==.(med)[paste('-',.(lwr))]^{paste('+',.(upr))})
  }
  if(!is.null(lim)){
    clpch=apply(chains,2, function(x) magclip(x, lim))
    clips=sapply(clpch, function(x) x$clip)
    chains=chains[apply(clips,1,all),]
  }
  op=par(mfrow=c(Npar,Npar), mai=rep(.03,4), oma=c(3.5,3.5,1.5,1.5), 
         family='serif')
  for (i in 1:Npar){
    for (j in 1:Npar){
      if (i < j){
        # upper panels
        plot.new(); next
      } 
      if (i == j){
        # diagonal panels
        lb=c(0,0,0,0)
        #if(j==1) lb=c(0,1,0,0)
        if(j==Npar) lb=c(1,0,0,0)
        xlab=ifelse(lb[1]==1, labels[j], '')
        bh=mcmc_hist(chains[,j], labels=lb, freq=F, lwd=1.5, xlab=xlab, 
                     minorn=minorn, family='serif', nbins=histbins)
        mtext(bqs[[j]], cex=.8)
        # displaying confidence inteval polygon
        fi=findInterval(qvec[,j], bh$breaks)
        polyx=c(rep(qvec[1,j],2),bh$breaks[(fi[1]+1):(fi[2])],rep(qvec[2,j],2))
        polyy=c(0,bh$counts[(fi[1]+1)],bh$counts[(fi[1]+1):(fi[2]+1)],0)
        polygon(polyx, polyy, col="#FF00004D", border = NA)
        abline(v=mdvec[j], col='red')
      }
      if (i > j){
        # lower panels
        lb=c(0,0,0,0)
        if(j==1) lb[2]=1
        if(i==Npar) lb[1]=1
        xlab=ifelse(lb[1]==1, labels[j], '')
        ylab=ifelse(lb[2]==1, labels[i], '')
        mcmc_plot(chains[,j],chains[,i], labels=lb, xlab=xlab, ylab=ylab, 
                  minorn=minorn, family='serif', gridsize=gridsize,
                  drawpoints=drawpoints, cont=cont*100)
        abline(v=mdvec[j], h=mdvec[i], col='red')
      }
    }
  }
  par(op)
}

# Reproducible example
set.seed(123)
ch1=rnorm(1e6)
ch2=rgamma(1e6,10)
ch3=ch1+runif(1e6)
chains=cbind(ch1,ch2,ch3)
labels=c(expression(alpha),expression(beta),expression(gamma))
cornerplot(chains=chains, labels=labels, samples=5e4, drawpoints=T)

