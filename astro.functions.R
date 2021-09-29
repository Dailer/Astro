

# Function for make fraction plots by bins
fracbin=function(variable, classe, nbin='hist', error='binomial', plot=T, erbar=T, ...){
  # nbin can be a integer number or 'hist'
  # error can be binomial, bootstrap or FALSE
  require(binom)
  require(Hmisc)
  require(magicaxis)
  require(boot)
  require(OneR)
  
  stopifnot(is.numeric(nbin) | nbin=='hist')
  df=na.omit(data.frame(var=variable, class=classe))

  # if(is.numeric(nbin)) cl=bin(df$var, nbin, 1:nbin, method)
  # if(nbin=='hist'){
  #   h=hist(df$var, 'fd', plot=F)
  #   b=h$breaks; bcen=h$mids; nbin=length(b)-1
  #   cl=cut(df$var, h$breaks, F)
  # }
  # lst=split(df, cl)
  # if(is.numeric(nbin)) bcen=unlist(lapply(sp, median))
  # nob=c(table(cl))

  if(is.numeric(nbin)){
    bins=seq(min(df$var), max(df$var), length.out = nbin+1)
    binss=seq(min(df$var), max(df$var), length.out = 2*nbin+1)
    bcen=binss[seq(2, length(binss), 2)]
  } else if(nbin=='hist'){
    hs=hist(df$var, 'fd', plot = FALSE)
    bins=hs$breaks
    bcen=hs$mids
    nbin=length(bins)-1
  } else {
    stop('nbin must be integer or "hist"')
  }
  lst=list()
  # lista de objetos de cada bin
  for(i in 1:nbin){
    dfm=subset(df, df$var > bins[i] & df$var <= bins[i+1])
    lst=c(lst, list(dfm))
  }
  cl=as.numeric(levels(as.factor(classe)))
  # numero de objetos por bin
  nob=unlist(lapply(lst, nrow))
  #print(lst)
  
    
  fr=list()
  # lista de fracciones de cada bin
  for(j in 1:nbin){
    ls=lst[[j]]
    ln=vector(length = length(cl))
    for(k in 1:length(cl)){
      ln[k]=length(which(ls$class == cl[k]))/nrow(ls)
    }
    fr=c(fr, list(ln))
  }
  
  # redimensionando lista 
  fin=list()
  for(g in 1:length(cl)){
    c1=lapply(fr, function(x) x[g])
    c1=unlist(c1)
    fin=c(fin, list(c1))
  }
  
  # creando data frame final
  dfin=as.data.frame(fin, col.names = F)
  dfin=data.frame(bcen, dfin)
  cln=paste0('c', as.character(cl))
  names(dfin)=c('bin', cln)
  
  lencl=length(cl)
  
  if(error!=FALSE){
    
    if(error == 'binomial'){
      # Finding binomial confidence intervals
      
      cflt=list()
      for(p in 1:nbin){
        lt=lst[[p]]
        bg=data.frame(lower=vector(), upper=vector())
        for(q in 1:lencl){
          bc=binom.confint(length(which(lt$class==cl[q])), nrow(lt), methods = 'exact')
          bi=bc[5:6]
          bg=rbind(bg, bi)
        }
        cflt=c(cflt, list(bg))
      }
      
      ls2=list()
      for(v in 1:lencl){
        lp=lapply(cflt, function(x) x[v,])
        dfm=data.frame(matrix(unlist(lp), ncol = 2, byrow = T))
        nms=c(paste0('c', cl[v], '.lower'), paste0('c', cl[v], '.upper'))
        colnames(dfm)=nms
        ls2=c(ls2, list(dfm))
      }
      
      dfint=do.call("cbind", ls2)
      dfin=cbind(dfin, dfint)
      
    }
    
    if(error == 'bootstrap'){
      # Finding Bootstrap mean and conf. intervals
      samplemean=function(x, d) {return(mean(x[d]))}
      fr2=list()
      for(j in 1:nbin){
        ls=lst[[j]]
        dci=data.frame()
        
        for(k in 1:length(cl)){
          sam=subset(ls, ls$class == k)[,1]
          
          if(length(sam) > 1){
            bres1=boot(sam, samplemean, R=1000)
            b.mean=mean(bres1$t)
            b.ci=boot.ci(bres1)$normal[1,2:3]
            bcm=c(b.mean, b.ci)
            #bcm=c(bcm[1], bcm[1]/bcm[2], bcm[3]/bcm[1])
            bcm=c(b.mean, sd(bres1$t), sd(bres1$t))
            dci=rbind(dci, bcm)
          } else if(length(sam) == 1){
            dci=rbind(dci, c(b.mean, 0, 0))
          } else {
            dci=rbind(dci, c(0,0,0))}
          names(dci)=c('mean', 'lower', 'upper')
          
        }
        
        fr2=c(fr2, list(dci))
      }
      
      ls3=list()
      for(v in 1:lencl){
        lp=lapply(fr2, function(x) x[v,])
        dfm=data.frame(matrix(unlist(lp), ncol = 3, byrow = T))
        nms=c(paste0('c', cl[v], '.mean'), paste0('c', cl[v], '.lower'), paste0('c', cl[v], '.upper'))
        colnames(dfm)=nms
        ls3=c(ls3, list(dfm))
      }
      
      dfint=do.call("cbind", ls3)
      dfin=cbind(dfin, dfint)
      
    }
  }
  
  dfin=data.frame(dfin, n.obj=nob)
  co=c('#3B4992','#EE0000','#008B45','#631879','#008280','#BB0021','#BB0021',
       '#5F559B','#808180','#1B1919')
  #co=rainbow(lencl+1)
  #co=sample(co, lencl)
  if(plot==TRUE){
    magplot(dfin$bin, rep(0, length(dfin$bin)), ylim = c(-0.05,1.05), pch='', ylab = 'Fraction', ...)
    abline(h=c(0,1), lty=2, col='gray')
    for(w in 1:lencl){
      points(dfin$bin, dfin[,w+1], pch=w-1, col=co[w], type='b')
        if(error!=FALSE){
          lwnm=paste0(cln[w], '.lower')
          lwind=grep(lwnm, colnames(dfin))
          uppind=lwind+1
          if(error=='binomial'){
            errbar(dfin$bin, dfin[,w+1], dfin[,lwind], dfin[,uppind], add=T, pch='', errbar.col=co[w])
          }
          if(error=='bootstrap'){
            errbar(dfin$bin, dfin[,w+1], dfin[,w+1]+dfin[,lwind], dfin[,w+1]-dfin[,uppind], 
                   add=T, pch='', errbar.col=w+1)
          }
        }
    }
    legend('top', legend = cl, pch=(1:lencl)-1, col=co[1:lencl], cex=.8)
  }
  return(dfin)
}

# Fraction of class y in bins of continuous x
fracbin=function(x, class, nbin='hist', method='length', error='binomial', plot=T, 
                 legend=T, add=F, cols=NULL, ...){
  # nbin can be an integer number or 'hist'
  # method can be 'length', 'content' or 'clusters', only applies if nbin is numeric
  # error can be 'binomial' (95% CI), 'bootstrap' or FALSE

  stopifnot(is.numeric(nbin) | nbin=='hist')
  df=na.omit(data.frame(var=x, class=factor(class)))
  cn=levels(df$class); ncl=length(cn)
  
  if(is.numeric(nbin)){
    require(OneR)
    nbin=as.integer(nbin)
    cl=bin(df$var, nbin, 1:nbin, method)
    n=c(table(cl))
  }
  if(nbin=='hist'){
    h=hist(df$var, plot=F)
    cl=cut(df$var, h$breaks, F)
    n=h$counts
    nbin=length(h$breaks)-1
  }
  l.var=split(df$var, cl)
  bc=unlist(lapply(l.var, median))
  sp.cl=split(df$class, cl)
  lt.cl=lapply(sp.cl, function(x) table(x))
  l.cl=lapply(sp.cl, function(x) round(table(x)/length(x), 5))
  f=do.call(rbind, l.cl)
  colnames(f)=paste0('f', cn)
  d=data.frame(bin=bc, f)
  
  if(error=='binomial'){
    require(binom)
    l.ic=lapply(lt.cl, function(x) binom.confint(x, sum(x), methods='exact')[,5:6])
    ics=round(do.call(rbind, lapply(l.ic, function(x) c(t(x)))), 5)
    colnames(ics)=paste0(c('lower','upper'), rep(cn, each=2))
    d=cbind(d, ics)
  }
  if(error=='bootstrap'){
    require(boot)
    for(i in 1:nbin){
    # continuar!  
    }
  }
  if(plot){
    if(is.null(cols)){
      cols=c('#3B4992','#008B45','#EE0000','#631879','#008280','#BB0021','#BB0021','#5F559B',
             '#808180','#1B1919')
    }
    pchs=c(16,0,17,15,4,1,8,2,3,7,5,9,6,10,11,18,12,13,14)
    require(magicaxis)
    if(!add){
      magplot(range(d$bin), c(0,1), pch='', ylab='Fraction', las=1, side=1:4, 
              labels=c(1,1,0,0), ...)
      grid(lty=1, col='#0000001A'); abline(h=0:1, lty=1, col='gray')
    }
    for(i in 1:ncl){
      j=cn[i]
      fname=paste0('f',j)
      lname=paste0('lower',j)
      uname=paste0('upper',j)
      points(d[,1], d[,fname], type='b', pch=pchs[i], col=cols[i])
      arrows(d[,1], d[,lname], d[,1], d[,uname], col=cols[i], length=0)
    }
    if(legend){
      legend('topright', col=cols[1:ncl], pch=pchs[1:ncl], legend = cn, inset=0.01,
             bty='n', y.intersp = .7)
    }
  }
  return(invisible(d))
}

# Median values of continuous y in bins of continuous x
medbin=function(x, y, nbin='hist', method='length', err.type='ci', min.cts=1, plot=T, lim=NULL, 
                xlab=NULL, ylab=NULL, poly=F, add=F, col=1, cex.axis=1, cex.lab=1, ...){
  # nbin can be an integer, a vector or 'hist'
  # method can be 'length', 'content' or 'clusters', only applies if nbin is numeric
  # err.type can be 'ci' (mean and 95% CI) or 'q' (median and quantile error)
  stopifnot(is.numeric(nbin) | nbin=='hist')
  df=na.omit(data.frame(x, y))
  
  if(is.numeric(nbin) & length(nbin)==1){
    require(OneR)
    nbin=as.integer(nbin)
    cl=bin(df$x, nbin, 1:nbin, method)
    n=c(table(cl))
    lx=split(df$x, cl)
    bk=c(unlist(lapply(lx, min)), max(df$x))
  }
  if(is.numeric(nbin) & length(nbin)>1){
    bk=nbin; nbin=length(bk)-1
    cl=cut(df$x, bk, F)
    n=c(table(cl))
    lx=split(df$x, cl)
  }
  if(nbin=='hist'){
    h=hist(df$x, plot=F)
    bk=h$breaks
    cl=cut(df$x, bk, F)
    n=h$counts
    nbin=length(h$breaks)-1
  }
  ly=split(df$y, cl)
  ncts=unlist(lapply(ly, length))
  bc=bk[-(nbin+1)]+diff(bk)/2
  
  if(err.type=='q'){
    qf=function(x) return(0.7415*(quantile(x,.75, na.rm=T)-quantile(x,.25, na.rm=T)))
    m=unlist(lapply(ly, median, na.rm=T))
    e=unlist(lapply(ly, qf))
    lw=m-e; up=m+e
    cname='median'
  }
  if(err.type=='ci'){
    require(Rmisc)
    cif=function(x) return(CI(c(na.omit(x)), .95))
    me=data.frame(do.call(rbind, lapply(ly, cif)))
    m=me$mean; lw=me$lower; up=me$upper
    cname='mean'
    #print(me)
  }
  if(length(m)!=nbin){
    bc=bc[as.numeric(names(n))]
    warning(nbin-length(n),' bins ommited')
  } 
  d=data.frame(bin=bc, mean=m, lower=lw, upper=up, N=ncts)
  colnames(d)[2]=cname
  d=na.omit(d)
  d=subset(d, !is.infinite(bin))
  d=subset(d, N>=min.cts)
  if(plot){
    require(magicaxis)
    #if(is.null(lim)) lim=range(c(d$lower,d$upper))
    if(is.null(lim)) lim=range(y, na.rm=T)
    if(!add){
      rg=max(d$bin)-min(d$bin)
      xrg=c(min(d$bin)-0.02*rg, max(d$bin)+0.02*rg)
      magplot(xrg, range(c(d$lower,d$upper)), pch='', side=1:4, labels=c(1,1,0,0),
              las=1, ylab=ylab, xlab=xlab, ylim=lim, family=par()$family,  
              cex.axis=cex.axis, cex.lab=cex.lab)
      #grid(lty=1, col='#0000001A')
    }
    if(poly){
      ss=smooth.spline(d$bin, d$mean, spar=.4)
      bin=ss$x
      s_mean=ss$y
      s_lower=smooth.spline(d$bin, d$lower, spar=.4)$y
      s_upper=smooth.spline(d$bin, d$upper, spar=.4)$y
      polygon(c(bin, rev(bin)), c(s_lower, rev(s_upper)), border=NA,col=col)
      lines(bin, s_mean, ...)
      d$mean=s_mean
      d$lower=s_lower
      d$upper=s_upper
    }else{
      points(d$bin, d[,2], type='b', col=col, ...)
      arrows(d$bin, d$lower, d$bin, d$upper, length=.05, code=3, angle=90, col=col, ...)
    }
  }
  return(invisible(round(d, 4)))
}

# Median values of continuous x in y classes
medclass=function(x, y, err.type='ci', lim=range(x, na.rm=T), lab='', plot=T, add=F, ...){
  # err.type can be 'ci' (mean and 95% CI) or 'q' (median and quantile error)
  cn=levels(factor(y))
  n=length(cn)
  sp=split(x, y)
  
  if(err.type=='q'){
    qf=function(x) 0.7415*diff(quantile(x, c(.25,.75), na.rm = T))
    m=unlist(lapply(sp, median, na.rm=T))
    e=unlist(lapply(sp, qf))
    lw=m-e; up=m+e
  }
  if(err.type=='ci'){
    require(Rmisc)
    cif=function(x) return(Rmisc::CI(c(na.omit(x)), .95))
    me=do.call(rbind, lapply(sp, cif))
    m=me[,2]; lw=me[,1]; up=me[,3]
  }
  if(plot){
    require(magicaxis)
    if(!add){
      magplot(c(0.8,n+0.2), range(x, na.rm=T), pch='', side=2, las=1, ylab=lab, ylim=lim)
      magaxis(1, majorn=n, minorn=0, xlab='class', labels = F)
      axis(1, 1:n, cn, tick = F, line=-0.5)
      grid(nx=NA, ny=NULL, lty=1, col='#0000001A') 
      abline(v=1:n, col="#0000001A", lty=1)
    }
    points(1:n, m, type='b', ...)
    arrows(1:n, lw, 1:n, up, length=.05, code=3, angle=90, ...)
  }
  data=round(data.frame(mean=m, lower=lw, upper=up),4)
  return(invisible(data))
}

# fractions of x classes in y classes with CI (95%) and test of proportion (if y has two classes)
propclass=function(x, y, type=1, plot=T){
  # type 1: twoSamplePermutationTestProportion
  # type 2: prop.test
  require(binom)
  tb=table(x, y)
  l=ncol(tb)
  n=nrow(tb)
  pt=round(prop.table(tb, 2), 3)
  trials=apply(tb, 2, sum)
  bconf=binom.confint(c(tb), rep(trials,each=n), 0.95, 'exact')
  bconf=round(bconf[,5:6], 3)
  bconf=split(bconf, rep(1:l, each=n))
  names(bconf)=names(trials)
  bconf=do.call(cbind, bconf)
  rownames(bconf)=rownames(tb)
  if(l>2) message('Warning: y must have two classes to perform the proportion test')
  if(l==2){
    if(type==1) require(EnvStats)
    k="Number Successes and Trials"
    pvs=numeric(n)
    for(i in 1:n){
      suc=tb[i,]
      if(type==1){
        pvs[i]=twoSamplePermutationTestProportion(suc, trials, k)$p.value
      }
      if(type==2){
        pvs[i]=prop.test(suc, trials)$p.value
      }
    }
    m=cbind(pt, p.value=round(pvs,3))
    m=data.frame(m, diff=pvs<0.05)
    cat('Fractions of x in y and proportion test p-value\n')
  }else{ 
    cat('Fractions of x in y\n')
    m=as.data.frame(cbind(pt))
  }
  if(plot){
    pch=c(16,0,17,3,15,4,1,8,2,7,5,9,6,10,11,18,12,13,14)
    col=c('#0F5AB0','#EE2111','#1A983D','#EC569A','#FB6B0C',
          '#1A9BED','#A1C721','#FEC10F','#1CA08C','#9A703E')
    la=seq(1, ncol(bconf), 2); ua=la+1
    if(is.numeric(y)) xn=as.numeric(names(trials))
    if(is.character(y)) xn=1:l
    for(i in 1:n){
      if(i==1){
        plot(xn, unlist(m[1,1:l]), pch=16, type='b', ylim=range(m), ylab='fraction',
             col=col[1], xlab='class', panel.first = grid(), axes=F)
        axis(1, xn, sort(unique(y[!is.na(y)]))); axis(2); box()
        arrows(xn, unlist(bconf[1,la]), xn, unlist(bconf[1,ua]), length = 0,
               col=col[1])
      }else{
        points(xn, unlist(m[i,1:l]), pch=pch[i], type='b', col=col[i])
        arrows(xn, unlist(bconf[i,la]), xn, unlist(bconf[i,ua]), length = 0,
               col=col[i])
      }
    }
    legend('right', pch=pch[1:n], lty=1, col=col[1:n], bty='n', y.intersp=.7, 
           legend=rownames(m))
  }
  return(list(fractions=m, CI=bconf))
}

# Function for download images from 6df web site
fits6df=function(input, select = 'id', delta = 2){
  
  #id=sample(gs$ID, 12)
  require(FITSio)
  githubfile='https://raw.githubusercontent.com/Dailer/Astro/master/6dfcatalog.txt'
  
  if(select == 'id'){
    id=as.character(input)
    lid=length(id)
    cat(paste0('\nGalaxies: ',lid,'\n'))
    n=1:lid
    lst=data.frame(n=n, id)
    cat('Reading 6dfGSzDR3 catalog\n')
    catal=read.table(githubfile)
    mcat=merge(lst, catal, by.x = 'id', by.y = 'ID')
    mcat=mcat[order(mcat$n),]
  }
  
  if(select == 'radec'){
    require(RANN)
    lid=nrow(input)
    cat(paste0('\nGalaxies: ',lid,'\n'))
    cat('Reading 6dfGSzDR3 catalog\n')
    catal=read.table(githubfile)
    cat('matching by coordinate distance\n')
    nn=nn2(catal[,2:3], input, k=1)
    mcat=data.frame(catal[nn$nn.idx,], dist=nn$nn.dists*3600)
    n=1:lid
    mcat=data.frame(n=n, mcat)
    colnames(mcat)[3:4]=c('nearest.ra','nearest.dec')
    wd=which(mcat$dist > delta)
    
    if(length(wd) != 0){
      mcat=mcat[-wd,]
      cat(paste0(length(wd),' galaxies rejected by distance (d > ',delta,' arcsec)\n'))
      lid=nrow(mcat)
      n=1:lid
    }
  }
  
  d1='http://www-wfau.roe.ac.uk/6dFGS/cgi-bin/show.cgi?release=dr3&targetname='
  
  lks=vector(length = lid)
  cat('Retrieving downloading links\n')
  for(i in 1:lid){
    d2=as.character(mcat$ID[i])
    dr=paste0(d1,d2)
    rl=readLines(dr, n=20)
    lk=rl[grep('Displaying', rl)]
    lk=strsplit(lk, '"')[[1]][2]
    lks[i]=lk
  }
  
  cat(paste0('Downloading ',lid,' FITS files\n'))
  pb=txtProgressBar(0,lid,0,'=',style=3)
  for(j in n){
    dest=tail(unlist(strsplit(lks[j], "/")),1)
    dest=paste0(mcat$n[j],'.fits')
    download.file(url = lks[j], destfile = dest, method = "wget", quiet = TRUE)
    setTxtProgressBar(pb, j)
  }
  
  z=zqual=rao=deco=vector(length = lid)
  for(i in n){
    #print(i)
    if(mcat$SPECID[i] == 1){
      rf=readFITS(paste0(mcat$n[i],'.fits'), hdu = 6)$hdr
      z[i]=as.numeric(rf[grep('Z_HELIO', rf)+1])
      zqual[i]=as.numeric(rf[grep('QUALITY', rf)+1])
    } else {
      z[i]=zqual[i]=NA
    }
  }
  
  mcat$SPECID=NULL
  outcat=data.frame(mcat, z=z, zqual=zqual)
  write.table(outcat, 'INFO.dat')
  
  return(outcat)
}

# Function for download images from GAMA web site
gamaim=function(input, filter='r', select = 'id', delta = 2){
  
  #id=sample(gs$ID, 12)
  require(FITSio)
  require(R.utils)
  githubfile='https://raw.githubusercontent.com/Dailer/Astro/master/GAMADR2.txt'
  
  if(select == 'id'){
    id=as.character(input)
    lid=length(id)
    cat(paste0('\nGalaxies: ',lid,'\n'))
    n=1:lid
    lst=data.frame(n=n, id)
    cat('Reading 2MASSDR2 catalog\n')
    if(!exists('catal')) catal=read.table(githubfile)
    mcat=merge(lst, catal, by.x = 'id', by.y = 'cataid')
    mcat=mcat[order(mcat$n),]
    if(nrow(mcat)==0){stop('Object not found in GAMA catalog')}
  }
  
  if(select == 'radec'){
    require(RANN)
    lid=nrow(input)
    cat(paste0('\nGalaxies: ',lid,'\n'))
    cat('Reading 2MASSDR2 catalog\n')
    catal=read.table(githubfile)
    cat('matching by coordinate distance\n')
    nn=nn2(catal[,2:3], input, k=1)
    mcat=data.frame(catal[nn$nn.idx,], dist=nn$nn.dists*3600)
    n=1:lid
    mcat=data.frame(n=n, mcat)
    colnames(mcat)[3:4]=c('nearest.ra','nearest.dec')
    wd=which(mcat$dist > delta)
    
    if(length(wd) != 0){
      mcat=mcat[-wd,]
      cat(paste0(length(wd),' galaxies rejected by distance (d > ',delta,' arcsec)\n'))
      lid=nrow(mcat)
      n=1:lid
    }
  }
  
  d1='http://www.gama-survey.org/dr2/data/imaging/gama/SersicPhotometry/v07/'
  
  downdir='./GAMA.IMAGES'
  if(!file.exists(downdir)){
    cat('Creating new folder: GAMA.IMAGES\n')
    dir.create(downdir)
  }
  
  setwd(downdir)
  sid=mcat$sigmaid
  name=vector(length = lid)
  cat(paste0('Downloading ',lid,' compressed files\n'))
  pb=txtProgressBar(0,lid,0,'=',style=3)
  for(i in 1:lid){
    nz=8-nchar(sid[i])
    d2=paste0('S', paste0(rep('0',nz), collapse = ''), sid[i],'/moddata.tar.gz')
    dr=paste0(d1,d2)
    dest=paste0('S', paste0(rep('0',nz), collapse = ''), sid[i], '.tar.gz')
    name[i]=dest
    download.file(url = dr, destfile = dest, method = "wget", quiet = TRUE)
    setTxtProgressBar(pb, i) 
  }
  
  cat('\n')
  mcat=data.frame(mcat, file=name)
  
  #name=list.files(downdir, pattern = 'gz')
  
  #lfiles=list.files(imagedir, pattern = '.fits')
  remfil=setdiff(c('u','g','r','i','z','H','J','K','Y'), filter)
  noim=c()
  cat("Extracting files...\n")
  pb=txtProgressBar(0,lid,0,'=',style=3)
  for(k in 1:lid){
    untar(name[k])
    unlink(remfil, recursive = TRUE)
    setwd(paste0('./',filter))
    untar(list.files(pattern = 'tar.gz'))
    fit=list.files(pattern = 'objim')
    if(length(fit)!=0){
      rf=readFITS(fit, hdu=2)
      rfp=rf$imDat; rfh=rf$header
      writeFITSim(rfp, paste0('../', mcat$n[k], '.fits'), header = rfh)
    } else {
      noim=c(noim, k)
    }
    setwd('../')
    unlink(filter, recursive = TRUE)
    setTxtProgressBar(pb, k)
  }
  
  if(length(noim)!=0){
    mcat=mcat[-noim,]
    cat('\n',length(noim),'galaxies images were not found at',filter,'band\n')
    }
  
  cat('Done. Info is in INFO.dat file\n')
  write.table(mcat, 'INFO.dat')
  setwd('../')
  
  return(mcat)
}

# Function for making match by distance
matchdist=function(radec1, radec2, maxdist=0.1, conv=3600){
  
  require(RANN)
  
  nn=nn2(radec2, radec1, k=1)
  idx=nn$nn.idx
  dist=nn$nn.dists
  
  df=data.frame(id=1:nrow(radec1), radec1, radec2[idx,], idx=idx, dist=dist*conv)
  df2=subset(df, dist <= maxdist)
  df3=merge(df[1:5], df2[c(1,6,7)], by='id', all.x = T)
  
  df4=data.frame(index.1=df3$id, index.2=df3$idx, distance=df3$dist)
  
  hist(df4$distance, breaks=round(nrow(df4)/100), main='Histogram of distance', xlab='distance (arcsec)')
  lrej=length(which(is.na(df4$distance)))
  nr=nrow(radec1)
  mtext(paste0('Distance < ', maxdist, ' arcsec'), line = -2, adj = .8)
  mtext(paste0('Rejected: ', lrej,' (', round(lrej*100/nr, 1), '%)'), line = -4, adj = .8)
  
  # cl=ifelse(df$dist > maxdist, 2, 1)
  # df=data.frame(df, col=cl)
  # plot(df[,2], df[,3], pch='.', xlab='X', ylab='Y')
  # points(df[,4], df[,5], pch=4, col=df$col, cex=.5)
  
  return(df4)
    
}

# Find galaxy coordinates from names
galname=function(name){
  # e.g. 'NGC-3640'
  require(celestial)
  require(XML)
  dirbase='http://simbad.u-strasbg.fr/simbad/sim-id?Ident='
  if(grepl('-', name) & grepl('NGC', name)) name=gsub('-','',name)
  galdir=paste0(dirbase,name)
  rl=readLines(galdir, n=500)
  if(any(grepl('Identifier not found in the database', rl))){
    ra=dec=z=NA
  }else{
    coo=rl[grep('ICRS', rl)[1]+11]
    z=rl[grep('Redshift', rl)+8]
    z=ifelse(length(z)!=0, as.numeric(strsplit(z, ' ')[[1]][2]), NA)
    if(grepl('<', coo)){
      coo=gsub('[[:alpha:]]','',coo)
      coo=XML::xpathApply(htmlParse(coo, asText=TRUE),"//body//text()", xmlValue)[[1]]
      coo=gsub('>>','',coo)
    }
    coo=unlist(strsplit(coo, ' '))
    cora=coo[1:3]; codec=coo[4:6]
    ra=paste(cora, collapse  = ':')
    dec=paste(codec, collapse = ':')
    
    if(length(coo)==0){
      ra=dec=NA
      message(paste('No entry could be found in Simbad for', name))
    }else{
      ra=hms2deg(ra)
      dec=dms2deg(dec)
    }
  }
  return(data.frame(name=name, ra=ra, dec=dec, z=z))
}

# Function for make composite image from SDSS DR12
rgbimage=function(radec, size=120, R='i', G='r', B='g', DR=12){
  
  start.time=Sys.time()
  nr=nrow(radec)
  n=seq_len(nr)
  #radec=cbind(id=n, radec)
  wna=c()
  rdcn=colnames(radec)
  
  if('name' %in% rdcn == TRUE){
    cat('Object name input, retrieving coordinates from http://simbad.u-strasbg.fr... \n')
    glname=as.character(radec$name)
    gn=galname(glname)
    ra=gn[1]; dec=gn[2]
    radec=data.frame(ra=ra, dec=dec)
    rdcn=colnames(radec)
  }
  
  if("run" %in% rdcn == FALSE || "camcol" %in% rdcn == FALSE || "field" %in% rdcn == FALSE){
    
    if(ncol(radec) < 2){
      stop('Error: Must provide ra and dec')
    }
    if("ra" %in% rdcn == FALSE || "dec" %in% rdcn == FALSE){
      stop('Error: ra and dec column names must be set')
    }
    if(is.factor(radec$ra) || is.factor(radec$dec)){
      dec=dms2deg(as.character((radec$dec)))
      ra=hms2deg(as.character((radec$ra)))
      radec$ra=ra; radec$dec=dec
    }
    
    # Extracting run, camcol and field info
    ra=radec$ra; dec=radec$dec
    
    basdir1.1='http://cas.sdss.org/'
    basdirdr=switch(DR, '12'='dr12', '13'='dr13', '14'='dr14')
    basdir1.2='/en/tools/search/x_results.aspx?searchtool='
    basdir2='Radial&uband=&gband=&rband=&iband=&zband=&jband=&hband=&kband=&'
    basdir3='TaskName=Skyserver.Search.Radial&ReturnHtml=true&whichphotometry'
    basdir4='=optical&coordtype=equatorial&ra='
    basdir5='&min_u=0&max_u=20&min_g=0&max_g=20&min_r=0&max_r=20&min_i=0&max_i=20&min_z=0&max_z=20&'
    basdir6='min_j=0&max_j=20&min_h=0&max_h=20&min_k=0&max_k=20&format=html&TableName=&limit=1'
      
    cat(paste0("Extracting RUN CAMCOL FIELD info from DR",DR,", please wait..."),"\n")
    dir=paste0(basdir1.1, basdirdr, basdir1.2, basdir2, basdir3, basdir4, sprintf("%.5f", ra), 
                 '&dec=', sprintf("%.5f", dec), '&radius=0.1',basdir5,basdir6)
      
    options(warn=-1)
    p=NULL; attempt=1
    while( is.null(p) && attempt <= 5 ) {
      attempt=attempt + 1
      tryCatch(try(p <- readLines(dir)), error = function(e) {})
    }
    options(warn=0)
    
    lgr=function(str, p){ l=length(grep(str, p)); return(l) }
    if (lgr('No objects',p) != 0 || lgr('fix errors',p) != 0 
        || lgr('No entries',p) != 0 || lgr('between 0',p) != 0){
      run=camcol=field=sdss=NA
      
     } else {
        
        g2=grep('BGCOLOR=#', p)
        
        us2=unlist(strsplit(p[g2], 'id='))[2]
        us4=unlist(strsplit(us2, '</font>'))
        us4.1=unlist(strsplit(us4[1], '>'))[1]
        us4=us4[-c(1,length(us4))]
        us5=strsplit(us4, '>')
        us6=unlist(lapply(us5, function(x) x[4]))
        us7=as.double(us6)[1:6]
        
        run=us7[1]
        camcol=us7[3]
        field=us7[4]
        sdss=us4.1
      }
      
      setTxtProgressBar(txtProgressBar(0,nr,0,'=',style=3), j)
      
      radec=data.frame(radec, run=run, camcol=camcol, field=field, sdssid=sdss)
      cat("\n")
        
      end.time=Sys.time()
      time.taken=round(as.numeric(end.time-start.time), 3)
        
      wna=which(is.na(run))
      if(length(wna)==nr){stop('No object found in DR12 Data Base')}
      l=length(wna)
      if(l >= 1){
        cat("Done, but info for", l, "object(s) have not been found:",wna, "  ", time.taken, 
            attr(unclass(end.time-start.time),'units'),"\n")
      } else {
        cat("Done, All info found!  ", time.taken, attr(unclass(end.time-start.time),'units'), "\n")
      }
      
      
    
  } 
  
  radec=data.frame(id=1:nr, radec)
  
  Link_list=c()
  basedir='https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/'
  #ft=c('u','r','z')
  ft=c(B, G, R)
  for(i in 1:3){
    
    link_DR12=paste0(basedir,as.character(run),'/',as.character(camcol),
                     '/frame-',ft[i],'-',sprintf("%06.f", run),'-',as.character(camcol),'-',
                     sprintf("%04.f", field),'.fits.bz2')
    
    Link_list[i]=link_DR12
  }
  
  # baixando as imagens
  dnames=vector(length = 3)
  start.time=Sys.time()
  cat("Downloading", length(Link_list), "image FITS files from DR12 Data Base, please wait...\n")
  pb=txtProgressBar(0,length(Link_list),0,'=',style=3)
  for(filename in Link_list){
    dest=tail(unlist(strsplit(filename, "/")),1)
    try(download.file(url = filename, destfile = dest, method = "wget", quiet = F))
    mi=match(filename, Link_list)
    setTxtProgressBar(pb, mi) 
    dnames[mi]=dest
  }
  end.time=Sys.time()
  time.taken=round(as.numeric(end.time-start.time), 3)
  cat('\nDone  ', time.taken, attr(unclass(end.time-start.time),'units'),'\n')
  
  require(R.utils)
  # extraindo as imagens
  cat("Extracting files...\n")
  start.time=Sys.time()
  pb=txtProgressBar(0,length(Link_list),0,'=',style=3)
  for(k in 1:length(Link_list)){
    nm=dnames[k]
    bunzip2(nm, unlist(strsplit(nm, '.bz2')), overwrite=TRUE)
    setTxtProgressBar(pb, k)
  }
  end.time=Sys.time()
  time.taken=round(as.numeric(end.time-start.time), 3)
  cat('\nDone  ', time.taken, attr(unclass(end.time-start.time),'units'),'\n')
  
  require(FITSio)
  require(magicaxis)
  fnames=vector(length = 3)
  cat("Extracting regions around target galaxies...\n")
  bim=list()
  pb=txtProgressBar(0,length(Link_list),0,'=',style=3)
  for(l in 1:length(Link_list)){
    name=unlist(strsplit(dnames[l], '.bz2'))
    fnames[l]=name
    fit=readFITS(name)
    cf=magcutoutWCS(fit, loc = c(ra,dec), box=c(300,300))$image
    # cf=cutfit(name, ra, dec, coord='radec',
    #           outname = paste0(ft[l],'-band.fits'), 
    #           write=TRUE, r=500)
    # writeFITSim(cf, paste0(ft[l],'-band.fits'), 
    #             crpixn = dim(cf)/2, crvaln = c(ra,dec), header = fit$header)
    bim=c(bim, list(cf))
    setTxtProgressBar(pb, l)
  }
  
  file.remove(fnames)
  
  B=magcutout(bim[[1]], box = c(size,size))$image
  G=magcutout(bim[[2]], box = c(size,size))$image
  R=magcutout(bim[[3]], box = c(size,size))$image
  
  writeFITSim(R, 'R.fits', crpixn = dim(R)/2, crvaln = c(ra,dec), header = fit$header)
  writeFITSim(G, 'G.fits', crpixn = dim(G)/2, crvaln = c(ra,dec), header = fit$header)
  writeFITSim(B, 'B.fits', crpixn = dim(B)/2, crvaln = c(ra,dec), header = fit$header)

  x = seq(0.5, dim(R)[1] - 0.5)
  y = seq(0.5, dim(R)[2] - 0.5)
  magimageRGB(x, y, R, G, B)
  title('RGB Composite Image')
  
}

# funcao para fazer short links em tinyurl.com
shortlink=function(url){
  require(RCurl)
  urle=curlEscape(url)
  d1='https://tinyurl.com/create.php?source=indexpage&url='
  d2='&submit=Make+TinyURL%21&alias='
  dd=paste0(d1, urle, d2)
  rl=readLines(dd)
  gp=grep('<div id=\"copyinfo\" data-clipboard-text=', rl)
  rlf=rl[gp]
  s1=unlist(strsplit(rlf, '"'))
  s2=s1[grep('http', s1)]
  return(s2)
}

# funcao para quantificar separacao entre dois picos de distribuicoes de densidade
sepind=function(x, y){
  require(DescTools)
  a=x; b=y
  if(mean(a)<mean(b)){
    x=a; y=b
  } else {
    x=b; y=a
  }
  dx=density(x, na.rm = T)
  dy=density(y, na.rm = T)
  # Encontrando separacao relativa
  mhx=max(dx$y); mhy=max(dy$y)
  mx=dx$x[which(dx$y==mhx)]
  my=dy$x[which(dy$y==mhy)]
  sep=abs(mx-my)
  nsep=sep/abs(min(dx$x)-max(dy$x))
  # Encontrando area de intersecao relativa
  dfx=data.frame(x1=dx$x, y1=dx$y)
  dfy=data.frame(x2=dy$x, y2=dy$y)
  dfx=subset(dfx, x1 >= mx)
  dfy=subset(dfy, x2 <= my)
  allseq=seq(mx, my, length.out = 500)
  apx=approx(dfx$x1, dfx$y1, allseq)
  apy=approx(dfy$x2, dfy$y2, allseq)
  df=data.frame(x=apx$x, dif=abs(apx$y-apy$y))
  df=na.omit(df)
  int.x=df[,1][which(df[,2]==min(df[,2]))]
  
  dfx=data.frame(x1=dx$x, y1=dx$y)
  dfy=data.frame(x2=dy$x, y2=dy$y)
  dfx1=subset(dfx, x1 >= int.x)
  dfy1=subset(dfy, x2 <= int.x)
  
  ax=AUC(dfx1$x1, dfx1$y1)
  ay=AUC(dfy1$x2, dfy1$y2)
  aint=ax+ay
  
  dfx2=subset(dfx, x1 <= int.x)
  dfy2=subset(dfy, x2 >= int.x)
  
  axt=AUC(dfx2$x1, dfx2$y1)
  ayt=AUC(dfy2$x2, dfy2$y2)
  aintt=axt+ayt
  arat=aint/aintt
  # Encontrando largura relativa dos picos
  hmx=mhx/2; hmy=mhy/2
  spx=splroot2(dx$x, dx$y-hmx)
  spy=splroot2(dy$x, dy$y-hmy)
  difx=spx[2]-spx[1]
  dify=spy[2]-spy[1]
  relx=difx/(max(dy$x)-min(dx$x))
  rely=dify/(max(dy$x)-min(dx$x))
  rel=mean(c(relx, rely))
  
  # Computando indice
  ind=nsep/(arat*rel)
  return(ind)
}

# function for Gradient Pattern Analysis INCOMPLETE!
gpa=function(image, plot=T){
  
  require(pracma)
  require(magicaxis)
  require(deldir)
  require(spatstat)
  require(raster)
  
  # mt=matrix(c(5,5,5,5,5,5,7,7,7,5,5,9,9,7,5,5,7,7,9,5,5,5,5,5,5), 5, 5, byrow = T)
  #mt=blur(as.im(image), sigma = 1)$v
  # r=raster(image)
  # projection(r)="+proj=longlat" # Any projection
  # r=sampleRegular(r, size=2e3, asRaster=TRUE)
  # mt=as.matrix(r)
  gft=pracma::gradient(mt)
  gnt=sqrt(gft$X^2+gft$Y^2)
  gnt=gnt/max(gnt)
  
  ###
  # gftx=as.matrix(imgradient(as.cimg(imt), axes = 'x', scheme = 3))
  # gfty=as.matrix(imgradient(as.cimg(imt), axes = 'y', scheme = 3))
  ###
  
  xy0=pracma::meshgrid(1:nrow(mt), 1:nrow(mt))
  x=as.vector(xy0$Y)-0.5
  y=as.vector(xy0$X)-0.5
  
  mv=max(unlist(gft))
  xn=as.vector(gft$Y)/mv
  yn=as.vector(gft$X)/mv
  arr=data.frame(x0=x, y0=y, x=x+xn, y=y+yn)
  
  imph=atan2(gft$X, gft$Y)
  imph=ifelse(imph<0, 2*pi+imph, imph)
  ph=atan2(yn, xn)
  ph=ifelse(ph<0, 2*pi+ph, ph)
  imf=data.frame(id=1:nrow(arr), arr, norm=sqrt((xn)^2+(yn)^2), phase=ph)
  
  ###
  
  ph1=imph*ifelse(imph < pi, pi+imph, imph-pi)
  ph2=.5*(gnt+ph1)*(gnt+ph1+1)+ph1
  
  np=data.frame(ph2=as.vector(ph2), imf)
  snp=split(np, np$ph2)
  tuni=sum(unlist(lapply(lapply(snp, nrow), function(x) which(x==1))))
  snp=snp[which(lapply(snp, function(x) nrow(x) > 1)==T)]
  snp=lapply(snp, function(x) cbind(x, ph3=ifelse(x$phase==min(x$phase), x$ph2, -1*x$ph2)))
  # snp=lapply(snp, function(x) x$ph2=ifelse(x$phase==min(x$phase), x$ph2, -1*x$ph2))
  snp=lapply(snp, function(x) abs(sum(x)/unique(abs(x))))
  tnp=sum(unlist(snp), na.rm =T)+tuni
  
  dnp=do.call('rbind', snp)
  rownames(dnp)=NULL
  dnp=dnp[order(dnp$id),]
  mnp=matrix(dnp$ph3, dim(ph2)[1], dim(ph2)[2])
  
  ps=np$phase; ps2=np$ph2
  i=1
  while(i != length(ps)){
    we=which(ps2 == ps2[i] & ps != ps[i])[1]
    if(length(we)!=0){
      ps=ps[-we]
      ps2=ps2[-we]
    }
    i=i+1
  }
  
  we=vector(length = nrow(np))
  for(i in 1:nrow(np)){
    we[i]=which(ps2 == ps2[i] & ps != ps[i])[1]
  }
  
  
  ###
  
  vu=data.frame(n=as.vector(gnt), v3)
  vu=unique(vu)
  
  td=arr[as.numeric(rownames(vu)),]
  td=subset(td, sqrt((td$x-td$x0)^2+(td$y-td$y0)^2) != 0)
  magimage(mt, magmap = F)
  arrows(td$x0, td$y0, td$x, td$y, col='red', length = 0.1, angle = 20)
  arrows(arr$x0, arr$y0, arr$x, arr$y, col='green', length = 0.1, angle = 20)
  
  ###

  sp=split(imf, imf$norm)
  sp=sp[which(lapply(sp, function(x) nrow(x) > 1)==T)]
  sp=lapply(sp, function(x) cbind(x, oph=ifelse(x$phase < pi, pi+x$phase, x$phase-pi)))
  sp=lapply(sp, function(x) x[which(x$phase %in% x$oph),])
  sp=sp[which(lapply(sp, function(x) nrow(x) > 1)==T)]
  sp=lapply(sp, function(x) cbind(x, php=x$phase*x$oph))
  sp2=unlist(lapply(sp, function(x) split(x, x$php)), recursive = F)
  sp2=lapply(sp2, function(x) split(x, x$phase))
  sp2=lapply(sp2, function(x) unlist(lapply(x, nrow)))
  
  dsp=do.call('rbind', sp)
  sd=split(dsp, dsp$phase)
  
  ###
  
  cb=data.frame(combinations(nrow(imf), 2), ndif=combn(imf$norm, 2, function(x) Reduce(`-`, x)), 
                ph.comb=combn(imf$phase, 2, function(x) Reduce(`-`, x)))
  scb=subset(cb, ndif == 0 & abs(ph.comb) == pi)
  
  ix=1:(nrow(mt)*ncol(mt))
  va1=setdiff(ix, unique(c(scb$X1, scb$X2)))
  va1=setdiff(va1, subset(imf, norm == 0)[,1])
  
  scb=scb[1:2]
  st=as.matrix(scb)
  i=1
  while(i <= nrow(st)){
    
    ids=as.vector(st[i,])
    wr1=which(st[,1]==ids[1] | st[,1]==ids[2])[-1]
    wr2=which(st[,2]==ids[2] | st[,2]==ids[1])[-1]
    st[wr1,1]=NA; st[wr2,2]=NA
    
    if(any(is.na(st))){
      st=st[-which(is.na(st), arr.ind = T)[,1],]
    }
    
    i=i+1
  }
  
  va2=setdiff(unique(c(scb$X1, scb$X2)), as.vector(st))
  va=sort(c(va1, va2))
  Va=length(va)
  
  # ga=imf[va,]
  # arrows(ga$x0, ga$y0, ga$x, ga$y, col='red', length = 0.1, angle = 20)
  
  imf$asy=FALSE
  imf$asy[va]=TRUE
  
  imfva=imf[va,]
  
  dd=deldir((imfva$x0+imfva$x)/2, (imfva$y0+imfva$y)/2, plotit = T)
  del=dd$delsgs
  
  Ta=nrow(del); g1=(Ta-Va)/Va
  
  #g2=Va*(2-)
  
  if(plot){
    magimage(mt, magmap = F)
    options(warn = -1)
    arrows(arr$x0, arr$y0, arr$x, arr$y, col='green', length = 0.1, angle = 20)
    segments(del$x1, del$y1, del$x2, del$y2, col=2)
    #options(warn = 0)
  }
  
  lst=list(Va=Va, data=imf, G1=g1, Del=del, G2image=gnt, G3image=imph)
  return(lst)
  
}

# functions for calculate photon asymmetry (Nurgaliev et al. 2013)
aphot_main=function(dat, sky, plot=T)
{
  
  phi=sort(dat$theta)/(2 * pi)
  n=N=length(phi)
  UN=vector(length = n-1)
  id=1:n
  
  # loop over all possible initial angles
  for(i in 1:(n-1)){
    if(i!=1){
      th=c(phi[i:n],phi[1:(i-1)])
    }else{
      th=phi
    }
    sum.terms=(th - mean(th) - (2 * id - 1)/(2 * n) + 0.5)^2
    u2=sum(sum.terms) + 1/(12 * n)
    # Watson, G. 1961, Biometrika, 48, 109, eq 26
    UN[i]=(u2 - 0.1/n + 0.1/(n^2)) * (1 + 0.8/n) 
    # 'WatsonTestRad' function in 'circular' package
    
    #UN[i]=1/(12*n)+sum((2*i+1)/(2*n)-th)^2-n*(.5-sum(th)/n)^2
    # Nurgaliev D., McDonald M., Benson B. A., Miller E. D., Stubbs C. W., Vikhlinin A., 2013, 
    # ApJ, 779, 112
  }
  
  WN=min(UN)
  dN=WN/N-1/(12*N)
  
  C=sum(as.numeric(dat$counts > dat$sky))
  dNC=N*(WN-1/12)/C^2
  
  dNres=c(dN, dNC, N, C)
  names(dNres)=c('dN', 'dNC', 'N', 'C')
  
  nmin=order(UN)[1]
  if(nmin!=1){
    thm=c(phi[nmin:n],phi[1:(nmin-1)])
  }else{
    thm=phi
  }
  
  if(plot){
    ec=Ecdf(thm*2*pi, col = 'deepskyblue', pl=F)
    plot(c(0,2*pi), c(0,1), pch='', xlab = 'Angle (rad)', ylab = 'Probability', mgp=c(1.8,.5,0))
    lines(ec$x, ec$y, col = 'deepskyblue', lwd = 2)
    abline(a=0, b=1/(2*pi))
    legend('topleft', legend = c('Uniform CDF', 'Empirical CDF'), col=c(1,'deepskyblue'),
           lty=1, cex=.7)
    mtext(text = paste('dN =', round(dN, 4)), side = 1, line = -2, adj = .9, cex=.85)
    mtext(text = paste('Annulus', unique(dat$ann)), side = 3, line = -4, adj = .1, cex=.7)
  }
  
  return(dNres)
}

mfrow=function(x){
  x=as.integer(x)
  if(x %% 2 == 1 & x != 1){x = x + 1}
  div=seq_len(abs(x))
  factors=div[x %% div == 0L]
  rf=as.integer(sqrt(x))
  difs=abs(factors-rf)
  wm=which(difs==min(difs))
  f1=factors[wm]
  f2=x/f1
  if(x==1){f1=1; f2=1}
  if(x==2){f1=2; f2=2}
  return(sort(c(f1,f2)))
}

aphot=function(image, sky, xo, yo, radius, annuli=4, seq=2, plot=T)
{
  
  require(magicaxis)
  require(Hmisc)
  require(IM)
  
  if(missing(xo) || missing(yo)){
    wm=which(image==max(image), arr.ind = T)
    xom=wm[1]; yom=wm[2]
    dim=dim(image)[1]; r=5
    yi=matrix(rep(1:dim,dim),ncol = dim,byrow = T)
    xi=matrix(rep(1:dim,dim),ncol = dim,byrow = F)
    x = (xi-xom); y = (yi-yom)
    circ = x^2/r^2 + y^2/r^2
    cirmask=ifelse(circ<=1,1,0)
    
    xyo=calcCentroid(image*cirmask)
    xo=xyo[1]; yo=xyo[2]
  }
  
  xyi=which(image==image, arr.ind = T)
  d=sqrt((xyi[,1]-xo)^2+(xyi[,2]-yo)^2)
  a=atan2(xyi[,2]-yo, xyi[,1]-xo)
  a=ifelse(a<0, 2*pi+a, a)
  
  df=data.frame(xyi, theta=a, dist=d, counts=as.vector(image), sky=as.vector(sky))
  #df=subset(df, counts > quantile(sky, probs = .5))
  df=subset.data.frame(df, counts > 0)
  rin=0.05*radius
  
  if(seq==1){
    ds=seq(rin^2, radius^2, length.out = annuli+1)[-1]
    ds=sqrt(ds)
  }else if(seq==2){
    ds=seq(rin, radius, length.out = annuli+1)[-1]
  }
  
  df=subset(df, dist <= max(ds) & dist > rin)
  sr=sapply(df$dist, function(x) which(ds >= x)[1])
  df$ann=sr
  spd=split(df, df$ann)
  
  op=par(mfrow=n2mfrow(annuli), mai=c(.5,.5,.2,.2))
  ldn=lapply(spd, aphot_main)
  par(op)
  ldn=do.call('rbind', ldn)
  
  dNC=ldn[,2]; C=ldn[,4]
  dNC[dNC < 0]=0
  aphval=100*sum(C*dNC)/sum(C)
  # names(aphval)='Aphot'
  pk=c(xo,yo); names(pk)=c('xo','yo')
  
  lout=list('Aphot'=aphval, 'Peak'=pk, 'Rings'=ldn)
  
  if(plot){
    magimage(image)
    cols=rainbow(n=length(unique(sr)))
    lp=lapply(spd, function(x) points(x[,1], x[,2], col=cols[unique(x$ann)], pch=20, cex=.5))
    points(xo, yo, col=3, pch=3)
    draw.circle(xo, yo, ds, border = 'white', lty = 2)
    npix=ldn[,3]
    legend(0,dim(image)[2] , col = cols, pch=20, cex=.6, legend = npix)
    mtext(text = paste('Aphot =', round(aphval, 4)), col = 'white', cex=.7, line=-1)
  }
  
  return(lout)
  
}

# function to convert coloured image file to grayscale
read.im2gray=function(path){
  require(imager)
  im=load.image(path)
  imv=as.vector(im[,,,1:3])
  xdim=dim(im)[1]; ydim=dim(im)[2]
  im=as.cimg(imv, x=xdim, y=ydim, cc=3)
  im=grayscale(im)
  im=as.matrix(im)
  return(im)
}

# funcao para plotar as medias de y en bins de x
varbin=function(x, y, plot=T, labs=c('x','y'), col=1){
  
  h=hist(x, plot = F)
  b=h$breaks; m=h$mids
  
  cl=apply(cbind(x), 1, function(x) which(b/x >= 1)[1])
  xc=data.frame(x, cl, y)
  
  l1=split(xc, xc$cl)
  ng=unlist(lapply(l1, nrow))
  mv=unlist(lapply(l1, function(x) mean(x$y)))
  mve=unlist(lapply(l1, function(x) std.error(x$y)))
  
  out=data.frame(x=m, mean.y=mv, mean.y.err=mve, n=ng)
  out=subset(out, n > 1)
  
  if(plot){
    require(astro)
    require(Hmisc)
    ly=c(min(mv-mve, na.rm = T), max(mv+mve, na.rm = T))
    aplot(m, mv, pch=20, xlim = range(x), ylim = ly, col=col, xlab = labs[1], ylab = labs[2])
    errbar(m, mv, mv+mve, mv-mve, add = T, errbar.col = col, pch='')
  }
  return(out)
}

# Function to convert arcsec to kpc
asec2kpc=function(asec, z, Om=0.27, Ol=1-Om, H0=71){
  require(cosmoFns)
  da=D.A(z, Om, Ol, H0)
  da2arcsec=asec*da/(180*3600/pi/1000)
  return(da2arcsec)
}

# Function to convert pixel to kpc
pix2kpc=function(pix=1, pixscale=1, z, omega.m=0.27, omega.lambda=0.73, H.0=71){
  da=cosmoFns::D.A(z, omega.m, omega.lambda, H.0)
  da2arcsec=da/(180*3600/pi/1000)
  phys.dist=pix*da2arcsec*pixscale
  return(phys.dist)
}

galdens=function(ra, dec, z){
  
  hz=hist(z, plot = FALSE)
  brk=hz$breaks
  hz=data.frame(counts=hz$counts, density=hz$density, mids=hz$mids)
  
  df=data.frame(id=1:length(ra), ra=ra, dec=dec, z=z)
  
  ct=cut(z, brk, labels = FALSE)
  zsp=split(df, ct)
  
  del=lapply(zsp, function(x) deldir::deldir(x$ra, x$dec))
  del.a=lapply(del, function(x) x$summary$del.area)
  da=do.call('c', del.a)
  df$A=1/da
  
  return(df)
  
}

# function to obtain local galaxy density in a cluster for flat universe
locdens=function(ra, dec, z, zref=median(z), radius=500, velocity=1000, plot=T){
  # zref: cluster redshift reference value
  # for isolated galaxies zref can be the proper z value
  # radius: radius limit in kpc
  # velocity: velocity limit in km/s
  require(cosmoFns)
  #zc=z*2.99792458e5 # in km/s
  zc=2.99792458e5*(z-median(z))/(1+median(z))
  # degree to kpc conversion factor
  dkpc=cosmoFns::D.A(zref, omega.m=0.315, omega.lambda=0.685, H.0=67.4)*1000*pi/180 
  if(length(dkpc)==1) dkpc=rep(dkpc, length(ra))
  # Area in Mp^2
  A=pi*(radius/1000)^2 
  # loop over all galaxies
  ng=length(ra)
  sig=vector(length = ng)
  for (i in 1:ng) {
    dist=sqrt((ra-ra[i])^2+(dec-dec[i])^2)
    zdiff=zc-zc[i]
    w.inn=which(dist <= radius/dkpc[i] & dist > 0 & abs(zdiff) <= velocity & zdiff != 0)
    nph=length(which(is.na(z[w.inn])))
    nvel=ng-nph
    ngvel=length(w.inn)+nph*length(w.inn)/nvel
    sig[i]=ngvel/A
  }
  if(plot){
    require(astro)
    op=par("mar"=c(5.1,4.1,4.1,5.1), "oma"=c(0,0,0,0))
    aplot(ra, dec, sig, pch=20, xlab = 'ra (deg)', ylab = 'dec (deg)', cb=T, 
          cbspan = 1, zlim = round(range(sig)), zlab = expression(Sigma~(Mpc^{-2})),
          bgcol="grey98")
    par(op)
  }
  return(sig)
}

nldens=function(ra, dec, zref, n=5, plot=T){
  require(RANN)
  require(cosmoFns)
  nn=nn2(cbind(ra,dec), k=n+1)
  dN=nn$nn.dists[,n+1]
  deg2Mpc=D.A(zref, 0.3, 0.7, 67.4)*pi/180 
  dN=deg2Mpc*dN
  s=n/(pi*dN^2)
  if(plot){
    require(astro)
    op=par(mar=c(5.1,4.1,4.1,5.1), oma=c(0,0,0,0))
    lgd=log10(s); zlim=round(range(lgd, na.rm=T),1)
    aplot(ra, dec, lgd, pch=20, xlab='RA', ylab='DEC', cb=T, zlim=zlim, las=1,
          cbspan=1, zlab=expression(log[10](Sigma/Mpc^{-2})), zline=3.9, asp=1)
    par(op)
  }
  return(s)
}

# function to obtain local galaxy density for flat universe via Voronoi tesselation
tessdens=function(ra, dec, z, rlim=3, vlim=1e3, nmax=20, Om=0.3, Ol=1-Om, H0=73, 
                  plot=T){
  require(cosmoFns)
  require(RANN)
  require(deldir)
  n=length(ra)
  v=3e5*z 
  da=D.A(z, Om, Ol, H0)
  #d=3000/(77106*z^3-95094*z^2+77665*z+0.1)
  #d=rlim/(da*pi/180)
  den=vector(length = n)
  pb=txtProgressBar(0,n,0,'=',style=3)
  for(i in 1:n){
    s=abs(v-v[i])<=1000
    if(sum(s)<2) {den[i]=-Inf; next;}
    ras=ra[s]
    decs=dec[s]
    k=min(nmax+1, sum(s))
    nn=nn2(cbind(ras,decs), cbind(ra[i],dec[i]), k=k)
    #idx=nn$nn.idx[nn$nn.dists<d[i]]
    #if(length(idx)<nmax+1) idx=nn$nn.idx
    idx=nn$nn.idx
    vor=deldir(ras[idx], decs[idx])
    # plot(vor, F, 'tess')
    sm=vor$summary
    a=sm$dir.area[1]
    den[i]=log10(1/(a*da[i]^2))
    setTxtProgressBar(pb, i)
  }
  if(plot){
    require(astro)
    op=par(mar=c(5.1,4.1,4.1,5.1), oma=c(0,0,0,0))
    zlim=round(range(den, na.rm=T), 1)
    aplot(ra, dec, den, pch=20, xlab='RA (deg)', ylab='Dec (deg)', cb=T, zline = 3.9, 
          zlim=zlim, las=1, cbspan = 1, zlab = expression(log[10](Sigma/Mpc^{-2})))
    par(op)
  }
  return(den)
}

# function to convert SDSS designation JHHMMSS.ss+DDMMSS.s to ra & dec in degrees
sdssid2radec=function(id){
  require(celestial)
  id=as.character(id)
  des2deg=function(id){
    sg=ifelse(grepl('\\+', id), '+', '-')
    ss=strsplit(id, "\\+|\\-")[[1]]
    ss=c(ss[1], paste0(ss[2], '0', collapse = ''))
    ss=gsub('[[:punct:]]', '', ss)
    ss=strsplit(ss, "")
    ss=lapply(ss, function(x) paste0(x[c(TRUE, FALSE)], x[c(FALSE, TRUE)]))
    ss=lapply(ss, function(x) c(x[1:2], paste0(x[3:4], collapse='.')))
    ss=lapply(ss, function(x) paste0(x, collapse = ':'))
    sra=ss[[1]]; sdec=paste0(sg, ss[[2]], collapse = '')
    ra=hms2deg(sra)
    dec=dms2deg(sdec)
    return(c(ra, dec))
  }
  coo=lapply(id, des2deg)
  coo=data.frame(do.call('rbind', coo))
  names(coo)=c('ra', 'dec')
  return(coo)
}

# Function to detach all packages but the base ones
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils",
                      "package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  # run detachAllPackages()
}

# Function to determine the peaks and pits of a vector
turnpts=function(x, n=512, plot=T, ...){
  require(pastecs)
  #require(ks); bw=hns(x)
  d=density(x, n=n, na.rm = T, ...)
  tp=pastecs::turnpoints(d$y)
  pks=d$x[tp$peaks]
  pits=d$x[tp$pits]
  if(plot){
    plot(d, main='')
    abline(v=pks, col=2)
    abline(v=pits, col=4)
  }
  return(list(peaks=pks, pits=pits))
}

# Function to find green valley region by fitting 2 gaussians with mclust
valley=function(x, y, n=10, delta.y=0.5, method = 'content', plot=T){
  require(OneR)
  require(mclust)
  xb=bin(x, n, labels = 1:n, method = method)
  sp=split(data.frame(x=x,y=y), xb)
  xc=unlist(lapply(sp, function(d) (max(d$x)+min(d$x))/2))
  yv=m.min=m.max=vector(length = n)
  for(i in 1:n){
    dens=densityMclust(sp[[i]]$y, G=2, verbose=FALSE)
    pars=dens$parameters
    #var=pars$variance$scale
    sd=cbind(sp[[i]], dens=dens$density)
    ss=subset(sd, y > pars$mean[1] & y < pars$mean[2])
    yv[i]=ss$y[order(ss$dens)[1]]
    m.min[i]=pars$mean[1]
    m.max[i]=pars$mean[2]
    # plot(dens, what = "density", data = sp[[i]]$y, xlab='y')
    # abline(v=yv[i], lty=2)
  }
  lr=lm(yv~xc); coefs=coef(lr)
  # reg.exp=paste(round(b, 3), '* x', round(a, 3))
  # l=list(model=coefs)
  if(plot){
    a=coefs[1]; b=coefs[2]
    plot(x, y, pch=20, cex=.7); abline(lr, col=2, lwd=1)
    points(xc, yv, pch=20, col=6, cex=1)
    ey=unlist(lapply(sp, function(d) plotrix::std.error(d$y)))
    arrows(xc, yv-ey, xc, yv+ey, 0.05, 90, 3, 6)
    abline(a+delta.y, b, col=2, lty=2, lwd=1)
    abline(a-delta.y, b, col=2, lty=2, lwd=1)
  }
  data=data.frame(x=xc, y=yv, mean1=m.max, mean2=m.min)
  return(list(regression=lr, data=data))
}

# Splining a polygon
# https://github.com/emudrak/BMGImagery/blob/master/spline.poly.R
spline.poly=function(xy, vertices, k=3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}

# probability of a (velocity) distribution be Gaussian
p.gaussian=function(vlos, nsamples=500, parallel=T, plot=T){
  if(parallel){
    require(doParallel)
    cl=makeCluster(detectCores())
    registerDoParallel(cl)
  }
  require(distrEx)
  require(mclust)
  require(robustbase)
  
  n=length(vlos)
  message('Running Mclust')
  ncomp=Mclust(vlos, verbose = F)$G
  if(ncomp == 1 & n >= 20){
    # hdt=function(x){
    #   1/(2.7186565328-0.2277443222*x+0.0549035863*x^2-0.0010958319*x^3
    #      +0.0000198378*x^4)^0.1192224458
    # }
    hdt=function(x){
      if(x <= 200){
        (3.729394-1.333472e-02*x+3.774593e-02*x^2-4.297352e-04*x^3+1.158357e-05*x^4)^-1.245888e-01
      }else{
        (-2206.678+34.14281*x-0.2150239*x^2+0.001007499*x^3)^-0.1457817
      }
    }
    th=hdt(n)
    #n.vlos=qnorm(rank(vlos)/(n+1), mean = 0, sd = 1)
    n.vlos=(vlos-mean(vlos))/sd(vlos)
    
    message('Finding Hellinger Distances on bootstraps samples')
    cls=foreach(i=1:nsamples, .packages = c('distrEx','robustbase'), .combine = c) %dopar% {
      s=sample(n.vlos, n, replace = T)
      hd=HellingerDist(Norm(), s)
      # hd=HellingerDist(Norm(), s , asis.smooth.discretize = "smooth", h.smooth=s_Qn(s)/2)
    }
    
    #prob=1-sum(cls > th)/500
    prob=pnorm(th, mean(cls), sd(cls))
    
    if(plot){
      lims=c(ifelse(min(cls)>0,min(cls)-0.05,0), ifelse(th>max(cls),th+0.1,max(cls)+0.05))
      d=density(cls, bw=ks::hns(cls))
      plot(d, xlim=lims, main='', lty=2, xlab='HD')
      abline(v=c(th,mean(cls)), col=2:1)
      polygon(c(0,th,th,0), c(0,0,max(d$y),max(d$y)), col="#FF00001A", border = NA)
      legend('topright', lty=c(2,1,1), col=c(1,1,2), bty='n', cex=.9, y.intersp = .8,
             legend=c('HD Distribution','Mean HD Dist.','HD Threshold'))
    }
    
  }else{
    message('Multimodality found by Mclust')
    prob=-1
  }
  if(parallel) stopCluster(cl)
  
  return(prob)
}

# getting data from PhotoObjAll table in CASJOBS
getPhotoSDSS=function(objID, DR='DR12'){
  require(SciServer)
  Name = 'dailer';
  Password = 0x540a438
  token1=Authentication.login(Name, Password);
  Config.DataRelease=DR
  ng=length(objID)
  if(ng==1){
    q=paste('select * from PhotoObjAll where objID =', objID)
  }else{
    q1='select * from PhotoObjAll where '
    q2=paste('objID =',objID, ' or ', collapse = '')
    q=paste(q1, q2)
    q=substr(q,1,nchar(q)-5)
  }
  qs=SkyServer.sqlSearch(q, dataRelease = NULL)
  return(qs)
}

# Getting absolute magnitude
absmag=function(z, mag=17.77, evo=1.6*(z-0.1), OM=0.315, OL=1-OM, H0=67.4){
  require(CRAC)
  cosmo=list(omegaM0=OM, omegaL0=OL, omegaK=0, h=H0/100)
  Rgal=vapply(z, distance.transverse, 1, cosmo=cosmo)
  DistMod=5*log10(Rgal*(1+z))+25-evo
  M=mag-DistMod
  return(M)
}

absmag2=function(z, mag=17.77, evo=1.6*(z-0.1), OM=0.27, OL=1-OM, H0=71){
  require(cosmoFns)
  DistMod=5*log10(D.L(z, OM, OL, H0))+25-evo
  return(mag-DistMod)
}

# Kruskal-Wallis multicomparison test, taken from the pgirmess package
kruskalmc=function(resp, data=NULL, probs=0.05, cont=NULL){
  mf=model.frame(resp,data)
  resp=mf[,1]
  categ=mf[,2]
  db=na.omit(data.frame(resp,categ))
  if(nrow(db)!=length(resp)) 
    warning(paste(length(resp)-nrow(db),"lines including NA have been omitted"))
  resp=db[,1]
  categ=db[,2]
  lst=split(rank(resp), categ)
  name=names(lst)
  R=sapply(lst, mean)
  n=sapply(lst, length)
  N=length(resp)
  dif=abs(outer(R, R, "-"))
  if(is.null(cont)) {
    difv=vname=indices=NULL
    for(i in 1:(length(name)-1)){
      for(j in (i+1):length(name)) {
        vname=c(vname, paste(name[i], "-", name[j], sep = ""))
        indices=rbind(indices, c(i, j))
        difv=c(difv,dif[i,j])
      }
    }
    names(difv)=vname
    z=qnorm(probs/(length(lst)*(length(lst)-1)), lower.tail=F)
    lims=z*sqrt((N*(N+1)/12)*(1/n[indices[1:length(vname),1]]+1/n[indices[1:length(vname),2]]))
    names(lims)=vname
    stat="Multiple comparison test after Kruskal-Wallis"
  }else{
    vname=indices=NULL
    for(j in 2:length(dif[1, ])){
      vname=c(vname, paste(name[1], "-", name[j], sep = ""))
      indices=rbind(indices, c(1, j))
    }
    dif=dif[1, 2:length(dif[1, ])]
    names(dif)=vname
    difv=dif
    choice=pmatch(cont, c("two-tailed","one-tailed"), nomatch=3)
    if(choice==1) {
      z=qnorm(probs/(2*(length(lst)-1)), lower.tail=F)
      lims=z*sqrt((N*(N + 1))/12*(1/n[indices[1:length(vname), 1]]+
                                    1/n[indices[1:length(vname), 2]]))
      names(lims)=vname
      stat="Multiple comparison test after Kruskal-Wallis, treatments vs control (two-tailed)"
    }
    if(choice==2){
      z=qnorm(probs/(length(lst)-1), lower.tail=F)
      lims=z*sqrt((N*(N+1)/12)*(1/n[indices[1:length(vname), 1]]+
                                  1/n[indices[1:length(vname), 2]]))
      names(lims)=vname
      stat="Multiple comparison test after Kruskal-Wallis, treatment vs control (one-tailed)"
    }
    if(choice==3) 
      stop("Values must be 'one-tailed' or 'two-tailed', partial matching accepted")
  }
  output=list(statistic=stat, signif.level=probs, 
              dif.com=data.frame(obs.dif=difv, critical.dif=lims, 
                                 difference=ifelse((difv-lims)>0, TRUE, FALSE)))
  class(output)=c("mc", "list")
  output
}

# Geometric histogram separation from Sutter & Barchi (2017)
ghs=function(x1, x2, plot=T){
  h1=hist(x1, plot = F)
  h2=hist(x2, plot = F)
  both=c(h1$breaks, h2$breaks)
  n=length(both)-1 # number of total bins
  rnge=range(both)
  breaks=seq(rnge[1], rnge[2], length.out = n+1)
  h1=hist(x1, breaks, plot = F)
  h2=hist(x2, breaks, plot = F)
  dx=diff(rnge)/n
  dy=apply(cbind(h1$density, h2$density), 1, min) # min of intersection
  ao=sum(dx*dy) # the relative area 
  a_height=max(h1$density)
  b_height=max(h2$density)
  c_height=max(dy)  
  bcl=(a_height+b_height-2*c_height)/(a_height+b_height)
  bca=1-ao/(2-ao)
  ghs=(sqrt(bca)+bcl)/2
  if(plot){
    ymax=max(h1$counts,h2$counts)
    plot(c(rnge), c(0, ymax), pch='', ylab='Frequency', xlab='x')
    abline(h=0, col='gray')
    lines(c(breaks[1],breaks), c(0,h1$counts,0), type='s', col=2) 
    lines(c(breaks[1],breaks), c(0,h2$counts,0), type='s', col=4) 
  }
  return(ghs)
}

# Kruskal-Wallis multicomparison test for multiple properties
kw=function(data, class, probs=0.05, plot=T, ...){
  #print(call)
  require(pgirmess)
  n=ncol(data)
  k=vector('list', length = n)
  for(i in 1:n){
    k[[i]]=kruskalmc(unlist(data[,..i])~class, probs=probs)$dif.com
  }
  rn=rownames(k[[1]])
  cn=colnames(data)
  kl=lapply(k, function(x) x[,3])
  ku=data.table(rn, do.call(cbind, kl))
  colnames(ku)=c('comp', cn)
  if(plot){
    op=par(mai=c(.9,.9,0.3,0.3))
    image(x=1:n, y=1:length(rn), t(as.matrix(ku[length(rn):1,-1])), col=c(1,'white'), 
          xlab = '', ylab = 'classes', axes=F, cex.lab=1, asp=1, mgp=c(2,1,0), ...)
    # ylab = 'Kruskal-Wallis intraclass comparison'
    # magicaxis::magimage(x=1:n, y=1:length(rn), t(as.matrix(ku[length(rn):1,-1])), col=c(1,'white'), 
    #                     xlab = 'data columns', ylab = 'intraclass comparison', axes=F, magmap = F)
    #axis(1, at=1:n, labels = cn, cex.axis=.8)
    
    #axis(1, at=1:n, labels = 1:n, mgp=c(2,.5,0), cex.axis=.8)
    axis(1, at=1:n, labels = F)
    text(1:n, par("usr")[3]-0.2, labels=colnames(data), srt=45, pos=1, xpd=T, cex=.8)
    axis(2, at=1:length(rn), labels = rev(rn), las=1, mgp=c(2,.5,0), cex.axis=.75)
    segments((0:n)+0.5,0.5,(0:n)+0.5,length(rn)+0.5, col='brown', lwd=2)
    segments(0.5,(0:length(rn))+0.5,n+0.5,(0:length(rn))+0.5, col='brown', lwd=2)
    #abline(v=(0:n)+0.5, h=(0:length(rn))+0.5)
    box()
    #title('Kruskal-Wallis multicomparison result')
    # mtext('White cell = TRUE', adj = 0, cex=.7)
    #mtext(paste('TRUE fraction =',round(length(which(ku[,-1]==T))/prod(dim(ku[,-1])), 3)), adj = 1)
    par(op)
  }
  return(ku)
}

medist=function(data, class, plot=T){
  # require ghs
  cl=sort(unique(class))
  cb=combn(cl, 2)
  nc1=ncol(cb); nc2=ncol(data)
  m=matrix(NA, nc1, nc2)
  rownames(m)=apply(cb, 2, paste, collapse='-')
  colnames(m)=colnames(data)
  for(i in 1:nc1){
    cp=cb[,i]
    id1=cp[1]; id2=cp[2]
    for(j in 1:ncol(data)){
      v=unlist(data[,..j])
      v1=v[class==id1]
      v2=v[class==id2]
      #d=abs(median(v1,T)-median(v2,T))
      #m[i,j]=d/(max(v, na.rm=T)-min(v, na.rm=T))
      m[i,j]=ghs(v1, v2, plot = F)
    }
  }
  if(plot){
    image(1:nc2, 1:nc1, t(m[nc1:1,]), col = hcl.colors(length(m)), 
          asp=1, xlab = '', ylab = '', axes=F)
    axis(1, 1:nc2, colnames(m), las=2)
    axis(2, 1:nc1, rev(rownames(m)), las=2)
    box()
  }
  return(round(m/max(m) ,3))
}

# multiple comparison using median and quantile error of the median
mc=function(data, class, plot=T){
  qf=function(x) return(0.7415*(quantile(x,.75, na.rm=T)-quantile(x,.25, na.rm=T)))
  sp=split(data, class)
  l=length(sp)
  med=lapply(sp, function(x) apply(x, 2, median, na.rm=T))
  err=lapply(sp, function(x) apply(x, 2, qf))
  med=data.frame(do.call(rbind, med))
  err=data.frame(do.call(rbind, err))
  colnames(err)=paste0(colnames(err),'_err')
  
  n=ncol(data)
  cb=combn(names(sp), 2)
  nms=apply(cb, 2, paste, collapse = '-')
  m=matrix(nrow = length(nms), ncol = n)
  rownames(m)=nms; colnames(m)=colnames(med)
  for(i in 1:ncol(cb)){
    m[i,]=abs(med[cb[1,i],]-med[cb[2,i],]) > err[cb[1,i],]+err[cb[2,i],]
  }
  if(length(nms)>1) m=m[order(nms),]
  rn=sort(nms)
  if(plot){
    z=t(as.matrix(m[length(nms):1,]))
    if(nrow(m)==1) z=t(z)
    image(x=1:n, y=1:length(nms), z, col=c(1,'white'), asp=1, mgp=c(2,1,0),
          xlab = 'data columns', ylab = '', axes=F)
    axis(1, at=1:n, labels = 1:n, mgp=c(2,.5,0), cex.axis=.7)
    axis(2, at=1:length(nms), labels = rev(rn), las=1, mgp=c(2,.5,0), cex.axis=.7)
    segments((0:n)+0.5,0.5,(0:n)+0.5,length(nms)+0.5, col='#17377A', lwd=2)
    segments(0.5,(0:length(nms))+0.5,n+0.5,(0:length(nms))+0.5, col='#17377A', lwd=2)
    box()
  }
  med=round(med, 3); err=round(err, 3)
  t1=t(med); colnames(t1)=paste0('median.',colnames(t1))
  t2=t(err); colnames(t2)=paste0('error.',colnames(t2))
  m=as.matrix(t(m)); colnames(m)=paste0('comp.',rn)
  tme=data.frame(t1, t2)
  st=data.frame(matrix(paste(t1, t2, sep='±'), nrow(t1), ncol(t1)))
  dimnames(st)=list(rownames(t1), rownames(med))
  return(list(values=tme, intraclass_comparison=m, string=st))
}

kwdm=function(data, class, probs=0.05, names=NULL, plot=T, sortnames=T, ...){
  #depend on kw(), medist(), mc() and ghs() functions
  require(pgirmess)
  require(plotrix)
  if(is.null(names)) names=colnames(data)
  k=kw(data, class, probs, plot=F)
  m=medist(data, class, plot=F)
  me=mc(data, class, plot=F)
  #cn=colnames(data)
  rn=gsub('-','×',k$comp)
  n=ncol(data); l=length(rn)
  tk=t(as.matrix(k[l:1,-1]))
  g=expand.grid(1:n,1:l)
  g$f=as.vector(t(m[l:1,]))
  #g$d=ifelse(as.vector(t(me[[3]][l:1,])), 1, 2)
  g$dd=ifelse(as.vector(me[[2]][,l:1]), 1, 2)
  g=g[c(tk),]
  if(plot){
    op=par(mai=c(1.2,1.2,.1,.1), ...)
    image(1:n, 1:l, tk, col=c(1,'white'), xlab='', ylab='', axes=F, asp=1)
    axis(1, at=1:n, labels = names, las=2, cex.axis=1); box()
    axis(2, at=1:l, labels = rev(rn), las=1, mgp=c(2,.5,0), cex.axis=1)
    segments((0:n)+0.5,0.5,(0:n)+0.5,l+0.5, lwd=1)
    segments(0.5,(0:l)+0.5,n+0.5,(0:l)+0.5, lwd=1)
    #points(g[,1:2], pch=1, cex=g[,3]*144/(n*l))
    clr=c('gray','white')
    dc=apply(g, 1, function(x) draw.circle(x[1], x[2], .45*x[3], border=1, col=clr[x[4]]))
    par(op)
  }
  lout=list(KW=k, relative_distances=m)
  return(invisible(c(lout,me)))
}

plot.kw=function(x, rows='all', sortnames=F){
  k=x$KW
  m=x$relative_distances
  me=x[3:5]
  if(is.numeric(rows)){
    k=k[rows,]
    m=m[rows,]
    me[[2]]=me[[2]][,rows]
  }
  rn=gsub('-','×',k$comp)
  n=ncol(k)-1 
  l=length(rn)
  tk=t(as.matrix(k[l:1,-1]))
  g=expand.grid(1:n,1:l)
  g$f=as.vector(t(m[l:1,]))
  g$dd=ifelse(as.vector(me[[2]][,l:1]), 1, 2)
  g=g[c(tk),]
  names=colnames(k)[-1]
  # Plotting...
  op=par(mai=c(1.2,1.2,.1,.1))
  image(1:n, 1:l, tk, col=c(1,'white'), xlab='', ylab='', axes=F, asp=1)
  axis(1, at=1:n, labels = names, las=2, cex.axis=1); box()
  axis(2, at=1:l, labels = rev(rn), las=1, mgp=c(2,.5,0), cex.axis=1)
  segments((0:n)+0.5,0.5,(0:n)+0.5,l+0.5, lwd=1)
  segments(0.5,(0:l)+0.5,n+0.5,(0:l)+0.5, lwd=1)
  clr=c('gray','white')
  dc=apply(g, 1, function(x) draw.circle(x[1], x[2], .45*x[3], border=1, col=clr[x[4]]))
  par(op)
  return(invisible(k))
}

# permutation test for multiple variables (and two classes)
perm=function(data, class, test=c('permTS','ks','t','wilcox')){
  require(perm)
  if(class(data)[1]!='data.frame')
    data=as.data.frame(data)
  sp=split(data, class)
  #nsp=names(sp)
  #nms=paste(nsp, collapse = '-')
  nc=ncol(data)
  p=rep(NA, nc)
  sp1=sp[[1]]; sp2=sp[[2]]
  pb=txtProgressBar(0, nc, style=3)
  for(i in 1:nc){
    v1=sp1[,i]
    v2=sp2[,i]
    if(test[1]=='permTS')
      p[i]=permTS(v1[!is.na(v1)], v2[!is.na(v2)])$p.value
    if(test[1]=='ks')
      p[i]=ks.test(v1[!is.na(v1)], v2[!is.na(v2)])$p.value
    if(test[1]=='t')
      p[i]=t.test(v1[!is.na(v1)], v2[!is.na(v2)])$p.value
    if(test[1]=='wilcox')
      p[i]=wilcox.test(v1[!is.na(v1)], v2[!is.na(v2)])$p.value
    setTxtProgressBar(pb, i)
  }
  d=data.frame(p.value=round(p, 4), signif=p<=0.05)
  rownames(d)=colnames(data)
  return(d)
}

multperm=function(data, class, plot=T){
  require(rcompanion)
  cn=colnames(data)
  nc=length(cn)
  cb=choose(length(unique(class)),2)
  ptm1=ptm2=matrix(NA, nc, cb)
  pb=txtProgressBar(0, nc, style=3)
  for(i in 1:nc){
    x=eval(parse(text=paste0("data$", cn[i])))
    PT=pairwisePermutationTest(x=x, g=class)
    ptm1[i,]=PT$p.adjust<0.05
    ptm2[i,]=PT$p.value<0.05
    setTxtProgressBar(pb, i)
  }
  tm1=ptm1[,cb:1]
  tm2=ptm2[,cb:1]
  comp=gsub(' = 0', '', PT$Comparison)
  comp=gsub(' - ', '×', comp)
  if(plot){
    op=par(mai=c(1.2,1.2,.1,.1))
    image(1:nc, 1:cb, tm2, col=c('gray95','white'), xlab='', ylab='', axes=F, asp=1)
    image(1:nc, 1:cb, tm1, col=c(1,'#FFFFFF00'), add=T)
    axis(1, at=1:nc, labels = cn, las=2, cex.axis=1); box()
    axis(2, at=1:cb, labels = rev(comp), las=1, mgp=c(2,.5,0), cex.axis=1)
    segments((0:nc)+0.5,0.5,(0:nc)+0.5,cb+0.5, lwd=1)
    segments(0.5,(0:cb)+0.5,nc+0.5,(0:cb)+0.5, lwd=1)
    par(op)
  }
}

# Kernel density based global two-sample comparison test (kde.test) by classes
kd=function(data, class, verbose=T, ...){
  require(ks)
  sp=split(data, class)
  n=length(sp)
  nsp=names(sp)
  cb=combn(1:n, 2)
  ncb=matrix(NA, nrow(cb), ncol(cb))
  ncb[]=nsp[cb]
  nms=apply(ncb, 2, paste, collapse = '-')
  k=vector(length = ncol(cb))
  for(i in 1:ncol(cb)){
    if(verbose) cat('comparing',nms[i],'\n')
    k[i]=kde.test(sp[[cb[1,i]]], sp[[cb[2,i]]], verbose=verbose, ...)$pvalue
  }
  if(verbose) cat('Done.\n')
  ks=ifelse(k < 0.001, '< 0.001', as.character(round(k,4)))
  difs=ifelse(k < 0.05, TRUE, FALSE)
  df=data.frame('p-value'=ks, 'difference'=difs)
  rownames(df)=nms
  return(df)
}

# Biweight location estimator
biwLoc=function(x){
  M=median(x, na.rm=T); d=x-M
  u=d/(6*mad(x, M, 2.21914, na.rm=T))
  uT=(1-u^2)^2; wu=which(abs(u)<1)
  top=sum((d*uT)[wu])
  bottom=sum(uT[wu])
  CBI=M+top/bottom
  if(is.nan(CBI)) CBI=x[1]
  return(CBI)
}

# Biweight scale estimator
biwScale=function(x){
  M=median(x, na.rm=T); d=x-M
  u=d/(9*mad(x, M, 2.21914, na.rm=T))
  uT=1-u^2; wu=which(abs(u)<1)
  top=sqrt(sum((d^2*uT^4)[wu]))
  bottom=abs(sum((uT*(1-5*u^2))[wu]))
  SBI=sqrt(length(wu))*top/bottom
  if(is.nan(SBI)) SBI=x[1]
  return(SBI)
}

# VDP using Bergond et al. (2006) prescription
vdp=function(dproj, vlos, Nb=1000, N=100, norm=F, plot=T, smooth=T, parallel='no', ...){
  # require functions biwLoc & biwScale
  # vlos must be the LOS velocity NORMALIZED to vdisp
  require(astro)
  require(boot)
  require(ks)
  ng=length(dproj)
  vlos=vlos[order(dproj)]; dproj=sort(dproj)
  r=seq(0, max(dproj), length.out = N)
  #r=seq(min(dproj), max(dproj), length.out = N)
  #hw=biwScale(dproj)
  hw=2*hns(dproj, 2)
  vmean=biwLoc(vlos)
  sigma_v=(vlos-vmean)^2
  bf=function(x, idx, w.) sqrt(sum(x[idx], na.rm=T)/w.)
  lc=lwr=upr=vector(length = N)
  pb=txtProgressBar(0, N, style=3)
  for(i in 1:N){
    wi=(1/hw)*exp(-(r[i]-dproj)^2/(2*hw^2))
    swi=sum(wi)
    w_sigma_v=wi*sigma_v
    b=boot(w_sigma_v, bf, R=Nb, w.=swi, parallel = parallel)
    lc[i]=b$t0
    err=biwScale(b$t)
    lwr[i]=b$t0-err
    upr[i]=b$t0+err
    setTxtProgressBar(pb, i)
  }
  if(smooth){
    lwr=smooth.spline(r, lwr)$y
    upr=smooth.spline(r, upr)$y
    lc=smooth.spline(r, lc)$y
  }
  if(plot){
    aplot(c(min(r), max(r)), c(min(lwr), max(upr)), pch='', xlab = expression(R/R[200]), las=1, 
          ylab = expression(sigma[p]/sigma[R[200]]), main='Velocity Dispersion Profile',...)
    polygon(c(r, rev(r)), c(lwr, rev(upr)), col='#0000001A', border = NA)
    lines(r, lc, col=1)
  }
  df=data.frame(r_norm=r, s_r=lc, lwr=lwr, upr=upr)
  return(invisible(df))
}

# VDP using Bergond et al. (2006) prescription
vdp_=function(dproj, vlos, Nb=1000, N=100, norm=F, plot=T){
  # require functions biwLoc & biwScale
  # vlos must be the LOS velocity NORMALIZED to vdisp
  require(astro)
  ng=length(dproj)
  vlos=vlos[order(dproj)]; dproj=sort(dproj)
  r=seq(0, max(dproj), length.out = N)
  #r=seq(min(dproj), max(dproj), length.out = N)
  hw=biwScale(dproj)
  #hw=2*ks::hns(dproj, 2)
  #hw=3*bw.nrd0(dproj)
  #print(c(hw))
  vmean=biwLoc(vlos)
  sigma_v=(vlos-vmean)^2
  sv=vector(length = Nb)
  sr=matrix(NA, N, Nb)
  pb=txtProgressBar(0, N, style=3)
  for(i in 1:N){
    wi=(1/hw)*exp(-(r[i]-dproj)^2/(2*hw^2))
    swi=sum(wi)
    w_sigma_v=wi*sigma_v
    for(j in 1:Nb){
      ws=sample(w_sigma_v, ng, replace=T)
      sv[j]=sqrt(sum(ws, na.rm=T)/swi)
    }
    sr[i,]=sv
    setTxtProgressBar(pb, i)
  }
  s=apply(sr, 1, biwLoc)
  if(norm){
    vcen=biwScale(vlos[dproj<1])
    s=s/vcen
  }
  sds=apply(sr, 1, biwScale)
  lwr=as.numeric(s-sds); upr=as.numeric(s+sds)
  lwr=smooth.spline(r, lwr)$y
  upr=smooth.spline(r, upr)$y
  ss=smooth.spline(r, s)$y
  if(plot){
    aplot(c(min(r), max(r)), c(min(lwr), max(upr)), pch='', xlab = expression(R/R[200]), las=1, 
          ylab = expression(sigma[p]/sigma[R[200]]), main='Velocity Dispersion Profile')
    polygon(c(r, rev(r)), c(lwr, rev(upr)), col='gray90', border = 'gray90')
    lines(r, ss, col=1)
  }
  df=data.frame(r_norm=r, s_r=ss, lwr=lwr, upr=upr)
  return(invisible(df))
}

# VDP using Bergond et al. (2006) prescription
vdp__=function(dproj, vlos, N=100, sigma=0.15, Nb=1000, plot=T){
  require(astro)
  ng=length(dproj)
  vlos=vlos[order(dproj)]; dproj=sort(dproj)
  r=seq(0, max(dproj), length.out = N)
  #r=seq(min(dproj), max(dproj), length.out = N)
  #hw=max(dproj)*sigma
  hw=biwScale(dproj)
  sigma_v=vlos^2
  sv=vector(length = N)
  s.r=vector("list", length = Nb)
  for (j in 1:Nb) {
    for (i in 1:N) {
      wi=(1/hw)*exp(-(r[i]-dproj)^2/(2*hw^2))
      if (j == 1) {
        ws=wi*sigma_v
      }else{
        ws=sample(wi*sigma_v, ng, replace = T)
      }
      sv[i]=sqrt(sum(ws)/sum(wi))
    }
    s.r[[j]]=sv
  }
  s.r=as.data.frame(do.call('rbind', s.r))
  #s.r=s.r/sigma200; 
  s=as.numeric(s.r[1,])
  #ci=function(x) qnorm(0.975)*sd(x)/sqrt(1000)
  qf=function(x) return(0.7415*(quantile(x,.75, na.rm=T)-quantile(x,.25, na.rm=T)))
  #ci=function(x) qt((1-0.95)/2, length(x)-1)*qf(x)
  sds=apply(s.r, 2, qf)#; sds=abs(sds)
  lwr=as.numeric(s-sds); upr=as.numeric(s+sds)
  lwr=smooth.spline(r, lwr)$y
  upr=smooth.spline(r, upr)$y
  rp=r #rp=r/r200
  if(plot){
    aplot(c(min(rp), max(rp)), c(min(lwr), max(upr)), pch='', xlab = expression(Radius/R[200]), 
          ylab = expression(sigma[p]/sigma[r[200]]), main='Velocity Dispersion Profile')
    polygon(c(rp, rev(rp)), c(lwr, rev(upr)), col='gray90', border = 'gray90')
    #points(rp, s, pch=20, type = 'b', cex=0.8, col=1)
    lines(rp, s, col=1)
  }
  return(data.frame(r_norm=rp, s_r=s, lwr=lwr, upr=upr))
}

# cumulative VDP
cvdp=function(dproj, vlos, Nb=1000, spar=.8, plot=T, smooth=T, ...){
  # require function biwScale
  # vlos must be the LOS velocity NORMALIZED to vdisp
  require(astro)
  ng=length(dproj)
  vlos=vlos[order(dproj)]; dproj=sort(dproj)
  mu=sig=vector(length = ng-9)
  bf=function(x, idx) biwScale(x[idx])
  #bf=function(x, idx) median(x[idx])
  pb=txtProgressBar(10, ng, style=3)
  for(i in 10:ng){
    v=vlos[1:i]
    b=boot(v, bf, R=Nb, ...)
    mu[i-9]=b$t0
    sig[i-9]=biwScale(b$t)
    setTxtProgressBar(pb, i)
  }
  lwr=mu-sig; upr=mu+sig
  r=dproj[10:ng]
  if(smooth){
    mu.s=smooth.spline(r, mu, spar=spar)$y
    lwr.s=smooth.spline(r, lwr, spar=spar)$y
    upr.s=smooth.spline(r, upr, spar=spar)$y
  }
  if(plot){
    tit='Cumulative Velocity Dispersion Profile'
    aplot(range(r), c(min(lwr), max(upr)), pch='', xlab = expression(Radius/R[200]), 
          ylab = expression(sigma[p]/sigma[r[200]]), main=tit)
    polygon(c(r, rev(r)), c(lwr, rev(upr)), col='gray90', border = NA)
    lines(r, mu, col='gray80', lty=1); box()
    if(smooth){
      lines(r, mu.s); lines(r, lwr.s, lty=2); lines(r, upr.s, lty=2)
    }
  }
  df=data.frame(rad=r, sig=mu, lwr=lwr, upr=upr)
  return(invisible(df))
}

# cumulative VDP
cvdp_=function(dproj, vlos, Nb=1000, spar=.8, plot=T){
  # require functions biwLoc & biwScale
  # vlos must be the LOS velocity NORMALIZED to vdisp
  require(astro)
  ng=length(dproj)
  vlos=vlos[order(dproj)]; dproj=sort(dproj)
  #vcen=biwScale(vlos[dproj<1])
  s=matrix(NA, Nb, ng-9)
  fb=function(x) biwScale(sample(x, length(x), T))
  pb=txtProgressBar(10, ng, style=3)
  for (i in 10:ng){
    v=vlos[1:i]
    for(j in 1:Nb){
      s[j,i-9]=biwScale(sample(v, i, T))
    }
    setTxtProgressBar(pb, i)
  }
  mu=apply(s, 2, biwLoc)#/vcen
  sig=apply(s, 2, biwScale)
  lwr=mu-sig; upr=mu+sig
  r=dproj[10:ng]; rp=r
  if(plot){
    tit='Cumulative Velocity Dispersion Profile'
    ss1=smooth.spline(rp, mu, spar=spar)
    ss2=smooth.spline(rp, lwr, spar=spar)
    ss3=smooth.spline(rp, upr, spar=spar)
    aplot(range(rp), c(min(ss2$y), max(ss3$y)), pch='', xlab = expression(Radius/R[200]), 
          ylab = expression(sigma[p]/sigma[r[200]]), main=tit)
    polygon(c(rp, rev(rp)), c(lwr, rev(upr)), col='gray90', border = NA)
    lines(rp, mu, col='gray80', lty=1); box()
    #lines(lowess(rp, mu, 0.2))
    lines(ss1); lines(ss2, lty=2); lines(ss3, lty=2)
  }
  return(data.frame(rad=rp, sig=mu, lwr=lwr, upr=upr))
}

# Velocity kurtosis profile using Bergond et al. (2006) prescription
kvp=function(dproj, vlos, Nb=1000, N=100, plot=T){
  # require functions biwLoc & biwScale
  # vlos must be the LOS velocity NORMALIZED to vdisp
  w.mean=function(x, mu) sum(mu*x)/sum(mu)
  w.sd=function(x, mu) ((sum(mu*x^2)/sum(mu))-w.mean(x, mu)^2)^0.5
  w.kurtosis=function(x, mu) ((sum(mu*(x-w.mean(x,mu))^4)/sum(mu))/w.sd(x,mu)^4)-3
  kurtosis=function(x) length(x)*sum((x-mean(x))^4)/(sum((x-mean(x))^2)^2)-3
  require(astro)
  ng=length(dproj)
  vlos=vlos[order(dproj)]; dproj=sort(dproj)
  k0=kurtosis(vlos[dproj<1])
  r=seq(0, max(dproj), length.out = N)
  hw=biwScale(dproj)
  sv=vector(length = Nb)
  sr=matrix(NA, N, Nb)
  pb=txtProgressBar(0, N, style=3)
  for(i in 1:N){
    wi=(1/hw)*exp(-(r[i]-dproj)^2/(2*hw^2))
    for(j in 1:Nb){
      si=sample(1:ng, ng, replace=T)
      sv[j]=w.kurtosis(vlos[si], wi[si])
    }
    sr[i,]=sv
    setTxtProgressBar(pb, i)
  }
  s=apply(sr, 1, biwLoc)/k0
  sds=apply(sr, 1, biwScale)
  lwr=as.numeric(s-sds); upr=as.numeric(s+sds)
  lwr=smooth.spline(r, lwr)$y
  upr=smooth.spline(r, upr)$y
  ss=smooth.spline(r, s)$y
  if(plot){
    aplot(c(min(r), max(r)), c(min(lwr), max(upr)), pch='', xlab = expression(R/R[200]), las=1, 
          ylab = expression(g[u]/g[u[R200]]), main='Velocity Kurtosis Profile')
    polygon(c(r, rev(r)), c(lwr, rev(upr)), col='gray90', border = 'gray90')
    lines(r, ss, col=1)
  }
  return(data.frame(r_norm=r, s_r=ss, lwr=lwr, upr=upr))
}

# Weighted median (continuous) profiles of y along with x
varprof=function(x, y, Nb=1000, N=100, err.type='q', norm=F, plot=T, smooth=T){
  # require functions biwLoc & biwScale
  require(astro)
  require(spatstat)
  n=length(x)
  y=y[order(x)]; x=sort(x)
  #r=seq(0, max(x), length.out = N)
  r=seq(min(x), max(x), length.out = N)
  hw=biwScale(x)
  sv=vector(length = Nb)
  sr=matrix(NA, N, Nb)
  pb=txtProgressBar(0, N, style=3)
  for(i in 1:N){
    wi=(1/hw)*exp(-(r[i]-x)^2/(2*hw^2))
    for(j in 1:Nb){
      si=sample(1:n, n, replace=T)
      sv[j]=weighted.median(y[si], wi[si])
    }
    sr[i,]=sv
    setTxtProgressBar(pb, i)
  }
  if(err.type=='q'){
    s=apply(sr, 1, biwLoc)
    sds=apply(sr, 1, biwScale)
  }
  if(err.type=='ci'){
    require(Rmisc)
    ci=apply(sr, 1, CI)
    s=apply(ci, 2, function(x) x[2])
    sds=apply(ci, 2, function(x) x[1]-x[2])
  }
  if(norm){
    ycen=biwScale(y[x<1])
    s=s/ycen
  }
  lwr=as.numeric(s-sds); upr=as.numeric(s+sds)
  if(smooth){
    lwr=smooth.spline(r, lwr)$y
    upr=smooth.spline(r, upr)$y
    s=smooth.spline(r, s)$y
  }
  if(plot){
    aplot(c(min(r), max(r)), c(min(lwr), max(upr)), pch='', xlab='x', las=1, ylab='y')
    polygon(c(r, rev(r)), c(lwr, rev(upr)), col='gray90', border = 'gray90')
    lines(r, s, col=1)
  }
  df=data.frame(x=r, y=s, lwr=lwr, upr=upr)
  return(invisible(df))
}

# perform nearest neighbor search in casjobs
CASNN=function(ra, dec, DR='DR12', tablename='temp', verbose=T){
  # maximum recommended: 20000
  require(data.table)
  require(SciServer)
  tbldrp=tryCatch(SkyQuery.dropTable(tablename), error=function(e){})
  Name = 'dailer';
  Password = 0x540a438
  Authentication.login(Name, Password);
  radec=data.frame(id=1:length(ra), ra=ra, dec=dec)
  if(verbose) message('Uploading data to MyDB')
  CasJobs.uploadDataFrameToTable(radec, tablename)
  q1='SELECT m.id, n.objid, n.distance, o.petroMag_u-o.extinction_u as u_pmag, 
  o.petroMag_g-o.extinction_g as g_pmag, o.petroMag_r-o.extinction_r as r_pmag, 
  o.petroMag_i-o.extinction_i as i_pmag, o.petroMag_z-o.extinction_z as z_pmag
  FROM '
  #E_bv=extinction_r/2.751 FROM ' #kcorrR01 as k01
  q2=' AS m OUTER APPLY dbo.fGetNearestObjEq(m.ra, m.dec, 0.5) AS n 
  LEFT JOIN PhotoObj AS o ON n.objid=o.objid JOIN Photoz as h 
  on h.objid=o.objid order by id'
  query=paste0(q1, 'MyDB.', tablename, q2)
  if(verbose) message('Executing nearest neighbor search')
  q=CasJobs.executeQuery(query, context=DR, format="csv")
  if(verbose) message('Done')
  q=fread(q)
  q[q=='null' | q==-9999]=NA
  tryCatch(SkyQuery.dropTable(tablename, datasetName="MyDB"), 
           error=function(e) message('Warning: could not drop table from MyDB'))
  return(q)
}

CASNNloop=function(ra, dec, n=20000, DR='DR12', tablename='temp'){
  data=data.table(ra,dec)
  nr=nrow(data)
  l=ceiling(nr/n)
  sp=split(data, rep(1:l, each=n, length.out=nr))
  nrs=unlist(lapply(sp, nrow))
  nrs=cumsum(nrs)-nrs
  dl=vector('list', length = l)
  pb=txtProgressBar(0, l, 0, '=', style=3)
  for(i in 1:l){
    d=sp[[i]]
    nn=CASNN(d$ra, d$dec, DR=DR, tablename=tablename, verbose = F)
    nn$id=nn$id+nrs[i]
    dl[[i]]=nn
    setTxtProgressBar(pb, i)
  }
  res=data.table(do.call(rbind, dl))
  return(res)
}

# search for some photometric variables from photo objid in casjobs
CASPH=function(objid, DR='DR12', tablename='temp'){
  require(data.table)
  require(SciServer)
  require(bit64, quietly = T)
  tbldrp=tryCatch(SkyQuery.dropTable(tablename), error=function(e){})
  Name = 'dailer';
  Password = 0x540a438
  Authentication.login(Name, Password);
  message('Uploading data to MyDB')
  objid[is.na(objid)]=1e18
  data=data.frame(id=1:length(objid), objid=as.integer64(objid))
  CasJobs.uploadDataFrameToTable(data, tablename)
  query='select
  id, t.objid, ra, dec, run, rerun, camcol, field, petroMag_u as upmag, 
  petroMagErr_u as upmag_err, petroMag_g as gpmag, petroMagErr_g as gpmag_err, 
  petroMag_r as rpmag, petroMagErr_r as rpmag_err, petroMag_i as ipmag, 
  petroMagErr_i as ipmag_err, petroMag_z as zpmag, petroMagErr_z as zpmag_err, 
  extinction_u as ext_u, extinction_g as ext_g, extinction_r as ext_r, 
  extinction_i as ext_i, extinction_z as ext_z
  from MyDB.temp as t 
  left join PhotoObjAll as p on p.objID=t.objid 
  order by id'
  message('Executing query')
  q=CasJobs.executeQuery(query, context=DR, format="csv")
  message('Done')
  q=fread(q)
  q[q=='null' | q==-9999 | q=='-9999']=NA
  q[,3:ncol(q)]=lapply(q[,3:ncol(q)], as.numeric)
  tryCatch(SkyQuery.dropTable(tablename), 
           error=function(e) message('Warning: could not drop table from MyDB'))
  return(q)
}

# search for some spectroscopic variables from photo objid in casjobs
CASPEC=function(objid, DR='DR12', tablename='temp', verbose=T){
  require(data.table)
  require(SciServer)
  require(bit64, quietly = T)
  tbldrp=tryCatch(SkyQuery.dropTable(tablename), error=function(e){})
  Name = 'dailer';
  Password = 0x540a438
  Authentication.login(Name, Password);
  objid[is.na(objid)]=1e18
  data=data.frame(id=1:length(objid), objid=as.integer64(objid))
  if(verbose) message('Uploading data to MyDB')
  CasJobs.uploadDataFrameToTable(data, tablename)
  dr=as.integer(strsplit(DR, 'DR')[[1]][2])
  if(dr>=12){
    query='select
    id, s.bestobjID as objid, s.specObjID as specobjid, s.z as zspec, s.zErr,
    zWarning as zW, sciencePrimary as prim, s.plate, s.mjd, s.fiberid, 
    snMedian as SN, bptclass as bpt, lgm_tot_p50 as smass, sfr_tot_p50 as sfr, 
    specsfr_tot_p50 as ssfr, oh_p50 as oh, d4000_n as d4000n, lick_hb as hb, 
    lick_hd_a as hd, lick_hg_a as hg, h_delta_eqw as Hd_EW, h_gamma_eqw as Hg_EW,
    h_beta_eqw as Hb_EW, h_alpha_eqw as Ha_EW, EW_Ha_6562, Flux_OIII_5006, 
    Flux_SII_6716, Flux_SII_6730, Flux_Ha_6562, Flux_Hb_4861
    from MyDB.temp as t 
    left join SpecObjAll as s on s.bestObjID=t.objid 
    left join galSpecExtra as x on s.specObjID=x.specObjID
    left join galSpecIndx as i on s.specObjID=i.specObjID
    left join galSpecLine as l on s.specObjID=l.specObjID
    left join emissionLinesPort as e on s.specObjID=e.specObjID
    order by id'
  }else{
    query='select
    id, s.bestobjID as objid, s.specObjID as specobjid, s.z as zspec, s.zErr,
    zWarning as zW, sciencePrimary as prim, s.plate, s.mjd, s.fiberid, s.eClass
    from MyDB.temp as t 
    left join SpecObjAll as s on s.bestObjID=t.objid 
    order by id'
  }
  
  if(verbose) message('Executing query')
  q=CasJobs.executeQuery(query, context=DR, format="csv")
  if(verbose) message('Done')
  q=fread(q)
  q[q=='null' | q==-9999]=NA
  q[,4:ncol(q)]=lapply(q[,4:ncol(q)], as.numeric)
  if(dr>=12){
    q=setorder(q, id, zErr, -bpt, -hb, na.last = T)
  }else{
    q=setorder(q, id, zErr, na.last = T)
  }
  wq=which(duplicated(q$id))
  if(length(wq)>0) q=q[-wq]
  q$objid=as.integer64(q$objid)
  q$specobjid=as.integer64(q$specobjid)
  #dt=merge(data.table(data), q[,-2], by='id', all.x=T, sort=F)
  tryCatch(SkyQuery.dropTable(tablename, datasetName="MyDB"), 
           error=function(e) message('Warning: could not drop table from MyDB'))
  return(q)
}

CASPECloop=function(objid, n=1000, DR='DR12', tablename='temp'){
  data=data.table(objid)
  nr=nrow(data)
  l=ceiling(nr/n)
  sp=split(data, rep(1:l, each=n, length.out=nr))
  nrs=unlist(lapply(sp, nrow))
  nrs=cumsum(nrs)-nrs
  dl=vector('list', length = l)
  pb=txtProgressBar(0, l, 0, '=', style=3)
  for(i in 1:l){
    d=sp[[i]]
    cp=CASPEC(d$objid, DR=DR, tablename=tablename, verbose = F)
    cp$id=cp$id+nrs[i]
    dl[[i]]=cp
    setTxtProgressBar(pb, i)
  }
  res=data.table(do.call(rbind, dl))
  return(res)
}

# search for some spectroscopic variables from spec objid in casjobs
CASPEC2=function(specid, DR='DR12', tablename='temp', verbose=T){
  require(data.table)
  require(SciServer)
  require(bit64, quietly = T)
  tbldrp=tryCatch(SkyQuery.dropTable(tablename), error=function(e){})
  Name = 'dailer';
  Password = 0x540a438
  Authentication.login(Name, Password);
  specid[is.na(specid)]=1e18
  data=data.frame(id=1:length(specid), specid=as.integer64(specid))
  if(verbose) message('Uploading data to MyDB')
  CasJobs.uploadDataFrameToTable(data, tablename)
  dr=as.integer(strsplit(DR, 'DR')[[1]][2])
  if(dr>=12){
    query='select
    id, s.bestobjID as objid, s.z as zspec, s.zErr,
    zWarning as zW, sciencePrimary as prim, s.plate, s.mjd, s.fiberid, 
    snMedian as SN, bptclass as bpt, lgm_tot_p50 as smass, sfr_tot_p50 as sfr, 
    specsfr_tot_p50 as ssfr, oh_p50 as oh, d4000_n as d4000n, lick_hb as hb, 
    lick_hd_a as hd, lick_hg_a as hg, h_delta_eqw as Hd_EW, h_gamma_eqw as Hg_EW,
    h_beta_eqw as Hb_EW, h_alpha_eqw as Ha_EW, EW_Ha_6562, Flux_OIII_5006, 
    Flux_SII_6716, Flux_SII_6730, Flux_Ha_6562, Flux_Hb_4861
    from MyDB.temp as t 
    left join SpecObjAll as s on s.specObjID=t.specid 
    left join galSpecExtra as x on t.specid=x.specObjID
    left join galSpecIndx as i on t.specid=i.specObjID
    left join galSpecLine as l on t.specid=l.specObjID
    left join emissionLinesPort as e on t.specid=e.specObjID
    order by id'
  }else{
    query='select
    id, s.bestobjID as objid, s.z as zspec, s.zErr,
    zWarning as zW, sciencePrimary as prim, s.plate, s.mjd, s.fiberid, s.eClass
    from MyDB.temp as t 
    left join SpecObjAll as s on s.specObjID=t.specid
    order by id'
  }
  
  if(verbose) message('Executing query')
  q=CasJobs.executeQuery(query, context=DR, format="csv")
  if(verbose) message('Done')
  q=fread(q)
  q[q=='null' | q==-9999]=NA
  q[,3:ncol(q)]=lapply(q[,3:ncol(q)], as.numeric)
  if(dr>=12){
    q=setorder(q, id, zErr, -bpt, -hb, na.last = T)
  }else{
    q=setorder(q, id, zErr, na.last = T)
  }
  wq=which(duplicated(q$id))
  if(length(wq)>0) q=q[-wq]
  q$objid=as.integer64(q$objid)
  tryCatch(SkyQuery.dropTable(tablename, datasetName="MyDB"), 
           error=function(e) message('Warning: could not drop table from MyDB'))
  return(q)
}

# make median values w/ errors per bin
binloc=function(x, y, nbin=10, fun='median', method='length', plot=T, ...){
  stopifnot(fun %in% c('median','mean','biw'))
  stopifnot(method %in% c('length','content','clusters'))
  require(plotrix)
  require(OneR)
  require(data.table)
  isfact=length(unique(x))<length(x)/100
  if(isfact){bin=x}else{bin=bin(x, nbin, labels=1:nbin, method = method)}
  xy=data.table(x=x, y=y)
  qf=function(x) return(0.7415*(quantile(x,.75, na.rm=T)-quantile(x,.25, na.rm=T)))
  floc=switch(fun, 'median'=function(x) median(x,T), 'mean'=function(x) mean(x,0,T), 'biw'=biwLoc)
  fsc=switch(fun, 'median'=qf, 'mean'=std.error, 'biw'=biwScale)
  xm=xy[, list(x=floc(x), y=floc(y), yerr=fsc(y)), bin]
  xm=setorder(xm, bin); xm=xm[yerr!=0]
  if(plot){
    plot(x, y, pch=20, cex=.7, col='gray', ...)
    points(xm$x, xm$y, pch=20, col=2, type='b')
    arrows(xm$x, xm$y-xm$yerr, xm$x, xm$y+xm$yerr, length=.03, code=3, col=2, angle=90)
  }
  return(xm)
}

# distance between two positions at same redshift, in kiloparsecs
dist_kpc=function(ra1, dec1, z1, ra2, dec2, OM=0.3, OL=1-OM, H0=70){
  require(cosmoFns)
  decm=0.5*(dec1+dec2)
  diffra=(ra2-ra1)*cos(decm*pi/180)
  #diffra=ra2-ra1
  diffdec=dec2-dec1
  dist_deg=sqrt(diffra^2+diffdec^2)
  daa=D.A(z1, OM, OL, H0)
  d=dist_deg*daa*pi/180
  return(d*1000)
}

# distance between two positions at different redshift, in Megaparsecs
dist_Mpc=function(ra1, dec1, z1, ra2, dec2, z2, Om=0.3, Ol=0.7, H0=70){
  require(cosmoFns)
  decm=0.5*(dec1+dec2)
  diffra=(ra2-ra1)*cos(decm*pi/180)
  diffdec=dec2-dec1
  dist_deg=sqrt(diffra^2+diffdec^2)
  da=D.A(z1, Om, Ol, H0)
  r=dist_deg*da*pi/180
  d1=D.M(z1, Om, Ol, H0)/(1+z1)
  d2=D.M(z2, Om, Ol, H0)/(1+z2)
  d=d2-d1*cos(asin(dist_deg/d1))
  R=sqrt(d^2+r^2)
  return(R)
}

# retrieve galactic extinction (from Schlegel)
extinction=function(ra, dec, use_ned=T, verbose=T, interpolate=F){
  n=length(ra)
  if(use_ned){
    require(RCurl)
    require(stringr)
    if(verbose) message('Retrieveng dust extinction(s)')
    form1='http://nedwww.ipac.caltech.edu/cgi-bin/nph-calc?'
    form2='in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&lon='
    form3='d&pa=0.0&out_csys=Equatorial&out_equinox=J2000.0'
    exl=vector('list', length = n)
    for(i in 1:n){
      if(verbose) print(i)
      q=paste0(form1, form2, ra[i], 'd&lat=', dec[i], form3)
      rl=readLines(q, 90)
      rl=rl[grep('SDSS    ',rl)]
      rl=gsub("\\s+"," ", rl)
      df=str_split_fixed(cbind(rl)[,1], " ", 8)
      #exl[[i]]=data.frame(band=df[,2], Schlafly=df[,4], Schlegel=df[,8])
      exl[[i]]=as.numeric(df[,8])
    }
    dc=data.frame(do.call(rbind, exl))
    colnames(dc)=paste0('ext_',df[,2])
    dc=data.frame(Ebv=round(dc$ext_r/2.751,3),dc)
    return(dc)
  }else{
    require(astrolibR)
    require(FITSio)
    galcoords=glactc(ra, dec, 2000, j=1, degree = T)
    if(!exists('ngp_dustmap') | !exists('sgp_dustmap')){
      if(verbose) message('Reading Schlegel dust maps')
      # maplk='https://svn.sdss.org/public/data/sdss/catalogs/dust/trunk/maps/'
      # ngp_dustmap=readFITS(paste0(maplk,'SFD_dust_4096_ngp.fits'))$imDat
      # sgp_dustmap=readFITS(paste0(maplk,'SFD_dust_4096_sgp.fits'))$imDat
      ngp_dustmap=readFITS('~/MEGAsync/Schlegel Dust Maps/SFD_dust_4096_ngp.fits')$imDat
      sgp_dustmap=readFITS('~/MEGAsync/Schlegel Dust Maps/SFD_dust_4096_sgp.fits')$imDat
    }
    assign(c('ngp_dustmap'), ngp_dustmap, envir=globalenv())
    assign(c('sgp_dustmap'), sgp_dustmap, envir=globalenv())
    k=ifelse(galcoords$gb > 0, 1, -1)
    x=2048.5+sqrt(1-k*sin(galcoords$gb*pi/180))*cos(galcoords$gl*pi/180)*2048
    y=2048.5-k*sqrt(1-k*sin(galcoords$gb*pi/180))*sin(galcoords$gl*pi/180)*2048
    if(interpolate){
      require(pracma)
      n_ebv=interp2(1:4096, 1:4096, ngp_dustmap, y[k==1], x[k==1], 'linear')
      s_ebv=interp2(1:4096, 1:4096, sgp_dustmap, y[k==-1], x[k==-1], 'linear')
    }else{
      n_ebv=ngp_dustmap[cbind(round(x[k==1]),round(y[k==1]))]
      s_ebv=sgp_dustmap[cbind(round(x[k==-1]),round(y[k==-1]))]
    }
    ebv=vector(length = n)
    ebv[which(k==1)]=n_ebv; ebv[which(k==-1)]=s_ebv
    dc=data.frame(Ebv=ebv, ext_u=5.155*ebv, ext_g=3.793*ebv, ext_r=2.751*ebv, 
                  ext_i=2.086*ebv, ext_z=1.479*ebv)
    return(dc)
  }
}

kpc2deg=function(d_kpc, z){
  m=77106.47*z^3-95094.4*z^2+77664.79*z+0.1006072
  return(d_kpc/m)
}

# getting membership for several PPS regions
pps.gal=function(rproj, absu, type=1, plot=T, ...){
  # type: 1 (Mahajan); 2 (Rhee); 3 (Pasquali); 4 (Oman & Hudson); 5 (< R200)
  if(type==1){
    res=numeric(length(rproj))
    res[absu >= 1 & absu <= 2 & rproj >= 1.5 & rproj <= 2]=3
    res[absu < 1 & rproj >= 1 & rproj <= 1.5]=2
    res[absu <= 1 & rproj <= 0.5]=1
    col=c(1,'#DE1F1A','#28BF11','#5649E3')
  }
  if(type==2){
    require(sp)
    x1=c(1.966,1.977,1.926,1.802,1.722,1.651,1.572,1.503,1.455,1.432,1.425,1.34,1.34,10,10,1.966)
    y1=c(0,0.148,0.293,0.446,0.6,0.694,0.911,1.064,1.202,1.344,1.503,1.758,3.5,3.5,-0.1,-0.1)
    x2=c(0,0,0.092,0.156,0.202,0.248,0.333,0.448,0.6,0.747,0.802,0.851,0.871,0.857,0.857,0)
    y2=c(3.5,2,1.8,1.67,1.526,1.38,1.214,1.048,0.952,1.06,1.214,1.457,1.653,2.161,10,10)
    x3=c(0.454,0.431,0.4,0.344,0.37,0.451,0.756,1.053,1.2,1.303,1.354,1.38,1.4,1.372,0.454)
    y3=c(0,0.16,0.348,0.533,0.6,0.673,0.633,0.651,0.574,0.444,0.362,0.3,0.143,-0.1,-0.1)
    x4=c(0,0.154,0.182,0.251,0.269,0.287,0.4,0.453,0,0)
    y4=c(1.837,1.536,1.365,1.056,0.921,0.755,0.332,-0.1,-0.1,1.837)
    p1=point.in.polygon(rproj, absu, x1, y1)*4
    p2=point.in.polygon(rproj, absu, x2, y2)*3
    p3=point.in.polygon(rproj, absu, x3, y3)*2
    p4=point.in.polygon(rproj, absu, x4, y4)
    res=p1+p2+p3+p4
    col=c(1,'#DE1F1A','#FFBF06','#28BF11','#5649E3')
  }
  if(type==3){
    a=c(0.011, 1.916, 3.061, 3.578, 3.599, 3.256, 2.681, 2.006)
    b=c(-4.81, -5.752, -6.326, -6.532, -6.37, -5.84, -4.942, -3.676)
    c=c(1.455, 2.38, 3.089, 3.582, 3.859, 3.92, 3.765, 3.394)
    res=(absu <= a[1]*rproj^2+b[1]*rproj+c[1])
    for(i in 2:8){
      w=(absu <= a[i]*rproj^2+b[i]*rproj+c[i] &  rproj <= 1 &
           absu > a[i-1]*rproj^2+b[i-1]*rproj+c[i-1])*i
      res=res+w
    }
    col=c(1,'#D12500','#DB6A03','#B2FB2E','#13AD00','#9CCADE','#065B95','#1D1F85','#600F52')
  }
  if(type==4){
    res=(absu < -4*sqrt(3)/3/1.1363*rproj+2*sqrt(3))*1
    col=c(1,'#DE1F1A')
  }
  if(type==5){
    res=(rproj<=1)*1
    col=c(1,'#DE1F1A')
  }
  if(plot){
    plot(rproj, absu, pch=20, col=col[res+1], xlab=expression(R/R[200]), 
         ylab=expression(group("|",v[los]/sigma,"|")), ...)
    legend('topright', bty='n', pch=20, col=col, legend=sort(unique(res)), 
           cex=.8, y.intersp = .7)
  }
  return(res)
}

linesRhee=function(...){
  #x1=c(1.966,1.977,1.926,1.802,1.722,1.651,1.572,1.503,1.455,1.432,1.425,1.34,1.34,1.34,1.34,1.34,1.34)
  #y1=c(0,0.148,0.293,0.446,0.6,0.694,0.911,1.064,1.202,1.344,1.503,1.758,2,2.5,2.7,2.9,3)
  x2=c(0,0.092,0.156,0.202,0.248,0.333,0.448,0.6,0.747,0.802,0.851,0.871,0.871,0.871)
  y2=c(2,1.8,1.67,1.526,1.38,1.214,1.048,0.952,1.06,1.214,1.457,1.653,2.161,3)
  x3=c(0.14,0.239,0.3,0.37,0.451,0.756,1.053,1.2,1.303,1.354,1.379,1.4,1.372)
  y3=c(0,0.293,0.411,0.6,0.673,0.633,0.651,0.574,0.444,0.362,0.3,0.143,0)
  #x4=c(0,0.154,0.182,0.251,0.269,0.287,0.4,0.453)
  #y4=c(1.837,1.536,1.365,1.056,0.921,0.755,0.332,0)
  x1=c(1.34,1.34, 1.425, 1.432, 1.455, 1.503, 1.572, 1.651, 1.722, 1.802, 1.926, 1.966, 1.977,2.02)
  y1=c(3,2.563, 1.622, 1.614, 1.397, 1.199, 0.951, 0.758, 0.605, 0.442, 0.194, 0.113, 0.091,0)
  x4=c(0, 0.154, 0.182, 0.251, 0.269, 0.287, 0.4, 0.453,0.46)
  y4=c(1.935, 1.373, 1.271, 0.998, 0.912, 0.833, 0.282, 0.024,0)
  #s1=supsmu(x1, y1)
  #s4=supsmu(x4, y4)
  lines(x1, y1, ...)
  lines(x2, y2, ...)
  lines(x3, y3, ...)
  lines(x4, y4, ...)
}

linesPasquali=function(...){
  a=c(0.011, 1.916, 3.061, 3.578, 3.599, 3.256, 2.681, 2.006)
  b=c(-4.81, -5.752, -6.326, -6.532, -6.37, -5.84, -4.942, -3.676)
  c=c(1.455, 2.38, 3.089, 3.582, 3.859, 3.92, 3.765, 3.394)
  f=function(x, a, b, c){a*x^2+b*x+c}
  for(i in 1:8) curve(f(x, a[i], b[i], c[i]), 0, 1, add=T,...)
}

drawppslines=function(){
  abline(1.5,-1.25, col=2)
  lines(c(0,.5,.5), c(1,1,0), col=4)
  lines(c(1,1,1.5,1.5), c(0,1,1,0), col=4)
  lines(c(1.5,1.5,2,2,1.5), c(1,2,2,1,1), col=4)
  linesRhee(col=3)
  linesPasquali(col='purple')
  mtext('Jaffé', 3, -2, adj=.95, col=2)
  mtext('Mahajan', 3, -2.8, adj=.95, col=4)
  mtext('Rhee', 3, -3.6, adj=.95, col=3)
  mtext('Pasquali', 3, -4.4, adj=.95, col='purple')
}

pchoilines=function(...) lines(c(1,1,2.6,3.5), c(0.5,0.3,-0.15,-0.15), ...)

# Dressler-Shectman test for detecting substructures in galaxy clusters
ds.test=function(ra, dec, z, Nboot=1000, plot=T){
  require(RANN)
  zcl=biwLoc(z)
  vi=299792*(z-zcl)/(1+zcl)
  v=biwLoc(vi)
  sig=biwScale(vi)
  N=length(vi)
  Nnn=floor(sqrt(N))
  k=(Nnn+1)/sig^2
  nn=nn2(cbind(ra,dec), k=Nnn+1)$nn.idx
  di2=matrix(NA, N, Nboot+1)
  pb=txtProgressBar(0, 1000, style=3)
  for(j in 1:(Nboot+1)){
    vis=sample(vi, N, T)
    if(j==1) vis=vi
    for(i in 1:N){
      vinn=vis[nn[i,]]
      vloc=biwLoc(vinn)
      sigloc=biwScale(vinn)
      di2[i,j]=k*((vloc-v)^2+(sigloc-sig)^2)
    }
    setTxtProgressBar(pb, j)
  }
  close(pb)
  di=sqrt(di2)
  Ddev=apply(di, 2, sum)/N
  prob=sum(Ddev[-1]<Ddev[1])/Nboot
  if(plot){
    op=par(mfrow=c(2,2), mai=c(0.65,0.65,0.2,0.2))
    zz=(di[,1]-min(di[,1]))/max(di[,1]-min(di[,1]))/20
    #zz=exp(di[,1]/2)/200
    plot(ra, dec, pch=20, cex=.5, asp=1, xlab='RA', ylab='Dec')
    circle.fun=function(x) draw.circle(x[1],x[2],x[3])
    ap=apply(data.frame(ra,dec,zz), 1, circle.fun)
    hist(Ddev[-1], 'fd', border='white', main='DS bootstrap distr.', col='gray',
         xlab='DS statistic', xlim=c(min(Ddev),max(Ddev[1],max(Ddev))))
    abline(v=median(Ddev), col='white', lty=2, lwd=2)
    abline(v=Ddev[1], lty=2, lwd=2)
    par(op)
  }
  return(c(DS=Ddev[1],Prob=prob))
}

# Dressler-Shectman test for detecting substructures in galaxy clusters
ds.test=function(ra, dec, z, zclus=NULL, N=10, plot=T){
  require(RANN)
  n=length(ra)
  #if(n<20) N=round(sqrt(n))
  N=sqrt(n)
  if(is.null(zclus)) zclus=biwLoc(z)
  vlos=3e5*(z-zclus)/(1+zclus)
  scl=biwScale(vlos)
  vcl=biwLoc(vlos)
  k=N+1; c=k/scl^2
  d=vector(length = n)
  for(i in 1:n){
    s=abs(vlos-vlos[i])<=1000
    k=ifelse(sum(s)<k, sum(s), k)
    nn=nn2(cbind(ra[s],dec[s]), cbind(ra[i],dec[i]), k=k)
    idx=nn$nn.idx[-1]
    vls=vlos[s][idx]
    vgal=biwLoc(vls)
    sgal=biwScale(vls)
    if(sgal>scl) next;
    d[i]=(vgal-vcl)^2+(sgal-scl)^2
  }
  di=sqrt(c*d); dst=sum(di); sub=dst/n>1
  #print(dst)
  if(plot){
    require(astro)
    op=par(mar=c(5.1,4.1,4.1,5.1), oma=c(0,0,0,0))
    zlim=round(range(di, na.rm=T), 1)
    aplot(ra, dec, di, pch=20, xlab = 'RA', ylab = 'DEC', cb=T, zlim=zlim, las=1,
          cbspan = 1, zlab = expression(delta[i]), zline = 3.9)
    title('Dressler-Shectman test')
    mtext(bquote(N[nn]==.(N)), adj=0); mtext(paste('substructures:',sub))
    par(op)
  }
  return(list(substructure=sub, di=di, dn=dst/n))
}

# hexagonal binning plot
hexplot=function(x, y, z, nbin=NULL, fun='median', mincts=5, cbspan=0.95, 
                 cbsep=0.01, cbwidth=0.03, zcex=1, pal=NULL, zrescale=F, 
                 zbreaks=5, logzlab=F, ...){
  # 'x', 'y', 'z' must be vectors
  # 'fun' should be one of 'mean', 'median', 'rmedian', 'mode', 'median.disp'
  require(hexbin)
  #zbreaks=5
  zround=ifelse(fun=='mode', 0, 2)
  d=data.frame(x=x, y=y, z=z)
  
  getmode=function(v){
    v=v[!is.na(v)]; uniqv=unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  qf=function(x) 0.7415*diff(quantile(x, c(.25,.75), na.rm=T))
  
  if(is.null(nbin)){
    nbin=max(c(length(hist(x, plot=F)$breaks), length(hist(y, plot=F)$breaks)))
    nz=F
    while(!nz){
      nbin=nbin-1
      hb=hexbin(x, y, nbin, IDs=T)
      sp=split(d, hb@cID)
      n=unlist(lapply(sp, function(w) length(w$z[!is.na(w$z)])))
      nz=length(which(n==0))==0
    }
    message('number of bins = ',nbin)
  }else{
    hb=hexbin(x, y, nbin, IDs=T)
    sp=split(d, hb@cID)
  }
  
  if(fun=='mean')
    nt=unlist(lapply(sp, function(w) mean(w$z, na.rm=T)))
  if(fun=='median')
    nt=unlist(lapply(sp, function(w) median(w$z, na.rm=T)))
  if(fun=='rmedian') # rounded median
    nt=round(unlist(lapply(sp, function(w) median(w$z, na.rm=T))))
  if(fun=='mode')
    nt=unlist(lapply(sp, function(w) getmode(w$z)))
  if(fun=='median.disp')
    nt=unlist(lapply(sp, function(w) qf(w$z)))
  
  hc=hcell2xy(hb)
  dx=getmode((hc$x[-1]-hc$x[-length(hc$x)])/2)
  dy=(hc$y[-1]-hc$y[-length(hc$y)])/3
  dy=getmode(dy[dy!=0])
  hexC=hexcoords(dx, dy, sep=NA)
  #cnt=cut(nt, hist(nt, plot=F)$breaks, labels=F)
  hex.polygon=function(x, y, hexC, dx, dy=NULL, col=rainbow(10), border=1){
    n=length(x); n7=rep.int(7, n)
    polygon(x=rep.int(hexC$x, n)+rep.int(x, n7), y=rep.int(hexC$y, n)+rep.int(y, n7), 
            col=col, border=border)
  }
  cts=hb@count
  c=cts>=mincts
  nt=nt[c]
  
  if(is.null(pal)){
    pal=c('#5E4FA2','#3288BD','#66C2A5','#ABDDA4','#E6F598','#FFFFBF', 
          '#FEE08B','#FDAE61','#F46D43','#D53E4F','#9E0142')
  }
  
  if(zrescale){
    snt=sort(unique(nt))
    lnt=length(snt)
    clpal=colorRampPalette(pal)(lnt)
    mnt=match(nt, snt)
    pal.=clpal[mnt]
    zseq=seq(min(snt), max(snt), length.out=length(snt))
    nx=approx(unique(nt), unique(mnt), xout = zseq)$y
    zbreaks=ifelse(zbreaks>lnt, lnt, zbreaks)
    if(is.null(zbreaks))
      zbreaks=length(hist(nx, plot = F)$mids)
    clpal=clpal[nx]
  }else{
    ntn=(nt-min(nt))/max(nt-min(nt))
    pal.=colorRamp(pal)(ntn)
    pal.=apply(pal., 1, function(x) rgb(x[1],x[2],x[3], maxColorValue=255))
    clpal=colorRampPalette(pal)(length(unique(pal.)))
  }

  xrange=c(min(hc$x[c])-dx/2,max(hc$x[c])+dx/2)
  dy2=sqrt(dx^2+dy^2)
  #dy2=dy2
  yrange=c(min(hc$y[c])-dy2-dy,max(hc$y[c])+dy2+dy)
  magplot(xrange, yrange, pch='', side=1:4, labels=c(1,1,0,0), grid=F, ...)
  hex.polygon(hc$x[c], hc$y[c], hexC, col=pal., border=NA)
  
  cusr=par('usr')
  xi=(1-cbspan)/2; xf=1-(1-cbspan)/2
  xl=xi*(cusr[2]-cusr[1])+cusr[1]
  xr=xf*(cusr[2]-cusr[1])+cusr[1]
  yb=(1+cbsep)*(cusr[4]-cusr[3])+cusr[3]
  yt=(1+cbsep+cbwidth)*(cusr[4]-cusr[3])+cusr[3]
  clab=round(seq(min(nt), max(nt), length.out=zbreaks), zround)
  if(logzlab)
    clab=round(10^clab) # in case of log(z)
  color.legend(xl, yb, xr, yt, clab, clpal, zcex)
  gr=na.omit(data.frame(x=hc$x[c],y=hc$y[c],z=nt,zcol=pal.,counts=cts[c]))
  
  return(invisible(list(class=hb@cID, grid=gr)))
}

# squared binning plot
mplot=function(x, y, z=NULL, nbin=100, fun='median', mincts=1, cbspan=0.95, 
               cbsep=0.01, cbwidth=0.03, zcex=1, pal=NULL, plot=T, magmap=F, 
               zrescale=F, zbreaks=5, ...){
  # 'x', 'y', 'z' must be vectors; if 'z' is null an histogram is plotted
  # 'fun' should be one of 'mean', 'median', 'rmedian', 'mode', 'median.disp'
  require(raster)
  #zbreaks=5
  zround=ifelse(fun=='mode', 0, 2)
  d=cbind(x, y)
  rast=raster()
  extent(rast)=extent(d)
  dim(rast)=c(nbin, nbin, 1)
  
  getmode=function(v){
    v=v[!is.na(v)]; uniqv=unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  qf=function(x) 0.7415*diff(quantile(x, c(.25,.75), na.rm=T))
  
  if(is.null(z)){
    mr=rasterize(d, rast, fun='count')
    m=as.matrix(mr)
    m[m < mincts]=NA
  }else{
    f=switch(fun, 'mean'=mean, 'median'=median, 
             'median.disp'=qf, 'mode'=getmode)
    mr=rasterize(d, rast, z, fun=function(x,...) f(x))
    mc=rasterize(d, rast, z, fun='count')
    m=as.matrix(mr)
    c=as.matrix(mc)
    m[c < mincts]=NA
  }
  xcut=seq(mr@extent@xmin, mr@extent@xmax, length.out = nbin+1)
  ycut=seq(mr@extent@ymin, mr@extent@ymax, length.out = nbin+1)
  midpoints=function(x) (x[-1]+x[-length(x)])/2
  xs=midpoints(xcut)
  ys=midpoints(ycut)
  m=t(m[nbin:1,])
  l=list(x=xs, y=ys, z=m)
  
  if(plot){
    require(magicaxis)
    if(is.null(pal)){
      pal=c('#5E4FA2','#3288BD','#66C2A5','#ABDDA4','#E6F598','#FFFFBF', 
            '#FEE08B','#FDAE61','#F46D43','#D53E4F','#9E0142')
    }
    cm=na.omit(c(m))
    um=unique(cm)
    lum=length(um)
    clpal=colorRampPalette(pal)(lum)
    zbreaks=ifelse(zbreaks>lum, lum, zbreaks)
    
    if(zrescale){
      sm=sort(um)
      m=matrix(match(c(m), sm), nbin, nbin)
      mcm=match(cm, sm)
      zseq=seq(min(sm), max(sm), length.out=length(sm))
      nx=approx(um, unique(mcm), xout = zseq)$y
      if(is.null(zbreaks))
        zbreaks=length(hist(nx, plot = F)$mids)
      clpal=clpal[nx]
      l$z=scale(m)
      l$z=m
    }
    
    dx=diff(range(xs))*0.07
    dy=diff(range(ys))*0.07
    xlim=c(min(xs)-dx, max(xs)+dx)
    ylim=c(min(ys)-dy, max(ys)+dy)
    magimage(xs, ys, m, col=pal, magmap=magmap, side=1:4, labels=c(1,1,0,0),
             asp=NA, useRaster=T, xlim=xlim, ylim=ylim, ...)
    
    cusr=par('usr')
    xi=(1-cbspan)/2; xf=1-(1-cbspan)/2
    xl=xi*(cusr[2]-cusr[1])+cusr[1]
    xr=xf*(cusr[2]-cusr[1])+cusr[1]
    yb=(1+cbsep)*(cusr[4]-cusr[3])+cusr[3]
    yt=(1+cbsep+cbwidth)*(cusr[4]-cusr[3])+cusr[3]
    clab=round(seq(min(cm), max(cm), length.out=zbreaks), zround)
    color.legend(xl, yb, xr, yt, clab, clpal, zcex)
  }
  invisible(l)
}

# get a,b coefficients from points coordinates
getab=function(x1,x2,y1,y2){
  b=(y2-y1)/(x2-x1); a=y1-b*x1
  ab=c(a,b); names(ab)=c('a','b')
  return(ab) #y=b*x+a
}

# get a,b coefficients from loc=locator() and draws the line
getab2=function(loc, draw=T){
  x=loc[[1]]; y=loc[[2]]
  b=(y[2]-y[1])/(x[2]-x[1]); a=y[1]-b*x[1]
  ab=c(a,b); names(ab)=c('a','b')
  if(draw) abline(a, b)
  return(ab)
}

# using sklearn.mixture.GaussianMixture from python
GMM=function(data, G=2:9, plot=T){
  require(reticulate)
  stopifnot(py_module_available('sklearn'))
  call=match.call()
  tol=0.001; reg_covar=1e-06
  max_iter=100L; n_init=1L
  init_params="kmeans"
  data=as.matrix(data)
  if(length(G)==1 & G[1]<2) 
    stop('No classification with the specified number of clusters')
  sk=import('sklearn')
  mix=sk$mixture
  GM=mix$GaussianMixture
  cv_types=c('spherical','tied','diag','full')
  iseq=as.integer(G)
  lG=length(G)
  bic=aic=matrix(NA, lG, 4, dimnames=list(G,cv_types))
  lowest_bic=Inf
  for(j in 1:4){
    cv_type=cv_types[j]
    for(i in iseq){
      gmm=GM(i, cv_type, tol=tol, reg_covar=reg_covar, max_iter=max_iter,
             n_init=n_init, init_params=init_params)
      gmm$fit(data)
      nbic=gmm$bic(data)
      bic[as.character(i),j]=nbic
      aic[as.character(i),j]=gmm$aic(data)
      if(nbic < lowest_bic){
        lowest_bic=nbic 
        bstm = gmm
      }
    }
  }
  n=nrow(data); d=ncol(data); ncomp=bstm$n_components
  cls=bstm$predict(data)+1; covs=bstm$covariances_ 
  aic.=bstm$aic(data); prob=bstm$predict_proba(data)
  dimnames(prob)=list(1:n); sc=bstm$score(data)
  scs=bstm$score_samples; wgh=bstm$weights_
  cvg=bstm$converged_; prec=bstm$precisions_
  cholprec=bstm$precisions_cholesky_
  pars=bstm$get_params()
  pars=c(list(mean=bstm$means_, covariances=covs), pars)
  #lwb=bstm$lower_bound_; smpl=bstm$sample()
  message(paste0('GMM model components: '),'(',pars$covariance_type,', ',ncomp,')')
  gmres=list(call=call, data=data, n=n, d=d, G=ncomp, cv_types=cv_types, BIC=bic, 
             bic=lowest_bic, AIC=aic, aic=aic., parameters=pars, probs=prob, 
             classification=cls, loglik=sc, logprobs=scs, weights=wgh, 
             converged=cvg, precisions=prec, precisions_cholesky=cholprec)
  if(plot) plot.GMM(gmres, what="all")
  return(invisible(gmres))
}

# plotting GMM object
plot.GMM=function(x, what=c("BIC", "classification", "uncertainty", "all"), add.ellipses=T){
  require(magicaxis)
  require(plotrix)
  require(RColorBrewer)
  stopifnot(what %in% c('BIC','classification','uncertainty','cl','c','bic',
                        'b','un','u','all'))
  bic=x$BIC; G=1:nrow(bic); cv_types=x$cv_types; ncomp=x$G
  data=x$data; cls=x$classification; means=x$parameters$mean
  covs=x$parameters$covariances; probs=x$probs; d=x$d
  if(d>1){
    ei=tryCatch(apply(covs, 1, eigen), error=function(e) NA)
    if(is.na(ei[1])){
      add.ellipses=F
      message("Warning: non-square matrix in 'eigen', add.ellipses set to FALSE")
    }else{
      ab=do.call(rbind,lapply(ei, function(x) 1.1*sqrt(2)*sqrt(x[[1]]))) #95%<-1.7; 68%<-1.1
      angles=unlist(lapply(ei, function(x) atan2(x[[2]][1,1], x[[2]][1,2])))
    }
  }
  if(what[1]=="all" & d<3) op=par(mfrow=c(2,2), mai=c(0.6,0.6,0.15,0.15))
  if(what[1] %in% c('BIC','bic','b') | what[1]=="all"){
    pch=c(17,2,16,0); col=c("#597DC3","#CC0000","#218B21","#508476")
    magplot(range(G), range(bic), pch='', ylab='BIC', grid=T, side=2, 
            grid.col="gray90")
    magaxis(1, max(G), 0, xlab='Number of components')
    for (i in 1:4) points(G, bic[,i], type='b', pch=pch[i], col=col[i])
    legend('topright', bty='n', legend=cv_types, pch=pch, col=col, cex=.9, 
           y.intersp = ifelse(what[1]=="all" & d<3, .5, .7))
  }
  col=c("#004586","#FF410F","#579D1C","#FFD321","#7E0C20","#83CAFF",
        "#AED003","#4B1F6F","#FF9511","#C51A0C","#0C83D1","#303F00")
  #col=sample(col)
  ecol=rgb(t(rbind(col2rgb(col)/255)), alpha=.3)
  pch=c(16, 0, 17, 3, 15, 4, 1, 8, 2, 7, 5, 9, 6, 10, 11, 18, 12, 13, 14)
  cps=1:ncomp; cols=col[cps]; ecols=ecol[cps]; pchs=pch[cps]
  diag.panel=function(w,...){
    usr=par("usr"); on.exit(par(usr))
    par(usr=c(usr[1:2], 0, 1.5))
    sp=split(w, cls)
    den=lapply(sp, density)
    for(i in 1:length(sp)){
      x=den[[i]]$x; y=den[[i]]$y
      polygon(x, y/max(y), col=ecols[i], border = NA)
    }
  }
  if(what[1] %in% c('classification','cl','c') | what[1]=="all"){
    if(d==1){
      name=as.character(x$call)[2]
      magplot(data, cls, majorn=ncomp, minorn = 0, side=2, ylim=c(.8,ncomp+.2), 
              ylab='Classification', pch=''); magaxis(1, xlab=name)
      arrows(data,cls-.15,data,cls+.15,length = 0, col=cols[cls])
    }
    if(d==2){
      dmnames=dimnames(data)[[2]]
      if(is.null(dmnames)) dmnames=cps
      magplot(data, col=cols[cls], pch=pchs[cls], xlab=dmnames[1], ylab=dmnames[2])
      if(add.ellipses){
        apply(means, 1, function(x) points(x[1], x[2], pch=10))
        draw.ellipse(means, a=ab[,1], b=ab[,2], angle=pi/2+angles, deg=F)
      }
    }
    if(d>2){
      pairs(data, col=cols[cls],pch=pchs[cls],diag.panel=diag.panel,upper.panel=NULL)
    }
  }
  if(what[1] %in% c('uncertainty','un','u') | what[1]=="all"){
    if(d==1){
      name=as.character(x$call)[2]
      err=1-apply(probs, 1, function(x) sort(x, decreasing=T)[1])
      magplot(data, err, pch='', xlab=name, ylab='Uncertainty')
      arrows(data, 0, data, err, length = 0, col=cols[cls])
      dd=density(data)
      lines(dd$x, dd$y/max(dd$y)*max(err), col='gray', lty=2)
    }
    if(d==2){
      dmnames=dimnames(data)[[2]]
      prb=1.6-apply(probs, 1, function(x) sort(x, decreasing=T)[1])^3
      magplot(data, col=cols[cls], pch=20, xlab=dmnames[1], ylab=dmnames[2], cex=prb)
      if(add.ellipses){
        apply(means, 1, function(x) points(x[1], x[2], pch=10))
        draw.ellipse(means, a=ab[,1], b=ab[,2], angle=pi/2+angles, deg=F)
      }
    }
    if(d>2){
      prb=1.6-apply(probs, 1, function(x) sort(x, decreasing=T)[1])^3
      pairs(data, col=cols[cls], pch=20, diag.panel = diag.panel, cex=prb,
            upper.panel = NULL)
    }
  }
  if(d==2 & what[1]=="all"){
    rg=apply(data, 2, range)
    xx=seq(rg[1,1], rg[2,1], length.out = 100)
    yy=seq(rg[1,2], rg[2,2], length.out = 100)
    xy=as.matrix(expand.grid(xx,yy))
    z=-x$logprobs(xy); z=matrix(z, length(xx), length(xx)); z[z<0]=NA
    magplot(data, pch='.', xlab=dmnames[1], ylab=dmnames[2], col='gray')
    pal=(brewer.pal(11,'Spectral'))
    contour(xx, yy, log10(z), add = T, nlevels = 10, drawlabels=F, col=pal, lwd=2)
  }
  if(what[1]=="all" & d<3) par(op)
}

# splitting dataframe in parts (by rows) with n rows
splitdf=function(data, n){
  nr=nrow(data)
  sp=split(data, rep(1:ceiling(nr/n), each=n, length.out=nr))
  return(sp)
}

# weighted kernel density estimation in 2D from python
kde_2d_weighted=function(x, y, w=NULL, n=100, h=c('factor','scott','silverman')){
  require(reticulate)
  ghdir='https://raw.githubusercontent.com/Dailer/Astro/master/kde_2d_weighted.py'
  if(!py_module_available('kde_2d_weighted')) 
    k=import_from_path('kde_2d_weighted', ghdir)
  xy=t(cbind(x,y))
  kg=k$gaussian_kde(xy, weights=w)
  if(h[1]=='factor')
    h=kg$factor
  kg$set_bandwidth(bw_method=h[1])
  xr=seq(min(x), max(x), length.out = n)
  yr=seq(min(y), max(y), length.out = n)
  xyr=expand.grid(xr, yr)
  zz=kg(list(xyr[,1], xyr[,2]))
  zim=matrix(zz, n, n)
  return(list(x=xr, y=yr, z=zim))
}

# plot lines from linear regression with confidence/prediction intervals as shaded polygon
lines.lm=function(x, y, interval='confidence', level=.95, n=100, plot.int=T, col.int='#BEBEBE80', 
                  col=1, lwd=1, lty=1, ...){
  md=lm(y~x)
  xr=range(x, na.rm=T)
  xoff=0.1*diff(xr)
  nwx=seq(xr[1]-xoff, xr[2]+xoff, length.out=n)
  prd=predict.lm(md, newdata=data.frame(x=nwx), interval=interval, level=level)
  lwr=prd[,2]; upr=prd[,3]
  abline(md, col=col, lwd=lwd, lty=lty)
  if(plot.int)
    polygon(c(nwx,rev(nwx)), c(lwr,rev(upr)), border=NA, col=col.int)
}

test.unimod=function(data){
  require(diptest)
  require(silvermantest)
  require(multimode)
  library(ks)
  #require(modes)
  require(pastecs)
  if(length(dim(data))!=2 & class(data)=='numeric') data=cbind(data)
  stopifnot(length(dim(data))==2)
  cat('Testing unimodality using the dip, Silverman and Cheng & Hall tests\n')
  n=ncol(data)
  dt=st=cht=tp=lm=numeric(n)
  tPb=txtProgressBar(0, n, 0, '=', style=3)
  for(i in 1:n){
    x=unlist(data[,i])
    x=c(na.omit(x))
    bw=hns(x, deriv.order = 2)
    dt[i]=dip.test(x)$p.value
    st[i]=silverman.test(x, k=1, R=500)@p_value
    cht[i]=modetest(x, method='CH', B=500)$p.value
    #nm[i]=nmodes(x, bw=bw)
    #am[i]=nrow(amps(x)$Peaks)
    tp[i]=sum(turnpoints(density(x, bw=bw)$y)$peaks)
    #lm[i]=length(locmodes(x, n=1e4)$locations)
    setTxtProgressBar(tPb, i)
  }
  cat('\n\n')
  t.res=data.frame(dip.pv=round(dt,3), silv.pv=st, cheng.pv=cht)
  res=apply(t.res<0.05, 1, sum)>1
  t.res$is.unimodal=!res
  #t.res$Nm.1=nm; t.res$Nm.2=am; 
  t.res$N.modes=tp
  #t.res$N.modes=lm
  rownames(t.res)=colnames(data)
  return(t.res)
}

# create id number from two numbers (pairing function)
id2=function(x, y){
  n=0.5*(x+y)*(x+y+1)+x
  return(n)
}

# create id number from three positive numbers
id3=function(x,y,z){
  require(bit64)
  n=as.integer64(choose(x,1)+choose(x+y+1,2)+choose(x+y+z+2,3))  
  return(n)
}

# rotate bivariate data (matrix, 1st column: x, 2nd column: y) by an angle
rotdata=function(data, angle, type='rad'){
  if(type=='deg') angle=angle*pi/180
  rot=matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), 2, 2) %*% t(data)
  return(t(rot))
}

rotdatalm=function(x, y, asp=1, plot=T){
  df=data.frame(x=x,y=y)
  reg=lm(y~x, df)
  angle=-coef(reg)[2]
  rot=t(matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), 2, 2) %*% t(as.matrix(df)))
  if(plot){
    op=par(mfrow=c(2,2), mai=c(.65,.65,.1,.1))
    plot(x, y, asp=asp)
    abline(reg, col=2)
    plot(rot[,1], rot[,2], asp=asp)
    par(op)
  }
  return(rot)
}

# 
showsdssimage=function(id, filter='r'){
  require(data.table)
  require(SciServer)
  require(FITSio)
  require(R.utils)
  require(magicaxis)
  require(celestial)
  Name = 'dailer';
  Password = 0x540a438
  Authentication.login(Name, Password);
  message('Searching data in CAS')
  q='SELECT run, camcol, field, ra, dec, petroRad_r from PhotoObjAll where objid='
  query=paste0(q, id)
  q=CasJobs.executeQuery(query, context='DR12', format="csv")
  d=fread(q)
  run=d$run; camcol=d$camcol; field=d$field; ra=d$ra; dec=d$dec; rp=c(d[,6])[[1]]
  basedir='https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/'
  filename=paste0('frame-',filter,'-',sprintf("%06.f", run),'-',as.character(camcol),'-',
                  sprintf("%04.f", field),'.fits.bz2')
  link=paste0(basedir,run,'/',camcol,'/',filename)
  message('Downloading image')
  system(paste('wget',link), intern = T)
  bunzip2(filename, overwrite=TRUE)
  ft=readFITS(gsub('.bz2','',filename))
  im=ft$imDat; hdr=ft$hdr
  boxsize=round(rep(3*rp,2))
  imcut=magcutoutWCS(im, hdr, c(ra, dec), box=boxsize)
  message('Plotting')
  magimage(imcut$image)
  return(invisible(list(data=data.table(id,d), image=imcut$image)))
}

# read all files in a directry with a name pattern
readall=function(path='.', pattern=NULL, ...){
  require(data.table)
  lf=list.files(path, pattern = pattern)
  lf=paste0(path,'\\', lf)
  n=length(lf)
  ls=list('vector', length=n)
  for(i in 1:n) ls[[i]]=fread(lf[i], ...)
  d=data.table(do.call(rbind, ls))
  message(n, ' items read')
  return(d)
}

# Function for clustering points with dbscan
cluspts=function(x, y, k=5, plot=F){
  require(RANN)
  require(dbscan)
  d=nn2(cbind(x,y), k=k+1)$nn.dists[,k+1]
  tp=turnpts(d, plot=F)
  d0=mean(c(tp$peaks[1], tp$pits[1]))
  if(is.na(d0)) d0=mean(c(tp$peaks[1], max(d)))
  db=dbscan(cbind(x,y), d0)
  cl=db$cluster
  id=as.numeric(names(sort(table(cl), T))[1])
  cl=ifelse(cl==id, 1, 0)
  if(plot) plot(x, y, col=cl+1, pch=20)
  return(cl)
}

# Two-Point Correlation Function
tpcf=function(x, y, nbin=6, type='DP', method='content', knn=5, full_region=T, plot=T){
  # require turnpts and cluspts function
  # type is the TPCF estimator: 
  #  'DP': Davis and Peebles (1983) 
  #  'H': Hamilton (1993)
  #  'LS': Landy and Szalay (1993)
  require(sp)
  require(flexclust)
  require(OneR)
  
  if(full_region){
    xy=cbind(x,y)
    pp=1
    mi.x=min(x); ma.x=max(x)
    mi.y=min(y); ma.y=max(y)
    pbuf=cbind(c(mi.x,ma.x,ma.x,mi.x,mi.x),c(mi.y,mi.y,ma.y,ma.y,mi.y))
    Ain=(max(x)-min(x))*(max(y)-min(y))
  }else{
    require(sf)
    require(pracma)
    cl=cluspts(x, y, knn)
    xi=x[cl==1]; yi=y[cl==1]
    xy=cbind(xi,yi)
    ch=chull(xi, yi)
    coords=xy[c(ch,ch[1]),]
    p=st_polygon(list(coords))
    pbuf=st_buffer(p, .01)[[1]]
    pp=point.in.polygon(x, y, pbuf[,1], pbuf[,2])
    Ain=abs(polyarea(pbuf[,1], pbuf[,2]))
  }
  
  Aout=(max(x)-min(x))*(max(y)-min(y))
  Ar=Aout/Ain
  Ar=2
  
  ND=nrow(xy)
  rxy=cbind(runif(round(ND*Ar),min(x),max(x)), runif(round(ND*Ar),min(y),max(y)))
  rp=point.in.polygon(rxy[,1], rxy[,2], pbuf[,1], pbuf[,2])
  ri=rxy[rp==1,]
  NR=nrow(ri)
  
  DD=as.numeric(dist(xy))
  RR=as.numeric(dist(ri))
  
  xyr=rbind(xy, ri)
  DR_all=as.numeric(dist(xyr))
  DR=DR1=setdiff(DR_all, c(DD,RR))
  #DDRR=c(DD, RR, DR)
  
  n_annuli=nbin+1
  #mn=min(min(DD), min(RR), min(DR))
  #mx=max(max(DD), max(RR), max(DR))
  #annuli=head(c(0.1, mn*10^(log10(mx/mn)*(1:n_annuli)/(n_annuli-1))), -1)
  #annuli=10^(seq(log10(0.2), log10(mx), length.out=nbin+1))
  DRS=log10(DD)
  DRS=DRS[DRS>log10(0.2)]
  bn=bin(DRS, nbin, labels=1:nbin, method=method)
  sp=split(DRS, bn)
  bks=do.call(rbind,lapply(sp, range))
  annuli=10^c(bks[,1], bks[nbin,2])
  bins=10^unlist(lapply(sp, median))
  Index.Annuli=1:n_annuli
  #bins=(head(annuli,-1)+annuli[-1])/2
  #print(bins)
  
  DD_sep=cut(DD, annuli, labels=F)
  RR_sep=cut(RR, annuli, labels=F)
  DR_sep=cut(DR, annuli, labels=F)
  
  DD_tab=data.frame(table(DD_sep))
  RR_tab=data.frame(table(RR_sep))
  DR_tab=data.frame(table(DR_sep))
  
  adf=data.frame(bin=1:nbin)
  adf=merge(adf, DD_tab, by.x='bin', by.y='DD_sep', all.x=T)
  adf=merge(adf, RR_tab, by.x='bin', by.y='RR_sep', all.x=T)
  adf=merge(adf, DR_tab, by.x='bin', by.y='DR_sep', all.x=T)
  adf[is.na(adf)]=0
  
  DD=adf[,2]; delta_DD=sqrt(DD)/1
  RR=adf[,3]; delta_RR=sqrt(RR)/1
  DR=adf[,4]; delta_DR=sqrt(DR)/1
  
  if(type=='LS'){
    ACF=NR*(NR-1)/(ND*(ND-1))*DD/RR-(NR-1)/ND*DR/RR+1 # Landy & Szalay TPCF
    Delta_ACF=sqrt(((NR*(NR-1)/(ND*(ND-1)))^2*((delta_DD/RR)^2+(DD*delta_RR/RR^2)^2))+
                     ((2*(NR-1)/ND)^2*((delta_DR/RR)^2+(DR*delta_RR/RR^2)^2)))
  }
  if(type=='H'){
    ACF=4*ND*NR/((ND-1)*(NR-1))*DD*RR/DR^2-1 # Hamilton TPCF
    Delta_ACF=sqrt((4*ND*NR/((ND-1)*(NR-1)))^2*((delta_DD*RR/DR^2)^2+(delta_RR*DD/DR^2)^2+
                                                  (delta_DR/DR^3*RR*DD)^2))
  }
  if(type=='DP'){
    ACF=2*NR/(ND-1)*DD/DR-1 # Davis & Peebles TPCF
    Delta_ACF=sqrt((2*NR/(ND-1))^2*((delta_DD/DR)^2+(DD*delta_DR/DR^2)^2))
  }
  
  idx=which(bins > 0 & bins < 3)
  lmfit_acf_h=lm(log10(ACF+1)[idx]~log10(bins)[idx])
  summ=summary(lmfit_acf_h)
  cf=coef(summ)
  interc=cf[1,1]; expont=cf[2,1]
  delta_interc=cf[,"Std. Error"][1] 
  delta_expont=cf[,"Std. Error"][2] 
  r_sqrd=summ$r.squared
  
  if(plot){
    #op=par(mai=c(0.65,0.65,0.1,0.1), mfrow=c(2,2))
    require(magicaxis)
    magplot(rxy, pch=4, col=ifelse(rp==1, 'gray50', 'gray90'), xlim=range(x), ylim=c(-3,3),
            xlab=expression(R/R[200]), ylab=expression(v[los]/sigma))
    points(x, y, pch=19, col=pp+1)
    lines(pbuf, col='red')
    box()
    
    # tcol=rgb(t(rbind(col2rgb(2:4)/255)), alpha=.4)
    # barplot(DR, col=tcol[1], names.arg=1:nbin) #red
    # barplot(RR, col=tcol[3], add=T) # blue
    # barplot(DD, col=tcol[2], add=T) # green
    
    magplot(bins, ACF+1, pch=19, log='xy', xlab="Separation r", ylab=expression(1+xi(r)), las=1,
            ylim=c(min(ACF+1-Delta_ACF), max(ACF+1+Delta_ACF)))
    arrows(bins, ACF+1-Delta_ACF, bins, ACF+1+Delta_ACF, length=0.05, angle=90, code=3)
    abline(lmfit_acf_h, col="blue")
    
    legend('bottomleft', legend=c(sprintf("Exponent = %+3.2f \u00B1 %3.2f", expont, delta_expont), 
                                sprintf("Intercept = %+3.2f \u00B1 %3.2f", interc, delta_interc), 
                                sprintf("R\UB2 = %3.2f", r_sqrd)), bty='n', cex=.8)
    # maghist(DR1, 'fd', log='x', col=tcol[1], verbose = F)
    # maghist(RR1, 'fd', log='x', col=tcol[3], add=T, verbose = F)
    # maghist(DD1, 'fd', log='x', col=tcol[2], add=T, verbose = F)
    #par(op)
  }
  param=round(c(gamma=expont, R2=r_sqrd),5)
  return(list(parameters=param, data=round(data.frame(xi=ACF,r=bins,xi.err=Delta_ACF),5)))
}

# Measuring intergalactic pair separation distribution following Capelato (1980)
intergal=function(ra, dec, z=NULL, r200=1, vdisp=1e3, Om=0.3, Ol=1-Om, H0=70, plot=T, ...){
  # if z is not provided, ra and dec are considered as projected normalized 
  # distances and velocities, respectively
  require(cosmoFns)
  require(minpack.lm)
  if(!is.null(z)){
    raclus=median(ra)
    declus=median(dec)
    zclus=median(z)
    # Obtaining projected velocities and distances
    vlos=299792*(z-zclus)/(1+zclus)/vdisp
    groupdist=D.A(zclus, Om, Ol, H0) # in Mpc
    zsep=sin(declus*pi/180)*sin(dec*pi/180)
    xysep=cos(declus*pi/180)*cos(dec*pi/180)*cos((raclus-ra)*pi/180)
    xyzsep=zsep+xysep
    xyzsep[xyzsep>1]=1
    angle=acos(xyzsep)
    dproj=groupdist*angle/r200 # in Mpc
  }else{
    dproj=ra; vproj=dec
  }
  # Obtaining pairwise distances
  dpr=c(dist(cbind(dproj,dproj)))
  dpv=c(dist(cbind(vproj,vproj)))
  # Binning distances using Freedman-Diaconis method
  hr=hist(dpr, 'fd', plot=F)
  hv=hist(dpv, 'fd', plot=F)
  xr=hr$mids; yr=hr$counts/length(dpr); yr=yr/max(yr)
  xv=hv$mids; yv=hv$counts/length(dpv); yv=yv/max(yv)
  # Fitting to Capelato's exponential model
  ftr=nlsLM(yr~pi0*exp(-(xr/s0)^beta), start=c(pi0=max(yr), s0=1, beta=1))
  ftv=nlsLM(yv~pi0*exp(-(xv/s0)^beta), start=c(pi0=max(yv), s0=1, beta=1))
  cfr=coef(summary(ftr))
  cfv=coef(summary(ftv))
  # Plotting histograms and best fit curves
  if(plot){
    require(magicaxis)
    par(mai=c(0.7,0.7,0.1,0.1), mfrow=c(2,2)) 
    newxr=seq(min(dpr), max(dpr), length.out = 100)
    newyr=predict(ftr, data.frame(xr=newxr))
    magplot(xr, yr, type='s', xlab='Radii separation', las=1, mtline=c(2,2.5),
            ylab=expression(P(s)/N[pairs]), side=1:4, labels=c(1,1,0,0),...)
    lines(newxr, newyr, col=2)
    spr1=sprintf('%3.2f \u00B1 %3.2f', cfr[1,1], cfr[1,2])
    spr2=sprintf('%3.2f \u00B1 %3.2f', cfr[2,1], cfr[2,2])
    spr3=sprintf('%3.2f \u00B1 %3.2f', cfr[3,1], cfr[3,2])
    mtext(bquote(pi[0]==~.(spr1)), 3, -2, adj=.95, cex=.8)
    mtext(bquote(s[0]==~.(spr2)), 3, -3, adj=.95, cex=.8)
    mtext(bquote(beta==~.(spr3)), 3, -4, adj=.95, cex=.8)
    
    newxv=seq(min(dpv), max(dpv), length.out = 100)
    newyv=predict(ftv, data.frame(xv=newxv))
    magplot(xv, yv, type='s', xlab='Velocities separation', las=1, mtline=c(2,2.5),
            ylab=expression(P(s)/N[pairs]), side=1:4, labels=c(1,1,0,0),...)
    lines(newxv, newyv, col=2)
    spv1=sprintf('%3.2f \u00B1 %3.2f', cfv[1,1], cfv[1,2])
    spv2=sprintf('%3.2f \u00B1 %3.2f', cfv[2,1], cfv[2,2])
    spv3=sprintf('%3.2f \u00B1 %3.2f', cfv[3,1], cfv[3,2])
    mtext(bquote(pi[0]==~.(spv1)), 3, -2, adj=.95, cex=.8)
    mtext(bquote(s[0]==~.(spv2)), 3, -3, adj=.95, cex=.8)
    mtext(bquote(beta==~.(spv3)), 3, -4, adj=.95, cex=.8)
  }
  l=list('radius_dist_fit'=cfr, 'radius_dist'=data.frame(x=xr,y=yr), 
         'vel_dist_fit'=cfv, 'vel_dist'=data.frame(x=xv,y=yv))
  return(l) 
}

# Optimized histogram
opthist=function(x, breaks='Knuth', Nmin=2, verbose=F, plot=T, ...){
  require(Rcpp)
  require(magicaxis)
  # require(hexbin)
  x=c(na.omit(x))
  sourceCpp('~/MEGAsync/PAPERS/Histogram Optimization/hist.cpp')
  N=NULL
  if(breaks %in% c('Hogg','H','h')){
    N=hhist(x, Nmin)
  }
  if(breaks %in% c('Knuth','K','k')){
    N=khist(x, Nmin)
  }
  if(breaks %in% c('Shimazaki','S','s')){
    N=sshist(x, Nmin)
  }
  if(!is.null(N)){
    if(N==100)  message('Warning: number of maximum bins (100) reached')
    breaks=seq(min(x), max(x), length.out = N+1)
  } 
  hst=maghist(x, breaks, verbose=verbose, plot=plot, ...)
  return(invisible(hst))
}

# 1D histogram bin optimization with Shimazaki, Knuth & Hogg methods
hst=function(x, breaks='Hogg', plot=T, col='lightgray', border=1, bwd=1, bty=1, 
             freq=T, ...){
  # 'x' must be a vector 
  # 'breaks' should be one of 'Shimazaki', 'Knuth', 'Hogg'
  x=c(na.omit(x))
  a=2
  b=100
  # Shimazaki method
  sshist=function(x){
    n=b-a+1
    rg=range(x)
    r=diff(rg)
    ntot=length(x)
    C=vector(length = n)
    for(i in 1:n){
      Ni=i+a-1
      D=r/Ni
      brk=seq(rg[1], rg[2], length.out = Ni+1)
      ki=hist(x, brk, plot=F)$counts
      k=ntot/Ni
      v=sum((ki-k)^2)/Ni
      C[i]=(2*k-v)/D^2
    }
    idx=which.min(C)
    optN=idx+a-1
    return(optN)
  }
  # Knuth method
  khist=function(x){
    n=b-a+1
    l=length(x)
    rg=range(x)
    logp=vector(length = n)
    for(i in 1:n){
      Ni=i+a-1
      brk=seq(rg[1], rg[2], length.out=Ni+1)
      cts=hist(x, brk, plot=F)$counts
      A=l*log(Ni)+lgamma(Ni/2)-lgamma(l+Ni/2)
      B=-Ni*lgamma(0.5)+sum(lgamma(cts+0.5))
      logp[i]=A+B
    }
    idx=which.max(logp)
    optN=idx+a-1
    return(optN)
  }
  # Hogg method
  hhist=function(x){
    alpha=1
    n=b-a+1
    rg=range(x)
    r=diff(rg)
    L=vector(length = n)
    for(i in 1:n){
      Nb=i+a-1
      D=r/Nb
      brk=seq(rg[1], rg[2], length.out = Nb+1)
      Ni=hist(x, brk, plot=F)$counts
      S=sum(Ni+alpha)
      Li=vector(length = Nb)
      for(j in 1:Nb){
        Li[j]=Ni[j]*log((Ni[j]+alpha-1)/(D*(S-1)))
      }
      L[i]=sum(Li, na.rm = T)
    }
    idx=which.max(L)
    optN=idx+a-1
    return(optN)
  }
  N=NULL
  brk=breaks[1]
  if(brk %in% c('Hogg','H','h')){
    N=hhist(x)
  }else if(brk %in% c('Knuth','K','k')){
    N=khist(x)
  }else if(brk %in% c('Shimazaki','S','s')){
    N=sshist(x)
  }else{
    stop("'breaks' should be one of 'Hogg', 'Knuth', 'Shimazaki'")
  }
  if(!is.null(N)){
    message('optimal number of bins: ', N)
    if(N==b)  message('Warning: number of maximum bins (',b,') reached')
    if(N>b)  message('Warning: number of maximum bins (',b,') excedeed')
  }
  if(length(breaks)==1 & is.null(N)) N=breaks
  breaks=seq(min(x), max(x), length.out = N+1)
  hst=hist(x, breaks, plot=plot, col=col, border=NA, freq=freq,...)
  if(plot){
    brk=hst$breaks
    if(freq) cts=hst$counts
    else cts=hst$density
    lines(c(brk[1],brk), c(0,cts,0), type='s', col=border, lwd=bwd, lty=bty)
  }
  return(invisible(hst))
}

# 2D histogram bin optimization with Shimazaki, Knuth & Hogg methods
hst2d=function(x, breaks='k', plot=F){
  # 'x' must be a two columns matrix
  x=na.omit(x)
  a=2
  b=100
  # Shimazaki method
  sshist_2d=function(x, a=2, b=100, plot=F){
    n=b-a+1
    rg1=diff(range(x[,1]))
    rg2=diff(range(x[,2]))
    C=numeric(n)
    ntot=nrow(x)
    for(i in 1:n){
      Ni=i+a-1
      Ni2=Ni^2
      ki=c(gplots::hist2d(x[,1], x[,2], Ni, show=F)$counts)
      D=rg1*rg2/Ni2
      k=ntot/Ni2
      v=sum((ki-k)^2)/Ni2
      C[i]=(2*k-v)/D^2
    }
    idx=which.min(C)
    optN=idx+a-1
    if(plot) plot(a:b, C, type='l')
    return(optN)
  }
  # Knuth method
  khist_2d=function(x, a=2, b=100, plot=F){
    n=b-a+1
    l=nrow(x)
    logp=numeric(n)
    for(i in 1:n){
      Ni=i+a-1
      Ni2=Ni^2
      cts=c(gplots::hist2d(x[,1], x[,2], Ni, show=F)$counts)
      A=l*log(Ni2)+lgamma(Ni2/2)-lgamma(l+Ni2/2)
      B=-Ni2*lgamma(0.5)+sum(lgamma(cts+0.5))
      logp[i]=A+B
    }
    idx=which.max(logp)
    optN=idx+a-1
    if(plot) plot(a:b, logp, type='l')
    return(optN)
  }
  # Hogg method
  hhist_2d=function(x, a=2, b=100, alpha=1, plot=F){
    n=b-a+1
    rg1=diff(range(x[,1]))
    rg2=diff(range(x[,2]))
    L=numeric(n)
    for(i in 1:n){
      Nb=i+a-1
      Nb2=Nb^2
      D=rg1*rg2/Nb2
      Ni=c(gplots::hist2d(x[,1], x[,2], Nb, show=F)$counts)
      S=sum(Ni+alpha)
      Li=numeric(Nb2)
      for(j in 1:Nb2){
        Li[j]=Ni[j]*log((Ni[j]+alpha-1)/(D*(S-1)))
      }
      L[i]=sum(Li, na.rm = T)
    }
    idx=which.max(L)
    optN=idx+a-1
    if(plot) plot(a:b, L, type='l')
    return(optN)
  }
  N=NULL
  if(breaks %in% c('Hogg','H','h')){
    N=hhist_2d(x, a, b)
  }else if(breaks %in% c('Knuth','K','k')){
    N=khist_2d(x, a, b)
  }else if(breaks %in% c('Shimazaki','S','s')){
    N=sshist_2d(x, a, b)
  }else{
    stop("'breaks' should be one of 'Hogg', 'Knuth', 'Shimazaki'")
  }
  if(!is.null(N)){
    message('optimal number of bins: ', N)
    if(N==b)  message('Warning: number of maximum bins (',b,') reached')
    if(N>b)  message('Warning: number of maximum bins (',b,') excedeed')
  }
  if(plot){
    col=c('#9E0142','#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF', 
          '#E6F598','#ABDDA4','#66C2A5','#3288BD','#5E4FA2')
    col=rev(col)
    hst=gplots::hist2d(x[,1], x[,2], N, col=col, show=T)
    box()
    return(invisible(hst))
  }
}

Matrix2DataFrame=function(mat,x=NULL,y=NULL,xlab="x",ylab="y",zlab="z"){
  
  # check if x and y are difined, if not, gives default values
  if(is.null(y)){
    y=seq(0,dim(mat)[2]-1)
  }
  if(is.null(x)){
    x=seq(0,dim(mat)[1]-1)
  }
  
  # if dimensions do not match, return error message
  dim=c(length(x),length(y))
  if(sum(dim==dim(mat))!=2)
  {
    cat("ERROR: mat must be of dimension x,y.")
    return(invisible())
  }
  
  # create output data frame
  df=expand.grid(x,y)
  # add z column to output data frame
  df$z=c(mat)
  
  # name columns of output data frame and return
  colnames(df)=c(xlab,ylab,zlab)
  return(df)
}

kde2dWeighted=function (x, y, w, h, n, lims = c(range(x), range(y)),proba.min=1E-6){
  if (missing(n)){
    n=c(length(unique(x)),
        length(unique(y)))
  }
  nx <- length(x)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  n<-rep(n, length.out = 2L)
  gx <- seq(lims[1], lims[2], length = n[1]) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n[2]) # gridpoints y
  
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w,n[1]), nrow=n[1], ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n[1], nx)) %*% 
    t(matrix(dnorm(ay), n[2], nx))/(sum(w) * h[1] * h[2]) # z is the density
  
  z[z<proba.min]=0
  return(Matrix2DataFrame(mat=z,x=gx,y=gy))
}

oner_breaks=function(x, n, method='content'){
  require(OneR)
  bn=bin(x, n, method = method, labels = 1:n)
  sp=split(x, bn)
  bks=unlist(lapply(sp, range))
  d=diff(bks)/2+bks[-2*n]
  ib=d[2*(1:(n-1))]
  b=c(min(bks), ib, max(bks))
  names(b)=NULL
  tb=median(table(bn))
  message('Median content: ',tb)
  return(b)
}

# Image interpolation using fields::interp.surface
interp.img=function(x, y, z, nbin=100, plot=T, ...){
  require(fields)
  require(reshape2)
  require(magicaxis)
  if(class(x)=='list' & length(x)==3){
    xyz=x
    x=xyz[[1]]; y=xyz[[2]]; z=xyz[[3]]
  }else{
    xyz=list(x=x, y=y, z=z)
  }
  rx=range(x, na.rm=T)
  ry=range(y, na.rm=T)
  xs=seq(rx[1], rx[2], length.out = nbin)
  ys=seq(ry[1], ry[2], length.out = nbin)
  eg=expand.grid(xs, ys)
  ip=interp.surface(xyz, eg)
  df=(data.frame(x=eg[,1], y=eg[,2], z=ip))
  ac=acast(df, y~x, value.var = 'z')
  m=t(ac)
  l=list(x=xs, y=ys, z=m)
  if(plot) magimage(l, magmap = F, ...)
  return(invisible(l))
}

# Vector field for a z matrix using the gradient
vector.field=function(x, y, z, add=F, alog=F, arrow.col='black', arrow.scale=1,
                      ...){
  require(pracma)
  if(class(x)=='list' & length(x)==3){
    xyz=x; x=xyz[[1]]; y=xyz[[2]]; z=xyz[[3]]
  }
  g=gradient(z)
  gn=sqrt(g$X^2+g$Y^2)
  if(alog) gn=log10(gn+1)
  gp=atan2(g$X, g$Y)
  w=which(!is.na(z), arr.ind = T)
  n=gn[w]; ph=gp[w]
  d=data.frame(x=x[w[,1]],y=y[w[,2]],norm=n,phase=ph)
  d$x1=d$x+d$norm*arrow.scale*cos(d$phase)
  d$y1=d$y+d$norm*arrow.scale*sin(d$phase)
  if(!add) image(x, y, z, ...)
  arrows(d$x, d$y, d$x1, d$y1, length = 0.05, col=arrow.col)
  invisible(na.omit(d))
}

density.plot=function(x, y=NULL, n=100, col=NULL, fillcol=NA, norm=F, contours=F, 
                      nlevels=10, add=F, ...){
  require(magicaxis)
  if(is.null(y)){
    den=density(x, na.rm = T)
    if(norm){
      den$y=den$y/max(den$y)
    }
    if(is.null(col)) col='black'
    if(add){
      polygon(den$x, den$y, col=fillcol, border=NA)
      lines(den, col=col[1], ...)
    }else{
      magplot(den, col=col[1], ylab='Density', ...)
    }
  }else{
    require(MASS)
    d=na.omit(data.frame(x=x, y=y))
    den=kde2d(d$x, d$y, n=n)
    if(is.null(col)){
      col=c("#5E4FA2","#3288BD","#66C2A5","#ABDDA4","#E6F598","#FFFFBF", 
            "#FEE08B","#FDAE61","#F46D43","#D53E4F","#9E0142")
    }
    if(!add) magimage(den, col=col, magmap=F, asp=NA, ...)
    if(contours) contour(den, nlevels=nlevels, drawlabels=F, add=T, col=col)
  }
  return(invisible(den))
}

# Get radius limit in degrees for rlim (in Mpc) and redshift differences for vlim (in km/s)
dvol=function(ra, dec, z, rlim=1, vlim=1000, Om=0.3, Ol=0.7, H0=72){
  require(cosmoFns)
  delta_z=vlim*(1+z)/299792
  delta_rad=rlim/(D.A(z, Om, Ol, H0)*pi/180) # rlim (in Mpc) to degree conversion
  v=c('delta_rad'=delta_rad, 'delta_z'=delta_z)
  return(v)
}

iter.lm=function(x, y, type=1, plot=T, ...){
  d=na.omit(data.frame(x=sort(x), y=y[order(x)]))
  xseq=seq(min(x), max(x), length.out=1e2)
  xdif=xseq[100]-xseq[1]
  n1=ct=0
  while(n1!=nrow(d)){
    n1=nrow(d)
    ft=lm(y~x, d)
    pr=predict(ft, data.frame(x=xseq), interval='prediction')
    b1=(pr[100,2]-pr[1,2])/xdif
    a1=pr[1,2]-b1*xseq[1]
    b2=(pr[100,3]-pr[1,3])/xdif
    a2=pr[1,3]-b2*xseq[1]
    if(type==1) w=which(d$y>b1*d$x+a1)
    if(type==2) w=which(d$y<b2*d$x+a2)
    if(type==3) w=which(d$y>b1*d$x+a1 & d$y<b2*d$x+a2)
    d=d[w,]
    ct=ct+1
  }
  message('iterations: ', ct)
  if(plot){
    plot(x, y, ...)
    abline(ft, col=2)
    lines(xseq, pr[,2], col=2, lty=2)
    lines(xseq, pr[,3], col=2, lty=2)
  }
  cf=coef(ft)
  df=data.frame(fit=c(cf), lwr=c(a1,b1), upr=c(a2,b2))
  rownames(df)=c('a','b')
  return(round(df,5))
}

# Aperture multipole moments of 2D coordinates within radii r1 and r2 in Mpc
# A&A 635, A195 (2020)
multip.mom=function(ra, dec, ra0, dec0, z0, r1=0, r2=1, m=0:10, N=1e3, 
                    Om=0.27, Ol=1-Om, H0=72, plot=T){
  require(cosmoFns)
  require(plotrix)
  #require(spatstat)
  da=D.A(z0, Om, Ol, H0) # Angular diameter distance in Mpc
  dz=sin(dec0*pi/180)*sin(dec*pi/180)
  dr=cos(dec0*pi/180)*cos(dec*pi/180)*cos((ra0-ra)*pi/180)
  R=da*acos(dz+dr)
  #R=(sqrt((ra-ra0)^2+(dec-dec0)^2)*pi/180)*da # Polar radii in Mpc
  phi=atan2(dec-dec0, ra-ra0) # Polar angles in radians
  dr=r2-r1
  rc=(r2+r1)/2
  w=(1/dr)*exp(-(rc-R)^2/(2*dr^2)) # Window function
  n=length(ra)
  nm=length(m)
  Qm=matrix(NA, N, nm)
  for(i in 1:N){
    idx=sample(1:n, n, replace = T) # Resampling
    for(j in 1:nm) Qm[i,j]=abs(sum(w[idx]*exp(1i*m[j]*phi[idx])))
  }
  qf=function(x) 0.7415*(quantile(x,.75)-quantile(x,.25))
  loc=apply(Qm, 2, median)
  scl=apply(Qm, 2, qf)
  qm=loc[m>0]
  w=(qm-min(qm))/(max(qm)-min(qm))
  Pm.med=spatstat::weighted.median(m[m>1], w[m>1]) # Median angular scale
  Pm.err=0.7415*diff(spatstat::weighted.quantile(m[m>1], w[m>1], c(.25,.75)))
  lwr=loc-scl
  upr=loc+scl
  if(plot){
    op=par(mai=c(0.55,0.55,0.1,0.1))
    layout(matrix(c(1,1,1,1,2,2), 3, 2, T), heights = c(1,1,1))
    plot(ra, dec, pch=20, cex=.7, xlab='RA', ylab='Dec', asp=1)
    draw.circle(ra0, dec0, r1/(da*pi/180), border=2)
    draw.circle(ra0, dec0, r2/(da*pi/180), border=2)
    w=which(m>0)
    plot(m[w], loc[w], type='l', ylim = c(min(lwr[w]), max(upr[w])), xlab='m', ylab='Qm')
    polygon(c(m[w],rev(m[w])), c(lwr[w],rev(upr[w])), col='#00000033', border = NA)
    layout(1)
    par(op)
  }
  df=data.frame(m=m, Qm=loc, Qm_err=scl)
  if(any(m==0)) df$Qm_Q0=loc/loc[m==0]
  names(Pm.err)=NULL
  pm=round(c('max'=m[m>1][which.max(loc[m>1])], 'median'=Pm.med, 'error'=Pm.err), 3)
  return(list('Multipole_moments'=df, 'Angular_scale'=pm))
}

# Sigma clipping
sigma.clip=function(x, nclip=3, N.max=5){
  mean=mean(x); sigma=sd(x)
  clip.lo=mean-(nclip*sigma)
  clip.up=mean+(nclip*sigma)
  x=x[x<clip.up & x>clip.lo]
  if(N.max>0){
    N.max=N.max-1
    x=Recall(x, nclip=nclip, N.max=N.max)
  }
  return(x)
}

# Analogous to sigma clipping but using robust stats
scale.clip=function(x, nclip=3, N.max=5){
  require(RobStatTM)
  lsm=locScaleM(x, 'bisquare', na.rm=T)
  mu=lsm$mu; disp=lsm$disper
  clip.lo=mu-(nclip*disp)
  clip.up=mu+(nclip*disp)
  x=x[x<clip.up & x>clip.lo]
  if(N.max>0){
    N.max=N.max-1
    x=Recall(x, nclip=nclip, N.max=N.max)
  }
  return(x)
}

# Binning function to use in shifting-gapper
binning=function(x, dx=0.6, nmin=15){
  # dx: minimum size of bins (Mpc)
  # nmin: minimum number of object per bin
  x.max=max(x)
  x=sort(x)
  n=length(x)
  xb=x[1]
  xlim=0
  while(n>=nmin){
    xlim=ifelse(x[nmin]<dx+xlim, dx+xlim, x[nmin])
    xb=c(xb, xlim)
    x=x[x>xlim]
    n=length(x)
  }
  xb[length(xb)]=x.max
  return(xb)
}

# Shifting-gapper technique following Lopes et al. (2009) recipe
shifting.gapper=function(x, y, dx=0.6, nmin=15, dxmax=0.6, vgap=300, vcenter=T, plot=T){
  # require binning()
  # dx: minimum size of bins (Mpc)
  # nmin: minimum number of object per bin
  # dxmax: maximum distance allowed between objects (Mpc)
  # vcenter: if TRUE re-defines the center of vlos inside 1 Mpc
  require(RobStatTM) # robust statistics
  N=length(x)
  wd=which(duplicated(y))
  if(length(wd)>0){
    #message('Warning: ', length(wd), ' duplicated velocities')
    dy=runif(length(wd), 1e-4, 2e-4)
    y[wd]=y[wd]+dy
  }
  if(vcenter) y=y-locScaleM(y[x<1], psi='bisquare')$mu
  x.=sort(x)
  y.=y[order(x)]
  w.inc=1:N
  wo=c()
  n=1
  # Loop until no more otliers are detected
  while(n!=0){
    x.=x.[w.inc]
    y.=y.[w.inc]
    br=binning(x., dx, nmin)
    Nb=length(br)-1
    fint=findInterval(x., br, F, T, T)
    y0=y.[x.<min(br[2],1)]
    #vlim=max(abs(range(sigma.clip(y0, 2.5))))
    #y0=y.[x.<1]
    lsM=locScaleM(y0, psi='bisquare')
    vlim=3*lsM$disper
    v.out=c()
    # Loop for every bin running the shifting-gapper
    for(i in 1:Nb){
      wi=fint==i
      xi=x.[wi]
      yi=y.[wi]
      dxi=diff(xi)
      wdx=which(dxi>=dxmax)
      if(length(wdx)>0){ # radial offset > 0.6 Mpc
        xoff=xi[wdx[1]]
        v.out=c(v.out, y[x>xoff])
        break;
      }
      w=which(abs(yi)>vlim)
      v.out=c(v.out, yi[w])
      vinf=sort(yi[yi<0])
      vsup=sort(yi[yi>0])
      
      facgap=max(vlim/5, vgap)
      Ni=sum(wi)
      dv=facgap*(1+exp(-(Ni-6)/33))
      
      winf=which(diff(vinf)>dv)
      wsup=which(diff(vsup)>dv)
      any.inf=length(winf)>0
      any.sup=length(wsup)>0
      if(any.inf) v.out=c(v.out, vinf[1:max(winf)])
      if(any.sup) v.out=c(v.out, vsup[(wsup[1]+1):length(vsup)])
      
      valid=!any.inf | (any.inf & !any.sup)
      if(valid){
        vlim.=max(abs(yi))
        if(any.inf & !any.sup & length(vsup)>0) vlim.=max(vsup)
        if(!any.inf & any.sup & length(vinf)>0) vlim.=max(abs(vinf))
        vlim=(vlim+vlim.)/2
      }
    }
    w.out=which(y %in% v.out)
    w.inc=which(!y. %in% v.out)
    wo=c(wo, w.out)
    n=length(v.out)
  }
  outl=rep(0, N)
  outl[wo]=1
  if(plot){
    plot(x, y, pch=c(20,1)[outl+1], ylim=c(-4000,4000), col=outl+1,
         xlab=expression(R[proj]~(Mpc)), ylab=expression(V[los]~(km/s)),
         panel.first=abline(h=0, col='#0000001A', lty=2))
  } 
  return(invisible(outl))
}

# vector as color for 3D scatterplot representation based on astro::aplot
ztocol=function(z, zcol=NULL){
  if(is.null(zcol)) zcol=hsv(seq(2/3, 0, len = 200))
  zlim=range(z, na.rm=T)
  base=(z-zlim[1])/(zlim[2]-zlim[1])
  col=zcol[(base*(length(zcol)-1))+1]
  return(col)
}

# convert keyword position to user coordinates
postocoord=function(pos='topleft', n, inset=0.1, cex=1, height.cex=0.05){
  usr = par("usr")
  xlo2 = usr[1]; xhi2 = usr[2]
  ylo2 = usr[3]; yhi2 = usr[4]
  xdim = par("pin")[1]; ydim = par("pin")[2]
  
  xlo = xlo2 + (((xhi2 - xlo2)/xdim) * inset)
  xhi = xhi2 - (((xhi2 - xlo2)/xdim) * inset)
  ylo = ylo2 + (((yhi2 - ylo2)/ydim) * inset)
  yhi = yhi2 - (((yhi2 - ylo2)/ydim) * inset)
  
  xmid = ((xhi2 - xlo2)/2) + xlo2
  ymid = ((yhi2 - ylo2)/2) + ylo2
  #textwidth = strwidth(1)*(4*n) * cex
  #textheight = strheight(1)*0.5 * cex
  textwidth = 0.15*(4*n) * cex
  textheight = height.cex * cex
  whitespace = 0.01
  rectwidth = textwidth + (((xhi2 - xlo2)/xdim) * whitespace)
  rectheight = textheight + (((yhi2 - ylo2)/ydim) * whitespace)
  if (pos == "topleft") {
    xleft = xlo
    xright = xlo + rectwidth
    ybottom = yhi - rectheight
    ytop = yhi
  }
  else if (pos == "top") {
    xleft = xmid - rectwidth/2
    xright = xmid + rectwidth/2
    ybottom = yhi - rectheight
    ytop = yhi
  }
  else if (pos == "topright") {
    xleft = xhi - rectwidth
    xright = xhi
    ybottom = yhi - rectheight
    ytop = yhi
  }
  else if (pos == "right") {
    xleft = xhi - rectwidth
    xright = xhi
    ybottom = ymid - rectheight/2
    ytop = ymid + rectheight/2
  }
  else if (pos == "bottomright") {
    xleft = xhi - rectwidth
    xright = xhi
    ybottom = ylo
    ytop = ylo + rectheight
  }
  else if (pos == "bottom") {
    xleft = xmid - rectwidth/2
    xright = xmid + rectwidth/2
    ybottom = ylo
    ytop = ylo + rectheight
  }
  else if (pos == "bottomleft") {
    xleft = xlo
    xright = xlo + rectwidth
    ybottom = ylo
    ytop = ylo + rectheight
  }
  else if (pos == "left") {
    xleft = xlo
    xright = xlo + rectwidth
    ybottom = ymid - rectheight/2
    ytop = ymid + rectheight/2
  }
  else if (pos == "centre") {
    xleft = xmid - rectwidth/2
    xright = xmid + rectwidth/2
    ybottom = ymid - rectheight/2
    ytop = ymid + rectheight/2	
  }
  return(c(xleft, ybottom, xright, ytop))
}

# Add horizontal colour bar to a plot
colorbar=function(pos='topleft', zmin, zmax, zcol=NULL, n=5, cex=1, inset=0.1, dig=2){
  # depends on postocoord()
  require(plotrix)
  if(is.null(zcol)) zcol=hsv(seq(2/3, 0, len = 200))
  leg=seq(zmin, zmax, length.out = n)
  leg=round(leg, digits = dig)
  rcol=colorRampPalette(zcol)(n)
  cd=postocoord(pos, n=n, cex=cex, inset=inset)
  align=ifelse(grepl('bottom', pos), 'lt', 'rb')
  color.legend(cd[1], cd[2], cd[3], cd[4], leg, rcol, cex, align = align)
}

# get nearest objects to data rows with columns ra & dec from a table to be matched 
# and returns the table indexes within a given radius in arcsec
getnearest=function(data, table, radius=0.5, plot=T){
  stopifnot(ncol(data)==2 & ncol(table)==2)
  require(RANN)
  nn=nn2(table, data, k=1)
  idx=nn$nn.idx[,1]
  dist=nn$nn.dists[,1]*3600
  idx[dist > radius]=NA
  perc=round(sum(!is.na(idx))*100/nrow(data), 1)
  message(paste0('Done. ',perc,'% matched'))
  if(plot){
    require(magicaxis)
    maghist(dist, 'fd', log='x', xlab='Distance (arcsec)', ylab='Counts', 
            col='gray', border=NA, side=1:4, labels=c(1,1,0,0), verbose=F)
    abline(v=radius, h=0, col=c('gray',2))
  }
  #dist[dist > radius]=NA
  #return(list(idx=idx,dist=dist))
  return(idx)
}

# Get projected radii and velocities from coordinates and redshift
pps=function(ra, dec, z, raclus, declus, zclus, Om=.3, Ol=1-Om, H0=67){
  require(cosmoFns)
  vlos=299792*(z-zclus)/(1+zclus)
  groupdist=D.A(zclus, Om, Ol, H0) # in Mpc
  zsep=sin(declus*pi/180)*sin(dec*pi/180)
  xysep=cos(declus*pi/180)*cos(dec*pi/180)*cos((raclus-ra)*pi/180)
  xyzsep=zsep+xysep
  xyzsep[xyzsep>1]=1
  angle=acos(xyzsep)
  dproj=groupdist*angle # in Mpc
  # see https://arxiv.org/pdf/1711.10018.pdf https://arxiv.org/pdf/1902.05276.pdf
  #     http://publications.lib.chalmers.se/records/fulltext/159140.pdf p.12
  return(data.frame(dproj, vlos))
}

# Limiting mass completeness calculation following Pozzetti et al. (2010)
masscomp=function(z, mag, mass, nbin=10, ...){
  df=na.omit(data.frame(z, mag, mass))
  breaks=seq(min(df$z), max(df$z), length.out=nbin)
  zc=(breaks[-1]+breaks[-length(breaks)])/2
  ct=cut(df$z, breaks, F)
  df$ct=ct
  dft=data.frame(table(ct))
  df=merge(df, dft, by='ct', sort = F)
  df$Freq=floor(df$Freq/5)
  sp=split(df, df$ct)
  sp=lapply(sp, function(x) x[order(x$mag, decreasing = T),])
  Mlim=unlist(lapply(sp, function(x) quantile(x[1:x$Freq[1],]$mass, .95)))
  model=nls(Mlim~a*zc^2+b*zc+c, start=list(a=1, b=1, c=0))
  p=coef(model)
  fitline=function(x, a, b, c){a*x^2+b*x+c}
  # plotting
  plot(z, mass, pch='.', ...)
  points(zc, Mlim, pch=15, col='red')
  curve(fitline(x, p[1], p[2], p[3]), min(df$z), max(df$z), add=T, col='red')
  return(model)
}

# Classifies points per squared bin in a 2D histogram
hist2Dclass=function(x, y, nbins=c(10,10), breaks=NULL){
  if(!is.null(breaks)){
    x.cuts=breaks[[1]]
    y.cuts=breaks[[2]]
  }else{
    if(length(nbins)==1) nbins=rep(nbins,2)
    x.cuts=seq(from = min(x), to = max(x), length = nbins[1] + 1)
    y.cuts=seq(from = min(y), to = max(y), length = nbins[2] + 1)
  }
  index.x=cut(x, x.cuts, include.lowest = TRUE)
  index.y=cut(y, y.cuts, include.lowest = TRUE)
  nx=as.numeric(index.x)
  ny=as.numeric(index.y)
  p=0.5*(nx+ny)*(nx+ny+1)+nx
  k=match(p, sort(unique(p)))
  #m=tapply(x, list(index.x, index.y), length)
  return(k)
}

# Finding which properties are more correlated with another one 
# following Blanton et al. (2005) approach
blanton.pred=function(x, n=5, method='length'){
  # require hist2Dclass()
  # x must be a matrix with first column being the variable of reference
  # n: number of bins
  # one-dimensional case
  require(OneR)
  require(plotrix)
  xnames=colnames(x)
  if(!is.matrix(x)) x=as.matrix(x)
  nc=ncol(x)-1
  var_Y=var(x[,1])
  f=function(x){sum((x-(sum(x)/length(x)))^2)}
  var_X=numeric(nc)
  brk=vector('list', length = nc)
  for(i in 1:nc){
    YX=na.omit(x[,c(1,i+1)])
    Y=YX[,1]
    X=YX[,2]
    bin_class=bin(X, n, 1:n, method = method)
    l=split(Y, bin_class)
    s=unlist(lapply(l, f))
    var_X[i]=sum(s)/(length(Y)-1)
    if(method=='content'){
      lx=split(X, bin_class)
      rg=unlist(lapply(lx, range))
      d=diff(rg)/2+rg[-2*n]
      ib=d[2*(1:(n-1))]
      b=c(min(rg), ib, max(rg))
      names(b)=NULL
      brk[[i]]=b
    }else{
      brk[[i]]=seq(min(X), max(X), length.out=n+1)
    }
  }
  vdif1=var_X-var_Y
  d1=data.frame(var_diff=round(vdif1,4), row.names = xnames[-1])
  # two-dimensional case
  comb=combn(nc, 2)
  n.comb=combn(xnames[-1],2)
  pn=apply(n.comb, 2, paste, collapse='-')
  ncb=ncol(comb)
  var_XX=hce=numeric(ncb)
  hcs=vector('list', length=ncb)
  for (i in 1:ncb){
    idx=comb[,i]+1
    YXX=na.omit(x[,c(1,idx)])
    Y=YXX[,1]
    if(method=='content'){
      hc=hist2Dclass(YXX[,2], YXX[,3], n, breaks=brk[idx-1])
    }else{
      hc=hist2Dclass(YXX[,2], YXX[,3], n)
    }
    l=split(Y, hc)
    s=unlist(lapply(l, f))
    var_XX[i]=sum(s)/(length(Y)-1)
    hct=as.numeric(table(hc))
    hce[i]=std.error(hct)
    hcs[[i]]=summary(hct)
  }
  hcs=data.frame(do.call(rbind,hcs))
  hcs=cbind(hcs,'std.error'=round(hce,2))
  rownames(hcs)=pn
  vdif2=var_XX-var_Y
  hcs$var_diff=round(vdif2,5)
  # creating output dataframe
  d2=data.frame(matrix(NA, nc, nc), row.names = xnames[-1])
  colnames(d2) = xnames[-1]
  d2[t(comb)]=round(vdif2,4)
  d=cbind(d1,d2)
  l1=xnames[-1][order(var_X)]
  l2=pn[order(var_XX)]
  l3=data.frame(do.call(rbind,brk))
  rownames(l3)=xnames[-1]
  l=list(d, 'best ordered single'=l1, 'best ordered pairs'=l2, 
         '1D & 2D breaks'=l3, '2D binning counts statistics'=hcs)
  names(l)[[1]]=paste('variables as predictors of',xnames[1])
  return(l)
}

# weighted variance
weighted.var=function(x, w, na.rm = FALSE){
  if(na.rm){
    w=w[i <- !is.na(x)]
    x=x[i]
  }
  sum.w=sum(w)
  sum.w2=sum(w^2)
  mean.w=sum(x*w)/sum(w)
  wvar=(sum.w/(sum.w^2-sum.w2))*sum(w*(x-mean.w)^2, na.rm=na.rm)
  return(wvar)
}

# Find redshift limit with maximum number of objects in a complete sample
zopt=function(z, Mr, OM=0.3, OL=1-OM, H0=70){
  # require absmag()
  d=cbind(z, Mr)
  f=function(z, d) sum(d[,1]<=z & d[,2]<=absmag(z, OM=OM, H0=H0))
  n=optimise(f, range(z), d, maximum = T)
  r=c('z'=n[[1]], 'mag'=absmag(n[[1]], OM=OM, H0=H0), 'N'=n[[2]])
  return(r)
}

# Geometric Histogram Separation for using in hd.test function
ghs.hd=function(x1, x2, plot=T){
  h1=hist(x1, 'fd', plot = F)
  h2=hist(x2, 'fd', plot = F)
  both=c(h1$breaks, h2$breaks)
  n=length(both)-1 # number of total bins
  rnge=range(both)
  breaks=seq(rnge[1], rnge[2], length.out = n+1)
  h1=hist(x1, breaks, plot = F)
  h2=hist(x2, breaks, plot = F)
  dx=diff(rnge)/n
  dy=apply(cbind(h1$density, h2$density), 1, min) # min of intersection
  ao=sum(dx*dy) # the relative area 
  a_height=max(h1$density)
  b_height=max(h2$density)
  c_height=max(dy)  
  bcl=(a_height+b_height-2*c_height)/(a_height+b_height)
  bca=1-ao/(2-ao)
  ghs=(sqrt(bca)+bcl)/2
  R2=cor(h1$counts,h2$counts)^2
  if(plot){
    ymax=max(h1$counts,h2$counts)
    plot(c(rnge), c(0, ymax), pch='', ylab='Frequency', xlab='log(HD)')
    abline(h=0, col='gray')
    lines(c(breaks[1],breaks), c(0,h1$counts,0), type='s', col='red') 
    lines(c(breaks[1],breaks), c(0,h2$counts,0), type='s', col='black') 
    legend('topleft', lty=1, col=c('red','black'), bty='n', x.intersp=.5, 
           legend=c('Simulation','Resampling'))
  }
  return(list('ghs'=ghs,'R2'=R2))
}

# Gaussianity indicator using the HD distance separation
# require the HD_simulations.RDS file!
hd.test=function(x, Nsamp=1000, bd.file='HD_simulations.RDS', plot=T){
  # x: vector of the LOS cluster velocities inside R200
  require(distrEx)
  db=readRDS(bd.file)
  n=length(x)
  if(length(x)<10) 
    stop('too few elements - at least 10 values are needed')
  sim=db[,n]
  hns=function(x) (4/(length(x)*3))^(1/5)*sd(x)
  x=(x-mean(x))/sd(x)
  hdr=vector(length=Nsamp)
  pb=txtProgressBar(1,Nsamp,1,'=',style=3)
  for(i in 1:Nsamp){
    samp=sample(x, n, replace=T)
    hs=hns(samp)
    hdr[i]=HellingerDist(Norm(), samp, asis.smooth.discretize="smooth", h.smooth=hs)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  comp=ghs.hd(log10(sim), log10(hdr), plot=plot)
  u=1-comp[[1]]
  R2=comp[[2]]
  #r=ghc(log10(sim), log10(hdr), plot=plot)
  if(plot){
    mtext(bquote(N==.(n)), line=.5, adj=.1)
    mtext(bquote(ghu==.(round(u,3))), line=.5, adj=.5)
    mtext(bquote(R^2==.(round(R2,3))), line=.5, adj=.9)
    legend('topleft', lty=1, col=c('red','black'), bty='n', x.intersp=.5, 
           legend=c('Simulation','Resampling'))
  } 
  res=list('ghu'=u, 'R2'=R2, 'HD_values'=hdr)
  return(invisible(res))
}
