#this is a comment in the first 50 lines for Vim's filetype detection.
library(gplots);
library(RSQLite);
library(akima);
library(plotrix);
library(fields);
library(ks);
library(boot);
library(Hmisc);
library(splancs);
library(ggplot2);
library(plyr);
library(parallel);
library(reshape);
source("~/R/colormaps.r");


readData2sqlite <- function(fn,tabName){
  con <- dbConnect(SQLite(),fn);
  q <- dbSendQuery(con,statement=sprintf("select * from %s",tabName));
  x <- fetch(q,-1);
  dbClearResult(q);
  dbDisconnect(con);
  x;
}

getDataConnection <- function(fn){
  con <- dbConnect(SQLite(),fn);
  f <- function(query){
    q <- dbSendQuery(con,statement=query);
    x <- fetch(q,-1);
    dbClearResult(q);
    x;
  }
  f;
}

plotDRMTheta <- function(r,ee,nbins=50,breaks=seq(0,180,length.out=nbins+1),mult=5,add=FALSE,...){
  hp <- hist(r$thp,breaks=breaks,plot=FALSE);
  ha <- hist(r$th,breaks=breaks,plot=FALSE);
  e <- ha$mids;
  a <- ha$counts/hp$counts*100*100/mult;
  if(add){
    lines(e,a,type='l',...);
  }else{
    plot(e,a,type='l',...); #,xlab=expression(theta*"  "*(degree)),ylab=expression(A[eff]*"  "*(cm^2)),main="effective area, 3 MeV photons [PRELIMINARY]");
  }
  mtext(ee,side=4,at=a[length(a)],cex=0.6);
}

plotDRMPhi <- function(r,ee,nbins=50,breaks=seq(-180,180,length.out=nbins+1),mult=5,add=FALSE,...){
  hp <- hist(r$php,breaks=breaks,plot=FALSE);
  ha <- hist(r$ph,breaks=breaks,plot=FALSE);
  e <- ha$mids;
  a <- ha$counts/hp$counts*100*100/mult;
  if(add){
    lines(e,a,type='l',...);
  }else{
    plot(e,a,type='l',...); #,xlab=expression(theta*"  "*(degree)),ylab=expression(A[eff]*"  "*(cm^2)),main="effective area, 3 MeV photons [PRELIMINARY]");
  }
  mtext(ee,side=4,at=a[length(a)],cex=0.6);
}

readDRM <- function(dbfn,fullEnergy=FALSE){
  db <- getDataConnection(dbfn);
  a <- db("select * from evtsBGO, pcles where evtsBGO.priIdx = pcles.idx and evtsBGO.totEDep>0.01");
  #a <- db("select * from evtsCZT, pcles where evtsCZT.priIdx = pcles.idx and evtsCZT.totEDep>0.01");
  p <- db("select * from pcles");

  e0 <- p$ee[1];
  sel <- if(fullEnergy){a$totEDep>0.95*e0;}else{a$totEDep>0;}

  thp <- atan2(sqrt(p$px**2+p$py**2),-p$pz)*180/pi;
  th <- atan2(sqrt(a$px[sel]**2+a$py[sel]**2),-a$pz[sel])*180/pi;
  php <- atan2(p$py,p$px)*180/pi;
  ph <- atan2(a$py[sel],a$px[sel])*180/pi;
  list(th=th,ph=ph,thp=thp,php=php);
}

plotAngleDRMs <- function(plotDRMFunc,fullEnergy=FALSE,breaks=seq(90,180,length.out=20),...){
  fnames <- c("test_22_1e+02.db","test_30.db","test_0.3.db","test_0.1.db","test_1.db","test_10.db","test_3.db");
  ees <- c(100,30,0.3,0.1,1,10,3);
  add<-FALSE;
  for(i in seq(length(fnames))){
    fn <- fnames[i];
    ee <- ees[i];
    r <- readDRM(fn,fullEnergy);
    plotDRMFunc(r,ee,xaxs='int',breaks=breaks,xlab=expression(theta*"  ("*degree*")"),ylab=expression("effective area (cm"^2*")"),yaxs="int",col=i,add=add,...);
    add<-TRUE;
  }
}

plotEnergyDRM <- function(fE){
  fnames <- c("test_22_1e+02.db","test_30.db","test_0.3.db","test_0.1.db","test_1.db","test_10.db","test_3.db");
  ees <- c(100,30,0.3,0.1,1,10,3);

  effarr <- function(dbfn,fullEnergy){
    r <- readDRM(dbfn,fullEnergy);
    sel <- r$th>160;
    selp <- r$thp>160;
    length(r$th[sel])/length(r$thp[selp])*100*100/5;
  }

  area <- mapply(effarr,fnames,fE);
  plot(ees,area);
  list(ee=ees,area=area);
}

measureEfficiency <- function(fn="filenameGoesHere",de=0,dp=0,mult=5,dist=-1){
  print(fn);
  if(dist==-1){
    # get distance from filename, like "eff_22_0.36_0.db" or "eff_22_0.36_0.2.db".
    dist = as.real(sub("eff_22_.*_(.*)\\.db","\\1",fn));
  }
  if(!is.data.frame(de)){
    db <- getDataConnection(fn);
    de <- db("select * from evtsBGO, pcles where evtsBGO.priIdx = pcles.idx");
  }
  if(!is.data.frame(dp)){
    db <- getDataConnection(fn);
    dp <- db("select * from pcles");
  }
  efficiency <- function(eth=-1){
    if(eth == -1){
      eth <- dp$ee[1]*0.97;
    }
    length(de$totEDep[de$totEDep>eth])/length(dp$ee)/mult;
  }

  #eths <- c(0.02,0.1,0.2,-1);
  #mapply(efficiency,eths);

  data.frame(ee=dp$ee[1],dist=dist,gt20=efficiency(0.02),gt100=efficiency(0.1),gt200=efficiency(0.2),fullEnergy=efficiency(-1));
}

efficiencyDataSet <- function(){
  rbind(measureEfficiency("eff_22_0.36_0.db"),
        measureEfficiency("eff_22_0.356_0.05.db"),
        measureEfficiency("eff_22_0.36_0.1.db"),
        measureEfficiency("eff_22_0.36_0.2.db"),
        measureEfficiency("eff_22_0.36_0.3.db"),
        measureEfficiency("eff_22_0.36_0.4.db"),
        measureEfficiency("eff_22_0.36_0.5.db"),

        measureEfficiency("eff_22_0.51_0.db"),
        measureEfficiency("eff_22_0.511_0.05.db"),
        measureEfficiency("eff_22_0.51_0.1.db"),
        measureEfficiency("eff_22_0.51_0.2.db"),
        measureEfficiency("eff_22_0.51_0.3.db"),
        measureEfficiency("eff_22_0.51_0.4.db"),
        measureEfficiency("eff_22_0.51_0.5.db"),

        measureEfficiency("eff_22_0.66_0.db"),
        measureEfficiency("eff_22_0.662_0.05.db"),
        measureEfficiency("eff_22_0.66_0.1.db"),
        measureEfficiency("eff_22_0.66_0.2.db"),
        measureEfficiency("eff_22_0.66_0.3.db"),
        measureEfficiency("eff_22_0.66_0.4.db"),
        measureEfficiency("eff_22_0.66_0.5.db"),

        measureEfficiency("eff_22_0.83_0.db"),
        measureEfficiency("eff_22_0.835_0.05.db"),
        measureEfficiency("eff_22_0.83_0.1.db"),
        measureEfficiency("eff_22_0.83_0.2.db"),
        measureEfficiency("eff_22_0.83_0.3.db"),
        measureEfficiency("eff_22_0.83_0.4.db"),
        measureEfficiency("eff_22_0.83_0.5.db"),

        measureEfficiency("eff_22_1.2_0.db"),
        measureEfficiency("eff_22_1.17_0.05.db"),
        measureEfficiency("eff_22_1.2_0.1.db"),
        measureEfficiency("eff_22_1.2_0.2.db"),
        measureEfficiency("eff_22_1.2_0.3.db"),
        measureEfficiency("eff_22_1.2_0.4.db"),
        measureEfficiency("eff_22_1.2_0.5.db"),

        measureEfficiency("eff_22_1.27_0.db"),
        measureEfficiency("eff_22_1.27_0.05.db"),
        measureEfficiency("eff_22_1.27_0.1.db"),
        measureEfficiency("eff_22_1.27_0.2.db"),
        measureEfficiency("eff_22_1.27_0.3.db"),
        measureEfficiency("eff_22_1.27_0.4.db"),
        measureEfficiency("eff_22_1.27_0.5.db"),

        measureEfficiency("eff_22_1.3_0.db"),
        measureEfficiency("eff_22_1.33_0.05.db"),
        measureEfficiency("eff_22_1.3_0.1.db"),
        measureEfficiency("eff_22_1.3_0.2.db"),
        measureEfficiency("eff_22_1.3_0.3.db"),
        measureEfficiency("eff_22_1.3_0.4.db"),
        measureEfficiency("eff_22_1.3_0.5.db"));
}

# i.e. plotEfficiencies(eds,eds$gt200)
plotEfficienciesVsDist <- function(eds,eff,title){
  ens <- unique(eds$ee);
  plot(0,0,xlim=c(0,50),ylim=c(0,100),xlab="distance (cm)",ylab="efficiency (%)",type='n',main=title);
  
  ctr <- 1;
  for(e in ens){
    ord = order(eds$dist[eds$ee==e]);
    lines(eds$dist[eds$ee==e][ord]*100,eff[eds$ee==e][ord]*100,col=colorBlindPalette[ctr],lwd=2);
    mtext(sprintf("%d keV",e*1000),4,at=eff[eds$ee==e][ord][length(eff[eds$ee==e])]*100,las=1,line=0.2,col=colorBlindPalette[ctr],cex=0.5);
    ctr <- ctr+1;
  }

}
plotEfficienciesVsEnergy <- function(eds,eff,title){
  dists <- unique(eds$dist);
  plot(0,0,xlim=c(0,1.4),ylim=c(0,100),xlab="energy (MeV)",ylab="efficiency (%)",type='n',main=title);
  
  ctr <- 1;
  for(d in dists){
    ord = order(eds$ee[eds$dist==d]);
    lines(eds$ee[eds$dist==d][ord],eff[eds$dist==d][ord]*100,col=colorBlindPalette[ctr],lwd=2);
    text(0.1,5+5*ctr,sprintf("%d cm",d*100),col=colorBlindPalette[ctr],cex=1.0,pos=4);
    ctr <- ctr+1;
  }

}

# read a BGO/CZT DRM pair
readDRMs <- function(fn){
  f <- file(fn,"rt");
  l <- readLines(f,5);
  close(f);
  priLine <- as.real(strsplit(l[3]," ")[[1]][-1:-3])
  outLine <- as.real(strsplit(l[4]," ")[[1]][-1:-5])
  print(l[1]);
  print(l[2]);

  a <- read.table(fn);
  bdf <- data.matrix(subset(a,a$V1=="BGO")[,-1]);
  cdf <- data.matrix(subset(a,a$V1=="CZT")[,-1]);

  om <- (outLine[-1] + outLine[-length(outLine)])/2;

  list(b=bdf,c=cdf,i=priLine,om=om,obks=outLine);
}

#vectorize over x.
binomCI <- function(x,n,conf.lev=0.6826895){
  #print(c(x,n,conf.lev));
  alpha <- 1-conf.lev;
  p.L <- function(x, alpha){
    y <- x; y[x==0] <- 1;
    ifelse(x == 0, 0, qbeta(alpha, y, n - y + 1));
  }
  p.U <- function(x, alpha){
    y <- x; y[x==n] <- 1;
    ifelse(x == n, 1, qbeta(1 - alpha, y + 1, n - y));
  }

  matrix(c(p.L(x, alpha), p.U(x, alpha)),ncol=2);
}

# effective areas measured in cm^2/keV
readDRMs_df <- function(fn,nPriPerE=1.0,rDisk1=1.0,rDisk0=0.0){
  f <- file(fn,"rt");
  l <- readLines(f,5);
  close(f);
  priLine <- as.real(strsplit(l[3]," ")[[1]][-1:-3])
  outLine <- as.real(strsplit(l[4]," ")[[1]][-1:-5])
  #print(l[1]);
  #print(l[2]);

  print("reading file, loading matrices...");
  a <- read.table(fn);
  bdf <- data.matrix(subset(a,a$V1=="BGO")[,-1]);
  cdf <- data.matrix(subset(a,a$V1=="CZT")[,-1]);

  om <- (outLine[-1] + outLine[-length(outLine)])/2;
  deo <- (outLine[-1] - outLine[-length(outLine)]);

  #list(b=bdf,c=cdf,i=priLine,om=om,obks=outLine);

  bdf <- data.frame(melt(t(matrix(bdf,ncol=length(om),dimnames=list(priLine,om)))));
  cdf <- data.frame(melt(t(matrix(cdf,ncol=length(om),dimnames=list(priLine,om)))));

  e1 <- outLine[findInterval(bdf$X1,outLine)];
  e2 <- outLine[findInterval(bdf$X1,outLine)+1];

  #print("calculating error bars...");
  #pb <- create_progress_bar("text"); pb$init(length(bdf$value)+length(cdf$value));
  #bcis <- mapply(function(suc,try){pb$step(); binom.test(suc,try,conf.level=0.6826895)$conf.int},bdf$value,nPriPerE);
  #ccis <- mapply(function(suc,try){pb$step(); binom.test(suc,try,conf.level=0.6826895)$conf.int},cdf$value,nPriPerE);
  #pb$term();

  ## multi-core apply, doesn't work, for now... results are the wrong shape, or something.
  #bcis <- mclapply(bdf$value,function(suc,try){pb$step(); binom.test(suc,try,conf.level=0.6826895)$conf.int},nPriPerE,mc.cores=4);
  #ccis <- mclapply(cdf$value,function(suc,try){pb$step(); binom.test(suc,try,conf.level=0.6826895)$conf.int},nPriPerE,mc.cores=4);

  bcis <- binomCI(bdf$value,nPriPerE);
  ccis <- binomCI(cdf$value,nPriPerE);

  # convert to cm^2/keV
  norm <- pi*(rDisk1**2-rDisk0**2)*100^2/((e2-e1)*1000.0)/nPriPerE;

  x <- data.frame(outE=bdf$X1
                  ,inE=bdf$X2
                  ,ctsB=bdf$value
                  ,cBmin=bcis[,1]*nPriPerE
                  ,cBmax=bcis[,2]*nPriPerE
                  ,areaB=bdf$value*norm
                  ,aBmin=bcis[,1]*nPriPerE*norm
                  ,aBmax=bcis[,2]*nPriPerE*norm
                  ,ctsC=cdf$value
                  ,cCmin=ccis[,1]*nPriPerE
                  ,cCmax=ccis[,2]*nPriPerE
                  ,areaC=cdf$value*norm
                  ,aCmin=ccis[,1]*nPriPerE*norm
                  ,aCmax=ccis[,2]*nPriPerE*norm
                  ,outEBinLow=e1
                  ,outEBinHigh=e2);
  attr(x,"nPriPerE") <- nPriPerE;
  attr(x,"rDisk") <- sqrt(rDisk1**2+rDisk0**2)
  attr(x,"rDisk1") <- rDisk1;
  attr(x,"rDisk0") <- rDisk0;
  x;
}

lalf <- function(x){y <- x; y[y==0]<-1; y <- log10(y); y[x==0] <- min(y); y};
imageDRM <- function(a,drm){
  p <- ggplot(data=a) + geom_tile(aes(x=inE,y=outE,
                                      #alpha=areaB));
                                      alpha=lalf(areaB)));
                                      #fill=lalf(areaB)));
  #p <- p + scale_fill_gradientn(colours=jet.colors(7));
  p <- p + scale_alpha_continuous(name=expression(log[10](cm^2/keV)));
  p <- p + scale_x_log10()+scale_y_log10();
  p <- p + ylab("energy deposited (MeV)") + xlab("input energy (MeV)");
  p <- p + theme_bw();
  print(p);
}

lineDRM <- function(a,e,justGeom=FALSE,col="black",...){
  ine <- as.real(levels(factor(a$inE)));
  e <- ine[findInterval(e,ine)];
  a <- subset(a,a$inE == e);

  if(justGeom){
    list(geom_line(data=a,aes(x=outE,y=areaB),colour=col)
         ,geom_ribbon(data=a,aes(x=outE,ymin=aBmin,ymax=aBmax),alpha=0.3,fill=col));
  }else{
    #p <- ggplot(data=a) + geom_line(aes(x=outE,y=areaB));
    p <- ggplot(data=a);

    #print(a);

    #p <- p + geom_ribbon(aes(x=outE,ymin=cBmin,ymax=cBmax),alpha=0.3,fill=col);
    #p <- p + geom_line(aes(x=outE,y=ctsB));
    p <- p + geom_ribbon(aes(x=outE,ymin=aBmin,ymax=aBmax),alpha=0.3,fill=col);
    p <- p + geom_line(aes(x=outE,y=areaB),colour=col);
    p <- p + geom_point(data=subset(a,a$areaB>0),aes(x=outE,y=areaB),colour=col);

    #p <- p + xlim(1,e*1.2);
    p <- p + scale_x_log10(limits=c(0.01,e*1.2));
    p <- p + scale_y_log10(limits=range(c(a$areaB[a$areaB>0],a$aBmax)));
    #p <- p + ylim(0,200);
    
    p <- p + ylab(expression("effective area (cm"^2/"keV)"));
    p <- p + xlab("energy deposited (MeV)");

    rest <- list(...);
    if(length(rest)>0){
      for(xx in rest){
        print(xx);
        p <- p + xx;
      }
    }
    p;
  }
}

lineDRM_multi <- function(a,e){
  ine <- as.real(levels(factor(a$inE)));

  e <- ine[findInterval(e,ine)];

  a <- subset(a,a$inE %in% e);

  p <- ggplot(data=a) + geom_line(aes(x=outE,y=ctsB,group=inE,color=sprintf("%.2f MeV",inE)));
  #p <- p + scale_color_brewer();
  p <- p + scale_color_discrete(name=expression(E[pri]));
  p <- p + theme_bw();
  p <- p + scale_x_log10();

  print(p);
}

integrateDRM <- function(a,eIn,eOutMin,eOutMax){
  a <- subset(a,a$inE==eIn);

  #f <- approxfun(a$outEBinLow,a$areaB,method="constant",f=0,yleft=0,yright=0);
  #integrate(f,eOutMin,eOutMax,rel.tol=0.001,abs.tol=0.001,subdivisions=1000);

  e1 <- a$outEBinLow;
  e2 <- a$outEBinHigh;

  de <- (e2-e1)*1000; # keV
  integral <- sum(de[a$outEBinLow > eOutMin]*a$areaB[a$outEBinLow > eOutMin]);
  bdySel <- a$outEBinLow < eOutMin & a$outEBinHigh > eOutMin;
  bdyFrac <- (a$outEBinHigh-eOutMin)/de;
  integral <- integral + sum(bdySel*de*a$areaB*bdyFrac);
  integral;
}

integrateCrosSecs <- function(a){
  e1 <- a$outEBinLow;
  e2 <- a$outEBinHigh;
  de <- (e2-e1)*1000; # keV

  eIn <- as.real(levels(factor(a$inE)));
  anyDepSig <- mapply(function(e){integrateDRM(a,e,0,100)},eIn);
  fullEDepSig <- mapply(function(e){integrateDRM(a,e,e*0.9,100)},eIn);

  sig <- data.frame(eIn=eIn,sig=c(anyDepSig,fullEDepSig)
                    ,type=c(rep("total",length(anyDepSig)),rep("full-energy",length(fullEDepSig))));

  p <- ggplot(data=sig) + geom_line(aes(x=eIn,y=sig,group=type,color=type));
  p <- p + scale_x_log10()+scale_y_log10();
  p <- p + ylab(expression("integrated cross section ("*cm^2*")")) + xlab("input energy (MeV)");
  p <- p + theme_bw();
  print(p);
}

subsetEADF <- function(df,sel){
  x <- subset(df,sel);
  attr(x,"nPriPerE") = attr(df,"nPriPerE");
  attr(x,"rDisk") = attr(df,"rDisk");
  x;
}

totalEffArea <- function(df){
  nPriPerE <- attr(df,"nPriPerE");
  #print(nPriPerE);
  sumB <- sum(df$ctsB);
  sumC <- sum(df$ctsC);
  bci <- binomCI(sumB,nPriPerE);
  cci <- binomCI(sumC,nPriPerE);

  x <- df[1,];
  x$outEBinHigh <- df$outEBinHigh[length(df$outEBinHigh)];
  x$outE <- (x$outEBinLow+x$outEBinHigh)/2;
  x$ctsB <- sumB;
  x$cBmin <- bci[1]*nPriPerE;
  x$cBmax <- bci[2]*nPriPerE;
  x$ctsC <- sumC;
  x$cCmin <- cci[1]*nPriPerE;
  x$cCmax <- cci[2]*nPriPerE;
  norm <- pi*attr(df,"rDisk")**2*100^2/nPriPerE;
  de <- 1000*(x$outEBinHigh-x$outEBinLow);
  x$areaC <- x$ctsC*norm/de
  x$areaB <- x$ctsB*norm/de
  x$areaCcmsq <- x$ctsC*norm;
  x$areaBcmsq <- x$ctsB*norm;
  x$aBmin=bci[1]*nPriPerE*norm/de;
  x$aBmax=bci[2]*nPriPerE*norm/de;
  x$aCmin=cci[1]*nPriPerE*norm/de;
  x$aCmax=cci[2]*nPriPerE*norm/de;
  x;
}

combineDRMOutBins <- function(a,ncomb){
  a$merger <- floor(ncomb:(length(a$outEBinHigh)+ncomb-1)/ncomb);
  b <- ddply(a,.(merger),totalEffArea,.progress="text");
  attr(b,"nPriPerE") <- attr(a,"nPriPerE");
  attr(b,"rDisk") <- attr(a,"rDisk");
  b;
}

averageDRMs <- function(drms){
  d <- do.call(rbind,drms);
  totnPriPerE = sum(sapply(drms,function(x){attr(x,"nPriPerE")}));
  avgdf <- function(df){
               x <- df[1,];
               bci <- binomCI(sum(df$ctsB),totnPriPerE);
               cci <- binomCI(sum(df$ctsC),totnPriPerE);
               x$ctsB <- sum(df$ctsB);
               x$cBmin <- bci[1]*totnPriPerE;
               x$cBmax <- bci[2]*totnPriPerE;
               x$ctsC <- sum(df$ctsC);
               x$cCmin <- cci[1]*totnPriPerE;
               x$cCmax <- cci[2]*totnPriPerE;
               norm <- pi*attr(df,"rDisk")**2*100^2/(x$outEBinHigh-x$outEBinLow)/totnPriPerE;
               x$areaC <- x$ctsC*norm;
               x$areaB <- x$ctsB*norm;
               x$aBmin=bci[1]*totnPriPerE*norm;
               x$aBmax=bci[2]*totnPriPerE*norm;
               x$aCmin=cci[1]*totnPriPerE*norm;
               x$aCmax=cci[2]*totnPriPerE*norm;
               x}
  d <- ddply(d,.(inE,outE),avgdf,.progress="text");
  attr(d,"nPriPerE") <- totnPriPerE;
  attr(d,"rDisk") <- attr(drms[[1]],"rDisk");
  d;
}


interpolateDRM <- function(a,epri){
  ine <- as.real(levels(factor(a$inE)));
  i <- findInterval(e,ine)
  e1 <- ine[i];
  e2 <- ine[i+1];

  d1 <- subset(a,a$inE == e1);
  d2 <- subset(a,a$inE == e2);

}

drmCDF <- function(drm){
  norm <- pi*attr(drm,"rDisk")**2*100.0^2/attr(drm,"nPriPerE");
  data.frame(outE = drm$outE, inE=drm$inE, aB=cumsum(drm$ctsB)*norm, aC=cumsum(drm$ctsC)*norm);
}

bgoMat <- function(drm){
  drm <- drm[order(drm$inE,drm$outE),];
  outE <- as.real(levels(factor(drm$outE)));
  inE <- as.real(levels(factor(drm$inE)));
  matrix(drm$ctsB,nrow=length(outE));
}
