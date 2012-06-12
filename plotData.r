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
source("~/R/utils.r");

plotAngleDRMs <- function(plotDRMFunc,fullEnergy=FALSE,breaks=seq(90,180,length.out=20),...){
  fnames <- c("test_22_1e+02.db","test_30.db","test_0.3.db","test_0.1.db","test_1.db","test_10.db","test_3.db");
  ees <- c(100,30,0.3,0.1,1,10,3);
  add<-FALSE;
  for(i in seq_along(fnames)){
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
  rbind(measureEfficiency("efficiency/eff_22_0.36_0.db"), measureEfficiency("efficiency/eff_22_0.356_0.05.db"),
        measureEfficiency("efficiency/eff_22_0.36_0.1.db"), measureEfficiency("efficiency/eff_22_0.36_0.2.db"),
        measureEfficiency("efficiency/eff_22_0.36_0.3.db"), measureEfficiency("efficiency/eff_22_0.36_0.4.db"),
        measureEfficiency("efficiency/eff_22_0.36_0.5.db"),

        measureEfficiency("efficiency/eff_22_0.51_0.db"), measureEfficiency("efficiency/eff_22_0.511_0.05.db"),
        measureEfficiency("efficiency/eff_22_0.51_0.1.db"), measureEfficiency("efficiency/eff_22_0.51_0.2.db"),
        measureEfficiency("efficiency/eff_22_0.51_0.3.db"), measureEfficiency("efficiency/eff_22_0.51_0.4.db"),
        measureEfficiency("efficiency/eff_22_0.51_0.5.db"),

        measureEfficiency("efficiency/eff_22_0.66_0.db"), measureEfficiency("efficiency/eff_22_0.662_0.05.db"),
        measureEfficiency("efficiency/eff_22_0.66_0.1.db"), measureEfficiency("efficiency/eff_22_0.66_0.2.db"),
        measureEfficiency("efficiency/eff_22_0.66_0.3.db"), measureEfficiency("efficiency/eff_22_0.66_0.4.db"),
        measureEfficiency("efficiency/eff_22_0.66_0.5.db"),

        measureEfficiency("efficiency/eff_22_0.83_0.db"), measureEfficiency("efficiency/eff_22_0.835_0.05.db"),
        measureEfficiency("efficiency/eff_22_0.83_0.1.db"), measureEfficiency("efficiency/eff_22_0.83_0.2.db"),
        measureEfficiency("efficiency/eff_22_0.83_0.3.db"), measureEfficiency("efficiency/eff_22_0.83_0.4.db"),
        measureEfficiency("efficiency/eff_22_0.83_0.5.db"),

        measureEfficiency("efficiency/eff_22_1.2_0.db"), measureEfficiency("efficiency/eff_22_1.17_0.05.db"),
        measureEfficiency("efficiency/eff_22_1.2_0.1.db"), measureEfficiency("efficiency/eff_22_1.2_0.2.db"),
        measureEfficiency("efficiency/eff_22_1.2_0.3.db"), measureEfficiency("efficiency/eff_22_1.2_0.4.db"),
        measureEfficiency("efficiency/eff_22_1.2_0.5.db"),

        measureEfficiency("efficiency/eff_22_1.27_0.db"), measureEfficiency("efficiency/eff_22_1.27_0.05.db"),
        measureEfficiency("efficiency/eff_22_1.27_0.1.db"), measureEfficiency("efficiency/eff_22_1.27_0.2.db"),
        measureEfficiency("efficiency/eff_22_1.27_0.3.db"), measureEfficiency("efficiency/eff_22_1.27_0.4.db"),
        measureEfficiency("efficiency/eff_22_1.27_0.5.db"),

        measureEfficiency("efficiency/eff_22_1.3_0.db"), measureEfficiency("efficiency/eff_22_1.33_0.05.db"),
        measureEfficiency("efficiency/eff_22_1.3_0.1.db"), measureEfficiency("efficiency/eff_22_1.3_0.2.db"),
        measureEfficiency("efficiency/eff_22_1.3_0.3.db"), measureEfficiency("efficiency/eff_22_1.3_0.4.db"),
        measureEfficiency("efficiency/eff_22_1.3_0.5.db"));
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

#vectorized over x.
# default confidence level gives 1 sigma error bar if normal approx holds.
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

sumRowGroupsMat <- function(nRows,nGrp){
  outer(seq(nRows/nGrp),seq(nRows),function(i,j){ifelse(j/nGrp-i<=0 & j/nGrp-i>-1,1,0)});
}

# effective areas measured in cm^2/keV
readDRMs_df <- function(fn,nPriPerE=1.0,rDisk1=1.0,rDisk0=0.0,combineOutBins=1){
  f <- file(fn,"rt");
  l <- readLines(f,5);
  close(f);
  priLine <- as.real(strsplit(l[3]," ")[[1]][-1:-3])
  outLine <- as.real(strsplit(l[4]," ")[[1]][-1:-5])
  #print(l[1]);
  #print(l[2]);

  print("reading file, loading matrices...");
  a <- read.table(fn);
  bdf <- t(data.matrix(subset(a,a$V1=="BGO")[,-1]));
  cdf <- t(data.matrix(subset(a,a$V1=="CZT")[,-1]));

  combMat <- sumRowGroupsMat(dim(bdf)[1],combineOutBins);

  bdf <- combMat %*% bdf;
  cdf <- combMat %*% cdf;

  outLine <- outLine[seq(1,length(outLine),by=combineOutBins)];

  om <- (outLine[-1] + outLine[-length(outLine)])/2;
  deo <- (outLine[-1] - outLine[-length(outLine)]);

  bdf <- data.frame(melt(matrix(bdf,nrow=length(om),dimnames=list(om,priLine))));
  cdf <- data.frame(melt(matrix(cdf,nrow=length(om),dimnames=list(om,priLine))));

  e1 <- outLine[findInterval(bdf$X1,outLine)];
  e2 <- outLine[findInterval(bdf$X1,outLine)+1];

  # convert to cm^2/keV
  norm <- pi*(rDisk1**2-rDisk0**2)*100^2/((e2-e1)*1000.0)/nPriPerE;

  x <- data.frame(outE=bdf$X1
                  ,inE=bdf$X2
                  ,ctsB=bdf$value
                  ,areaB=bdf$value*norm
                  ,ctsC=cdf$value
                  ,areaC=cdf$value*norm
                  ,outEBinLow=e1
                  ,outEBinHigh=e2);
  attr(x,"nPriPerE") <- nPriPerE;
  attr(x,"rDisk") <- sqrt(rDisk1**2+rDisk0**2)
  attr(x,"rDisk1") <- rDisk1;
  attr(x,"rDisk0") <- rDisk0;
  attr(x,"norm") <- norm;
  attr(x,"outEBins") <- outLine;
  x;
}

addErrorBars_df <- function(df){
  nPriPerE <- attr(df,"nPriPerE");
  norm <- attr(df,"norm");

  print("calculating BGO error bars...");
  bcis <- binomCI(df$ctsB,nPriPerE);
  print("calculating CZT error bars...");
  ccis <- binomCI(df$ctsC,nPriPerE);
  print("done");

  df$cBmin <- bcis[,1]*nPriPerE;
  df$cBmax <- bcis[,2]*nPriPerE;
  df$aBmin <- bcis[,1]*nPriPerE*norm;
  df$aBmax <- bcis[,2]*nPriPerE*norm;

  df$cCmin <- ccis[,1]*nPriPerE;
  df$cCmax <- ccis[,2]*nPriPerE;
  df$aCmin <- ccis[,1]*nPriPerE*norm;
  df$aCmax <- ccis[,2]*nPriPerE*norm;

  df;
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
    p <- ggplot(data=a) + theme_bw();

    #print(a);

    #p <- p + geom_ribbon(aes(x=outE,ymin=cBmin,ymax=cBmax),alpha=0.3,fill=col);
    #p <- p + geom_line(aes(x=outE,y=ctsB));
    #p <- p + geom_ribbon(aes(x=outE,ymin=aBmin,ymax=aBmax),alpha=0.3,fill=col);
    p <- p + geom_line(aes(x=outE,y=areaB),colour=col);
    p <- p + geom_point(data=subset(a,a$areaB>0),aes(x=outE,y=areaB),colour=col);

    #p <- p + xlim(1,e*1.2);
    p <- p + scale_x_log10(limits=c(0.01,e*1.2));
    p <- p + scale_y_log10();
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

lineDRM_multi <- function(a,e,geomOnly=FALSE,col='black'){
  ine <- as.real(levels(factor(a$inE)));

  e <- ine[findInterval(e,ine)];

  a <- subset(a,a$inE %in% e);
  a$estr <- sprintf("%.2f MeV",a$inE);
  a$estr <- ordered(factor(a$estr),levels=sprintf("%.2f MeV",sort(unique(a$inE))))

  g <- geom_line(data=a,aes(x=outE,y=ctsB,group=inE,color=estr));

  p <- ggplot();
  #p <- p + scale_color_brewer();
  p <- p + scale_color_manual(values=cbpr,name=expression(E[pri]));
  p <- p + theme_bw();
  p <- p + scale_x_log10();
  #p <- p + xlim(limits=c(0,5));
  p <- p + scale_y_log10();

  p <- p + xlab("Energy deposited (MeV)");

  p <- p + ylab("counts per bin");

  if(geomOnly){
    g;
  }else{
    p+g;
  }
}

# assumes d1 and d2 are subsets of DRM data frames describing same primary energy.
compareDRMs <- function(d1,d2){
  range <- c(0.01,1.1*max(d1$outE[d1$ctsB>0],d2$outE[d2$ctsB>0],na.rm=TRUE));
  #range[1] <- range[2]-0.1;
  p <- ggplot() + theme_bw();
  p <- p + geom_line(data=d1,aes(x=outE,y=ctsB),col='black');
  p <- p + geom_line(data=d2,aes(x=outE,y=ctsB),col='blue');
  p <- p + scale_x_log10(limits=range);
  p <- p + scale_y_log10();
  p <- p + xlab("Energy deposited (MeV)");
  p;
}

testInterpBG <- function(df,e){
  ntrp <- interpolateDRMbg(df,e);
  e <- ntrp$e;
  dfs <- subset(df,df$inE==e);
  p1 <- compareDRMs(dfs,interpDRMtoDF(ntrp));
  p2 <- drmDiffPlot(dfs,interpDRMtoDF(ntrp));
  multiplot(p1,p2,cols=1);
  print(e);
}

testInterp <- function(df,e){
  eIn <- unique(df$inE);
  ei <- findInterval(e,eIn);
  e1 <- eIn[ei-1];
  e2 <- eIn[ei+1];
  e <- eIn[ei];
  print(c(e1,e2,e));
  ntrp <- interpolateDRMs_givenE(df,e1,e2,e);
  dfs <- subset(df,df$inE==e);
  p1 <- compareDRMs(dfs,ntrp);
  p2 <- drmDiffPlot(dfs,ntrp);
  multiplot(p1,p2,cols=1);
  print(e);
}

drmDiffPlot <- function(d1,d2){
  diff <- data.frame(outE=d1$outE,dc=d1$ctsB-d2$ctsB,type="delta(ctsB)");
  diffsig <- data.frame(outE=d1$outE,dc=(d1$ctsB-d2$ctsB)/sqrt(d1$ctsB+1),type="delta(ctsB)/sqrt(ctsB+1)");
  diff <- rbind(diff,diffsig);

  range <- c(0.01,1.1*max(d1$outE[d1$ctsB>0],d2$outE[d2$ctsB>0],na.rm=TRUE));
  #range[1] <- range[2]-0.1;
  p <- ggplot() + theme_bw();
  p <- p + geom_line(data=diff,aes(x=outE,y=dc),col='black');
  p <- p + scale_x_log10(limits=range);
  p <- p + facet_wrap(~ type, ncol=1,scale = "free_y");
  p <- p + xlab("Energy deposited (MeV)");
  p <- p + ylab("Difference in ctsB");
  p;
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
  #bci <- binomCI(sumB,nPriPerE);
  #cci <- binomCI(sumC,nPriPerE);

  x <- df[1,];
  x$outEBinHigh <- df$outEBinHigh[length(df$outEBinHigh)];
  x$outE <- (x$outEBinLow+x$outEBinHigh)/2;
  x$ctsB <- sumB;
  #x$cBmin <- bci[1]*nPriPerE;
  #x$cBmax <- bci[2]*nPriPerE;
  x$ctsC <- sumC;
  #x$cCmin <- cci[1]*nPriPerE;
  #x$cCmax <- cci[2]*nPriPerE;
  norm <- pi*attr(df,"rDisk")**2*100^2/nPriPerE;
  de <- 1000*(x$outEBinHigh-x$outEBinLow);
  x$areaC <- x$ctsC*norm/de
  x$areaB <- x$ctsB*norm/de
  x$areaCcmsq <- x$ctsC*norm;
  x$areaBcmsq <- x$ctsB*norm;
  #x$aBmin=bci[1]*nPriPerE*norm/de;
  #x$aBmax=bci[2]*nPriPerE*norm/de;
  #x$aCmin=cci[1]*nPriPerE*norm/de;
  #x$aCmax=cci[2]*nPriPerE*norm/de;
  x;
}

combineDRMs <- function(dfs){

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

interpolateDRMbg <- function(df,e){
  inEs <- unique(df$inE);
  e <- inEs[findInterval(e,inEs)];
  dfa <- attributes(df);
  df <- subset(df, df$inE == e); attributes(df) <- dfa;
  outE <- df$outE;
  outEb <- attr(df,"outEBins");
  cts <- df$ctsB;
  eBinMid <- outE[findInterval(e,outEb)];
  eBinMin <- outEb[findInterval(e,outEb)];

  me <- 0.510998903;
  moveLines <- c(e,e-me,e-2*me);
  statLines <- c(me,2*me);

  lines <- if(e>2*me){c(moveLines,statLines);}else{c(e);}

  mvidxs <- if(e>2*me){c(1,2,3)}else{c(1)};
  stidxs <- if(e>2*me){c(4,5)}else{c()};

  toDrop <- findInterval(lines,outEb);
  outEd <- outE[-toDrop];
  ctsd <- cts[-toDrop];
  outEdd <- outE[toDrop];
  ctsdd <- cts[toDrop];

  ctsd <- ctsd[outEd<=eBinMid];
  outEd <- outEd[outEd<=eBinMid];

  lEmax <- max(outEd);
  f <- loess(ctsd~outEd,weights=1/(ctsd+1),span=11/length(ctsd),control=loess.control(surface="direct"));
  #f <- loess(ctsd~outEd,span=11/length(ctsd),control=loess.control(surface="direct"));


  bg <- function(oE){
    ans <- ifelse(oE>eBinMin,0,predict(f,oE));
    ans[ans<0 | is.na(ans)] <- 0;
    ans;
  }
  peakCts <- ctsdd - bg(lines);

  list(e=e,bg=bg,sl=statLines,slc=peakCts[stidxs],ml=moveLines,mlc=peakCts[mvidxs],outE=outE,outEb=outEb);
}

addLines <- function(outEb,lines,lcts,cts){
  ctr <- 1;
  for(ee in lines){
    i <- findInterval(ee,outEb);
    cts[i] <- cts[i] + lcts[ctr];
    ctr <- ctr + 1;
  }
  cts
}

interpDRMtoDF <- function(ntrp){
  bg <- ntrp$bg(ntrp$outE);
  ctsB <- addLines(ntrp$outEb,ntrp$sl,ntrp$slc,bg);
  ctsB <- addLines(ntrp$outEb,ntrp$ml,ntrp$mlc,ctsB);

  data.frame(outE=ntrp$outE,ctsB=ctsB,inE=ntrp$e);
}

interpolateH_bg <- function(i1,i2,e){
  f1 <- i1$bg;
  f2 <- i2$bg;
  oE <- i1$outE;

  pf1 <- (i2$e-e)/(i2$e-i1$e);
  pf2 <- 1-pf1;

  x <- oE/e;
  pf1*f1(x*i1$e)+pf2*f2(x*i2$e);
}

interpolateDRMs_givenE <- function(df,e1,e2,e){
  i1 <- interpolateDRMbg(df,e1);
  i2 <- interpolateDRMbg(df,e2);

  pf1 <- (e2-e)/(e2-e1);
  pf2 <- 1-pf1;

  cbg <- interpolateH_bg(i1,i2,e);

  outE <- i1$outE;

  outEb <- attr(df,"outEBins");
  
  cbg <- addLines(outEb,i1$sl,i1$slc*pf1,cbg);
  cbg <- addLines(outEb,i2$sl,i2$slc*pf2,cbg);

  n <- max(length(i1$ml),length(i2$ml));
  mlines <- pf1*c(i1$ml,rep(0,n-length(i1$ml)))+pf2*c(i2$ml,rep(0,n-length(i2$ml)));
  mlcts <- pf1*c(i1$mlc,rep(0,n-length(i1$mlc)))+pf2*c(i2$mlc,rep(0,n-length(i2$mlc)));

  cbg <- addLines(outEb,mlines,mlcts,cbg)

  dfa <- attributes(df);
  d <- data.frame(outE=outE,
                  inE=e,
                  ctsB=cbg);
  #attributes(d) <- dfa;
  d;
}

interpolateDRMs <- function(df,e){
  es <- unique(df$inE);
  if(e<min(es) || e>max(es)){
    print("out of range!");
    return(0);
  }
  e1 <- es[findInterval(e,es)]
  e2 <- es[findInterval(e,es)+1]
  interpolateDRMs_givenE(df,e1,e2,e);
}

drmConvolver <- function(df){
  inEs <- unique(df$inE);
  ntrps <- lapply(inEs,function(e){interpolateDRMbg(df,e)})
  outE <- ntrps[[1]]$outE;
  outEb <- ntrps[[1]]$outEb;
  bgs <- lapply(ntrps,function(ntrp){ntrp$bg(outE)});

  linBSL <- function(e){
    j1 <- findInterval(e,inEs);
    j2 <- j1+1;
    bg <- interpolateH_bg(ntrps[[j1]],ntrps[[j2]],e);
    lns <- ntrps[[j2]]$sl;
    if(length(ntrps[[j2]]$slc)>0){
      p1 <- (inEs[j2]-e)/(inEs[j2]-inEs[j1]); p2 <- 1-p1;
      lcts <- p1*c(ntrps[[j1]]$slc,rep(0,2-length(ntrps[[j1]]$slc))) + p2*ntrps[[j2]]$slc;
      addLines(outEb,lns,lcts,bg);
    }else{
      bg;
    }
  }

  bslConv <- function(e1,e2,sf){ #background and static line convolution
    oEs <- outE[outE>=e1 & outE<=e2];
    bwds <- diff(outEb)[outE>=e1 & outE<=e2];
    wts <- sf(oEs)*bwds;
    wts <- wts/sum(wts); 
    bg <- Reduce(function(v,i){v + linBSL(oEs[i])*wts[i]}, seq_along(oEs),rep(0,length(outE)));
  }

  mlConv <- function(e1,e2,sf){ #spectral line contribution to convolution
    poEs <- outE[outE>=e1 & outE<=e2];
    pbwds <- diff(outEb)[outE>=e1 & outE<=e2];
    oEs <- outE;
    bwds <- diff(outEb);
    wts <- sf(oEs);
    wts <- wts/sum(wts); 

    #sapply here simplifies down to a matrix, i.e. mlm[,1]= ntrps[[1]]$ml
    i1 <- findInterval(e1,inEs);
    i2 <- findInterval(e2,inEs)+1;
    mlm <- sapply(ntrps[i1:i2],function(xx)c(xx$ml,rep(0,3-length(xx$ml))));
    mlcm <- sapply(ntrps[i1:i2],function(xx)c(xx$mlc,rep(0,3-length(xx$mlc))));

    fs <- lapply(1:nrow(mlm),function(i){approxfun(mlm[i,],mlcm[i,],yleft=0,yright=0)});

    me <- 0.510998903;
    f1 <- fs[[1]](oEs)*(oEs>e1 & oEs<e2)*sf(oEs);
    f2 <- fs[[2]](oEs)*(oEs+me>e1 & oEs+me<e2)*sf(oEs+me);
    f3 <- fs[[3]](oEs)*(oEs+2*me>e1 & oEs+2*me<e2)*sf(oEs+2*me);

    n1 <- sum(sf(oEs)*(oEs>e1 & oEs<e2));
    n2 <- sum(sf(oEs+me)*(oEs+me>e1 & oEs+me<e2));
    n3 <- sum(sf(oEs+2*me)*(oEs+2*me>e1 & oEs+2*me<e2));

    print(c(n1,n2,n3));
    print(c(max(f1),max(f2),max(f3)));

    
    a1 <- if(n1>0){f1/n1}else{0};
    a2 <- if(n2>0){f2/n2}else{0};
    a3 <- if(n3>0){f3/n3}else{0};

    # TODO: doublecheck normalization.
    a1+a2+a3;
  }

  list(f=function(e1,e2,sf){ bslConv(e1,e2,sf)+mlConv(e1,e2,sf); },
       fb=function(e1,e2,sf)bslConv(e1,e2,sf),
       fm=function(e1,e2,sf)mlConv(e1,e2,sf),
       inE = inEs,
       outE = outE,
       outEb = outEb,
       attrs = attributes(df));
}

rebinDRMMat <- function(eInb,eOutb){
  eInbl <- head(eInb,-1);
  eInbu <- tail(eInb,-1);
  eOutbl <- head(eOutb,-1);
  eOutbu <- tail(eOutb,-1);
  outer(seq_along(eOutbl),seq_along(eInbl),
        function(i,j)pmax(0,(pmin(eInbu[j],eOutbu[i])-pmax(eInbl[j],eOutbl[i]))/(eInbu[j]-eInbl[j])));
}

makeDRM <- function(dCr,spect,inEBins,outEBins){
  elow <- head(inEBins,-1);
  ehigh <- tail(inEBins,-1);
  norm <- pi*(dCr$attrs$rDisk1**2-dCr$attrs$rDisk0**2)*100^2/dCr$attrs$nPriPerE
  mat <- mapply(dCr$f,elow,ehigh,MoreArgs=list(spect));
  rebinMat <- rebinDRMMat(dCr$outEb,outEBins);
  norm * rebinMat %*% mat;
}

applyDrmConvolution_makeDF <- function(dCr,e1,e2,spect){
  data.frame(ctsB=dCr$f(e1,e2,spect),
             outE=dCr$outE,
             inE=e1);
}

