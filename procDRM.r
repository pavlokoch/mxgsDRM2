# This file includes code to process histograms from GEANT simulations into detector response matrices.
# 
# Basic usage from the R prompt:
# a <- readDRMs_df("../results/mxgsDRM_1/mats_22_500000_0.60_0.00_0.00_30.00_0.01_1e+02_41.txt",combineOutBins=2,nPriPerE=500000,rDisk1=0.6,rDisk0=0.0);
# f <- drmConvolver(a);
# bins <- 10**seq(-1,2,length.out=40);
# primarySpectrum <- function(e){1/e};
# drm <- makeDRM(f,primarySpectrum,bins,bins);
# image.plot(bins,bins,drm);

# libraries necessary
library(ggplot2); # plotting library
library(reshape); # utilities for rearranging vectors and matrices.
source("~/R/utils.r"); # color maps, multiplot function.

# Calculate a confidence interval for probability p in binomial given
# observation of x successes out of n trials, vectorized over x.  default
# confidence level gives 1 sigma error bar if normal approx holds.
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

# construct a matrix that when multiplied by another matrix, gives a matrix
# with row groups summed together.
sumRowGroupsMat <- function(nRows,nGrp){
  outer(seq(nRows/nGrp),seq(nRows),function(i,j){ifelse(j/nGrp-i<=0 & j/nGrp-i>-1,1,0)});
}

# read GEANT output.
# effective areas measured in cm^2/keV
# returns a data table with columns for input and output energies, counts and
# effective areas for BGO and CZT layers.  Pay the most attention to the ctsB
# variable, as it is the easiest to understand.  The attributes store relevant
# parameters for later calculation.
readDRMs_df <- function(fn,nPriPerE=1.0,rDisk1=1.0,rDisk0=0.0,combineOutBins=1){
  f <- file(fn,"rt");
  l <- readLines(f,5);
  close(f);
  priLine <- as.real(strsplit(l[3]," ")[[1]][-1:-3])
  outLine <- as.real(strsplit(l[4]," ")[[1]][-1:-5])

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

# add columns to DRM data table describing error bars.
# This calculation may take a long time and/or use up all the memory on the
# computer.  Use with caution.
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

# plot the DRM simulation results at the given energies.
# Rounds desired energies down to the nearest value that was simulated.
# geomOnly is useful for those familiar with ggplot for stacking multiple
# plots.
#
# example:
# a <- readDRMs_df(...);
# lineDRM(a,c(1,10,100);
lineDRM <- function(a,e,geomOnly=FALSE){
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

plotDRM <- function(inBins,outBins,drm){
  image.plot(outBins,inBins,drm,log='xy',xlab="deposited energy (MeV)",ylab="primary energy (MeV)",legend.lab="effective area (cm^2)",legend.mar=4);
}

# compares two DRM data tables.
# Assumes d1 and d2 are subsets of DRM data frames describing same primary energy.
# d1 is shown in black, d2 is shown in blue.
compareDRMs <- function(d1,d2){
  range <- c(0.01,1.1*max(d1$outE[d1$ctsB>0],d2$outE[d2$ctsB>0],na.rm=TRUE));
  p <- ggplot() + theme_bw();
  p <- p + geom_line(data=d1,aes(x=outE,y=ctsB),col='black');
  p <- p + geom_line(data=d2,aes(x=outE,y=ctsB),col='blue');
  p <- p + scale_x_log10(limits=range);
  p <- p + scale_y_log10();
  p <- p + xlab("Energy deposited (MeV)");
  p;
}
  
# take the subset of a DRM data table, preserving attributes.
subsetEADF <- function(df,sel){
  x <- subset(df,sel);
  attributes(x) <- attributes(df);
  x;
}

# interpolate a DRM simulation, separating the continuum (bg) from the spectral
# lines.
interpolateDRMbg <- function(df,e){
  inEs <- unique(df$inE);
  e <- inEs[findInterval(e,inEs)]; # round down to nearest simulated E.
  dfa <- attributes(df);
  df <- subset(df, df$inE == e); attributes(df) <- dfa;
  outE <- df$outE;
  outEb <- attr(df,"outEBins");
  cts <- df$ctsB;
  eBinMid <- outE[findInterval(e,outEb)];
  eBinMin <- outEb[findInterval(e,outEb)];

  me <- 0.510998903;
  moveLines <- c(e,e-me,e-2*me); # moving lines: full-energy, one-escape, and two-escape peaks.
  statLines <- c(me,2*me); # static lines: annihilation, and twice-annihilation.

  # combine lines, ignoring annihilation and escape if primary energy too low.
  lines <- if(e>2*me){c(moveLines,statLines);}else{c(e);}

  # indices in lines array of static and moving lines.
  mvidxs <- if(e>2*me){c(1,2,3)}else{c(1)};
  stidxs <- if(e>2*me){c(4,5)}else{c()};

  # drop bins containing lines from DRM simulation results
  toDrop <- findInterval(lines,outEb);
  outEd <- outE[-toDrop];
  ctsd <- cts[-toDrop];
  outEdd <- outE[toDrop];
  ctsdd <- cts[toDrop];

  # drop everything at or above the full-energy peak.
  ctsd <- ctsd[outEd<=eBinMid];
  outEd <- outEd[outEd<=eBinMid];

  # do weighted LOESS smoothing.  span is the fraction of the dataset to use,
  # and the control variable allows for extrapolation (not really used, but may
  # prevent NA from popping up.  Weights are calculated as in poisson error
  # bars, to prevent bias.
  f <- loess(ctsd~outEd,weights=1/(ctsd+1),span=11/length(ctsd),control=loess.control(surface="direct"));

  # function to evaluate loess smoothed continuum, stripping out NA's and negatives.
  bg <- function(oE){
    ans <- ifelse(oE>eBinMin,0,predict(f,oE));
    ans[ans<0 | is.na(ans)] <- 0;
    ans;
  }

  # subtract continuum contribution from spectral lines.
  peakCts <- ctsdd - bg(lines);

  list(e=e,bg=bg,sl=statLines,slc=peakCts[stidxs],ml=moveLines,mlc=peakCts[mvidxs],outE=outE,outEb=outEb);
}

# add spectral line contributions (lines, lcts) to a histogram containing the
# continuum (outEb, cts).
addLines <- function(outEb,lines,lcts,cts){
  ctr <- 1;
  for(ee in lines){
    i <- findInterval(ee,outEb);
    cts[i] <- cts[i] + lcts[ctr];
    ctr <- ctr + 1;
  }
  cts
}

# convert the results of DRM bg interpolation to a DRM data table.
interpDRMtoDF <- function(ntrp){
  bg <- ntrp$bg(ntrp$outE);
  ctsB <- addLines(ntrp$outEb,ntrp$sl,ntrp$slc,bg);
  ctsB <- addLines(ntrp$outEb,ntrp$ml,ntrp$mlc,ctsB);

  data.frame(outE=ntrp$outE,ctsB=ctsB,inE=ntrp$e);
}

# interpolation helper function, takes two bg interpolations and an energy, and
# makes a linear interpolation along lines radiating out from the origin.
interpolateH_bg <- function(i1,i2,e){
  f1 <- i1$bg;
  f2 <- i2$bg;
  oE <- i1$outE;

  pf1 <- (i2$e-e)/(i2$e-i1$e);
  pf2 <- 1-pf1;

  x <- oE/e;
  pf1*f1(x*i1$e)+pf2*f2(x*i2$e);
}

# given a DRM data table and two energies to use as interpolation base points,
# interpolate in between.
interpolateDRMs_givenE <- function(df,e1,e2,e){
  # do background/line interpolation at input energies.
  i1 <- interpolateDRMbg(df,e1);
  i2 <- interpolateDRMbg(df,e2);

  # fractional contributions
  pf1 <- (e2-e)/(e2-e1);
  pf2 <- 1-pf1;

  # do background interpolation between energies.
  cbg <- interpolateH_bg(i1,i2,e);

  outE <- i1$outE;
  outEb <- attr(df,"outEBins");
  
  # add stationary spectral lines to background.
  cbg <- addLines(outEb,i1$sl,i1$slc*pf1,cbg);
  cbg <- addLines(outEb,i2$sl,i2$slc*pf2,cbg);

  # interpolate moving lines to final position, intensity.
  n <- max(length(i1$ml),length(i2$ml));
  mlines <- pf1*c(i1$ml,rep(0,n-length(i1$ml)))+pf2*c(i2$ml,rep(0,n-length(i2$ml)));
  mlcts <- pf1*c(i1$mlc,rep(0,n-length(i1$mlc)))+pf2*c(i2$mlc,rep(0,n-length(i2$mlc)));

  # add moving lines to background.
  cbg <- addLines(outEb,mlines,mlcts,cbg)

  dfa <- attributes(df);
  dfa$names <- c("outE","inE","ctsB");
  dfa$row.names = seq_along(outE);
  d <- data.frame(outE=outE,
                  inE=e,
                  ctsB=cbg);
  attributes(d) <- dfa;
  d;
}

# wrapper for interpolateDRMs_givenE that determines energies automatically.
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

# construct a function that convolves a DRM with a given primary spectrum
# between two given primary energies.
drmConvolver <- function(df){
  inEs <- unique(df$inE);
  ntrps <- lapply(inEs,function(e){interpolateDRMbg(df,e)})
  outE <- ntrps[[1]]$outE;
  outEb <- ntrps[[1]]$outEb;
  bgs <- lapply(ntrps,function(ntrp){ntrp$bg(outE)});

  # interpolation of background and static lines
  linBSL <- function(e){
    j1 <- findInterval(e,inEs);
    j2 <- j1+1;
    bg <- interpolateH_bg(ntrps[[j1]],ntrps[[j2]],e); # background interp.
    lns <- ntrps[[j2]]$sl;
    if(length(ntrps[[j2]]$slc)>0){ # static line interp.
      p1 <- (inEs[j2]-e)/(inEs[j2]-inEs[j1]); p2 <- 1-p1;
      lcts <- p1*c(ntrps[[j1]]$slc,rep(0,2-length(ntrps[[j1]]$slc))) + p2*ntrps[[j2]]$slc;
      addLines(outEb,lns,lcts,bg);
    }else{
      bg;
    }
  }

  #background and static line convolution
  bslConv <- function(e1,e2,sf){
    oEs <- outE[outE>=e1 & outE<=e2];
    bwds <- diff(outEb)[outE>=e1 & outE<=e2];
    wts <- sf(oEs)*bwds;
    wts <- wts/sum(wts); 
    # sum over primary energies given by output bins between energy limits.
    bg <- Reduce(function(v,i){v + linBSL(oEs[i])*wts[i]}, seq_along(oEs),rep(0,length(outE)));
  }

  # moving line convolution
  mlConv <- function(e1,e2,sf){
    print(c(e1,e2));
    poEs <- outE[outE>=e1 & outE<=e2];
    pbwds <- diff(outEb)[outE>=e1 & outE<=e2];
    oEs <- outE;
    bwds <- diff(outEb);

    #sapply here simplifies down to a matrix, i.e. mlm[,1]= ntrps[[1]]$ml
    i1 <- findInterval(e1,inEs);
    i2 <- findInterval(e2,inEs)+1;
    mlm <- sapply(ntrps[i1:i2],function(xx)c(xx$ml,rep(0,3-length(xx$ml))));
    mlcm <- sapply(ntrps[i1:i2],function(xx)c(xx$mlc,rep(0,3-length(xx$mlc))));

    # make functions interpolating linearly in spectral line positions, cts.
    # note that there are always 3 moving lines, so nrow(mlm) = 3.
    fs <- lapply(1:nrow(mlm),function(i){approxfun(mlm[i,],mlcm[i,],yleft=0,yright=0)});

    # use functions to determine counts in each output bin, weighted by primary
    # spectrum and bin width.  Note that primary spectrum is evaluated at the
    # primary energy necessary to put the line in question at the output bin in
    # question.
    me <- 0.510998903;
    f1 <- fs[[1]](oEs)*(oEs>e1 & oEs<e2)*sf(oEs)*bwds;
    f2 <- fs[[2]](oEs)*(oEs+me>e1 & oEs+me<e2)*sf(oEs+me)*bwds;
    f3 <- fs[[3]](oEs)*(oEs+2*me>e1 & oEs+2*me<e2)*sf(oEs+2*me)*bwds;

    n1 <- sum(sf(oEs)*bwds*(oEs>e1 & oEs<e2));
    n2 <- sum(sf(oEs+me)*bwds*(oEs+me>e1 & oEs+me<e2));
    n3 <- sum(sf(oEs+2*me)*bwds*(oEs+2*me>e1 & oEs+2*me<e2));

    a1 <- if(n1>0){f1/n1}else{0};
    a2 <- if(n2>0){f2/n2}else{0};
    a3 <- if(n3>0){f3/n3}else{0};

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
