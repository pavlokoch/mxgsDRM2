
# This file includes code to process histograms from GEANT simulations
# into detector response matrices.
# 
# WARNING: this file can be manually edited for use directly in R, but
# can also be automatically generated from the appendix in notes/docs.org.
# If you have emacs and org-mode, it is better to edit in emacs and use
# org-babel-tangle to generate the file.  That way the file and the document
# will stay in sync.
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
library(fields); # image.plot
#source("~/R/utils.r"); # color maps, multiplot function.

cbpr <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# from Robin Evans, http://cran.r-project.org/web/packages/rje/index.html
cubeHelix <- function (n, start = 0.5, r = -1.5, hue = 1, gamma = 1, reverse=FALSE) 
{
    M = matrix(c(-0.14861, -0.29227, 1.97294, 1.78277, -0.90649, 
        0), ncol = 2)
    lambda = seq(0, 1, length.out = n)
    l = rep(lambda^gamma, each = 3)
    phi = 2 * pi * (start/3 + r * lambda)
    t = rbind(cos(phi), sin(phi))
    out = l + hue * l * (1 - l)/2 * (M %*% t)
    out = pmin(pmax(out, 0), 1)
    out = apply(out, 2, function(x) rgb(x[1], x[2], x[3]))
    if(reverse){
        rev(out);
    }else{
        out;
    }
}


# Calculate a confidence interval for probability p in binomial given
# observation of x successes out of n trials, vectorized over x.  default
# confidence level gives 1 sigma error bar if normal approx holds.
# Used for adding error bars to a DRM.
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
# used for combining simulation output bins if desired (since they're really small by default).
sumRowGroupsMat <- function(nRows,nGrp){
  outer(seq(nRows/nGrp),seq(nRows),function(i,j){ifelse(j/nGrp-i<=0 & j/nGrp-i>-1,1,0)});
}

# read GEANT output.
# effective areas measured in cm^2/keV
# returns a data table with columns for input and output energies, counts and
# effective areas for BGO and CZT layers.  Pay the most attention to the ctsB
# variable, as it is the easiest to understand.  The attributes store relevant
# parameters for later calculation.
readDRMs_df <- function(fn,nPriPerE=1.0,rDisk1=0.6,rDisk0=0.0,combineOutBins=1,defaultDet='bgo'){
  f <- file(fn,"rt");
  l <- readLines(f,5);
  close(f);
  priLine <- as.double(strsplit(l[3]," ")[[1]][-1:-3])
  outLine <- as.double(strsplit(l[4]," ")[[1]][-1:-5])

  print("reading file, loading matrices...");
  a <- read.table(fn);
  bdf <- t(data.matrix(subset(a,a$V1=="BGO")[,-1]));
  cdf <- t(data.matrix(subset(a,a$V1=="CZT")[,-1]));

  combMat <- sumRowGroupsMat(dim(bdf)[1],combineOutBins);

  bdf <- combMat %*% bdf;
  cdf <- combMat %*% cdf;

  outLine <- outLine[seq(1,length(outLine),by=combineOutBins)];

  # centers and widths of output energy bins
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

  x$cts <- if(defaultDet=='bgo'){x$ctsB}else{x$ctsC}
  x$area <- if(defaultDet=='bgo'){x$areaB}else{x$areaC}
  
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
  ine <- as.double(levels(factor(a$inE)));

  e <- ine[findInterval(e,ine)];

  a <- subset(a,a$inE %in% e);
  a$estr <- sprintf("%.2f MeV",a$inE);
  a$estr <- ordered(factor(a$estr),levels=sprintf("%.2f MeV",sort(unique(a$inE))))

  g <- geom_line(data=a,aes(x=outE,y=cts,group=inE,color=estr));

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

#suggested additional argumens: log="xy", zlim=c(0,500)
plotDRM <- function(inBins,outBins,drm,...){
    image.plot(outBins,inBins,drm,xlab="deposited energy (MeV)",ylab="primary energy (MeV)",legend.lab="effective area (cm^2)",legend.mar=4,col=cubeHelix(1000,reverse=TRUE),...);
}

# makes plot in matrix sense, i.e. with origin at upper left instead of lower left.
plotDRMAsMat <- function(drm,...){
    image.plot(t(drm)[,ncol(drm):1],legend.lab="effective area (cm^2)",legend.mar=4,col=cubeHelix(1000,reverse=TRUE),...);
}

# compares two DRM data tables.
# Assumes d1 and d2 are subsets of DRM data frames describing same primary energy.
# d1 is shown in black, d2 is shown in blue.
compareDRMs <- function(d1,d2){
  p <- ggplot() + theme_bw();
  range <- c(0.01,1.1*max(d1$outE[d1$cts>0],d2$outE[d2$cts>0],na.rm=TRUE));
  p <- p + geom_line(data=d1,aes(x=outE,y=cts),col='black');
  p <- p + geom_line(data=d2,aes(x=outE,y=cts),col='blue');
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

# interpolate a DRM simulation over bins, separating the continuum (bg) from the spectral
# lines.  i.e. smooth a single DRM to limit noise.
# interpolation BETWEEN multiple DRMs is handled elsewhere.
interpolateDRMbg <- function(df,e){
  inEs <- unique(df$inE);
  e <- inEs[findInterval(e,inEs)]; # round down to nearest simulated E.
  dfa <- attributes(df);
  df <- subset(df, df$inE == e); attributes(df) <- dfa;
  outE <- df$outE;
  outEb <- attr(df,"outEBins");
  
  cts <- df$cts;
  
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

# convert the results of DRM bg interpolation/smoothing back to a DRM data table (histogram).
interpDRMtoDF <- function(ntrp){
  bg <- ntrp$bg(ntrp$outE);
  cts <- addLines(ntrp$outEb,ntrp$sl,ntrp$slc,bg);
  cts <- addLines(ntrp$outEb,ntrp$ml,ntrp$mlc,cts);

  data.frame(outE=ntrp$outE,cts=cts,inE=ntrp$e);
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
  dfa$names <- c("outE","inE","cts");
  dfa$row.names <- seq_along(outE);
  d <- data.frame(outE=outE, inE=e, cts=cbg)
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
  # returns a histogram of interpolated counts containing background + static line contributions
  linBSL <- function(e){
    j1 <- findInterval(e,inEs);
    j2 <- j1+1;
    bg <- interpolateH_bg(ntrps[[j1]],ntrps[[j2]],e); # background interp.
    lns <- ntrps[[j2]]$sl;
    if(length(ntrps[[j2]]$slc)>0){ # static line interp.
      p1 <- (inEs[j2]-e)/(inEs[j2]-inEs[j1]);
      p2 <- 1-p1;
      lcts <- p1*c(ntrps[[j1]]$slc,rep(0,2-length(ntrps[[j1]]$slc))) + p2*ntrps[[j2]]$slc;
      addLines(outEb,lns,lcts,bg);
    }else{
      bg;
    }
  }

  #background and static line convolution
  # returns a histogram summed from many results from linBSL
  # taken from many primary energies, weighted by the given $sf \propto dN/dE$ spectrum
  # grid in input energies is taken from output energies.  reasonable, but not required(?)
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
    #poEs <- outE[outE>=e1 & outE<=e2];  # unused?
    #pbwds <- diff(outEb)[outE>=e1 & outE<=e2];  # unused?
    oEs <- outE;
    bwds <- diff(outEb);

    # extract info from the interpolations for the range needed
    #sapply here simplifies down to a matrix, i.e. mlm[,1]= ntrps[[1]]$ml
    i1 <- findInterval(e1,inEs);
    i2 <- findInterval(e2,inEs)+1;
    mlm <- sapply(ntrps[i1:i2],function(xx)c(xx$ml,rep(0,3-length(xx$ml))));
    mlcm <- sapply(ntrps[i1:i2],function(xx)c(xx$mlc,rep(0,3-length(xx$mlc))));

    # make functions interpolating linearly in spectral line positions, cts.
    # note that there are always 3 moving lines (since the arrays above are padded),
    # so nrow(mlm) = 3.
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

# convolution with detector resolution
detResConvMat <- function(eOutb,relErrAt1MeV=0.05){
  el <- head(eOutb,-1);
  eu <- tail(eOutb,-1);
  em <- (el+eu)/2.0;
  w <- sqrt(em)*relErrAt1MeV;
  outer(seq_along(el),seq_along(el), function(i,j)(pnorm(eu[i],em[j],w[j])-pnorm(el[i],em[j],w[j])));
}

# like sumRowGroupsMat, but properly accounting for splitting of bins.
# used to decrease resolution of output bins to match desired outputs.
# assumes bins that split have uniform density, which isn't a terrible approximation
# since before splitting the bins are really small.
rebinDRMMat <- function(eInb,eOutb){
  eInbl <- head(eInb,-1);
  eInbu <- tail(eInb,-1);
  eOutbl <- head(eOutb,-1);
  eOutbu <- tail(eOutb,-1);
  outer(seq_along(eOutbl),seq_along(eInbl),
        function(i,j)pmax(0,(pmin(eInbu[j],eOutbu[i])-pmax(eInbl[j],eOutbl[i]))/(eInbu[j]-eInbl[j])));
}

makeDRM <- function(dCr,spect,inEBins,outEBins,relErrAt1MeV=0.05){
  elow <- head(inEBins,-1);
  ehigh <- tail(inEBins,-1);
  norm <- pi*(dCr$attrs$rDisk1**2-dCr$attrs$rDisk0**2)*100^2/dCr$attrs$nPriPerE
  mat <- mapply(dCr$f,elow,ehigh,MoreArgs=list(spect));
  rebinMat <- rebinDRMMat(dCr$outEb,outEBins);
  if(relErrAt1MeV>0){
      resConvMat <- detResConvMat(dCr$outEb,relErrAt1MeV);
      norm * rebinMat %*% resConvMat %*% mat;
  }else{
      norm * rebinMat %*% mat;
  }
}

applyDrmConvolution_makeDF <- function(dCr,e1,e2,spect){
  data.frame(cts=dCr$f(e1,e2,spect), outE=dCr$outE, inE=e1);
}

drmDiffPlot <- function(d1,d2){
  diff <- data.frame(outE=d1$outE,dc=d1$cts-d2$cts,type="delta(cts)");
  diffsig <- data.frame(outE=d1$outE,dc=(d1$cts-d2$cts)/sqrt(d1$cts+1),type="delta(cts)/sqrt(cts+1)");
  diff <- rbind(diff,diffsig);

  range <- c(0.01,1.1*max(d1$outE[d1$cts>0],d2$outE[d2$cts>0],na.rm=TRUE));
  #range[1] <- range[2]-0.1;
  p <- ggplot() + theme_bw();
  p <- p + geom_line(data=diff,aes(x=outE,y=dc),col='black');
  p <- p + scale_x_log10(limits=range);
  p <- p + facet_wrap(~ type, ncol=1,scale = "free_y");
  p <- p + xlab("Energy deposited (MeV)");
  p <- p + ylab("Difference in cts");
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

# for example, running R in the directory with the simulation results files:
# makeDRM_batch(fns=list.files("resultsDir",pattern="^mats_22_.*.txt",include.dirs=FALSE,full.name=TRUE),
#               nPriPerE=500000,
#               rDisk=0.6,
#               bins=10**seq(-1,2,length.out=40),
#               primarySpect=function(e){1/e})
makeDRM_batch <- function(fns,nPriPerE,rDisk1=0.6,rDisk0=0.0,bins,primarySpect,det='bgo'){
  for(fn in fns){
    print("*****************************");
    print(sprintf("working on %s...",fn));
    print("reading simulation results...");
    a <- readDRMs_df(fn,combineOutBins=2,nPriPerE=nPriPerE,rDisk1=rDisk1,rDisk0=rDisk0,defaultDet=det);
    print("making DRM convolver...");
    f <- drmConvolver(a);
    print("making DRM...");
    drm <- makeDRM(f,primarySpect,bins,bins);
    
    print("writing output...");
    dn = dirname(fn);
    fn = basename(fn);
    outDRMFn <- sprintf("%s/drm_%s",dn,fn);
    write.table(drm,outDRMFn,row.names=FALSE,col.names=FALSE)

    print("plotting...");
    plotFn <- sprintf("%s/drmPlot_logScale500Max_%s.pdf",dn,fn);
    pdf(plotFn,width=8,height=8);
    plotDRM(bins,bins,drm,log='xy',zlim=c(0,500),main=fn);
    dev.off();
    plotFn <- sprintf("%s/drmPlot_logScale_%s.pdf",dn,fn);
    pdf(plotFn,width=8,height=8);
    plotDRM(bins,bins,drm,log='xy',main=fn);
    dev.off();
    plotFn <- sprintf("%s/drmPlot_linScale_%s.pdf",dn,fn);
    pdf(plotFn,width=8,height=8);
    plotDRM(bins,bins,drm,main=fn);
    dev.off();
    
    print("done!");
  }
} 


# needs some additional setup to be run in the REPL...
# library(runParallel)
# cl <- makeCluster(4) # number of CPUs
# registerDoParallel(cl)
# then makeDRM_batch_parallel(...)
makeDRM_batch_parallel <- function(fns,nPriPerE,rDisk1=0.6,rDisk0=0.0,bins,primarySpect,det='bgo'){
    foreach(fn=fns) %dopar% {
        source("../procDRM.r");
        makeDRM_batch(c(fn),nPriPerE,rDisk1,rDisk0,bins,primarySpect,det)
    }
}

# note dropping 0-20 keV bin since it's not relevant to include all the way down to zero.
# note also division by 1000 to convert to MeV.
rhessiBins <- c(20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 100.0,
                120.0, 140.0, 160.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0,
                700.0, 800.0, 1000.0, 1200.0, 1400.0, 1700.0, 2000.0, 2500.0, 3000.0,
                3500.0, 4000.0, 5000.0, 6000.0, 7000.0, 9000.0, 11000.0, 13000.0,
                15000.0, 17500.0, 19999.0, 23000.0, 30000.0, 40000.0) / 1000.0;

logBins <- 10**seq(-2,2,length.out=41);
