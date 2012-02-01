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

  a <- read.table(fn);
  bdf <- data.matrix(subset(a,a$V1=="BGO")[,-1]);
  cdf <- data.matrix(subset(a,a$V1=="CZT")[,-1]);
  list(b=bdf,c=cdf,i=priLine,om=(outLine[-1] + outLine[-length(outLine)])/2,obks=outLine);
}
