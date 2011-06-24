#this is a comment in the first 50 lines for Vim's filetype detection.
library(rgl);
library(gplots);
library(RSQLite);
library(akima);
library(plotrix);
library(fields);
library(ks);
library(boot);
library(Hmisc);
library(splancs)
source("~/utils/colormaps.r");


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
  #a <- db("select * from evtsBGO, pcles where evtsBGO.priIdx = pcles.idx and evtsBGO.totEDep>0.01");
  a <- db("select * from evtsCZT, pcles where evtsCZT.priIdx = pcles.idx and evtsCZT.totEDep>0.01");
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
