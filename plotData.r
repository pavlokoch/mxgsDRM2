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

plotDRM <- function(r=c(),dbfn="",mult=0){
  if(length(r)==0){
    r <- list();
    db <- getDataConnection(dbfn);
    #r$a <- db("select * from evtsBGO, pcles where evtsBGO.priIdx = pcles.idx");
    r$a <- db("select * from evtsCZT, pcles where evtsCZT.priIdx = pcles.idx");
    r$p <- db("select * from pcles");
  }
  print(names(r));
  r$thp <- atan2(sqrt(r$p$px**2+r$p$py**2),-r$p$pz)*180/pi;
  r$th <- atan2(sqrt(r$a$px**2+r$a$py**2),-r$a$pz)*180/pi;
  nbins <- 20;
  r$hp <- hist(r$thp,breaks=seq(0,90,length.out=nbins));
  r$ha <- hist(r$th,breaks=seq(0,90,length.out=nbins));
  print(r$ha);
  plot(r$ha$mids,r$ha$counts/r$hp$counts*100*100/mult,type='l',xlab=expression(theta*"  "*(degree)),ylab=expression(A[eff]*"  "*(cm^2)),main="effective area, 3 MeV photons [PRELIMINARY]");
  list(a=r$a,p=r$p,th=r$th,thp=r$thp,ha=r$ha,hp=r$hp);
}
