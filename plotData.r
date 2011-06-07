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

plotDRM <- function(dbfn,mult){
  db <- getDataConnection(dbfn);
  a <- db("select * from evtsBGO, pcles where evtsBGO.priIdx = pcles.idx");
  p <- db("select * from pcles");
  thp <- atan2(sqrt(p$px**2+p$py**2),-p$pz)*180/pi;
  th <- atan2(sqrt(a$px**2+a$py**2),-a$pz)*180/pi;
  hp <- hist(thp,breaks=seq(0,90,length.out=40));
  ha <- hist(th,breaks=seq(0,90,length.out=40));
  plot(ha$mids,ha$counts/hp$counts*100*100/mult,type='l',ylim=c(0,800),xlim=c(0,90),xlab=expression(theta*"  "*(degree)),ylab=expression(A[eff]*"  "*(cm^2)),main="BGO effective area, 3 MeV photons [PRELIMINARY]");
  list(a=a,p=p,th=th,thp=thp,ha=ha,hp=hp);
}


