import os
import math
import random
import sys
import mri
# import pp
from multiprocessing import Pool
import numpy as np
import sqlite3

random.seed()

def runGEANT(dbf,nPri):
  #cmd = "/home/brant/geant4/bin/Linux-g++/mxgsDRM 0 %s %d %d 2>&1 | ./logring %s %s 500"%(dbf,nPri,nPri,"deleteme.log.crap","deleteme.log.crap2")
  cmd = "/home/brant/geant4/bin/Linux-g++/mxgsDRM 0 %s %d %d"%(dbf,nPri,nPri)
  print cmd
  os.system(cmd)
  return 0

def runSet(prefix):
  dbf = "%s.db"%prefix
  nPri = 1000000
  mri.makeGEANTInputDB(dbf,nPri,1,2,22,0,0,(0,0,0),3)
  runGEANT(dbf,nPri)

prefix = "test"
#nPriFile = 10000
#nPri = 10000

def rPS_h((theta,phi,eMeV)):
  name = "%s_%.2f_%.2f_%.2g"%(prefix,theta,phi,eMeV)
  runSet(nPriFile,nPri,theta,phi,eMeV,name,maxAlt,satAlt,igrfYear)
  return 0

runs = [(lat,31.930,ang*math.pi/180.0,run) for (lat) in [-13.070] for run in [0,1,2] for ang in [5,10,15,20,25,30,40,50]]

#print runs
#if __name__ == '__main__':
#  p = Pool(6)
#  result = p.map(rPS_h,runs)
#  print result
##rPS_h((4.555,-82,0))

runSet("test")
