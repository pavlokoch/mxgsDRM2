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

def runGEANT(dbf,nSim):
  #cmd = "/home/brant/geant4/bin/Linux-g++/mxgsDRM 0 %s %d %d 2>&1 | ./logring %s %s 500"%(dbf,nPri,nPri,"deleteme.log.crap","deleteme.log.crap2")
  cmd = "/home/brant/geant4/bin/Linux-g++/mxgsDRM 0 %s %d %d | ./logring %s %s 500"%(dbf,nSim,nSim,dbf+"_log.txt",dbf+"_log2.txt")
  print cmd
  os.system(cmd)
  return 0

#x="0" y="2829.7 + 145.8/2 + 863/2" z="-1399.6 + 1130/2 - 863 + 550/2"
mxgs_orig_in_columbus = (0, 2.8297 + 0.1458/2.0 + 0.8630/2.0, -1.3996 + 1.1300/2.0 - 0.8630 + 0.5500/2.0)
mxgs_orig_in_asim = (0, 0.200, -0.125)

def runSet(prefix,eMeV):
  dbf = "%s.db"%prefix
  nPri = 1000000
  nSim = 10000000
  mri.makeGEANTInputDB(dbf,nPri,1,20,22,0,0,mxgs_orig_in_columbus,eMeV)
  runGEANT(dbf,nSim)


prefix = "test"
#nPriFile = 10000
#nPri = 10000

def rPS_h(eMeV):
  name = "/mnt/sdb1/mxgsDRM/test_%.2g"%(eMeV)
  runSet(name,eMeV)
  return 0

runs = [0.001,0.003,0.01,0.03,0.1,0.3,1,3]
#runs = [1]

print runs
if __name__ == '__main__':
  p = Pool(2)
  result = p.map(rPS_h,runs)
  print result

#rPS_h(0.1)
