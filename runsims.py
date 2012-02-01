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

def runSet(prefix,eMeV,alt_m,pdgID):
  dbf = "%s.db"%prefix
  nPri = 100000
  nSim = 500000
  mri.makeGEANTInputDB2(dbf,nPri,pdgID,alt_m,eMeV)
  runGEANT(dbf,nSim)


#prefix = "test"
#nPriFile = 10000
#nPri = 10000

def rPS_h((eMeV,pdgID,alt_m)):
  name = "eff_%d_%.3g_%.2g"%(pdgID,eMeV,alt_m)
  runSet(name,eMeV,alt_m,pdgID)
  return 0

ba133 = 0.356
na22_1 = 0.511
na22_2 = 1.274
cs137 = 0.662
co60_1 = 1.173
co60_2 = 1.333
mn54 = 0.835
energies = [ba133,na22_1,na22_2,cs137,mn54,co60_1,co60_2]
#energies = [cs137]
#dists = [0.0,0.1,0.2,0.3,0.4,0.5]
dists = [0.05]
runs = [(e,22,d) for e in energies for d in dists]

print runs
if __name__ == '__main__':
  p = Pool(2)
  result = p.map(rPS_h,runs)
  print result

#rPS_h(runs[0])
