from math import *

workingDir = "~/sim-build";
resultsDir = "~/results";

outMinE = 0.001
outMaxE = 110.0
outNumE = 5000

def linRange(min,max,length):
  if(length<=1):
    return [min]
  else:
    return [min + float(n)/(length-1)*(max-min) for n in range(length)]
def logRange(min,max,length):
  if(length<=1):
    return [min]
  else:
    return [10**(log10(min) + float(n)/(length-1)*(log10(max)-log10(min))) for n in range(length)]

def mxgsDRMCmd(pdgID,nPriPerE,rad,theta,phi,emin,emax,numE):
  outfn = "%s/mats_%d_%d_%.2f_%.2f_%.2f_%.2g_%.2g_%d.txt"%(resultsDir,pdgID,nPriPerE,rad,theta,phi,emin,emax,numE)
  comment = outfn
  return "./mxgsDRM %d %d %f %f %f %g %g %d %g %g %d %s %s"%(pdgID,nPriPerE,rad,theta,phi,emin,emax,numE,outMinE,outMaxE,outNumE,outfn,comment)

def mpirunCmd(pdgID,nPriPerE,rad,(th0,th1,nth),(ph0,ph1,nph),(e0,e1,ne)):
  return "mpirun " + " : ".join(["-np 1 "+mxgsDRMCmd(pdgID,nPriPerE,rad,th,ph,e0,e1,ne) 
    for th in linRange(th0,th1,nth) for ph in linRange(ph0,ph1,nph)])

def commands(name,pdgID,nPriPerE,rad,theta,phi,energy):
  return ["module unload pgi","module load gcc"
      ,"cd ~/geant/geant4.9.5-install/bin/"
      ,"source geant4.sh"
      ,"cd ~/sim-build"
      ,"mkdir -p %s/%s"%(resultsDir,name)
      ,mpirunCmd(pdgID,nPriPerE,rad,theta,phi,energy)]

def printPBS(name,n,cmds,walltime):
  numNodes = int(n/2) + n%2
  print("#! /bin/sh -")
  print("#PBS -S /bin/sh")
  print("#PBS -N \"%s\""%name)
  print("#PBS -A fysisk")
  print("#PBS -l walltime=%s,nodes=%d:ppn=2,walltime=%s"%(wallitme,numNodes))
  print("#PBS -l pmem=500mb")
  print("#PBS -m abe");
  print("#PBS -M brant.carlson@ift.uib.no");
  print("#PBS -o %s.stdout.txt"%name)
  print("#PBS -e %s.stderr.txt"%name)
  print("")
  for cmd in cmds:
    print(cmd)

name = "testJob"
printPBS(name,4,commands(name,22,100,0.5,(0,1,2),(1,2,2),(0.1,1,10)),"00:10:00");
