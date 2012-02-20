from math import *
import sys

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

def mxgsDRMCmd(name,pdgID,nPriPerE,rad,rad0,theta,phi,emin,emax,numE):
  outfn = "%s/%s/mats_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2g_%.2g_%d.txt"%(resultsDir,name,pdgID,nPriPerE,rad,rad0,theta,phi,emin,emax,numE)
  comment = outfn
  return "./mxgsDRM %d %d %f %f %f %f %g %g %d %g %g %d %s %s"%(pdgID,nPriPerE,rad,rad0,theta,phi,emin,emax,numE,outMinE,outMaxE,outNumE,outfn,comment)

def mpirunCmd(name,pdgID,nPriPerE,rad,rad0,(th0,th1,nth),(ph0,ph1,nph),(e0,e1,ne)):
  #  for th in linRange(th0,th1,nth) for ph in linRange(ph0,ph1,nph)])
  cmds = []
  ctr=1
  for (ph,th) in [(ph,th) for th in linRange(th0,th1,nth) for ph in linRange(ph0,ph1,nph)]:
    cmds.append(("-H `head -%d $PBS_NODEFILE | tail -1` -np 1 "%ctr)+mxgsDRMCmd(name,pdgID,nPriPerE,rad,rad0,th,ph,e0,e1,ne))
    ctr = ctr+1
    
  return "mpirun " + " : ".join(cmds)

def commands(name,pdgID,nPriPerE,rad,rad0,theta,phi,energy):
  return ["module unload pgi","module load gcc"
      ,"cd ~/geant/geant4.9.5-install/bin/"
      ,"source geant4.sh"
      ,"cd ~/sim-build"
      ,"mkdir -p %s/%s"%(resultsDir,name)
      ,mpirunCmd(name,pdgID,nPriPerE,rad,rad0,theta,phi,energy)]

def printPBS(outf,name,n,cmds,walltime):
  numNodes = int(n/2) + n%2
  outf.write("#! /bin/sh -\n")
  outf.write("#PBS -S /bin/sh\n")
  outf.write("#PBS -N \"%s\"\n"%name)
  outf.write("#PBS -A fysisk\n")
  outf.write("#PBS -l walltime=%s,nodes=%d:ppn=2\n"%(walltime,numNodes))
  outf.write("#PBS -l pmem=500mb\n")
  outf.write("#PBS -m abe\n");
  outf.write("#PBS -M brant.carlson@ift.uib.no\n");
  outf.write("#PBS -o %s.stdout.txt\n"%name)
  outf.write("#PBS -e %s.stderr.txt\n"%name)
  outf.write("")
  for cmd in cmds:
    outf.write(cmd)
    outf.write("\n")


def walltimeStr(nPriPerE,nPriE,rate=250):
  sTot = nPriPerE*nPriE/rate
  h = int(sTot/3600)
  m = int((sTot-h*3600)/60)
  s = int((sTot-h*3600-m*60))
  return "%02d:%02d:%02d"%(h,m,s)

name = "testJob"
ntheta = 1
nphi = 1
ne = 7
nPriPerE = 100000

#printPBS(sys.stdout,name,ntheta*nphi,commands(name,22,nPriPerE,0.6,0.0,(0,45,ntheta),(0,270,nphi),(0.1,100,ne)),walltimeStr(nPriPerE,ne,250));

name = sys.argv[1]
rad = float(sys.argv[2])
printPBS(sys.stdout,name,ntheta*nphi,commands(name,22,nPriPerE,rad+1.0,rad,(0,45,ntheta),(0,270,nphi),(0.1,100,ne)),walltimeStr(nPriPerE,ne,250));


