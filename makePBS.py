from math import *
import sys

workingDir = "~/sim-build"
resultsDir = "~/results"

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

def blocks(lst,n):
  if(len(lst) % n != 0):
    print("WARNING: grid didn't divide evenly")
  d = len(lst)/n;
  return [lst[i:(i+d)] for i in range(0,len(lst),d)]

def mxgsDRMCmd(name,pdgID,nPriPerE,rad,rad0,theta,phi,emin,emax,numE):
  outfn = "%s/%s/mats_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2g_%.2g_%d.txt"%(resultsDir,name,pdgID,nPriPerE,rad,rad0,theta,phi,emin,emax,numE)
  comment = outfn
  return "./mxgsDRM %d %d %d %f %f %f %f %g %g %d %g %g %d %s %s"%(0,pdgID,nPriPerE,rad,rad0,theta,phi,emin,emax,numE,outMinE,outMaxE,outNumE,outfn,comment)

def nthLineCmd(filename,n):
  #return("head -%d %s | tail -1"%(n,filename))
  return("sed -n '%d{p;q;}' %s"%(n,filename))

def mpirunCmds(name,pdgID,nPriPerE,rad,rad0,(th0,th1,nth),(ph0,ph1,nph),(e0,e1,ne)):
  #  for th in linRange(th0,th1,nth) for ph in linRange(ph0,ph1,nph)])
  cmds = []

  for (ph,th) in [(ph,th) for th in linRange(th0,th1,nth) for ph in linRange(ph0,ph1,nph)]:
    cmds.append(mxgsDRMCmd(name,pdgID,nPriPerE,rad,rad0,th,ph,e0,e1,ne))
  return "\n".join(cmds)


def commands(name,pdgID,nPriPerE,rad,rad0,theta,phi,energy):
  return ["module unload pgi","module load gcc"
      ,"cd ~/geant/geant4.9.5-install/bin/"
      ,"source geant4.sh"
      ,"cd ~/sim-build"
      ,"mkdir -p %s/%s"%(resultsDir,name)
      ,mpirunCmds(name,pdgID,nPriPerE,rad,rad0,theta,phi,energy)]

def writePBS(outfn,name,n,cmds,walltime):
  outf = open(outfn,"w")
  numNodes = 1 #int(n/2) + n%2
  outf.write("#! /bin/sh -\n")
  outf.write("#PBS -S /bin/sh\n")
  outf.write("#PBS -N \"%s\"\n"%name)
  outf.write("#PBS -A fysisk\n")
  outf.write("#PBS -l walltime=%s,nodes=%d:ppn=1\n"%(walltime,numNodes))
  outf.write("#PBS -l pmem=500mb\n")
  outf.write("#PBS -m abe\n")
  outf.write("#PBS -M brant.carlson@ift.uib.no\n")
  outf.write("#PBS -o %s.stdout.txt\n"%name)
  outf.write("#PBS -e %s.stderr.txt\n"%name)
  outf.write("")
  for cmd in cmds:
    outf.write(cmd)
    outf.write("\n")
  outf.close()


# test runs suggest rate = 450primaries/second, spinup ~ 10 s (for full
# columbus, direction where paricles either hit ASIM or columbus).
# these estimates are padded in case particles shot in from another direction take longer.
def walltimeStr(nPriPerE,nPriE,nth,nph,rate=350,spinup=12):
  sTot = nth*nph*(spinup + nPriPerE*nPriE/rate)
  h = int(sTot/3600)
  m = int((sTot-h*3600)/60)
  s = int((sTot-h*3600-m*60))
  return "%02d:%02d:%02d"%(h,m,s)

def writeJobScript(name,(thmin,thmax,ntheta),(phmin,phmax,nphi),(emin,emax,ne,nPriPerE)):
  pdgID = 22
  writePBS("%s.pbs"%name,name,ntheta*nphi,
      commands(name,pdgID,nPriPerE,0.6,0.0,(thmin,thmax,ntheta),(phmin,phmax,nphi),(emin,emax,ne)),
      walltimeStr(nPriPerE,ne,ntheta,nphi,250))

def writeJobGrid(baseName,(thmin,thmax,nth,nthBlks),(phmin,phmax,nph,nphBlks),(emin,emax,ne,neBlks)):
  thetas = linRange(thmin,thmax,nth)
  phis = linRange(phmin,phmax,nph)
  es = logRange(emin,emax,ne)
  thblks = blocks(thetas,nthBlks)
  phblks = blocks(phis,nphBlks)
  eblks = blocks(es,neBlks)

  nPriPerE = 1000000

  for (thblk,phblk,eblk) in [((min(ths),max(ths),len(ths)),(min(phs),max(phs),len(phs)),(min(es),max(es),len(es),nPriPerE)) for ths in thblks for phs in phblks for es in eblks]:
    name = "%s_%.0f_%.0f_%.0f_%.0f_%.1g_%.1g"%(baseName,thblk[0],thblk[1],phblk[0],phblk[1],eblk[0],eblk[1])
    print(name)
    writeJobScript(name,thblk,phblk,eblk)

#writeJobScript("testJob",0,40,5,0,90,10,10000)

