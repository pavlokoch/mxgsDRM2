import random
from math import *
import sqlite3

# energies in MeV
# angles in radians
# positions in m
# times in ns

def transformationFunc_ZAxisToThetaPhiDir((theta,phi)):
  def xform(x0,y0,z0):
    return transformation_ZAxisToThetaPhiDir((theta,phi),x0,y0,z0)
  return xform

def transformation_ZAxisToThetaPhiDir((theta,phi),x0,y0,z0):
  (x1,y1,z1) = (x0*cos(theta)+z0*sin(theta),y0,z0*cos(theta)-x0*sin(theta))
  (x2,y2,z2) = (x1*cos(phi)-y1*sin(phi),y1*cos(phi)+x1*sin(phi),z1)
  return (x2,y2,z2)

def randomParticle_plane(planeWid,planeDist,pdgCode,origin,eMeV):
  (x0,y0,z0) = origin
  def randomPcle():
    theta = pi - pi/2*random.random() # PI OVER 2
    phi =2*pi*random.random()
    x = random.random()*planeWid-planeWid/2.0
    y = random.random()*planeWid-planeWid/2.0
    z = planeDist;
    (x,y,z) = transformation_ZAxisToThetaPhiDir((theta,phi),x,y,z)
    px = -sin(theta)*cos(phi);
    py = -sin(theta)*sin(phi);
    pz = -cos(theta);
    return (pdgCode,eMeV,x+x0,y+y0,z+z0,px,py,pz,0.0)
  return randomPcle

def randomParticle_planeThetaPhi(planeWid,planeDist,pdgCode,theta,phi,origin,eMeV):
  tr = transformationFunc_ZAxisToThetaPhiDir((theta,phi))
  (x0,y0,z0) = origin
  def randomPcle():
    x = random.random()*planeWid-planeWid/2.0
    y = random.random()*planeWid-planeWid/2.0
    z = planeDist;
    (x,y,z) = tr(x,y,z)
    px = -sin(theta)*cos(phi);
    py = -sin(theta)*sin(phi);
    pz = -cos(theta);
    return (pdgCode,eMeV,x,y,z,px,py,pz,0.0)
  return randomPcle

def resetDB(c):
  c.execute("DROP TABLE IF EXISTS pcles;")
  c.execute("DROP TABLE IF EXISTS evtsCZT;")
  c.execute("DROP TABLE IF EXISTS evtsBGO;")
  c.execute("DROP TABLE IF EXISTS hitsCZT;")
  c.execute("DROP TABLE IF EXISTS hitsBGO;")

  c.execute("CREATE TABLE pcles (idx INTEGER, oIdx INTEGER, pIdx INTEGER, gen INTEGER, pdgID INTEGER, ee REAL, x REAL, y REAL, z REAL, px REAL, py REAL, pz REAL, t REAL, wt REAL);")
  c.execute("CREATE TABLE evtsCZT (idx INTEGER, priIdx INTEGER, totEDep REAL);")
  c.execute("CREATE TABLE hitsCZT (evtIdx INTEGER, edep REAL, segIdx INTEGER, x REAL, y REAL, z REAL, t REAL);")
  c.execute("CREATE TABLE evtsBGO (idx INTEGER, priIdx INTEGER, totEDep REAL);")
  c.execute("CREATE TABLE hitsBGO (evtIdx INTEGER, edep REAL, segIdx INTEGER, x REAL, y REAL, z REAL, t REAL);")

def makeGEANTInputDB(dbfn,nPri,planeWid,planeDist,pdgCode,theta,phi,origin,eMeV):
  conn = sqlite3.connect(dbfn)
  c = conn.cursor()
  resetDB(c)
  ctr = 0
  #rpfn = randomParticle_planeThetaPhi(nPri,planeWid,planeDist,pdgCode,theta,phi,origin,eMeV)
  rpfn = randomParticle_plane(planeWid,planeDist,pdgCode,origin,eMeV)
  while(ctr<nPri):
    (pdgCode,eMeV,x,y,z,px,py,pz,t) = rpfn()
    c.execute("INSERT INTO pcles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);",(ctr,ctr,-1,0,pdgCode,eMeV,x,y,z,px,py,pz,t,1.0));
    ctr = ctr+1
  conn.commit()
  c.close()
  conn.close()

def makeGEANTInputDB2(dbfn,nPri,pdgCode,alt_m,eMeV):
  conn = sqlite3.connect(dbfn)
  c = conn.cursor()
  resetDB(c)
  ctr = 0
  (x0,y0,z0) = (0.0,-0.0560,0.0160) # center of top of BGO bar.  add 5 to get to top of housing.  from gdml, but in meters.
  dy = 0.1520/2.0; dx = 0.0520/2.0; # width and length of bgo bar from gdml, but in meters.
  rpfn = randomParticle_atSquare(pdgCode,eMeV,x0,y0,z0+0.005+alt_m, x0-dx,y0-dy,z0,x0+dx,y0+dy,z0)
  while(ctr<nPri):
    (pdgCode,eMeV,x,y,z,px,py,pz,t) = rpfn()
    c.execute("INSERT INTO pcles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);",(ctr,ctr,-1,0,pdgCode,eMeV,x,y,z,px,py,pz,t,1.0));
    ctr = ctr+1
  conn.commit()
  c.close()
  conn.close()

# generates a particle launched from x0,y0,z0 going toward the square defined by
# xs1,ys1,zs1 and xs2,ys2,zs2, assuming zs1==zs2 for now.
def randomParticle_atSquare(pdgCode,eMeV,x0,y0,z0,xs1,ys1,zs1,xs2,ys2,zs2):
  # zs1 assumed to == zs2.
  (mx,my,mz) = ((xs1+xs2)/2.0,(ys1+ys2)/2.0,(zs1+zs2)/2.0)
  if((mx-x0)**2+(my-y0)**2 == 0):
    if(z0>mz):
      theta0 = pi
    else:
      theta0 = 0.0
  else:
    theta0 = atan2(sqrt((mx-x0)**2+(my-y0)**2),(mz-z0))
  phi0 = atan2(my-y0,mx-x0)
  def findMaxAlongSeg(x1,y1,z1,x2,y2,z2):
    t = 0
    thmx = 0
    while(t<=1):
      xx = x1+(x2-x1)*t
      yy = y1+(y2-y1)*t
      zz = z1+(z2-z1)*t
      thmx = max(thmx,acos(((xx-x0)*(mx-x0)+(yy-y0)*(my-y0)+(zz-z0)*(mz-z0))/(sqrt((mx-x0)**2+(my-y0)**2+(mz-z0)**2)*sqrt((xx-x0)**2+(yy-y0)**2+(zz-z0)**2))))
      t += 0.05
    return thmx
  thetaMax = max(findMaxAlongSeg(xs1,ys1,zs1,xs2,ys1,zs1),findMaxAlongSeg(xs1,ys1,zs1,xs1,ys2,zs1),
      findMaxAlongSeg(xs1,ys2,zs2,xs2,ys2,zs2),findMaxAlongSeg(xs2,ys1,zs2,xs2,ys2,zs2))
  def randomPcle():
    xp=yp=1.0e100
    while(xp<xs1 or xp>xs2 or yp<ys1 or yp>ys2):
      phi = 2*pi*random.random()
      theta = acos(cos(thetaMax) + (1-cos(thetaMax))*random.random())
      (x,y,z) = (sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
      (xt,yt,zt) = transformation_ZAxisToThetaPhiDir((theta0,phi0),x,y,z)
      dr = abs((z0-mz)/zt)
      xp = x0+xt*dr
      yp = y0+yt*dr
    return (pdgCode,eMeV,x0,y0,z0,xt,yt,zt,0.0)
  return randomPcle
