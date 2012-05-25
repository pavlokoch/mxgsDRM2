w=43.0        # width of hole
dx=dy=46.0    # grid spacing
sep=dx-w      # separation between holes
x0=0.0        # origin
y0=0.0        # origin
z=0.0         # origin
thk=1.0       # thickness of mask
cmwid = 742   # width of overall mask

repx = 1      # number of repeats of mask pattern in x direction
repy = 1      # number of repeats of mask pattern in y direction
repsep_dx = 0 # separation between repetitions of mask pattern
repsep_dy = 0 # separation between repetitions of mask pattern


def name(i,j,suffix=""):
  return "cm_%02d_%02d%s"%(i,j,suffix)

def printPos(open,i,j,xoff,yoff):
  x=x0+i*dx
  if(open==1):
    y=y0-j*dy+dy/2-sep/2
  else:
    y=y0-j*dy
  print "<position name=\"%s\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(name(i,j,"_pos"),x-xoff,y+yoff,z)

def printVol(open,i,j):
  if(open==1):
    print "<volume name=\"%s\"> <materialref ref=\"W\"/> <solidref ref=\"cmopen\"/> </volume>"%(name(i,j,"_v"))
  else:
    print "<volume name=\"%s\"> <materialref ref=\"W\"/> <solidref ref=\"cmclosed\"/> </volume>"%(name(i,j,"_v"))

def printPhysVol(open,i,j):
  print "<physvol> <volumeref ref=\"%s\"/> <positionref ref=\"%s\"/> <rotationref ref=\"id\"/> </physvol>"%(name(i,j,"_v"),name(i,j,"_pos"))

def printPositions(pattern):
  xoff = len(pattern[0])*dx/2-w/2
  yoff = len(pattern)*dy/2-w/2-sep/2
  for j in range(len(pattern[0])):
    print "<position name=\"cm_ysep_%02d_pos\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(j,j*dx+dx/2-xoff,-len(pattern)/2.0*dy+dy/2+yoff,z)
  for j in range(len(pattern)):
    for i in range(len(pattern[j])):
      printPos(pattern[j][i],i,j,xoff,yoff)

def printRepSepPositions(nx,ny,pattern,xsep,ysep):
  xwid = len(pattern[1])*dx + xsep
  ywid = len(pattern)*dy + ysep
  xmin = -(float(nx)/2 - 0.5)*xwid
  ymin = -(float(ny)/2 - 0.5)*ywid
  print "<position name=\"xminsep_pos\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%((xsep+sep)/2.0-xwid*(nx/2.0)-sep/2.0,0,z)
  print "<position name=\"yminsep_pos\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(sep/2,(ysep+sep)/2.0-ywid*(ny/2.0)-sep/2.0,z)

  xmw = len(pattern[1])*dx*nx + (nx-1)*xsep + sep
  ymw = len(pattern)*dy*ny + (ny-1)*ysep + sep
  xbwid = (cmwid - xmw)/2.0
  ybwid = (cmwid - ymw)/2.0
  print "<position name=\"ymb_pos\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(0,-cmwid/2+ybwid/2,z)
  print "<position name=\"ypb_pos\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(0,cmwid/2-ybwid/2,z)
  print "<position name=\"xmb_pos\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(-cmwid/2+xbwid/2,0,z)
  print "<position name=\"xpb_pos\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(cmwid/2-xbwid/2,0,z)
  for i in range(nx):
    for j in range(ny):
      print "<position name=\"%s\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(name(i,j,"_reppos"),xmin + i*xwid+sep/2, ymin+j*ywid+sep/2,z)
      if(i < nx-1 and repsep_dx > 0):
        print "<position name=\"%s\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(name(i,j,"_xrepseppos"),xmin + i*xwid + xwid/2+sep/2, ymin+j*ywid+sep/2,z)
      if(j < ny-1 and repsep_dy > 0):
        print "<position name=\"%s\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(name(i,j,"_yrepseppos"),xmin + i*xwid+sep/2, ymin+j*ywid + ywid/2+sep/2,z)
      if(j < ny-1 and i<nx-1 and repsep_dx > 0 and repsep_dy > 0):
        print "<position name=\"%s\" unit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(name(i,j,"_xyrepseppos"),xmin + i*xwid + xwid/2+sep/2, ymin+j*ywid + ywid/2+sep/2,z)


def printBoxes(pattern):
  print "<box name=\"cmopen\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(w,sep,thk)
  print "<box name=\"cmclosed\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(w,dy,thk)
  print "<box name=\"cmysep\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(sep,len(pattern)*dy,thk)

def printRepSepBoxes(nx,ny,pattern,xsep,ysep):
  xwid = len(pattern[1])*dx
  ywid = len(pattern)*dy
  xmw = len(pattern[1])*dx*nx + (nx-1)*xsep + sep
  ymw = len(pattern)*dy*ny + (ny-1)*ysep + sep
  xbwid = (cmwid - xmw)/2.0
  ybwid = (cmwid - ymw)/2.0
  print "<box name=\"xrepsep\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(xsep,ywid,thk)
  print "<box name=\"yrepsep\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(xwid,ysep,thk)
  print "<box name=\"xyrepsep\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(xsep,ysep,thk)
  print "<box name=\"xminsep\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(sep,ywid*ny+(ny-1)*ysep+sep,thk)
  print "<box name=\"yminsep\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(xwid*nx+(nx-1)*xsep,sep,thk)

  print "<box name=\"yb\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(cmwid,ybwid,thk)
  print "<box name=\"xb\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(xbwid,ymw,thk)

def printVolumes(pattern):
  for j in range(len(pattern[0])):
    print "<volume name=\"cm_ysep_%02d_v\"> <materialref ref=\"W\"/> <solidref ref=\"cmysep\"/> </volume>"%(j)
  for j in range(len(pattern)):
    for i in range(len(pattern[j])):
      printVol(pattern[j][i],i,j)

def printRepSepVolumes(nx,ny):
  if(repsep_dy > 0):
    print "<volume name=\"yrepsep_v\"> <materialref ref=\"W\"/> <solidref ref=\"yrepsep\"/> </volume>"
  if(repsep_dx > 0):
    print "<volume name=\"xrepsep_v\"> <materialref ref=\"W\"/> <solidref ref=\"xrepsep\"/> </volume>"
  if(repsep_dy > 0 and repsep_dx > 0):
    print "<volume name=\"xyrepsep_v\"> <materialref ref=\"W\"/> <solidref ref=\"xyrepsep\"/> </volume>"

  print "<volume name=\"xminsep_v\"> <materialref ref=\"W\"/> <solidref ref=\"xminsep\"/> </volume>"
  print "<volume name=\"yminsep_v\"> <materialref ref=\"W\"/> <solidref ref=\"yminsep\"/> </volume>"
  print "<volume name=\"xb_v\"> <materialref ref=\"W\"/> <solidref ref=\"xb\"/> </volume>"
  print "<volume name=\"yb_v\"> <materialref ref=\"W\"/> <solidref ref=\"yb\"/> </volume>"

def printPhysVolumes(pattern):
  for j in range(len(pattern[0])):
    print "<physvol> <volumeref ref=\"cm_ysep_%02d_v\"/> <positionref ref=\"cm_ysep_%02d_pos\"/> <rotationref ref=\"id\"/> </physvol>"%(j,j)
  for j in range(len(pattern)):
    for i in range(len(pattern[j])):
      printPhysVol(pattern[j][i],i,j)

def printVolume(pattern):
  print "<volume name=\"codedMask_rep\">"
  print "<materialref ref=\"Vacuum\"/>"
  print "<solidref ref=\"codedMask_rep_mv\"/>"
  printPhysVolumes(pattern)
  print "<!--"
  print "<auxiliary auxtype=\"visibility\" auxvalue=\"invisible\"/>"
  print "-->"
  print "</volume>"

def printRepPhysVols(nx,ny):
  print "<physvol> <volumeref ref=\"xminsep_v\"/> <positionref ref=\"xminsep_pos\"/> <rotationref ref=\"id\"/> </physvol>"
  print "<physvol> <volumeref ref=\"yminsep_v\"/> <positionref ref=\"yminsep_pos\"/> <rotationref ref=\"id\"/> </physvol>"
  print "<physvol> <volumeref ref=\"xb_v\"/> <positionref ref=\"xmb_pos\"/> <rotationref ref=\"id\"/> </physvol>"
  print "<physvol> <volumeref ref=\"xb_v\"/> <positionref ref=\"xpb_pos\"/> <rotationref ref=\"id\"/> </physvol>"
  print "<physvol> <volumeref ref=\"yb_v\"/> <positionref ref=\"ymb_pos\"/> <rotationref ref=\"id\"/> </physvol>"
  print "<physvol> <volumeref ref=\"yb_v\"/> <positionref ref=\"ypb_pos\"/> <rotationref ref=\"id\"/> </physvol>"
  for i in range(nx):
    for j in range(ny):
      print "<physvol> <volumeref ref=\"codedMask_rep\"/> <positionref ref=\"%s\"/> <rotationref ref=\"id\"/> </physvol>"%name(i,j,"_reppos")

      # separators of repetitions of coded mask pattern.
      if(i < nx-1 and repsep_dx>0):
        print "<physvol> <volumeref ref=\"xrepsep_v\"/> <positionref ref=\"%s\"/> <rotationref ref=\"id\"/> </physvol>"%name(i,j,"_xrepseppos")
      if(j < ny-1 and repsep_dy>0):
        print "<physvol> <volumeref ref=\"yrepsep_v\"/> <positionref ref=\"%s\"/> <rotationref ref=\"id\"/> </physvol>"%name(i,j,"_yrepseppos")
      if(j < ny-1 and i<nx-1 and repsep_dx>0 and repsep_dy>0):
        print "<physvol> <volumeref ref=\"xyrepsep_v\"/> <positionref ref=\"%s\"/> <rotationref ref=\"id\"/> </physvol>"%name(i,j,"_xyrepseppos")

def printGDML(pattern):
  print "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>"
  print "<gdml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd\">"
  print ""
  print "<define>"
  print "<rotation name=\"id\"/>"
  printPositions(pattern)
  printRepSepPositions(repx,repy,pattern,repsep_dx,repsep_dy)
  print "</define>"
  print "<materials>"
  print "<material name=\"W\" Z=\"74\"> <D value=\"19.25\" unit=\"g/cm3\"/> <atom value=\"183.84\"/> </material>"
  print "<material name=\"Vacuum\" Z=\"1.0\"> <D value=\"1.0e-25\"/> <atom value=\"1.00794\"/> </material>"
  print "</materials>"
  print "<solids>"
  print "<box name=\"codedMask_rep_mv\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(len(pattern[0])*dx,len(pattern)*dy,1)
  print "<box name=\"codedMask_mv\" lunit=\"mm\" x=\"%f\" y=\"%f\" z=\"%f\"/>"%(cmwid,cmwid,1)
  printBoxes(pattern)
  printRepSepBoxes(repx,repy,pattern,repsep_dx,repsep_dy)
  print "</solids>"
  print "<structure>"
  printVolumes(pattern)
  printRepSepVolumes(repx,repy)
  print ""
  printVolume(pattern)
  print "<volume name=\"codedMask\">"
  print "<materialref ref=\"Vacuum\"/>"
  print "<solidref ref=\"codedMask_mv\"/>"
  printRepPhysVols(repx,repy)
  print ""
  print "<!--"
  print "<auxiliary auxtype=\"visibility\" auxvalue=\"invisible\"/>"
  print "-->"
  print "</volume>"
  print "</structure>"
  print "<setup name=\"Default\" version=\"1.0\">"
  print "<world ref=\"codedMask\"/>"
  print "</setup>"
  print "</gdml>"

pattern = [
    [1,1,1,1,0,1,1,1,0,1,1],
    [0,0,0,1,1,1,1,0,0,1,0],
    [0,1,1,0,1,1,0,0,1,1,1],
    [1,1,1,1,1,0,1,1,1,1,0],
    [1,0,1,1,0,0,1,0,0,1,1],
    [1,1,1,0,0,1,1,1,1,0,1],
    [1,1,1,1,0,1,1,0,0,1,1],
    [0,0,0,1,1,0,1,1,0,1,0],
    [1,1,1,0,1,1,0,1,1,1,1],
    [1,1,0,1,1,1,1,1,1,1,0],
    [1,0,0,1,0,1,1,0,0,1,1],
    [1,0,1,1,1,1,0,0,1,0,1],
    [1,0,1,1,0,1,1,1,1,1,1],
    [0,1,0,1,1,1,1,0,1,1,0],
    [1,0,1,0,1,1,0,0,1,0,0],
    [1,1,0,1,1,0,1,1,1,1,1],
    [1,0,0,1,0,0,1,0,1,1,0],
    [0,1,1,1,1,1,1,1,1,0,0],
    [1,1,1,1,1,1,1,0,0,1,1],
    [1,0,1,1,0,0,1,1,0,1,0],
    [1,1,1,0,0,1,0,1,1,1,1],
    [0,1,1,0,1,1,0,1,0,0,1],
    [1,0,0,1,0,1,1,0,1,0,0]]

pattern2 = [
    [0,0,0,0,0,1,1,1],
    [1,0,1,1,1,0,0,1],
    [1,0,1,0,1,1,0,1],
    [1,0,1,0,1,0,0,0],
    [0,1,0,1,0,0,1,0],
    [1,1,1,0,1,1,0,0],
    [1,1,1,1,1,0,0,0],
    [1,1,1,1,1,1,0,1]]

pattern3 = [
    [0,0,0,0,0,1,1,1, 0,0,0,0,0,1,1,1],
    [1,0,1,1,1,0,0,1, 1,0,1,1,1,0,0,1],
    [1,0,1,0,1,1,0,1, 1,0,1,0,1,1,0,1],
    [1,0,1,0,1,0,0,0, 1,0,1,0,1,0,0,0],
    [0,1,0,1,0,0,1,0, 0,1,0,1,0,0,1,0],
    [1,1,1,0,1,1,0,0, 1,1,1,0,1,1,0,0],
    [1,1,1,1,1,0,0,0, 0,1,1,1,1,0,0,0],
    [1,1,1,1,1,1,0,0, 0,0,1,1,1,1,0,1],

    [0,0,0,0,0,1,0,0, 0,0,0,0,0,1,1,1],
    [1,0,1,1,1,0,0,0, 0,0,1,1,1,0,0,1],
    [1,0,1,0,1,1,0,1, 1,0,1,0,1,1,0,1],
    [1,0,1,0,1,0,0,0, 1,0,1,0,1,0,0,0],
    [0,1,0,1,0,0,1,0, 0,1,0,1,0,0,1,0],
    [1,1,1,0,1,1,0,0, 1,1,1,0,1,1,0,0],
    [1,1,1,1,1,0,0,0, 1,1,1,1,1,0,0,0],
    [1,1,1,1,1,1,0,1, 1,1,1,1,1,1,0,1]]

printGDML(pattern3)
