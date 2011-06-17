#include "SamplingGeometry.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "CZTSD.hh"
#include "BGOSD.hh"
#include "G4SDManager.hh"

#include "HitRecorder.hh"

#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4PropagatorInField.hh"

#include "G4UniformMagField.hh"
#include "globals.hh"
#include <fstream>
#include "G4VisAttributes.hh"

#include "G4GDMLParser.hh"

using namespace std;

SamplingGeometry::SamplingGeometry(double x0, double y0, double z0)
  :G4VUserParallelWorld("SamplingSpheres")
{
  x=x0; y=y0; z=z0;
}

SamplingGeometry::~SamplingGeometry()
{;}

void SamplingGeometry::Construct()
{
  G4VPhysicalVolume *world = GetWorld();
  G4LogicalVolume *worldL = world->GetLogicalVolume();

  double r=0;
  char name[100];
  for(r=1.0*m; r<50.0*m; r*=3.0){
    sprintf(name,"sph_%.1f",r);
    G4Orb *sph = new G4Orb(name,r);
    sprintf(name,"sph_%.1fL",r);
    G4LogicalVolume *sphL = new G4LogicalVolume(sph,0,name,0,0,0);
    sprintf(name,"sph_%.1fP",r);
    new G4PVPlacement(0,G4ThreeVector(x,y,z),sphL,name,worldL,0,0);
  }
}
