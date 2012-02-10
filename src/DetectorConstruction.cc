#include "DetectorConstruction.hh"

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
#include "gsl/gsl_histogram.h"

using namespace std;

DetectorConstruction::DetectorConstruction(gsl_histogram *hb, gsl_histogram *hc){
  hBGO = hb; hCZT = hc;
}

DetectorConstruction::~DetectorConstruction()
{;}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  char fn[1000];
  //sprintf(fn,"%s/columbus.gdml",getenv("MY_GDML_PATH"));
  sprintf(fn,"columbus.gdml");
  //// NOTE: using asim.gdml, need to use asim_world volume instead of asim volume for world ref in asim.gdml.
  //// i.e. <world ref="asim_world"> vs <world ref="asim"> at end of asim.gdml.
  //// also have to adjust the target position.
  //
  G4GDMLParser parser;
  parser.Read(fn,false);
  G4VPhysicalVolume *W = parser.GetWorldVolume();
  //W->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);

  const G4GDMLAuxMapType* auxmap = parser.GetAuxMap();
  std::cout << "Found " << auxmap->size()
    << " volume(s) with auxiliary information."
    << G4endl ;
  char sdname[100];
  int bgoctr=0;
  int cztctr=2;
  for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin(); iter!=auxmap->end(); iter++) {
    for(G4GDMLAuxListType::const_iterator vit=(*iter).second.begin(); vit!=(*iter).second.end();vit++) {
      if ((*vit).type=="visibility" && (*vit).value=="invisible") {
        G4cout << "rendering logical volume " << (*iter).first->GetName()
          << " invisible."<< G4endl;
        G4LogicalVolume* myvol = (*iter).first;
        myvol->SetVisAttributes(G4VisAttributes::Invisible);
      }
      if ((*vit).type=="sensitive" && (*vit).value=="czt") {
        G4cout << "making logical volume " << (*iter).first->GetName() << " CZT sensitive."<< G4endl;
        G4LogicalVolume* myvol = (*iter).first;
        sprintf(sdname,"bgo_%d",bgoctr); ++bgoctr;
        CZTSD *czt = new CZTSD(sdname,hCZT);
        myvol->SetSensitiveDetector(czt);
        G4SDManager::GetSDMpointer()->AddNewDetector(czt);
      }
      if ((*vit).type=="sensitive" && (*vit).value=="bgo") {
        G4cout << "making logical volume " << (*iter).first->GetName() << " BGO sensitive."<< G4endl;
        G4LogicalVolume* myvol = (*iter).first;
        sprintf(sdname,"czt_%d",cztctr); ++cztctr;
        BGOSD *bgo = new BGOSD(sdname,hBGO);
        myvol->SetSensitiveDetector(bgo);
        G4SDManager::GetSDMpointer()->AddNewDetector(bgo);
      }
    }
  }

  return W;
}
