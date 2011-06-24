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

using namespace std;

DetectorConstruction::DetectorConstruction(sqlite3 *dbin){
  db = dbin;
}

DetectorConstruction::~DetectorConstruction()
{;}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4GDMLParser parser;
  parser.Read("columbus.gdml");
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
        CZTSD *czt = new CZTSD(sdname,db);
        myvol->SetSensitiveDetector(czt);
        G4SDManager::GetSDMpointer()->AddNewDetector(czt);
      }
      if ((*vit).type=="sensitive" && (*vit).value=="bgo") {
        G4cout << "making logical volume " << (*iter).first->GetName() << " BGO sensitive."<< G4endl;
        G4LogicalVolume* myvol = (*iter).first;
        sprintf(sdname,"czt_%d",cztctr); ++cztctr;
        BGOSD *bgo = new BGOSD(sdname,db);
        myvol->SetSensitiveDetector(bgo);
        G4SDManager::GetSDMpointer()->AddNewDetector(bgo);
      }
    }
  }

  return W;
}
