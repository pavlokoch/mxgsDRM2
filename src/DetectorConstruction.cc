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

#include "LayerSD.hh"
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

DetectorConstruction::DetectorConstruction(HitRecorder *hrec)
 :  experimentalHall_log(0), 
    experimentalHall_phys(0),
    targetdepth(targetdepth)
{pHRec = hrec;}

DetectorConstruction::~DetectorConstruction()
{;}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4GDMLParser parser;
  //parser.Read("iss_C3I2.gdml");
  parser.Read("test.gdml");
  G4VPhysicalVolume *W = parser.GetWorldVolume();
  //W->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);

  const G4GDMLAuxMapType* auxmap = parser.GetAuxMap();
  std::cout << "Found " << auxmap->size()
    << " volume(s) with auxiliary information."
    << G4endl ;
  for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin(); iter!=auxmap->end(); iter++) {
    for(G4GDMLAuxListType::const_iterator vit=(*iter).second.begin(); vit!=(*iter).second.end();vit++) {
      if ((*vit).type=="visibility") {
        G4cout << "rendering logical volume " << (*iter).first->GetName()
          << " invisible."<< G4endl;
        G4LogicalVolume* myvol = (*iter).first;
        myvol->SetVisAttributes(G4VisAttributes::Invisible);
      }
    }
  }

  return W;
}
