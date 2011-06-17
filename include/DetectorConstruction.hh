#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4LogicalVolume;
class G4PVPlacement;
//class G4VPhysicalVolume;

#include "HitRecorder.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UniformMagField.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"
#include <sqlite3.h>

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction(sqlite3 *dbin);
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

  private:
    sqlite3 *db;
};

#endif

