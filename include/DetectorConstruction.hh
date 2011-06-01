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

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction(sqlite3 *dbin);
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

  private:
    sqlite3 *db;
    
    // Logical volumes
    //
    G4LogicalVolume* experimentalHall_log;
    G4LogicalVolume** airLayer_log;

    // Physical volumes
    //
    G4PVPlacement* experimentalHall_phys;
    G4PVPlacement** airLayer_phys;

    void makeDetectorLayer(double xcenter, double ycenter, double zcenter, double
        xwid, double ywid, double zwid, G4Material *mat, G4LogicalVolume *parent,
        HitRecorder *pHRec, int layer);
    void makeHalfCollimatorCell(double altitude, double x, double y, double z,
        double depth, G4Material *mat, G4LogicalVolume *parent, G4RotationMatrix
        *zRot1, G4RotationMatrix *zRot2); 
    void makeCollimator( double xcorner, double ycorner, double zcorner,
        double collmod_width, double collmod_length,
        G4Material *mat, G4LogicalVolume *experimentalHall_log);
    
    // Magnetic field
    G4UniformMagField *bfield;

};

#endif

