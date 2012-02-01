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
#include "gsl/gsl_histogram.h"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction(gsl_histogram *hb, gsl_histogram *hc);
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

  private:
    gsl_histogram *hBGO, *hCZT;
};

#endif

