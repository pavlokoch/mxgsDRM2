#ifndef SamplingGeometry_H
#define SamplingGeometry_H 1

class G4LogicalVolume;
class G4PVPlacement;
//class G4VPhysicalVolume;

#include "HitRecorder.hh"
#include "G4VUserParallelWorld.hh"
#include "G4UniformMagField.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"
#include <sqlite3.h>

class SamplingGeometry : public G4VUserParallelWorld
{
  public:

    SamplingGeometry(double x0, double y0, double z0);
    ~SamplingGeometry();

    void Construct();

  private:
    double x,y,z;
};

#endif

