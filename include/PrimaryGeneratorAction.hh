#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "HitRecorder.hh"
#include "G4ParticleTable.hh"
#include <fstream>
#include <sqlite3.h>

using namespace std;

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(int pdgID, double startDiskRad, double startDiskRad0,
        double theta_deg, double phi_deg, 
        double targetX, double targetY, double targetZ, double targetDisp);
    ~PrimaryGeneratorAction();
    void setPriEn(double en);

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    void drawPrimary();
    G4ParticleGun* particleGun;
    G4ParticleTable *pTable;
    int pdgencoding,ctr;
    double enMeV,th,ph,vx,vy,vz,x,y,z,x0,y0,z0;
    double rDisk, rDisk0, disp;
};

#endif


