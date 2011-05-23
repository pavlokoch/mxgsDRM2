#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "HitRecorder.hh"
#include "G4ParticleTable.hh"
#include <fstream>

using namespace std;

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(char *infilename, HitRecorder *hrec,
        double width, double rad,double zoff);
    ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    void NextLine();

  private:
    HitRecorder *pHRec;
    G4ParticleGun* particleGun;
    G4ParticleTable *pTable;
    ifstream infile;
    int pdgencoding,ctr;
    double kinetic_energy_MeV,theta,vx,vy,vz,x,y,z;

    double planewidth,planerad,zoffset;

    void SetRandomInitialConditions();
};

#endif


