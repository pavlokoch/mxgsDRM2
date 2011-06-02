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
    PrimaryGeneratorAction(sqlite3 *dbin);
    ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    void NextGen(){genIn=genOut; ++genOut;}

  private:
    void readFromDB();
    G4ParticleGun* particleGun;
    G4ParticleTable *pTable;
    int pdgencoding,ctr;
    double enMeV,theta,vx,vy,vz,x,y,z;

    int genIn,genOut;
    sqlite3 *db;
    sqlite3_stmt* stmt;
    double currentWeight;
    double ptime;
    int currentSrcLineNumber;

    void SetRandomInitialConditions();
};

#endif


