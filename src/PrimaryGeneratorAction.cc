#include "PrimaryGeneratorAction.hh"

#include "HitRecorder.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"
#include <sqlite3.h>

using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction(sqlite3* dbin)
{
  ctr=0;
  particleGun = new G4ParticleGun(1);
  pTable = G4ParticleTable::GetParticleTable();
  genIn=0; genOut=1;
  stmt=0;
  currentWeight=1;
  ptime=0;
  db=dbin;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  sqlite3_finalize(stmt);
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  readFromDB();
  if(ctr%1000==0){
    cout << ctr << ": " << enMeV << ' ' 
      << vx << ' ' << vy << ' ' << vz << ' ' 
      << x << ' ' << y << ' ' << z <<  ' '
      << pTable->FindParticle(pdgencoding)->GetPDGMass()/kg << endl;
  }
  ++ctr;

  particleGun->SetParticleDefinition(pTable->FindParticle(pdgencoding));
  particleGun->SetParticlePosition(G4ThreeVector(x*m, y*m, z*m));
  G4ThreeVector v(vx,vy,vz);
  particleGun->SetParticleMomentumDirection(v);
  particleGun->SetParticleEnergy(enMeV*MeV);
  particleGun->SetParticleTime(ptime*s);

  particleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::readFromDB(){
  while(!stmt || sqlite3_step(stmt)!=SQLITE_ROW){
    cout << "initializing sqlite for source particle reads\n";
    sqlite3_finalize(stmt);
    sqlite3_prepare_v2(db,"SELECT olIDX,pdgID,ee,x,y,z,px,py,pz,t,wt FROM pcles WHERE gen=?",-1,&stmt,0);
    sqlite3_bind_int(stmt,1,genIn);
  }

  currentSrcLineNumber = sqlite3_column_int(stmt,0);
  pdgencoding = sqlite3_column_int(stmt,1);
  enMeV = sqlite3_column_double(stmt,2);
  x = sqlite3_column_double(stmt,3);
  y = sqlite3_column_double(stmt,4);
  z = sqlite3_column_double(stmt,5);
  vx = sqlite3_column_double(stmt,6);
  vy = sqlite3_column_double(stmt,7);
  vz = sqlite3_column_double(stmt,8);
  ptime = sqlite3_column_double(stmt,9);
  currentWeight = sqlite3_column_double(stmt,10);

  return;
}

// need to add member variables:
// currentWeight, ptime
// kinetic_energy_MeV -> enMeV
// source line number somehow.
