#include "PrimaryGeneratorAction.hh"

#include "HitRecorder.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"
#include <sqlite3.h>
#include "EventInfo.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction(int pdgID, double startDiskRad, 
    double startDiskRad0, double theta_deg, double phi_deg, 
    double targetX, double targetY, double targetZ, double targetDisp)
{
  ctr=0;
  particleGun = new G4ParticleGun(1);
  pTable = G4ParticleTable::GetParticleTable();
  pdgencoding = pdgID;
  th = theta_deg*CLHEP::pi/180.0; ph = phi_deg*CLHEP::pi/180.0;
  x0 = targetX; y0 = targetY; z0 = targetZ;
  rDisk = startDiskRad;
  rDisk0 = startDiskRad0;
  disp = targetDisp;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  drawPrimary();
  if(ctr%1000==0){
    cout << ctr << ": " << enMeV << ' ' 
      << vx << ' ' << vy << ' ' << vz << ' ' 
      << x << ' ' << y << ' ' << z <<  ' '
      << pTable->FindParticle(pdgencoding)->GetPDGMass()/CLHEP::kg << endl;
  }
  ++ctr;

  particleGun->SetParticleDefinition(pTable->FindParticle(pdgencoding));
  particleGun->SetParticlePosition(G4ThreeVector(x*CLHEP::m, y*CLHEP::m, z*CLHEP::m));
  G4ThreeVector v(vx,vy,vz);
  particleGun->SetParticleMomentumDirection(v);
  particleGun->SetParticleEnergy(enMeV*CLHEP::MeV);

  particleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::drawPrimary(){
  double zhx = 0;
  double zhy = 0;
  double zhz = 1;

  // rotate zh about y axis by theta
  double tx = zhx*cos(th) + zhz*sin(th);
  double ty = zhy;
  double tz = zhz*cos(th) - zhx*sin(th);

  // rotate zh about z axis by phi
  zhx = tx*cos(ph) - ty*sin(ph);
  zhy = ty*cos(ph) + tx*sin(ph);
  zhz = tz;

  // do likewise for an initial point
  //double r = rDisk*sqrt(CLHEP::RandFlat::shoot());
  double r = sqrt((rDisk*rDisk - rDisk0*rDisk0)*CLHEP::RandFlat::shoot() + rDisk0*rDisk0);
  //double r = rDisk*sqrt(0.75*CLHEP::RandFlat::shoot() + 0.25);

  double diskth = CLHEP::RandFlat::shoot(2*CLHEP::pi);
  x = r*cos(diskth);
  y = r*sin(diskth);
  z = 0;

  tx = x*cos(th) + z*sin(th);
  ty = y;
  tz = z*cos(th) - x*sin(th);

  x = tx*cos(ph) - ty*sin(ph);
  y = ty*cos(ph) + tx*sin(ph);
  z = tz;

  // displace initial point to target and away to source zone.
  x += x0 - zhx*disp;
  y += y0 - zhy*disp;
  z += z0 - zhz*disp;

  // set initial velocity back towards target.
  vx = zhx;
  vy = zhy;
  vz = zhz;

  //cerr << "INPUTTEST: " << enMeV << ' ' << x << ' ' << y << ' ' << z << ' ' << vx << ' ' << vy << ' ' << vz << endl;

  return;
}

void PrimaryGeneratorAction::setPriEn(double en){
  enMeV = en;
}

