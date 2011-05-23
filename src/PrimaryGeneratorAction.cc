#include "PrimaryGeneratorAction.hh"

#include "HitRecorder.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"

using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction(char *infilename, HitRecorder* hrec, double width, double rad, double zoff)
{
  ctr=0;
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  pHRec = hrec;
  planewidth=width; planerad=rad; zoffset=zoff;

  infile.open(infilename);

  pTable = G4ParticleTable::GetParticleTable();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  infile.close();
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //ReadFromInfile(); //sets private variables vx,vy,vz, etc.
  
  SetRandomInitialConditions();
  if(ctr%1000==0){
    cout << ctr << ": " << kinetic_energy_MeV << ' ' 
      << vx << ' ' << vy << ' ' << vz << ' ' 
      << x << ' ' << y << ' ' << z <<  ' '
      << pTable->FindParticle(pdgencoding)->GetPDGMass()/kg << endl;
  }
  ++ctr;

  particleGun->SetParticleDefinition(pTable->FindParticle(pdgencoding));

  pHRec->NewPrimary(ctr,pdgencoding,kinetic_energy_MeV,vx,vy,vz,x*m,y*m,z*m);

  G4ThreeVector v(vx,vy,vz);
  particleGun->SetParticleMomentumDirection(v);
  particleGun->SetParticleEnergy(kinetic_energy_MeV*MeV);

  particleGun->SetParticlePosition(G4ThreeVector(x*m, y*m, z*m));
  particleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::NextLine(){
  infile >> pdgencoding >> kinetic_energy_MeV >> theta;
}

void PrimaryGeneratorAction::SetRandomInitialConditions(){
  double x0=0,y0=0,z0=0;
  x0 = planewidth*(G4UniformRand()-0.5);
  y0 = planewidth*(G4UniformRand()-0.5);

  double phi=2*3.14159*(G4UniformRand()-0.5);

  // rotate by theta about y axis
  double x1,y1,z1;
  x1 = x0*cos(theta)+z0*sin(theta);
  y1 = y0;
  z1 = z0*cos(theta)-x0*sin(theta);

  // rotate by phi about z axis
  double x2,y2,z2;
  x2 = x1*cos(phi)-y1*sin(phi);
  y2 = y1*cos(phi)+x1*sin(phi);
  z2 = z1;

  // translate away from origin, offset to center of satellite.
  double x3,y3,z3;
  x3 = x2+planerad*sin(theta)*cos(phi);
  y3 = y2+planerad*sin(theta)*sin(phi);
  z3 = z2+planerad*cos(theta)+zoffset;

  // "output" step, set member variables.
  x=x3; y=y3; z=z3;

  vx=-sin(theta)*cos(phi);
  vy=-sin(theta)*sin(phi);
  vz=-cos(theta);

}
