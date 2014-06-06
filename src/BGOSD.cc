#include "BGOSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4StepPoint.hh"
#include "EventInfo.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "gsl/gsl_histogram.h"
#include "CLHEP/Units/SystemOfUnits.h"

BGOSD::BGOSD(G4String name, gsl_histogram *hb)
:G4VSensitiveDetector(name)
{
  h = hb;
  evtctr=0;
}

BGOSD::~BGOSD(){;}

void BGOSD::Initialize(G4HCofThisEvent*HCE){
  edep=0;
}

G4bool BGOSD::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
  G4double thisEdep = aStep->GetTotalEnergyDeposit()/CLHEP::MeV;
  if(thisEdep>0){
    edep += thisEdep;
  }
  return true;
}

void BGOSD::EndOfEvent(G4HCofThisEvent* /*HCE*/){
  if(edep>0){
    gsl_histogram_increment(h,edep);
  }
  ++evtctr;
}

