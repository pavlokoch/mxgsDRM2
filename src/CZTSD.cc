#include "CZTSD.hh"
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

CZTSD::CZTSD(G4String name, gsl_histogram *hc)
:G4VSensitiveDetector(name)
{
  h = hc;
  evtctr=0;
}

CZTSD::~CZTSD(){;}

void CZTSD::Initialize(G4HCofThisEvent*HCE){
  edep=0;
}

G4bool CZTSD::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
  G4double thisEdep = aStep->GetTotalEnergyDeposit()/CLHEP::MeV;
  if(thisEdep>0){
    edep += thisEdep;
  }
  return true;
}

void CZTSD::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
  if(edep>0){
    gsl_histogram_increment(h,edep);
  }
  ++evtctr;
}

