#include "LayerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

LayerSD::LayerSD(G4String name, HitRecorder *hrec, int lay, int wi)
:G4VSensitiveDetector(name)
{
  layer=lay;
  wire=wi;
  pHRec = hrec;
}

LayerSD::~LayerSD(){;}

void LayerSD::Initialize(G4HCofThisEvent*HCE) {;}

G4bool LayerSD::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  //G4double t = aStep->GetTrack()->GetGlobalTime();
  if(edep==0.) return true;

  pHRec->RecordHit(layer,wire,edep);

  return true;
}

void LayerSD::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{;}

