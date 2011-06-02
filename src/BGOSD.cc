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

BGOSD::BGOSD(G4String name, sqlite3 *dbin)
:G4VSensitiveDetector(name)
{
  db=dbin;
  stmt=0;
  evtctr=0;
}

BGOSD::~BGOSD(){;}

void BGOSD::Initialize(G4HCofThisEvent*HCE){
  priIdx = ((EventInfo*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation())->GetIdx();
  edep=0;
}

G4bool BGOSD::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{
  G4double thisEdep = aStep->GetTotalEnergyDeposit();
  if(thisEdep>0){
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle tch = preStepPoint->GetTouchableHandle();

    double x=preStepPoint->GetPosition().x();
    double y=preStepPoint->GetPosition().y();
    double z=preStepPoint->GetPosition().z();
    double t=preStepPoint->GetGlobalTime();

    int i3 = tch->GetCopyNumber(1);
    int i6 = tch->GetCopyNumber(3);
    int i12 = preStepPoint->GetPosition().y()>0 ? 1 : 0;

    int idx = i3+3*i6+6*i12;

    sqlite3_prepare_v2(db,"INSERT INTO hitsBGO VALUES (?,?,?,?,?,?,?)",-1,&stmt,0);
    sqlite3_bind_int(stmt,1,evtctr);
    sqlite3_bind_double(stmt,2,thisEdep/MeV);
    sqlite3_bind_int(stmt,3,idx);
    sqlite3_bind_double(stmt,4,x/m);
    sqlite3_bind_double(stmt,5,y/m);
    sqlite3_bind_double(stmt,6,z/m);
    sqlite3_bind_double(stmt,7,t/ns);
    while(sqlite3_step(stmt) != SQLITE_DONE){
      cout << "sqlite3_step failed? retrying...\n";
    };

    edep += thisEdep;
  }
  return true;
}

void BGOSD::EndOfEvent(G4HCofThisEvent* /*HCE*/){
  if(edep>0){
    sqlite3_prepare_v2(db,"INSERT INTO evtsBGO VALUES (?,?,?)",-1,&stmt,0);
    sqlite3_bind_int(stmt,1,evtctr);
    sqlite3_bind_int(stmt,2,priIdx);
    sqlite3_bind_double(stmt,3,edep/MeV);
    while(sqlite3_step(stmt) != SQLITE_DONE){
      cout << "sqlite3_step failed? retrying...\n";
    };
    sqlite3_finalize(stmt);
  }
  ++evtctr;
}

