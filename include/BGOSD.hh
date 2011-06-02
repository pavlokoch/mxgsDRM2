#ifndef BGOSD_h
#define BGOSD_h 1

#include "G4VSensitiveDetector.hh"
#include <fstream>
#include <sqlite3.h>

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class BGOSD : public G4VSensitiveDetector
{

  public:
      BGOSD(G4String name, sqlite3 *dbin);
      virtual ~BGOSD();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      sqlite3 *db;
      double edep;
      sqlite3_stmt *stmt;
      int evtctr;
      int priIdx;
};




#endif

