#ifndef CZTSD_h
#define CZTSD_h 1

#include "G4VSensitiveDetector.hh"
#include <fstream>
#include <sqlite3.h>

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class CZTSD : public G4VSensitiveDetector
{

  public:
      CZTSD(G4String name, sqlite3 *dbin);
      virtual ~CZTSD();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      int priIdx;
      sqlite3 *db;
      double edep;
      sqlite3_stmt *stmt;
      int evtctr;
      int npix;
};




#endif

