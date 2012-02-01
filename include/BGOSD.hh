#ifndef BGOSD_h
#define BGOSD_h 1

#include "G4VSensitiveDetector.hh"
#include <fstream>
#include "gsl/gsl_histogram.h"

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class BGOSD : public G4VSensitiveDetector
{

  public:
      BGOSD(G4String name, gsl_histogram *hb);
      virtual ~BGOSD();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      gsl_histogram *h;
      double edep;
      int evtctr;
};




#endif

