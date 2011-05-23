#ifndef LayerSD_h
#define LayerSD_h 1

#include "HitRecorder.hh"
#include "G4VSensitiveDetector.hh"
#include <fstream>

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class LayerSD : public G4VSensitiveDetector
{

  public:
      LayerSD(G4String name, HitRecorder *hrec, int lay, int wi);
      virtual ~LayerSD();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      HitRecorder *pHRec;
      int layer;
      int wire;
};




#endif

