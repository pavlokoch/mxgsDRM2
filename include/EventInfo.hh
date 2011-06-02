#ifndef EventInfo_h
#define EventInfo_h 1

#include "G4VUserEventInformation.hh"

using namespace std;

class EventInfo : public G4VUserEventInformation
{
  private:
      int pIdx;

  public:
      EventInfo(int pIdxIn){pIdx=pIdxIn;}
      int GetIdx(){return pIdx;}
      inline void Print()const{};
};

#endif

