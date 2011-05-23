#ifndef HitRecorder_H
#define HitRecorder_H 1
// non-geant approach to hit counting.
#include <fstream>
#include <iostream>
#include "globals.hh"
using namespace std;

const int PERIM_WIRE_NUM=-1;

class HitRecorder{
  private:
    ofstream outfile;
    int pri_pdgenc, pri_idx;  // primary particle info.
    double pri_kinetic_energy_MeV;
    double pri_vx,pri_vy,pri_vz,pri_x,pri_y,pri_z;
    double edep00,edep01;
    double edep10,edep11;
    double edep_perim;
  public:
    HitRecorder(char outfilename[]){
      pri_vx=pri_vy=pri_vz=pri_x=pri_y=pri_z=
        edep00=edep01=edep10=edep11=
        edep_perim=pri_pdgenc=pri_idx=0;
      outfile.open(outfilename);
    }

    void CutAndOutput(){
      if(edep_perim==0){ // perimeter veto got nothing
        // and only one wire on the inside got anything...
        if((edep00==0 && edep01==0 && edep10==0 && edep11>0)
            || (edep00==0 && edep01==0 && edep10>0 && edep11==0)
            || (edep00==0 && edep01>0 && edep10==0 && edep11==0)
            || (edep00>0 && edep01==0 && edep10==0 && edep11==0)){
          outfile << pri_idx << ' '
            << pri_pdgenc << ' '
            << pri_kinetic_energy_MeV << ' '
            << pri_vx << ' '
            << pri_vy << ' ' 
            << pri_vz << ' '
            << pri_x/m << ' '
            << pri_y/m << ' '
            << pri_z/m << "   "
            << edep00/keV << ' '
            << edep01/keV << ' '
            << edep10/keV << ' '
            << edep11/keV << endl;
        }
      }
    }

    ~HitRecorder(){
      CutAndOutput();
      outfile.close();
    }

    void NewPrimary(int idx, int pdgencoding, double kinetic_energy_MeV, 
        double vx, double vy, double vz, double x, double y, double z){
      CutAndOutput();
      pri_vx=vx;pri_vy=vy;pri_vz=vz;pri_x=x;pri_y=y;pri_z=z;
      pri_pdgenc=pdgencoding;pri_idx=idx;
      pri_kinetic_energy_MeV=kinetic_energy_MeV;
      edep00=edep01=edep10=edep11=edep_perim=0;
    }

    void RecordHit(int layer, int wire, double edep){
      //cout << layer << ' ' << wire << ' ' << edep << endl;
      if(wire==PERIM_WIRE_NUM)
        edep_perim+=edep;

      if(layer==0 && wire==0)
        edep00+=edep;
      if(layer==0 && wire==1)
        edep01+=edep;
      if(layer==1 && wire==0)
        edep10+=edep;
      if(layer==1 && wire==1)
        edep11+=edep;
    }
};


#endif
