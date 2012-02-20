#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "HitRecorder.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4VisExecutive.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"
#include "G4UItcsh.hh"
#include "LHEP.hh"
#include "QGSP_BERT_HP.hh"
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <fcntl.h>
#include "gsl/gsl_histogram.h"

int main(int argc, char** argv){
  int i,j,k;
  int nargs = 14;
  if(argc < nargs){
    cout << "usage:\n" 
      << argv[0] 
      << " priPDGID(22,11,-11,...) nPriPerE priStartDiskRad(m) priStartDiskRad0(m) theta(deg) phi(deg) Emin(MeV) Emax(MeV) numEnergies outEMin outEMax outNumE outputfileName  ...rest of arguments written as comment to output file...\n";
    return 1;
  }

  int priPDGID; sscanf(argv[1],"%d",&priPDGID);
  int nPriPerE; sscanf(argv[2],"%d",&nPriPerE);
  double priStartDiskRad; sscanf(argv[3],"%lf",&priStartDiskRad);
  double priStartDiskRad0; sscanf(argv[4],"%lf",&priStartDiskRad0);
  double priTheta; sscanf(argv[5],"%lf",&priTheta);
  double priPhi; sscanf(argv[6],"%lf",&priPhi);
  double priEMin; sscanf(argv[7],"%lf",&priEMin);
  double priEMax; sscanf(argv[8],"%lf",&priEMax);
  int priNumE; sscanf(argv[9],"%d",&priNumE);
  double outEMin; sscanf(argv[10],"%lf",&outEMin);
  double outEMax; sscanf(argv[11],"%lf",&outEMax);
  int outNumE; sscanf(argv[12],"%d",&outNumE);

  double *Epri = new double[priNumE];
  for(i=0; i<priNumE; ++i){
    if(i==0){
      Epri[i] = priEMin;
    }else{
      Epri[i] = pow(10,log10(priEMin) + ((double)i)/(priNumE-1)*(log10(priEMax)-log10(priEMin)));
    }
  }
  double *outBins = new double[outNumE+1];
  for(i=0; i<outNumE+1; ++i){
    outBins[i] = pow(10,log10(outEMin) + ((double)i)/outNumE*(log10(outEMax)-log10(outEMin)));
  }
  gsl_histogram *hBGO = gsl_histogram_alloc(outNumE);
  gsl_histogram_set_ranges(hBGO,outBins,outNumE+1);
  gsl_histogram *hCZT = gsl_histogram_alloc(outNumE);
  gsl_histogram_set_ranges(hCZT,outBins,outNumE+1);

  // start output process
  std::ofstream out(argv[13]);
  out << "#";
  for(i=nargs; i<argc; ++i){ out << ' ' << argv[i]; }
  out << endl;
  out << "#";
  for(i=0; i<nargs; ++i){ out << ' ' << argv[i]; }
  out << endl;
  out << "# primary energies:";
  for(i=0; i<priNumE; ++i){ out << ' ' << Epri[i]; }
  out << endl;
  out << "# output histogram bin boundaries:";
  for(i=0; i<outNumE+1; ++i){ out << ' ' << outBins[i]; }
  out << endl;
  out << "# One row per primary energy, one column per output histogram bin.\n";


  ///random number seeding code stolen from LISARunAction.cc
  long seeds[2];
  int rfp = open("/dev/random",O_RDONLY);
  long seed1,seed2;
  int x1 = read(rfp,&seed1,sizeof(seed1));
  x1 = read(rfp,&seed2,sizeof(seed2));
  close(rfp);
  seeds[0] = seed1;
  seeds[1] = seed2;
  G4cout << "seed1: " << seeds[0] << "; seed2: " << seeds[1] << G4endl;
  CLHEP::HepRandom::setTheSeeds(seeds);
  CLHEP::HepRandom::showEngineStatus();


  // construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  DetectorConstruction *world = new DetectorConstruction(hBGO,hCZT);
  runManager->SetUserInitialization(world);
  runManager->SetUserInitialization(new LHEP());

  // set mandatory user action class
  // targetY/Y/Z are the coordinates (m) of the center of MXGS in the colombus.gdml reference frame.
  ////// target for columbus.gdml.
  double targetX = 0.0;
  double targetY = 2.8297 + 0.1458/2.0 + 0.8630/2.0;
  double targetZ = -1.3996 + 1.1300/2.0 - 0.8630 + 0.5500/2.0;
  ////// target for asim.gdml.
  //double targetX = 0.0;
  //double targetY = 0.2;
  //double targetZ = -0.125;

  double targetDisp = 100;

  PrimaryGeneratorAction *pgen = new PrimaryGeneratorAction(priPDGID,
      priStartDiskRad,priStartDiskRad0,
      priTheta,priPhi,
      targetX,targetY,targetZ,targetDisp);
  runManager->SetUserAction(pgen);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // initialize G4 kernel
  runManager->Initialize();

  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4String command = "/run/beamOn 1";


  for(j=0; j<priNumE; ++j){ // loop over primary energies.
    // reset the histograms.
    gsl_histogram_reset(hBGO);
    gsl_histogram_reset(hCZT);

    // set primary energy
    pgen->setPriEn(Epri[j]);

    //// use for interactive session to doublecheck things via command line.
    //G4UIsession* session = new G4UIterminal(new G4UItcsh);
    //session->SessionStart();
    //delete session;

    // main event loop for this primary energy
    for(k=0; k<nPriPerE; ++k){
      UI->ApplyCommand(command); 
    }

    // write histogram line(s) to matrices.
    out << "BGO";
    for(k=0; k<outNumE; ++k){
      out << ' ' << gsl_histogram_get(hBGO,k);
    } out << endl;
    out << "CZT";
    for(k=0; k<outNumE; ++k){
      out << ' ' << gsl_histogram_get(hCZT,k);
    } out << endl;
  }


  // close matrix file.
  out.close();

  // job termination
  delete[] Epri;
  delete[] outBins;
  gsl_histogram_free(hBGO);
  gsl_histogram_free(hCZT);
  delete runManager;
  delete visManager;
  return 0;
}

/*
/vis/open OGL
/vis/scene/create
/vis/scene/add/volume
/vis/sceneHandler/attach
/vis/viewer/set/projection perspective 20 deg
/vis/viewer/set/viewpointThetaPhi 90 0
/vis/viewer/panTo 1.4226 3.7341
/vis/scene/add/trajectories
/vis/scene/add/hits
/run/beamOn 1

/vis/viewer/set/viewpointThetaPhi 90 0
/vis/viewer/panTo 1.4226 3.7341

/vis/open OGLIX
/vis/scene/create
/vis/scene/add/volume
/vis/sceneHandler/attach
/vis/viewer/flush
/vis/viewer/set/viewpointThetaPhi 90 0
/vis/viewer/zoom 30
/vis/scene/add/trajectories
/vis/scene/add/hits
/run/beamOn 1
*/
