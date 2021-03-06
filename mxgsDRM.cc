#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "HitRecorder.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4VisExecutive.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"
#include "G4UItcsh.hh"
//#include "LHEP.hh" // gone in 9.10?
//#include "QGSP_BERT_HP.hh" // takes a _long_ time to load the HP part, doesn't include low-energy EM stuff.
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <fcntl.h>
#include "gsl/gsl_histogram.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc, char** argv){
  int i,j,k;
  int nargs = 15;
  if(argc < nargs){
    cout << "usage:\n"
      << argv[0]
      << " interactive(0|1) priPDGID(22,11,-11,...) nPriPerE priStartDiskRad(m) priStartDiskRad0(m) theta(deg) phi(deg) Emin(MeV) Emax(MeV) numEnergies outEMin outEMax outNumE outputfileName  ...rest of arguments written as comment to output file...\n";
    return 1;
  }

  time_t t = time(0);
  struct tm tm = *localtime(&t);
  cout << "execution started: " << tm.tm_year << '/' << tm.tm_mon+1  << '/' << tm.tm_mday << ' '
    << tm.tm_hour << ':' << tm.tm_min << ':' << tm.tm_sec << endl;

  int interactive; sscanf(argv[1],"%d",&interactive);
  int priPDGID; sscanf(argv[2],"%d",&priPDGID);
  int nPriPerE; sscanf(argv[3],"%d",&nPriPerE);
  double priStartDiskRad; sscanf(argv[4],"%lf",&priStartDiskRad);
  double priStartDiskRad0; sscanf(argv[5],"%lf",&priStartDiskRad0);
  double priTheta; sscanf(argv[6],"%lf",&priTheta);
  double priPhi; sscanf(argv[7],"%lf",&priPhi);
  double priEMin; sscanf(argv[8],"%lf",&priEMin);
  double priEMax; sscanf(argv[9],"%lf",&priEMax);
  int priNumE; sscanf(argv[10],"%d",&priNumE);
  double outEMin; sscanf(argv[11],"%lf",&outEMin);
  double outEMax; sscanf(argv[12],"%lf",&outEMax);
  int outNumE; sscanf(argv[13],"%d",&outNumE);

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
  std::ofstream out(argv[14]);
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

  G4PhysListFactory factory;
  G4VModularPhysicsList* physlist = factory.GetReferencePhysList("QGSP_BERT_LIV");
  //physlist->RegisterPhysics(new G4OpticalPhysics());  // may be useful later? copy-pasted from qweak.jlab.org
  runManager->SetUserInitialization(physlist);

  // or, instead of using the factory,
  //runManager->SetUserInitialization(new QGSP_BERT_HP());

  // set mandatory user action class
  // targetX/Y/Z are the coordinates (m) of the center of MXGS in the colombus.gdml reference frame.
  ////// target for columbus.gdml.
  double targetX = 0.0 - 0.207;
  double targetY = 2.8297 + 0.1458/2.0 + 0.8630/2.0 - 0.046;
  double targetZ = -1.3996 + 1.1300/2.0 - 0.8630 + 0.5500/2.0;
  ////// target for asim.gdml.
  //double targetX = 0.0;
  //double targetY = 0.2;
  //double targetZ = -0.125;

  double targetDisp = 20;

  PrimaryGeneratorAction *pgen = new PrimaryGeneratorAction(priPDGID,
      priStartDiskRad,priStartDiskRad0,
      priTheta,priPhi,
      targetX,targetY,targetZ,targetDisp);
  runManager->SetUserAction(pgen);

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive("Errors"); //only init and error messages
  visManager->Initialize();
#endif

  // initialize G4 kernel
  runManager->Initialize();

  G4UImanager* UI=0;

  UI = G4UImanager::GetUIpointer();
  //  G4String command = "/run/beamOn 1";
  //
  std::ostringstream command;	//command

  for(j=0; j<priNumE; ++j){ // loop over primary energies.
    // reset the histograms.
    gsl_histogram_reset(hBGO);
    gsl_histogram_reset(hCZT);

    // set primary energy
    pgen->setPriEn(Epri[j]);

    if(interactive){
// Brant's original ui
//      G4UIsession* session = new G4UIterminal(new G4UItcsh);
//      session->SessionStart();
//      delete session;
// more advanced QT ui
#ifdef G4VIS_USE
		G4UIExecutive* ui = new G4UIExecutive(argc, argv);
		UI->ApplyCommand("/control/execute vis.mac"); 
		ui->SessionStart();   
#endif 
    }

    // main event loop for this primary energy
    command.str(""); //empty command
    command << "/run/beamOn " << nPriPerE << std::flush;
    UI->ApplyCommand(command.str());
    //    for(k=0; k<nPriPerE; ++k){
    //  UI->ApplyCommand(command);
    // }

    // write histogram line(s) to matrices.
    out << "BGO";
    for(k=0; k<outNumE; ++k){
      out << ' ' << gsl_histogram_get(hBGO,k);
    } out << endl << flush;
    out << "CZT";
    for(k=0; k<outNumE; ++k){
      out << ' ' << gsl_histogram_get(hCZT,k);
    } out << endl << flush;
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

  t = time(0);
  tm = *localtime(&t);
  cout << "execution finished: " << tm.tm_year+1900 << '/' << tm.tm_mon+1  << '/' << tm.tm_mday << ' '
    << tm.tm_hour << ':' << tm.tm_min << ':' << tm.tm_sec << endl;

  return 0;
}

/*
/vis/open OGL
/vis/scene/create
/vis/scene/add/volume
/vis/sceneHandler/attach
/vis/viewer/set/projection o
/vis/viewer/set/viewpointThetaPhi 90 0
/vis/viewer/panTo 1.4226 3.7341
/vis/viewer/set/viewpointThetaPhi 140 40
/vis/scene/add/trajectories
/vis/scene/add/hits
/vis/viewer/zoom 2
/vis/viewer/zoom 2
/vis/viewer/zoom 2



/vis/geometry/set/visibility asim -1 0
/vis/geometry/set/visibility czt -1 1
/vis/geometry/set/colour cztKapton -1 red


/vis/viewer/set/projection perspective 20 deg
/vis/viewer/set/viewpointThetaPhi 90 0
/vis/viewer/panTo 1.4226 3.7341
/vis/viewer/set/viewpointThetaPhi 140 40
/vis/scene/add/trajectories
/vis/scene/add/hits
/vis/viewer/zoom 2
/vis/viewer/zoom 2
/vis/viewer/zoom 2

/vis/geometry/set/visibility asim 1 0
/vis/geometry/set/visibility mxgs_world 1 0
/vis/geometry/set/visibility codedMask 1 0

/vis/geometry/set/colour codedMask -1 red

/vis/viewer/set/style surface
/vis/viewer/set/style wireframe

/run/beamOn 1


#set to aluminum but with much lower density than normal aluminum
/vis/geometry/set/colour Columbus3EPF_8_Logical_147892368 -1 blue


/vis/viewer/set/background white

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

DUMP VRML FILE:
/vis/sceneHandler/create VRML2FILE
/vis/viewer/create
/vis/drawVolume
/vis/viewer/flush

/vis/viewer/set/style surface

CHECK CODED MASK
/vis/open OGLIX
/vis/scene/create
/vis/scene/add/volume
/vis/sceneHandler/attach
/vis/viewer/flush

*/
