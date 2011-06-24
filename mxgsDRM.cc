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
// #include "SamplingGeometry.hh"

#include <sqlite3.h>

// ideas:
//
int main(int argc, char** argv)
{
    int nEvt=0, targetN=0;
    int action=0;
    sqlite3 *db = 0;
    int errcode=0;
    
    if(argc < 2){
        cout << "usage:\n" 
            << argv[0] << "action(0|1|2|3|4|5)\n";
        return 1;
    }
    sscanf(argv[1],"%d",&action);

    if((action==0 || action==1) && argc != 5){ // interactive session
      printf("usage for action==0,1: %s (0|1) dbfn nEvents targetN\n",argv[0]);
      return 1;
    }else{
      errcode = sqlite3_open(argv[2],&db);
      if(errcode != SQLITE_OK){
        printf("wtf! %s\n",sqlite3_errmsg(db));
        return 0;
      }
      sqlite3_exec(db,"pragma synchronous = 0;",0,0,0);
      sscanf(argv[3],"%d",&nEvt);
      sscanf(argv[4],"%d",&targetN);
    }

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
    DetectorConstruction *world = new DetectorConstruction(db);
    //world->RegisterParallelWorld(new SamplingGeometry(0,0,0));
    runManager->SetUserInitialization(world);
    runManager->SetUserInitialization(new LHEP());

    // set mandatory user action class
    PrimaryGeneratorAction *pgen = new PrimaryGeneratorAction(db);
    runManager->SetUserAction(pgen);

    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    // initialize G4 kernel
    runManager->Initialize();


    if(action==0){
      //// Use for automatic run.
      G4UImanager* UI = G4UImanager::GetUIpointer();
      G4String command = "/run/beamOn 1";
      for(int i=0; i<nEvt; ++i){
        UI->ApplyCommand(command); 
      }
    }

    if(action==1){
      // use for interactive session
      //bash$ export G4VIS_USE_OPENGLX=1
      G4UIsession* session = new G4UIterminal(new G4UItcsh);
      session->SessionStart();
      delete session;
      // see sample ui commands below
    }

    if(action==2){
      //// dump DAWN file of geometry.
      G4UImanager* UI = G4UImanager::GetUIpointer();
      G4String command = "/vis/open DAWNFILE";
      UI->ApplyCommand(command); 
      command = "/vis/drawVolume";
      UI->ApplyCommand(command); 
      command = "/vis/viewer/flush";
      UI->ApplyCommand(command); 
      ////////
    }

    if(action==3){
      //// dump VRML of geometry.
      G4UImanager* UI = G4UImanager::GetUIpointer();
      G4String command = "/vis/sceneHandler/create VRML2FILE";
      UI->ApplyCommand(command); 
      command = "/vis/viewer/create";
      UI->ApplyCommand(command); 
      //command = "/vis/viewer/set/style surface";
      //UI->ApplyCommand(command); 
      command = "/vis/drawVolume";
      UI->ApplyCommand(command); 
      command = "/vis/viewer/flush";
      UI->ApplyCommand(command); 
      ////////
    }

    if(action==4){
      //// dump Raytraced .jpeg of geometry.
      G4UImanager* UI = G4UImanager::GetUIpointer();
      G4String command = "/vis/open RayTracer";
      UI->ApplyCommand(command); 
      command = "/vis/rayTracer/lightDirection 0.3 0.3 -1";
      UI->ApplyCommand(command); 
      command = "/vis/viewer/set/viewpointThetaPhi 70 30";
      UI->ApplyCommand(command); 
      command = "/vis/drawVolume";
      UI->ApplyCommand(command); 
      command = "/vis/viewer/flush";
      UI->ApplyCommand(command); 
      ////////
    }

    if(action==5){
      //// check geometry
      G4UImanager* UI = G4UImanager::GetUIpointer();
      G4String command = "/geometry/test/grid_test true";
      UI->ApplyCommand(command); 
      ////////
    }
    
    // job termination
    sqlite3_close(db);
    delete runManager;
    delete visManager;
    return 0;
}


/*
/vis/open OGL
/vis/scene/create
/vis/scene/add/volume
/vis/sceneHandler/attach
/vis/viewer/set/viewpointThetaPhi 90 0
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
