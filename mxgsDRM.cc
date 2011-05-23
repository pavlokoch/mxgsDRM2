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

int main(int argc, char** argv)
{
    char infile[500],hitsfilename[500];
    double planewidth,planerad;
    int nlines=0;
    if(argc != 6){
        cout << "usage:\n" 
            << argv[0] << " inputfilename hitsoutfilename nlines planewidth(m) planerad(m)\n";
        return 1;
    }else{
        strcpy(infile,argv[1]);
        strcpy(hitsfilename,argv[2]);
        sscanf(argv[3],"%d",&nlines);
        sscanf(argv[4],"%lf",&planewidth);
        sscanf(argv[5],"%lf",&planerad);
    }

    cout << "called with parameters:\n"
      << "    inputfilename: " << infile << endl
      << "    hitsoutfilename: " << hitsfilename<< endl
      << "    nlines: " << nlines << endl
      << "    planewidth(m): " << planewidth << endl
      << "    planerad(m): " << planerad << endl;

    ///random number seeding code stolen from LISARunAction.cc
    long seeds[2];
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    // G4cout << "seed1: " << seeds[0] << "; seed2: " << seeds[1] << G4endl;
    CLHEP::HepRandom::setTheSeeds(seeds);
    CLHEP::HepRandom::showEngineStatus();

    HitRecorder hrec(hitsfilename);

    // construct the default run manager
    G4RunManager* runManager = new G4RunManager;

    // set mandatory initialization classes
    runManager->SetUserInitialization(new DetectorConstruction(&hrec));
    runManager->SetUserInitialization(new LHEP);

    // set mandatory user action class ; -0.5 is the approximate zoffset to the "center" of the detector.
    PrimaryGeneratorAction *pgen = new PrimaryGeneratorAction(infile,&hrec,planewidth,planerad,-0.5);
    runManager->SetUserAction(pgen);

    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    // initialize G4 kernel
    runManager->Initialize();


    ////// Use for automatic run.
    //G4UImanager* UI = G4UImanager::GetUIpointer();
    //G4String command = "/run/beamOn 100000";
    //for(int i=0; i<nlines; ++i){
    //  pgen->NextLine();
    //  UI->ApplyCommand(command); 
    //}
    //////

    //// use for interactive session
    ////bash$ export G4VIS_USE_OPENGLX=1
    //pgen->NextLine();
    //G4UIsession* session = new G4UIterminal(new G4UItcsh);
    //session->SessionStart();
    //delete session;
    //// see sample ui commands below

    //// use dump VRML of geometry.
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

    ////// check geometry
    //G4UImanager* UI = G4UImanager::GetUIpointer();
    //G4String command = "/geometry/test/grid_test true";
    //UI->ApplyCommand(command); 
    //////////
    
    // job termination
    delete runManager;
    delete visManager;
    return 0;
}


/*
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
