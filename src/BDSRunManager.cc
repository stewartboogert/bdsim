/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/
// G.A.Blair, Royal Holloway Univ of London. 1.07.02
// This class alters one part of the default G4RunManager to
// avoid closing the geometry after each event - this was a significant
// time overhead. 
// >>>Whether there are other effects of this needs to be checked <<<
//
// SPM: Altered BeamOn function to account for Placet synchronisation
//

#include "BDSGlobalConstants.hh"
#include "G4Timer.hh"

#include "BDSRunManager.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"


BDSRunManager* BDSRunManager::fRunManager = 0;

BDSRunManager* BDSRunManager::GetRunManager()
{ return fRunManager; }

BDSRunManager::BDSRunManager(){ fRunManager = this;}

BDSRunManager::~BDSRunManager(){ G4cout << "BDSRunManager deleting..." << G4endl; }

void BDSRunManager::BeamOn(G4int n_event,const char* macroFile,G4int n_select)
{
  G4bool cond = ConfirmBeamOnCondition();
  G4StackManager* SM;
  SM = G4EventManager::GetEventManager()->GetStackManager();
  if(cond)
  {
    numberOfEventToBeProcessed = n_event;
    RunInitialization();
    if(n_event>0) DoEventLoop(n_event,macroFile,n_select);
    RunTermination();
    while(!BDSGlobals->holdingQueue.empty()){
      BDSGlobals->setReadFromStack(true);
      SM->ClearPostponeStack();

      RunInitialization();
//      DoEventLoop(n_event,macroFile,n_select);
      DoEventLoop(BDSGlobals->holdingQueue.size(),macroFile,n_select);
      RunTermination();

      BDSGlobals->setReadFromStack(false);
    }
  }
}

void BDSRunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  //G4StateManager* stateManager = G4StateManager::GetStateManager();

  if(verboseLevel>0) 
  { timer->Start(); }

  G4String msg;
  if(macroFile!=0)
  { 
    if(n_select<0) n_select = n_event;
    msg = "/control/execute ";
    msg += macroFile;
  }
  else
  { n_select = -1; }

  G4int i_event;
  for( i_event=0; i_event<n_event; i_event++ )
  {


    //    stateManager->SetNewState(EventProc);

    currentEvent = GenerateEvent(i_event);

    if(currentEvent == NULL) G4cerr<<__FILE__<<" : "<<__LINE__<<"Event generation failed "<<G4endl;

    eventManager->ProcessOneEvent(currentEvent);

    AnalyzeEvent(currentEvent);

    if(i_event<n_select) G4UImanager::GetUIpointer()->ApplyCommand(msg);

    // gab: the following commented out - this is the only
    // difference between this and the normal G4 class
    // (significant time saving...)
    //stateManager->SetNewState(GeomClosed);

    StackPreviousEvent(currentEvent);
    currentEvent = 0;
    if(runAborted) break;
  }

  if(verboseLevel>0)
  {
    timer->Stop();
    G4cout << "Run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
    { G4cout << "  Run Aborted after " << i_event << " events processed." << G4endl; }
    else
    { G4cout << "  Number of events processed : " << n_event << G4endl; }
    G4cout << "  "  << *timer << G4endl;
  }
}


