/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/

//==========================================================
//==========================================================

#include "BDSExecOptions.hh"
#include "BDSGlobalConstants.hh" 
#include "BDSRunAction.hh"
#include "BDSRunManager.hh"
#include "BDSPhotonCounter.hh"

#include "G4Run.hh"
//#include "G4UImanager.hh"
//#include "G4VVisManager.hh"
#include "G4ios.hh"

#include <fstream>

#include "BDSTrackingFIFO.hh"

#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
#include <unistd.h> // for sleep, will not work on windows?
#endif


//==========================================================

BDSRunAction::BDSRunAction()
{
}

//==========================================================

BDSRunAction::~BDSRunAction()
{}

//==========================================================

void BDSRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //  if (G4VVisManager::GetConcreteInstance())
  //    {
      //      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      //  UI->ApplyCommand("/vis/scene/notifyHandlers");
  //    } 


}

//==========================================================

void BDSRunAction::EndOfRunAction(const G4Run* aRun)
{
  //Do the fifo at the end of the run.
  BDSTrackingFIFO* fifo = new BDSTrackingFIFO();
  fifo->doFifo();
  G4cout << "### Run " << aRun->GetRunID() << " end." << G4endl;
  G4cout << "Number of optical photons produced in run = " << BDSPhotonCounter::Instance()->nPhotons() << G4endl;
  G4cout << "Total energy of optical photons produced in run = " << BDSPhotonCounter::Instance()->energy()/CLHEP::GeV << " GeV" << G4endl;
}
//==========================================================
