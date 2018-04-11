/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2018.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "BDSPhysicsEMDissociation.hh"
#include "G4EMDissociation.hh"
#include "G4EMDissociationCrossSection.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4IonConstructor.hh"
#include "G4ProcessManager.hh"

BDSPhysicsEMDissociation::BDSPhysicsEMDissociation():
  G4VPhysicsConstructor("BDSPhysicsEMDissociation")
{;}

void BDSPhysicsEMDissociation::ConstructParticle()
{
  G4Gamma::Gamma();
  G4GenericIon::GenericIon();

  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void BDSPhysicsEMDissociation::ConstructProcess()
{
  if (Activated())
    {return;}

  G4HadronInelasticProcess* inelProcIon = new G4HadronInelasticProcess("ionInelastic", G4GenericIon::GenericIon());

  G4EMDissociationCrossSection* crossSectionData = new G4EMDissociationCrossSection();
  inelProcIon->AddDataSet(crossSectionData);

  G4EMDissociation* emdModel = new G4EMDissociation();
  inelProcIon->RegisterMe(emdModel);

  G4ProcessManager* pmanager = G4GenericIon::GenericIon()->GetProcessManager();
  pmanager->AddDiscreteProcess(inelProcIon);
  
  SetActivated();
}