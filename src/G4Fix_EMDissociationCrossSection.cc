/* 
Beam Delivery Simulation (BDSIM) Copyright (C) BDSIM Collaboration, 2001 - 2026.

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
#include "G4Fix_EMDissociationCrossSection.hh"

#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

G4double G4Fix_EMDissociationCrossSection::GetElementCrossSection(const G4DynamicParticle* dp,
                                                                  G4int                    Z,
                                                                  const G4Material*        mat)
{
  const G4double ekin = dp->GetKineticEnergy();
  const G4double mass = dp->GetMass();

  // exact beta from the kinetic energy: beta = sqrt(t(t+2))/(t+1), with t = Ekin/m
  const G4double t = ekin / mass;
  const G4double b = std::sqrt(t*(t + 2.0)) / (t + 1.0);

  // as in G4EMDissociationCrossSection::GetElementCrossSection, with the exact beta
  const G4ParticleDefinition* def = dp->GetDefinition();
  const G4double AP = def->GetBaryonNumber();
  const G4double ZP = def->GetPDGCharge();
  const G4double AT = G4NistManager::Instance()->GetAtomicMassAmu(Z);
  const G4double ZT = (G4double)Z;
  const G4double bmin = spectrum.GetClosestApproach(AP, ZP, AT, ZT, b);
  G4PhysicsFreeVector* xsP = GetCrossSectionForProjectile(AP, ZP, AT, ZT, b, bmin);
  G4PhysicsFreeVector* xsT = GetCrossSectionForTarget(AP, ZP, AT, ZT, b, bmin);
  const G4double xs = (*xsP)[0] + (*xsP)[1] + (*xsT)[0] + (*xsT)[1];
  delete xsP;
  delete xsT;
  return xs;
}
