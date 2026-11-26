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
#ifndef G4FIX_EMDISSOCIATIONCROSSSECTION_H
#define G4FIX_EMDISSOCIATIONCROSSSECTION_H

#include "G4EMDissociationCrossSection.hh"
#include "G4EMDissociationSpectrum.hh"
#include "globals.hh"

class G4DynamicParticle;
class G4Material;

/**
 * @brief Patched G4EMDissociationCrossSection that works at high Lorentz factors.
 *
 * G4FIX: temporary patch of a Geant4 bug - remove this class (and use
 * G4EMDissociationCrossSection again) once the fix is in a Geant4 release.
 *
 * G4EMDissociationCrossSection::GetElementCrossSection takes beta from
 * G4DynamicParticle::GetBeta(), which returns exactly 1 when the Lorentz factor
 * is above ~1000 (e.g. above ~180 TeV kinetic energy for Os-192, or for Pb-208
 * at LHC energies). The Lorentz factor 1/sqrt(1 - beta^2) is then infinite and
 * the cross section NaN, so EMD never happens. Here, beta is computed exactly
 * from the kinetic energy and the mass when the Lorentz factor is above 1000,
 * and the cross section is otherwise computed in the same way as in Geant4.
 * Below that, the Geant4 implementation is used unchanged.
 *
 * Only the cross section is affected: G4EMDissociation (the model) takes beta
 * from the four-momentum and works correctly.
 */

class G4Fix_EMDissociationCrossSection: public G4EMDissociationCrossSection
{
public:
  G4Fix_EMDissociationCrossSection() = default;
  virtual ~G4Fix_EMDissociationCrossSection() = default;

  virtual G4double GetElementCrossSection(const G4DynamicParticle* dp,
                                          G4int                    Z,
                                          const G4Material*        mat) override;

private:
  G4EMDissociationSpectrum spectrum;
};

#endif
