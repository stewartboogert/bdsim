/* 
G4ChannelingFastSimModelBDSim is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License
http://www.gnu.org/licenses/

G4ChannelingFastSimModelBDSim is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*/

//Author: Alexei Sytov alexei.sytov@cern.ch

#ifndef CRYSTALDEFLECTORCONSTRUCTOR_H
#define CRYSTALDEFLECTORCONSTRUCTOR_H

#include "BDSComponentConstructor.hh"

class CrystalDeflectorConstructor: public BDSComponentConstructor
{
public:
  CrystalDeflectorConstructor(){;}
  virtual ~CrystalDeflectorConstructor(){;}


	virtual BDSAcceleratorComponent* Construct(GMAD::Element const* elementIn,
						   GMAD::Element const* prevElementIn,
						   GMAD::Element const* nextElementIn,
						   const BDSBeamlineIntegral& integral);

  virtual BDSAcceleratorComponent* Construct(GMAD::Element const* elementIn,
					     GMAD::Element const* prevElementIn,
					     GMAD::Element const* nextElementIn,
					     G4double currentArcLengthIn,
					     G4double brhoIn,
					     G4double beta0In);
};

#endif
