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

#include "CrystalDeflector.hh"
#include "CrystalDeflectorConstructor.hh"

#include "BDSAcceleratorComponent.hh"
#include "BDSBeamlineIntegral.hh"
#include "parser/element.h" // for GMAD::Element
#include "BDSUtilities.hh"
#include "BDSGlobalConstants.hh"
#include "BDSException.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Types.hh" // for G4String etc

BDSAcceleratorComponent* CrystalDeflectorConstructor::Construct(GMAD::Element const* element,
                                                                 GMAD::Element const* prevElement,
                                                                 GMAD::Element const* nextElement,
                                                                 const BDSBeamlineIntegral& integral)
{
	// Pull needed values from the integral
	G4double currentArcLength = integral.arcLength;
	G4double brho = integral.designParticle.BRho();
	G4double beta0 = integral.designParticle.Beta();

	// Call the actual implementation
	return Construct(element, prevElement, nextElement, currentArcLength, brho, beta0);
}

BDSAcceleratorComponent* CrystalDeflectorConstructor::Construct(GMAD::Element const* element,
						       GMAD::Element const* /*prevElement*/,
						       GMAD::Element const* /*nextElement*/,
						       G4double /*currentArcLength*/,
						       G4double /*brhoIn*/,
						       G4double /*beta0In*/)
{

    G4double horizontalWidth = element->horizontalWidth;
    if(horizontalWidth<0.001){horizontalWidth = 1;}

    G4String vacuumMaterial = element->vacuumMaterial;
    if(vacuumMaterial == "")
    {vacuumMaterial=BDSGlobalConstants::Instance()->WorldMaterial();}

    G4String material = element->material;
    if(material == ""){material="G4_Si";}

    // the map is supplied at run time so we must check it actually has the key
    // to avoid throwing an exception.
  G4String params = G4String(element->userParameters);
    // use function from BDSUtilities to process user params string into
    // map of strings and values.
    std::map<G4String, G4String> map = BDS::GetUserParametersMap(params);

    G4String region;
    auto regionSearch = map.find("crystalRegion");
    if (regionSearch != map.end())
    {region = regionSearch->second;}
    else
    {region = "crystal";}

    G4double bendingAngle=0;
    auto radiusSearch = map.find("crystalBendingAngle");
    if (radiusSearch != map.end())
    {bendingAngle = std::stod(radiusSearch->second);}

    G4String crystalLattice;
    auto latticeSearch = map.find("crystalLattice");
    if (latticeSearch != map.end())
    {crystalLattice = latticeSearch->second;}
    else
    {crystalLattice = "(111)";}

    G4String colour;
    auto colourSearch = map.find("colour");
    if (colourSearch != map.end())
    {colour = colourSearch->second;}
    else
    {colour = "magenta";}

    G4bool radiationModel=false;
    auto radiationSearch = map.find("radiationModel");
    if (radiationSearch != map.end())
    {if(radiationSearch->second=="true"){radiationModel=true;};}

  // Here we pull out the information we want from the parser element and construct
  // one of our components and return it. BDSIM will delete this at the end.
  // If you don't need to use all the parameters, use an /*inline*/ comment to
  // avoid a compiler warning about an unused variable.

  BDSAcceleratorComponent* crystal = new CrystalDeflector(element->name,
                                                          element->l*CLHEP::m,
                                                          element->xsize*CLHEP::m,
                                                          element->ysize*CLHEP::m,
                                                          element->materialThickness*CLHEP::m,
                                                          element->axisX,
                                                          element->axisY,
                                                          element->axisZ,
                                                          horizontalWidth*CLHEP::m,
                                                          bendingAngle,
                                                          material,
                                                          vacuumMaterial,
                                                          crystalLattice,
                                                          region,
                                                          colour,
                                                          radiationModel,
                                                          element->aper1,
                                                          element->aper2,
                                                          element->apertureType);
  return crystal;
}
