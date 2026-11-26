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

#ifndef CRYSTALDEFLECTOR_H
#define CRYSTALDEFLECTOR_H

#include "BDSAcceleratorComponent.hh"

#include "globals.hh" // geant4 types / globals

// forward declare things that are only needed as pointers
class G4LogicalVolume;
class G4Material;

/**
 * @brief Custom dipole component.
 * 
 * A pipe then dipole yoke (straight) with custom vacuum
 * chamber.
 * 
 * Note all comments are in Doxygen style for later use.
 *
 * @author Someone
 */

class CrystalDeflector: public BDSAcceleratorComponent
{
public:
  /// Construct as we wish.
  CrystalDeflector(G4String name,
                   G4double l,
                   G4double xsize,
                   G4double ysize,
                   G4double materialThickness,
                   G4double axisX,
                   G4double axisY,
                   G4double axisZ,
                   G4double horizontalWidth,
                   G4double bendingAngle,
                   G4String material,
                   G4String vacuumMaterial,
                   G4String crystalLattice,
                   G4String region,
                   G4String colour,
                   G4bool radiationModel,
                   G4double aper1,
                   G4double aper2,
                   G4String apertureType);
  
  virtual ~CrystalDeflector(); ///< Remember destructor has to be virtual.

protected:
  /// The minimum bit you should implement to build some geometry. This
  /// should reassign the member containerLogicalVolume and containerSolid
  /// and should be a shape with flat incoming / outgoing ends that's big
  /// enough to encompass the whole object.
  virtual void BuildContainerLogicalVolume();

  /// The more detailed build method. The default implementation in BDSIM
  /// calls BuildContainerLogicalVolume() and then sets some visualisation
  /// attributes for the container. We override it here to a) call the base
  /// class versioni then complete our own more detailed construction of
  /// the required geometry. Here, this function calls a series of functions
  /// to do small bits of construction.
  virtual void Build();

private:
  // Remove the default constructor as we only ever want to use our one.
  CrystalDeflector() = delete;

  // Split the construction into steps for our own sanity.
  void BuildCrystal();

  // Parameters that determine how CrystalDeflector is built that need to get
  // from the constructor to Build which is called later on.
  //G4double horizontalWidth;

  // We cache but don't own the materials


  G4double crystalXSize;
  G4double crystalYSize;
  G4double crystalMaterialThickness;
  G4double crystalAxisX;
  G4double crystalAxisY;
  G4double crystalAxisZ;
  G4double containerHorizontalWidth;
  G4double crystalBendingAngle;
  G4String crystalLattice;
  G4String crystalRegionName;

  G4Material* vacuum;
  G4Material* crystalMaterial;
  G4String  crystalcolour;
  G4bool crystalRadiationModel;
  G4double aper1;
  G4double aper2;
  G4String apertureType;
};

#endif
