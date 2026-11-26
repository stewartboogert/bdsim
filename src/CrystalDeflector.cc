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

#include "BDSAcceleratorComponent.hh"
#include "BDSBeamPipe.hh"
#include "BDSBeamPipeFactory.hh"
#include "BDSBeamPipeInfo.hh"
#include "BDSBeamPipeType.hh"
#include "BDSColours.hh"
#include "BDSExtent.hh"
#include "BDSFieldBuilder.hh"
#include "BDSFieldInfo.hh"
#include "BDSFieldType.hh"
#include "BDSIntegratorType.hh"
#include "BDSMagnetGeometryType.hh"
#include "BDSMagnetOuter.hh"
#include "BDSMagnetOuterFactory.hh"
#include "BDSMagnetOuterInfo.hh"
#include "BDSMagnetStrength.hh"
#include "BDSMagnetType.hh"
#include "BDSMaterials.hh"
#include "BDSSamplerPlane.hh"
#include "BDSSamplerRegistry.hh"
#include "BDSSDSampler.hh" // so we can convert BDSSamplerSD* to G4VSensitiveDetector*
#include "BDSSDManager.hh"
#include "BDSSimpleComponent.hh"
#include "BDSUtilities.hh"

#include "globals.hh" // geant4 globals / types
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

#include "G4ChannelingFastSimModel.hh"

CrystalDeflector::CrystalDeflector(G4String name,
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
                                   G4double aper1In,
                                   G4double aper2In,
                                   G4String apertureTypeIn):
  BDSAcceleratorComponent(name,l, 0, "crystaldeflector"),
  crystalXSize(xsize),
  crystalYSize(ysize),
  crystalMaterialThickness(materialThickness),
  crystalAxisX(axisX),
  crystalAxisY(axisY),
  crystalAxisZ(axisZ),
  containerHorizontalWidth(horizontalWidth),
  crystalBendingAngle(bendingAngle),
  crystalLattice(crystalLattice),
  crystalRegionName(region),
  crystalcolour(colour),
  crystalRadiationModel(radiationModel),
  aper1(aper1In),
  aper2(aper2In),
  apertureType(apertureTypeIn)
{

  G4cout  << "name = "<< name << G4endl ;
  G4cout  << "l = "<< l << G4endl ;
  G4cout  << "xsize = "<< xsize << G4endl ;
  G4cout  << "ysize = "<< ysize << G4endl ;
  G4cout  << "materialThickness = "<< materialThickness << G4endl ;
  G4cout  << "axisX = "<< axisX << G4endl ;
  G4cout  << "axisY = "<< axisY << G4endl ;
  G4cout  << "axisZ = "<< axisZ << G4endl ;
  G4cout  << "horizontalWidth = "<< horizontalWidth << G4endl ;
  G4cout  << "bendingAngle = "<< bendingAngle << G4endl ;
  G4cout  << "material = "<< material << G4endl ;
  G4cout  << "vacuumMaterial = "<< vacuumMaterial << G4endl ;
  G4cout  << "crystalLattice = "<< crystalLattice << G4endl ;
  G4cout  << "region = "<< region << G4endl ;
  G4cout  << "colour = "<< colour << G4endl ;
  G4cout  << "radiationModel = "<< radiationModel << G4endl ;
  G4cout  << "aper1 = "<< aper1 << G4endl ;
  G4cout  << "aper2 = "<< aper1 << G4endl ;
  G4cout  << "apertureType = "<< apertureType << G4endl ;

  crystalMaterial   = BDSMaterials::Instance()->GetMaterial(material);
  vacuum = BDSMaterials::Instance()->GetMaterial(vacuumMaterial);

}

CrystalDeflector::~CrystalDeflector()
{

}

void CrystalDeflector::BuildContainerLogicalVolume()
{
  // containerSolid is a member of BDSAcceleratorComponent we should fill in
  containerSolid = new G4Box(name + "_container_solid",
                 containerHorizontalWidth * 0.5,
                 containerHorizontalWidth * 0.5,
			     chordLength * 0.5);


  G4cout  << "chordLength = "<< chordLength << G4endl ;

  // Set the extent of the component to its container. Without it the extent
  // stays zero, and the BDSIM link then sizes the output sampler from the
  // placement offset alone: particles exiting outside that (small) radius are
  // not recorded, and hence counted as lost.
  SetExtent(BDSExtent(containerHorizontalWidth * 0.5,
                      containerHorizontalWidth * 0.5,
                      chordLength * 0.5));

  // containerLogicalVolume is a member of BDSAcceleratorComponent we should fill in
  containerLogicalVolume = new G4LogicalVolume(containerSolid,
					       vacuum,
					       name + "_container_lv");
}

void CrystalDeflector::Build()
{
  // Call the base class implementation of this function. This builds the container
  // logical volume then sets some visualisation attributes for it.
  BDSAcceleratorComponent::Build();

  BuildCrystal();
}

void CrystalDeflector::BuildCrystal()
{
    //build crystal solid
    G4Box* crystalSolid = new G4Box(name + "_crystal_solid",
                   crystalXSize * 0.5,
                   crystalYSize * 0.5,
                   crystalMaterialThickness * 0.5);
    RegisterSolid(crystalSolid); // for deletion by bdsim

    // make a logical volume for the crystal
    G4LogicalVolume* crystalLV = new G4LogicalVolume(crystalSolid,
                            crystalMaterial,
                            name + "_crystal_lv");

    // visualisation attributes - make it nicely visible
    G4VisAttributes* crystalVis = new G4VisAttributes(*BDSColours::Instance()->GetColour(crystalcolour));
    crystalVis->SetVisibility(true);
    crystalLV->SetVisAttributes(crystalVis);
    RegisterVisAttributes(crystalVis);// for deletion by bdsim
    RegisterLogicalVolume(crystalLV); // for deletion by bdsim

    G4RotationMatrix* crystalRM = new G4RotationMatrix();
    crystalRM->rotateX(crystalAxisX);
    crystalRM->rotateY(crystalAxisY);
    crystalRM->rotateZ(crystalAxisZ);
    RegisterRotationMatrix(crystalRM); // for deletion by bdsim
    G4double crystalZPos = 0*CLHEP::cm;
    G4ThreeVector crystalPos = G4ThreeVector(0,0, crystalZPos);

    //physical volume
    auto crystalPV = new G4PVPlacement(crystalRM,
                      crystalPos,
                      crystalLV,
                      name + "_crystal",
                      containerLogicalVolume,
                      false,
                      0,
                      checkOverlaps);
    RegisterPhysicalVolume(crystalPV); // for deletion by bdsim

    G4Region* crystalRegion = new G4Region(crystalRegionName);
    crystalRegion->AddRootLogicalVolume(crystalLV);

    //create the channeling model for this region
    G4ChannelingFastSimModel* ChannelingModel = new G4ChannelingFastSimModel("ChannelingModel", crystalRegion);
    //activate the channeling model
    ChannelingModel->Input(crystalMaterial, crystalLattice);
    //setting bending angle of the crystal planes (default is 0)
    ChannelingModel->GetCrystalData()->SetBendingAngle(crystalBendingAngle,crystalLV);

    if(crystalRadiationModel){
        ChannelingModel->RadiationModelActivate();
        ChannelingModel->GetRadiationModel()->SetRadiationAngleFactor(4.);
        ChannelingModel->GetRadiationModel()->SetMinPhotonEnergy(1*CLHEP::MeV);
    }
}
