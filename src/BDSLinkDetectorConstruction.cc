/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2020.

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
#include "BDSAcceleratorModel.hh"
#include "BDSBeamline.hh"
#include "BDSBeamlineElement.hh"
#include "BDSCollimatorJaw.hh"
#include "BDSComponentFactory.hh"
#include "BDSCrystalInfo.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSExtent.hh"
#include "BDSExtentGlobal.hh"
#include "BDSGlobalConstants.hh"
#include "BDSLinkComponent.hh"
#include "BDSLinkDetectorConstruction.hh"
#include "BDSLinkOpaqueBox.hh"
#include "BDSLinkPrimaryGeneratorAction.hh"
#include "BDSLinkRegistry.hh"
#include "BDSMaterials.hh"
#include "BDSParser.hh"
#include "BDSSDManager.hh"
#include "BDSTiltOffset.hh"

#include "parser/element.h"
#include "parser/elementtype.h"

#include "G4Box.hh"
#include "G4ChannelingOptrMultiParticleChangeCrossSection.hh"
#include "G4PVPlacement.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include "G4VisAttributes.hh"

#include <set>
#include <vector>
#include <include/BDSCrystalInfo.hh>

class BDSParticleDefinition;

BDSLinkDetectorConstruction::BDSLinkDetectorConstruction():
  worldSolid(nullptr),
  worldPV(nullptr),
  linkBeamline(nullptr),
  linkRegistry(nullptr),
  primaryGeneratorAction(nullptr),
  designParticle(nullptr),
  crystalBiasing(nullptr)
{
  linkRegistry = new BDSLinkRegistry();
  BDSSDManager::Instance()->SetLinkRegistry(linkRegistry);
}

BDSLinkDetectorConstruction::~BDSLinkDetectorConstruction()
{
  delete linkBeamline;
  delete linkRegistry;
  delete crystalBiasing;
}

G4VPhysicalVolume* BDSLinkDetectorConstruction::Construct()
{
  BDSGlobalConstants* globalConstants = BDSGlobalConstants::Instance();

  auto componentFactory = new BDSComponentFactory(designParticle, nullptr, false);
  auto beamline = BDSParser::Instance()->GetBeamline();

  std::vector<BDSLinkOpaqueBox*> opaqueBoxes = {};
  linkBeamline = new BDSBeamline();

  for (auto elementIt = beamline.begin(); elementIt != beamline.end(); ++elementIt)
    {
      GMAD::ElementType eType = elementIt->type;

      if (eType == GMAD::ElementType::_LINE || eType == GMAD::ElementType::_REV_LINE)
        {continue;}
      
      std::set<GMAD::ElementType> acceptedTypes = {GMAD::ElementType::_ECOL,
						   GMAD::ElementType::_RCOL,
						   GMAD::ElementType::_JCOL,
						   GMAD::ElementType::_CRYSTALCOL,
						   GMAD::ElementType::_ELEMENT};
      auto search = acceptedTypes.find(eType);
      if (search == acceptedTypes.end())
	{throw BDSException(G4String("Unsupported element type for link = " + GMAD::typestr(eType)));}

      // Only need first argument, the rest pertain to beamlines.
      BDSAcceleratorComponent* component = componentFactory->CreateComponent(&(*elementIt),
									     nullptr,
									     nullptr,
									     0);

      BDSTiltOffset* to = new BDSTiltOffset(elementIt->offsetX * CLHEP::m,
                                            elementIt->offsetY * CLHEP::m,
                                            elementIt->tilt * CLHEP::rad);
      BDSLinkOpaqueBox* opaqueBox = new BDSLinkOpaqueBox(component,
                                                         to,
                                                         component->GetExtent().MaximumAbsTransverse());

      opaqueBoxes.push_back(opaqueBox);

      BDSLinkComponent* comp = new BDSLinkComponent(opaqueBox->GetName(),
							opaqueBox,
							opaqueBox->GetExtent().DZ());

      nameToElementIndex[elementIt->name] = (G4int)linkBeamline->size();
      linkBeamline->AddComponent(comp);
    }

  // update world extents and world solid
  UpdateWorldSolid();
  
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldSolid,
						 BDSMaterials::Instance()->GetMaterial("G4_Galactic"),
						 "world_lv");

  G4VisAttributes* debugWorldVis = new G4VisAttributes(*(BDSGlobalConstants::Instance()->ContainerVisAttr()));
  debugWorldVis->SetForceWireframe(true);//just wireframe so we can see inside it
  worldLV->SetVisAttributes(debugWorldVis);
  worldLV->SetUserLimits(globalConstants->DefaultUserLimits());

  worldPV = new G4PVPlacement(nullptr,
				   G4ThreeVector(),
				   worldLV,
				   "world_pv",
				   nullptr,
				   false,
				   0,
				   true);

  // place any defined link elements in input
  for (auto element : *linkBeamline)
    {PlaceOneComponent(element);}

  delete componentFactory;

  return worldPV;
}

void BDSLinkDetectorConstruction::AddLinkCollimatorJaw(const std::string& collimatorName,
						                                           const std::string& materialName,
                                                       G4double length,
                                                       G4double halfApertureLeft,
                                                       G4double halfApertureRight,
                                                       G4double rotation,
                                                       G4double xOffset,
                                                       G4double yOffset,
                                                       G4bool   buildLeftJaw,
                                                       G4bool   buildRightJaw,
                                                       G4bool   /*isACrystal*/,
                                                       G4double crystalAngle)
{
  auto componentFactory = new BDSComponentFactory(designParticle, nullptr, false);

  std::map<std::string, std::string> collimatorToCrystal =
    {
     {"cry.mio.b1", "stf75"},   //b1 h
     {"cry.mio.b2", "tcp76"},   //b2 h
     {"tcpv.a6l7.b1", "qmp34"}, // b1 v
     {"tcpv.a6r7.b2", "qmp53"}  // b2 v
    };
  auto searchC = collimatorToCrystal.find(collimatorName); // will use later if needed
  G4bool isACrystal = searchC != collimatorToCrystal.end();

  std::map<std::string, std::string> sixtrackToBDSIM =
      {
          {"CU", "Cu"},
          {"W",  "W"},
          {"C",  "G4_GRAPHITE_POROUS"},
          {"Si", "Si"},
          {"SI", "Si"}
      };
  std::string g4material;
  auto search = sixtrackToBDSIM.find(materialName);
  if (search != sixtrackToBDSIM.end())
    {g4material = search->second;}
  else
    {g4material = materialName;}

  // build component
  GMAD::Element el = GMAD::Element();
  el.type     = GMAD::ElementType::_JCOL;
  el.name     = collimatorName;
  el.material = g4material;
  el.l        = length / CLHEP::m;
  el.xsizeLeft  = halfApertureLeft / CLHEP::m;
  el.xsizeRight = halfApertureRight / CLHEP::m;
  el.ysize    = 0.2; // half size
  el.tilt     = rotation / CLHEP::rad;
  el.offsetX  = xOffset / CLHEP::m;
  el.offsetY  = yOffset / CLHEP::m;
  el.horizontalWidth = 2.0; // m

  // if we don't want to build a jaw, then we set it to outside the width.
  if (!buildLeftJaw)
    {el.xsizeLeft = el.horizontalWidth * 1.2;}
  if (!buildRightJaw)
    {el.xsizeRight = el.horizontalWidth * 1.2;}

  if (isACrystal)
    {
      // find the bending angle of this particular crystal
      // so we can add half of that on for the BDSIM convention of the 0 about the centre for crystals
      BDSCrystalInfo* ci = componentFactory->PrepareCrystalInfo(searchC->second);
      crystalAngle *= CLHEP::rad;
      crystalAngle += 0.5 * ci->bendingAngleYAxis;
      delete ci; // no longer needed

      el.type = GMAD::ElementType::_CRYSTALCOL;
      el.apertureType = "circularvacuum";
      // need a small margin in length as crystal may have angled face and be rotated
      el.l += 10e-6; // TODO - confirm margin with sixtrack interface backtracking on input side
      if (collimatorName.find("2") != std::string::npos) // b2
        {
          el.crystalLeft = collimatorToCrystal[collimatorName];
          el.crystalAngleYAxisLeft = crystalAngle + 0.5 * ci->bendingAngleYAxis;
        }
      else
        {
          el.crystalRight = collimatorToCrystal[collimatorName];
          el.crystalAngleYAxisRight = crystalAngle - 0.5 * ci->bendingAngleYAxis;
        }
    }
  else
    {el.region = "r1";} // stricter range cuts for default collimators

  BDSAcceleratorComponent* component = componentFactory->CreateComponent(&el,
									 nullptr,
									 nullptr,
									 0); 

  // wrap in box
  BDSTiltOffset* to = new BDSTiltOffset(el.offsetX * CLHEP::m,
                                        el.offsetY * CLHEP::m,
                                        el.tilt * CLHEP::rad);
  BDSLinkOpaqueBox* opaqueBox = new BDSLinkOpaqueBox(component,
						     to,
						     component->GetExtent().MaximumAbsTransverse());

  // add to beam line
  BDSLinkComponent* comp = new BDSLinkComponent(opaqueBox->GetName(),
						opaqueBox,
						opaqueBox->GetExtent().DZ());
  nameToElementIndex[collimatorName] = (G4int)linkBeamline->size();
  linkBeamline->AddComponent(comp);

  // update world extents and world solid
  UpdateWorldSolid();

  // place that one element
  PlaceOneComponent(linkBeamline->back());

  // update crystal biasing
  BuildPhysicsBias();
}

void BDSLinkDetectorConstruction::UpdateWorldSolid()
{
  BDSExtentGlobal we = linkBeamline->GetExtentGlobal();
  we = we.ExpandToEncompass(BDSExtentGlobal(BDSExtent(10*CLHEP::m, 10*CLHEP::m, 10*CLHEP::m))); // minimum size
  G4ThreeVector worldExtentAbs = we.GetMaximumExtentAbsolute();
  worldExtentAbs *= 1.2;

  if (!worldSolid)
    {
      worldSolid = new G4Box("world_solid",
                             worldExtentAbs.x(),
                             worldExtentAbs.y(),
                             worldExtentAbs.z());
    }
  else
    {
      worldSolid->SetXHalfLength(worldExtentAbs.x());
      worldSolid->SetYHalfLength(worldExtentAbs.y());
      worldSolid->SetZHalfLength(worldExtentAbs.z());
    }
  worldExtent = BDSExtent(worldExtentAbs);
  if (primaryGeneratorAction)
    {primaryGeneratorAction->SetWorldExtent(worldExtent);}
}

void BDSLinkDetectorConstruction::PlaceOneComponent(const BDSBeamlineElement* element)
{
  G4String placementName = element->GetPlacementName() + "_pv";
  G4Transform3D* placementTransform = element->GetPlacementTransform();
  G4int copyNumber = element->GetCopyNo();
  // auto pv =
  new G4PVPlacement(*placementTransform,
                    placementName,
                    element->GetContainerLogicalVolume(),
                    worldPV,
                    false,
                    copyNumber,
                    true);

  auto lc = dynamic_cast<BDSLinkComponent*>(element->GetAcceleratorComponent());
  if (!lc)
    {return;}
  BDSLinkOpaqueBox* el = lc->Component();
  G4Transform3D elCentreToStart = el->TransformToStart();
  G4Transform3D globalToStart = elCentreToStart * (*placementTransform);
  linkRegistry->Register(el, globalToStart);
}

void BDSLinkDetectorConstruction::BuildPhysicsBias()
{
  if (!crystalBiasing) // cache it because we may have to dynamically add later
    {crystalBiasing = new G4ChannelingOptrMultiParticleChangeCrossSection();}

  // crystal biasing necessary for implementation of variable density
  std::set<G4LogicalVolume*>* crystals = BDSAcceleratorModel::Instance()->VolumeSet("crystals");
  if (!crystals->empty())
    {
      G4cout << __METHOD_NAME__ << "Using crystal biasing: true" << G4endl; // to match print out style further down
      for (auto crystal : *crystals)
        {
          // if it hasn't been added already - g4 will whine for double registration
          if (!crystalBiasing->GetBiasingOperator(crystal))
            {crystalBiasing->AttachTo(crystal);}
        }
    }
}
