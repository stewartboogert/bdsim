/* BDSIM code.    Version 1.0
   Author: John C. Carter, Royal Holloway, Univ. of London.
   Last modified 02.12.2004
   Copyright (c) 2004 by J.C.Carter.  ALL RIGHTS RESERVED. 
*/

#include "BDSGlobalConstants.hh" // must be first in include list
#include "BDSElement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "BDSHelixStepper.hh"
#include "BDSQuadStepper.hh"
#include "BDSOctStepper.hh"
#include "BDSSextStepper.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixMixedStepper.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixHeum.hh"
#include "G4SimpleRunge.hh"
#include "G4ExactHelixStepper.hh"
#include "BDSAcceleratorComponent.hh"
#include "G4CashKarpRKF45.hh"

// geometry drivers
#include "ggmad.hh"
#include "BDSGeometrySQL.hh"

#ifdef USE_LCDD
#include "BDSGeometryLCDD.hh"
#endif

#ifdef USE_GDML
#include "BDSGeometryGDML.hh"
#endif

#include <map>

using namespace std;

//============================================================

typedef std::map<G4String,int> LogVolCountMap;
extern LogVolCountMap* LogVolCount;

typedef std::map<G4String,G4LogicalVolume*> LogVolMap;
extern LogVolMap* LogVol;

extern BDSMaterials* theMaterials;
extern G4RotationMatrix* RotY90;

//============================================================

BDSElement::BDSElement(G4String aName, G4String geometry, G4String bmap,
		       G4double aLength, G4double bpRad, G4double outR, G4String aTunnelMaterial, G4double aTunnelRadius, G4double aTunnelOffsetX, G4String aTunnelCavityMaterial):
  BDSAcceleratorComponent(
			  aName,
			  aLength,bpRad,0,0,
			  SetVisAttributes(), aTunnelMaterial, "", 0., 0., 0., 0., aTunnelRadius*m, aTunnelOffsetX*m, aTunnelCavityMaterial),
  itsField(NULL), itsMagField(NULL)
{
  itsFieldVolName="";
  itsFieldIsUniform=false;
  itsOuterR = outR;
  SetType(_ELEMENT);

  //Set marker volume lengths
  CalculateLengths();

  // WARNING: ALign in and out will only apply to the first instance of the
  //          element. Subsequent copies will have no alignement set.
  align_in_volume = NULL;
  align_out_volume = NULL;

  if(!(*LogVolCount)[itsName])
    {
#ifdef DEBUG 
      G4cout<<"BDSElement : starting build logical volume "<<
        itsName<<G4endl;
#endif
      BuildGeometry(); // build element box
      
#ifdef DEBUG 
      G4cout<<"BDSElement : end build logical volume "<<
        itsName<<G4endl;
#endif

      PlaceComponents(geometry,bmap); // place components (from file) and build filed maps
    }
  else
    {
      (*LogVolCount)[itsName]++;
      
      itsMarkerLogicalVolume=(*LogVol)[itsName];
    }
}

void BDSElement::BuildGeometry()
{
  // multiple element instances not implemented yet!!!

  // Build the logical volume 

#ifdef DEBUG 
  G4cout<<"BDSElement : creating logical volume"<<G4endl;
#endif
  G4double elementSizeX=itsOuterR+BDSGlobals->GetLengthSafety()/2, elementSizeY = itsOuterR+BDSGlobals->GetLengthSafety()/2;
  
  
  elementSizeX = std::max(elementSizeX, this->GetTunnelRadius()+2*std::abs(this->GetTunnelOffsetX()) + BDSGlobals->GetTunnelThickness()+BDSGlobals->GetTunnelSoilThickness() + 4*BDSGlobals->GetLengthSafety() );   
  elementSizeY = std::max(elementSizeY, this->GetTunnelRadius()+2*std::abs(BDSGlobals->GetTunnelOffsetY()) + BDSGlobals->GetTunnelThickness()+BDSGlobals->GetTunnelSoilThickness()+4*BDSGlobals->GetLengthSafety() );

  G4double elementSize=std::max(elementSizeX, elementSizeY); 
  
  itsMarkerLogicalVolume = 
    new G4LogicalVolume(new G4Box(itsName+"generic_element",
                                  elementSize,
                                  elementSize,   
				  itsLength/2),
			theMaterials->GetMaterial(BDSGlobals->GetVacuumMaterial()),
			itsName);
  
  (*LogVolCount)[itsName] = 1;
  (*LogVol)[itsName] = itsMarkerLogicalVolume;
#ifndef NOUSERLIMITS
  itsOuterUserLimits = new G4UserLimits();
  G4double stepfactor=1e-2;
  itsOuterUserLimits->SetMaxAllowedStep(itsLength*stepfactor);
  itsOuterUserLimits->SetUserMaxTime(BDSGlobals->GetMaxTime());
  if(BDSGlobals->GetThresholdCutCharged()>0){
    itsOuterUserLimits->SetUserMinEkine(BDSGlobals->GetThresholdCutCharged());
  }
  itsMarkerLogicalVolume->SetUserLimits(itsOuterUserLimits);
#endif

  //Build the tunnel
  if(BDSGlobals->GetBuildTunnel()){
    BuildTunnel();
  }
}

// place components 
void BDSElement::PlaceComponents(G4String geometry, G4String bmap)
{

  G4String gFormat="", bFormat="";
  G4String gFile="", bFile="";

  // get geometry format and file
  if(geometry != ""){
    G4int pos = geometry.find(":");
    gFormat="none";
    if(pos<0) { 
      G4cerr<<"WARNING: invalid geometry reference format : "<<geometry<<endl;
    }
    else {
      gFormat = geometry.substr(0,pos);
      gFile = geometry.substr(pos+1,geometry.length() - pos); 
    }
  }

  // get fieldmap format and file
  bFormat="none";
  if(bmap != ""){
    G4int pos = bmap.find(":");
    if(pos<0) {
      G4cerr<<"WARNING: invalid B map reference format : "<<bmap<<endl; 
    }
    else {
      bFormat = bmap.substr(0,pos);
      bFile = bmap.substr(pos+1,bmap.length() - pos); 
    }
  }

  G4cout<<"placing components:\n geometry format - "<<gFormat<<G4endl<<
    "file - "<<gFile<<G4endl;

  G4cout<<"bmap format - "<<bFormat<<G4endl<<
    "file - "<<bFile<<G4endl;


  // get the geometry for the driver
  // different drivers may interpret the fieldmap differently
  // so a field map without geometry is not allowed

  GGmadDriver *ggmad;
  BDSGeometrySQL *Mokka;

#ifdef USE_LCDD
  BDSGeometryLCDD *LCDD;
#endif
#ifdef USE_GDML
  BDSGeometryGDML *GDML;
#endif

  if(gFormat=="gmad") {

    ggmad = new GGmadDriver(gFile);
    ggmad->Construct(itsMarkerLogicalVolume);

    // set sensitive volumes
    //vector<G4LogicalVolume*> SensComps = ggmad->SensitiveComponents;
    //for(G4int id=0; id<SensComps.size(); id++)
    //  SetMultipleSensitiveVolumes(SensComps[id]);

    
    SetMultipleSensitiveVolumes(itsMarkerLogicalVolume);

    // attach magnetic field if present

    if(bFormat=="XY")
      {
	itsMagField = new BDSXYMagField(bFile);

	// build the magnetic field manager and transportation
	BuildMagField();
      }
  }
  else if(gFormat=="lcdd") {
#ifdef USE_LCDD
    LCDD = new BDSGeometryLCDD(gFile);
    //Make marker visible (temp debug)
    G4VisAttributes* VisAttLCDD = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    VisAttLCDD->SetForceSolid(true);  
    VisAttLCDD->SetVisibility(true);
    itsMarkerLogicalVolume->SetVisAttributes(VisAttLCDD);

    LCDD->Construct(itsMarkerLogicalVolume);
    SetMultipleSensitiveVolumes(itsMarkerLogicalVolume);
    if(bFormat=="XY"){
      itsMagField = new BDSXYMagField(bFile);
      // build the magnetic field manager and transportation
      BuildMagField(true);
    } else if ( bFormat == "mokka" ){
      G4ThreeVector uniField = G4ThreeVector(0,3.5*tesla,0);
      std::vector<G4ThreeVector> uniFieldVect;
      uniFieldVect.push_back(uniField);
      std::vector<G4double> nulld;
      std::vector<G4String> nulls;
      //      itsMagField = new BDSMagFieldSQL(bFile,itsLength,nulls, nulld, nulls, nulld, nulls, nulld, LCDD->FieldVol, uniFieldVect);
    } else if ( bFormat == "test" ){
    }
    else if ( bFormat == "none" ){
      itsFieldIsUniform=LCDD->GetFieldIsUniform();
      if(itsFieldIsUniform){
	G4cout << "BDSElement> using LCDD format uniform field..." << G4endl;
	itsUniformMagField=LCDD->GetUniformField();
      }else{
	itsMagField=LCDD->GetField();
      }
      itsFieldVolName=LCDD->GetFieldVolName();
      BuildMagField(true);
    }

    vector<G4LogicalVolume*> SensComps = LCDD->SensitiveComponents;
    for(G4int id=0; id<(G4int)SensComps.size(); id++)
      SetMultipleSensitiveVolumes(SensComps[id]);

#else
    G4cout << "LCDD support not selected during BDSIM configuration" << G4endl;
    G4Exception("Please re-compile BDSIM with USE_LCDD flag in Makefile");
#endif
  }
  else if(gFormat=="mokka") {

    Mokka = new BDSGeometrySQL(gFile,itsLength);
    Mokka->Construct(itsMarkerLogicalVolume);
    vector<G4LogicalVolume*> SensComps = Mokka->SensitiveComponents;
    for(unsigned int i=0; i<Mokka->GetMultiplePhysicalVolumes().size(); i++){
      SetMultiplePhysicalVolumes(Mokka->GetMultiplePhysicalVolumes().at(i));
    }
    for(G4int id=0; id<(G4int)SensComps.size(); id++)
      SetMultipleSensitiveVolumes(SensComps[id]);
    align_in_volume = Mokka->align_in_volume;
    align_out_volume = Mokka->align_out_volume;
    // attach magnetic field if present

    if(bFormat=="mokka" || bFormat=="none")
      {
	if(Mokka->HasFields || bFile!="")
	  // Check for field file or volumes with fields
	  // as there may be cases where there are no bFormats given
	  // in gmad file but fields might be set to volumes in SQL files
	  {
	    itsMagField = new BDSMagFieldSQL(bFile,itsLength,
					  Mokka->QuadVolBgrad,
					  Mokka->SextVolBgrad,
					  Mokka->OctVolBgrad,
					  Mokka->UniformFieldVolField,
					  Mokka->nPoleField,
					  Mokka->HasUniformField);


	    // build the magnetic field manager and transportation
	    BuildMagField(true);
	  }
      }
  }
  else if(gFormat=="gdml") {
#ifdef USE_GDML
    GDML = new BDSGeometryGDML(gFile);
    GDML->Construct(itsMarkerLogicalVolume);
#else
    G4cout << "GDML support not selected during BDSIM configuration" << G4endl;
    G4Exception("Please re-compile BDSIM with USE_GDML flag in Makefile");
#endif
  }
  else {
    G4cerr<<"geometry format "<<gFormat<<" not supported"<<G4endl;
  }




  //   // test - dump field values
//     G4cout<<"dumping field values..."<<G4endl;
//     G4double Point[4];
//     G4double BField[6];
//     ofstream testf("fields.dat");

//     for(G4double x=-1*m;x<1*m;x+=1*cm)
//        for(G4double y=-1*m;y<1*m;y+=1*cm)
// 	 for(G4double z=-1*m;z<1*m;z+=1*cm)
//       {
// 	Point[0] = x;
// 	Point[1] = y;
// 	Point[2] = z;
// 	Point[3] = 0;
// 	itsMagField->GetFieldValue(Point,BField);
// 	testf<<Point[0]<<" "<<Point[1]<<" "<<Point[2]<<" "<<
// 	  BField[0]<<" "<<BField[1]<<" "<<BField[2]<<G4endl;
//       }

//     testf.close();
//     G4cout<<"done"<<G4endl;


}



G4VisAttributes* BDSElement::SetVisAttributes()
{
  itsVisAttributes=new G4VisAttributes(G4Colour(0.5,0.5,1));
  return itsVisAttributes;
}


void BDSElement::BuildMagField(G4bool forceToAllDaughters)
{
  G4cout << "BDSElement.cc::BuildMagField Building magnetic field...." << G4endl;
  // create a field manager
  G4FieldManager* fieldManager = new G4FieldManager();


  if(!itsFieldIsUniform){
    itsEqRhs = new G4Mag_UsualEqRhs(itsMagField);
    if( (itsMagField->GetHasUniformField())&!(itsMagField->GetHasNPoleFields() || itsMagField->GetHasFieldMap())){
      //    itsFStepper = new G4ExactHelixStepper(itsEqRhs); 
      itsFStepper = new G4HelixMixedStepper(itsEqRhs,6); 
    } else {
      //        itsFStepper = new G4HelixMixedStepper(itsEqRhs,6); 
      itsFStepper = new G4HelixImplicitEuler(itsEqRhs); 
      //    itsFStepper = new G4HelixSimpleRunge(itsEqRhs);
      //    itsFStepper = new G4HelixHeum(itsEqRhs);
      //    itsFStepper = new G4ClassicalRK4(itsEqRhs);
    }
    //    //  itsFStepper = new G4HelixSimpleRunge(itsEqRhs); 
    //  itsFStepper = new G4HelixExplicitEuler(itsEqRhs); 
    fieldManager->SetDetectorField(itsMagField );
  } else {
    itsEqRhs = new G4Mag_UsualEqRhs(itsUniformMagField);
    //    itsFStepper = new G4HelixImplicitEuler(itsEqRhs); 
    //  itsFStepper = new G4CashKarpRKF45(itsEqRhs); 
    itsFStepper = new G4HelixMixedStepper(itsEqRhs, 6); 
    fieldManager->SetDetectorField(itsUniformMagField );
  }

  if(BDSGlobals->GetDeltaOneStep()>0){
    fieldManager->SetDeltaOneStep(BDSGlobals->GetDeltaOneStep());
  }
  if(BDSGlobals->GetMaximumEpsilonStep()>0){
  fieldManager->SetMaximumEpsilonStep(BDSGlobals->GetMaximumEpsilonStep());
  }
  if(BDSGlobals->GetMinimumEpsilonStep()>=0){
  fieldManager->SetMinimumEpsilonStep(BDSGlobals->GetMinimumEpsilonStep());
  }
  if(BDSGlobals->GetDeltaIntersection()>0){
  fieldManager->SetDeltaIntersection(BDSGlobals->GetDeltaIntersection());
  }

  G4MagInt_Driver* fIntgrDriver = new G4MagInt_Driver(BDSGlobals->GetChordStepMinimum(),
							itsFStepper, 
							itsFStepper->GetNumberOfVariables() );
  fChordFinder = new G4ChordFinder(fIntgrDriver);
  
  fChordFinder->SetDeltaChord(BDSGlobals->GetDeltaChord());
  
  fieldManager->SetChordFinder( fChordFinder ); 

  /*
    if(itsFieldVolName != ""){
    itsFieldVolName=itsFieldVolName+"_PhysiComp";
    G4cout << "Printing daughters... itsFieldVolName=" << itsFieldVolName << G4endl;
    for(int i=0;itsMarkerLogicalVolume->IsAncestor(itsMarkerLogicalVolume->GetDaughter(i)); i++){
    G4cout << itsMarkerLogicalVolume->GetDaughter(i)->GetName() << endl;
    if( itsMarkerLogicalVolume->GetDaughter(i)->GetName() == itsFieldVolName){
    itsMarkerLogicalVolume->GetDaughter(i)->GetLogicalVolume()->SetFieldManager(fieldManager,forceToAllDaughters);
    G4cout << "Field volume set..." << G4endl;
    break;
    }
    }
    } else {
    }
  */
  itsMarkerLogicalVolume->SetFieldManager(fieldManager,forceToAllDaughters);
}

// creates a field mesh in the reference frame of a physical volume
// from  b-field map value list 
// has to be called after the component is placed in the geometry
    void BDSElement::PrepareField(G4VPhysicalVolume *referenceVolume)
{
  if(!itsMagField) return;
  itsMagField->Prepare(referenceVolume);
}

// Rotates and positions the marker volume before it is placed in
// BDSDetectorConstruction. It aligns the marker volume so that the
// the beamline goes through the specified daugther volume (e.g. for mokka)
void BDSElement::AlignComponent(G4ThreeVector& TargetPos, 
				G4RotationMatrix *TargetRot, 
				G4RotationMatrix& globalRotation,
				G4ThreeVector& rtot,
				G4ThreeVector& rlast,
				G4ThreeVector& localX,
				G4ThreeVector& localY,
				G4ThreeVector& localZ)
{
  
  
  if(align_in_volume == NULL)
    {
      if(align_out_volume == NULL)
	{
	  // advance co-ords in usual way if no alignment volumes found
	  
	  rtot = rlast + localZ*(itsLength/2);
	  rlast = rtot + localZ*(itsLength/2);
	  return;
	}
      else 
	{
	  G4cout << "BDSElement : Aligning outgoing to SQL element " 
		 << align_out_volume->GetName() << G4endl;
	  G4RotationMatrix Trot = *TargetRot;
	  G4RotationMatrix trackedRot;
	  G4RotationMatrix outRot = *(align_out_volume->GetFrameRotation());
	  trackedRot.transform(outRot.inverse());
	  trackedRot.transform(Trot.inverse());
	  globalRotation = trackedRot;

	  G4ThreeVector outPos = align_out_volume->GetFrameTranslation();
	  G4ThreeVector diff = outPos;

	  G4ThreeVector zHalfAngle = G4ThreeVector(0.,0.,1.);

	  zHalfAngle.transform(globalRotation);

	  //moving positioning to outgoing alignment volume
	  rlast = TargetPos - ((outPos.unit()).transform(Trot.inverse()) )*diff.mag();
	  localX.transform(outRot.inverse());
	  localY.transform(outRot.inverse());
	  localZ.transform(outRot.inverse());

	  localX.transform(Trot.inverse());
	  localY.transform(Trot.inverse());
	  localZ.transform(Trot.inverse());

	  //moving position in Z be at least itsLength/2 away
	  rlast +=zHalfAngle*(itsLength/2 + diff.z());
	  return;
	}
    }

  if(align_in_volume != NULL)
    {
      G4cout << "BDSElement : Aligning incoming to SQL element " 
      	     << align_in_volume->GetName() << G4endl;
      
      const G4RotationMatrix* inRot = align_in_volume->GetFrameRotation();
      TargetRot->transform((*inRot).inverse());
      
      G4ThreeVector inPos = align_in_volume->GetFrameTranslation();
      inPos.transform((*TargetRot).inverse());
      TargetPos+=G4ThreeVector(inPos.x(), inPos.y(), 0.0);
      
      if(align_out_volume == NULL)
	{
	  // align outgoing (i.e. next component) to Marker Volume
	  
	  G4RotationMatrix Trot = *TargetRot;
	  globalRotation.transform(Trot.inverse());
	  
	  G4ThreeVector zHalfAngle = G4ThreeVector(0.,0.,1.);
	  zHalfAngle.transform(Trot.inverse());
	  
	  rlast = TargetPos + zHalfAngle*(itsLength/2);
	  localX.transform(Trot.inverse());
	  localY.transform(Trot.inverse());
	  localZ.transform(Trot.inverse());
	  return;
	}

      else
	{
	  G4cout << "BDSElement : Aligning outgoing to SQL element " 
		 << align_out_volume->GetName() << G4endl;
	  G4RotationMatrix Trot = *TargetRot;
	  G4RotationMatrix trackedRot;
	  G4RotationMatrix outRot = *(align_out_volume->GetFrameRotation());
	  trackedRot.transform(outRot.inverse());
	  trackedRot.transform(Trot.inverse());
	  globalRotation = trackedRot;

	  G4ThreeVector outPos = align_out_volume->GetFrameTranslation();
	  G4ThreeVector diff = outPos;

	  G4ThreeVector zHalfAngle = G4ThreeVector(0.,0.,1.);

	  zHalfAngle.transform(globalRotation);

	  //moving positioning to outgoing alignment volume
	  rlast = TargetPos - ((outPos.unit()).transform(Trot.inverse()) )*diff.mag();
	  localX.transform(outRot.inverse());
	  localY.transform(outRot.inverse());
	  localZ.transform(outRot.inverse());

	  localX.transform(Trot.inverse());
	  localY.transform(Trot.inverse());
	  localZ.transform(Trot.inverse());

	  //moving position in Z be at least itsLength/2 away
	  rlast +=zHalfAngle*(itsLength/2 + diff.z());
	  return;
	}
    }
  
}

BDSElement::~BDSElement()
{

  delete itsVisAttributes;
  delete itsMarkerLogicalVolume;
  delete fChordFinder;
  delete itsFStepper;
  delete itsFEquation;
}
