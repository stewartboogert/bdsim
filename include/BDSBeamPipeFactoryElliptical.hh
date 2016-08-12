#ifndef BDSBEAMPIPEFACTORYELLIPTICAL_H
#define BDSBEAMPIPEFACTORYELLIPTICAL_H

#include "BDSBeamPipeFactoryBase.hh"
#include "BDSBeamPipe.hh"

/**
 * @brief Factory for elliptical beam pipes
 * 
 * singleton pattern
 * 
 * @author Laurie Nevay <laurie.nevay@rhul.ac.uk>
 */

class BDSBeamPipeFactoryElliptical: public BDSBeamPipeFactoryBase
{
public:
  static BDSBeamPipeFactoryElliptical* Instance(); /// singleton accessor
  
  virtual ~BDSBeamPipeFactoryElliptical();

  virtual BDSBeamPipe* CreateBeamPipe(G4String    nameIn,                     // name
				      G4double    lengthIn,                   // length [mm]
				      G4double    aper1 = 0,                  // aperture parameter 1
				      G4double    aper2 = 0,                  // aperture parameter 2
				      G4double    aper3 = 0,                  // aperture parameter 3
				      G4double    aper4 = 0,                  // aperture parameter 4
				      G4Material* vacuumMaterialIn = nullptr,    // vacuum material
				      G4double    beamPipeThicknessIn = 0,    // beampipe thickness [mm]
				      G4Material* beamPipeMaterialIn = nullptr   // beampipe material
				      );
  
  virtual BDSBeamPipe* CreateBeamPipeAngledInOut(G4String    nameIn,
						 G4double    lengthIn,
						 G4double    angleInIn,  // the normal angle of the input face
						 G4double    angleOutIn, // the normal angle of the output face
						 G4double    aper1 = 0,
						 G4double    aper2 = 0,
						 G4double    aper3 = 0,
						 G4double    aper4 = 0,
						 G4Material* vacuumMaterialIn = nullptr,
						 G4double    beamPipeThicknessIn = 0,
						 G4Material* beamPipeMaterialIn = nullptr
						 );

private:
  BDSBeamPipeFactoryElliptical(); /// private default constructor - singelton pattern
  static BDSBeamPipeFactoryElliptical* _instance;

  //abstract common build features to one function
  //use member variables unique to this factory to pass them around

  /// only the solids are unique, once we have those, the logical volumes and placement in the
  /// container are the same.  group all this functionality together
  BDSBeamPipe* CommonFinalConstruction(G4String    nameIn,
				       G4Material* vacuumMaterialIn,
				       G4Material* beamPipeMaterialIn,
				       G4double    lengthIn,
				       G4double    aper1In,
				       G4double    aper2In,
				       G4double    beamPipeThicknessIn);

  /// The angled ones have degeneracy in the geant4 solids they used so we can
  /// avoid code duplication by grouping common construction tasks.
  void CreateGeneralAngledSolids(G4String      nameIn,
				 G4double      lengthIn,
				 G4double      aper1In,
				 G4double      aper2In,
				 G4double      beamPipeThicknessIn,
				 G4ThreeVector inputfaceIn,
				 G4ThreeVector outputfaceIn);
};
  
#endif
