#include "BDSIntegratorDipole2.hh"

#include "globals.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

#include "CLHEP/Units/SystemOfUnits.h"

BDSIntegratorDipole2::BDSIntegratorDipole2(G4Mag_EqRhs* eqOfMIn,
					   G4double     minimumRadiusOfCurvatureIn):
  G4MagHelicalStepper(eqOfMIn),
  minimumRadiusOfCurvature(minimumRadiusOfCurvatureIn),
  backupStepper(nullptr)
{
  //backupStepper = new G4ClassicalRK4(eqOfMIn, 6);
}

BDSIntegratorDipole2::~BDSIntegratorDipole2()
{
  delete backupStepper;
}

void BDSIntegratorDipole2::DumbStepper(const G4double yIn[],
				       G4ThreeVector  field,
				       G4double       stepLength,
				       G4double       yOut[])
{
  AdvanceHelix(yIn, field, stepLength, yOut);
}

void BDSIntegratorDipole2::Stepper(const G4double yIn[],
				   const G4double dydx[],
				   G4double       h,
				   G4double       yOut[],
				   G4double       yErr[])
{
  // Extra storage arrays.
  G4double yTemp[7], yTemp2[7];
  
  // Arrays for field querying (g4 interface)
  G4double bO[4], bM[4]; // original and mid point
  GetEquationOfMotion()->GetFieldValue(yIn, bO);
  G4ThreeVector bOriginal = G4ThreeVector(bO[0],bO[1],bO[2]);

  // Do a full step - the result we use
  AdvanceHelix(yIn, bOriginal, h, yOut);

  // Now we have the radius of curvature to check on.
  // If it's smaller than limit, we artificially advance the particle
  // along its helix axis (parallel to the field) so it'll hit something
  // and finish in a timely manner.
  const G4double radiusOfCurvature = GetRadHelix();
  if (radiusOfCurvature < minimumRadiusOfCurvature)
    {
      AdvanceHelixForSpiralling(yIn, bOriginal, h, yOut, yErr);
      return;
    }

  // error estimation - do two half steps and compare difference to
  // the result from one full step
  AdvanceHelix(yIn, bOriginal, h*0.5, yTemp); // first step

  // resample field at midway point (although if pure dipole, this is
  // unnecessary) - could go outside the range of the field though
  GetEquationOfMotion()->GetFieldValue(yTemp, bM);
  G4ThreeVector bMid = G4ThreeVector(bM[0],bM[1],bM[2]);

  AdvanceHelix(yTemp, bMid, h*0.5, yTemp2); // second step
  
  // Error estimation
  for(G4int i = 0; i < 6; i++)
    {yErr[i] = yOut[i] - yTemp2[i];}
}

void BDSIntegratorDipole2::AdvanceHelixForSpiralling(const G4double yIn[],
						     G4ThreeVector  field,
						     G4double       h,
						     G4double       yOut[],
						     G4double       yErr[])
{
  AdvanceHelix(yIn, field, h, yOut);

  G4ThreeVector unitMomentum    = G4ThreeVector(yIn[3], yIn[4], yIn[5]);
  G4ThreeVector unitField       = field.unit();
  //G4ThreeVector unitSideways    = unitField.cross(unitMomentum.unit());
  G4ThreeVector correctPosition = G4ThreeVector(yOut[0], yOut[1], yOut[2]);
  G4ThreeVector delta = h*unitField;
  G4ThreeVector newPosition     = correctPosition + delta;
  //G4cout << newPosition << delta << G4endl;
  //SetAngCurve(7);
  //SetRadHelix(0);
  yOut[0] = newPosition[0];
  yOut[1] = newPosition[1];
  yOut[2] = newPosition[2];

  //SetAngCurve(7);
  //SetRadHelix(backupStepper->DistChord()*0.5);

  // empirically chosen error values that give the fewest warnings
  // for spiralling particles in showers in strong dipole fields
  for(G4int i = 0; i < 3; i++)
    {yErr[i] = 1e-20;}
  for(G4int i = 3; i < 6; i++)
    {yErr[i] = 1e-40;}
}
