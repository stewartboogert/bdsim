/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/

#ifndef BDSDrift_h
#define BDSDrift_h 1

#include "globals.hh"
#include "BDSMultipole.hh"

class BDSDrift :public BDSMultipole
{
public:
  BDSDrift(G4String aName, 
	   G4double aLength,
	   G4double startAper, 
	   G4double endAper = BDSGlobalConstants::Instance()->GetBeampipeRadius(), 
	   G4String aTunnelMaterial = "", 
	   G4bool aperset = false, 
	   G4double aper = 0, 
	   G4double tunnelOffsetX = BDSGlobalConstants::Instance()->GetTunnelOffsetX());

  BDSDrift(G4String aName, 
	   G4double aLength,
           std::list<G4double> blmLocZ, 
	   std::list<G4double> blmLocTheta, 
	   G4double startAper, 
	   G4double endAper = BDSGlobalConstants::Instance()->GetBeampipeRadius(), 
	   G4String aTunnelMaterial = "", 
	   G4bool aperset = false, 
	   G4double aper = 0, 
	   G4double tunnelOffsetX = BDSGlobalConstants::Instance()->GetTunnelOffsetX(), 
	   G4double phiAngleIn = 0, 
	   G4double phiAngleOut = 0);
  ~BDSDrift();

protected:
  virtual void Build();

private:
  virtual void BuildBeampipe(G4String materialName = "");
  virtual void BuildBPFieldAndStepper();
  virtual void BuildBLMs();
  G4double itsStartOuterR;
  G4double itsEndOuterR;
  G4bool itsAperset;
  void init();
};

#endif
