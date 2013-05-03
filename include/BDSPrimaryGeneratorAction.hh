/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 

   Modified 22.03.05 by J.C.Carter, Royal Holloway, Univ. of London.
   Added GABs SynchGen code
*/


#ifndef BDSPrimaryGeneratorAction_h
#define BDSPrimaryGeneratorAction_h 

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <fstream>
#include "BDSSynchrotronRadiation.hh"

#include "BDSGen5Shell.hh"

using std::ifstream;
using std::ostream;

class G4ParticleGun;
class G4Event;
class BDSDetectorConstruction;


class BDSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  BDSPrimaryGeneratorAction(BDSDetectorConstruction*);    
  ~BDSPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  
private:
  G4ParticleGun*              particleGun;	  //pointer a to G4 service class
  BDSDetectorConstruction*    BDSDetector;  //pointer to the geometry
  
  G4double weight;
 
  // beam data:
  G4double beta_x,sig_z,KineticEnergy,
    gamma,beta_y,sig_dp,emit_x,emit_y,charge;

  G4double sig_x,sig_xp,sig_y,sig_yp,sig_t;
  
  
  ifstream InputBunchFile;
  ifstream ExtractBunchFile;

  BDSGen5Shell* itsBDSGen5Shell;

  G4double logXfactor;
  G4double logYfactor;

  G4double itsSynchCritEng;
  BDSSynchrotronRadiation* itsBDSSynchrotronRadiation;

};

#endif


