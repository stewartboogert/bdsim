/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/

#ifndef BDSGlobalConstants_h
#define BDSGlobalConstants_h 1

#include "G4Timer.hh" // This must be the first in the include list

#include <fstream>
#include <string>
#include <set>
#include <list>

// gab tmp
#include "G4ios.hh"

#include "G4ThreeVector.hh"
#include "globals.hh"
#include "BDSAcceleratorType.hh"
#include "G4FieldManager.hh"

using std::ofstream;
using std::istream;
using std::ifstream;
using std::ostream;
using std::list;

class BDSGlobalConstants 
{
  public:
  BDSGlobalConstants (const G4String& CardsFileName);
  ~BDSGlobalConstants();
  G4int ReadCard(G4String& AccName, G4String& BunchType);
  void SetLogFile(ofstream & os);
  void StripHeader(istream& is);
  
  BDSAcceleratorType* GetAcceleratorType();
  G4double GetBackgroundScaleFactor();

  G4String GetOutputNtupleFileName();

  void SetWorldSizeZ(G4double WorldSizeZ);
  G4double GetWorldSizeZ();

  void SetTotalS(G4double TotalS);
  G4double GetTotalS();  

  G4double GetComponentBoxSize();
  G4double GetMagnetPoleSize();
  G4double GetMagnetPoleRadius();

  G4double GetTunnelRadius(); 
  G4double GetHorizontalBeamlineOffset(); 
  G4double GetVerticalBeamlineOffset(); 

  G4double GetBeampipeRadius(); 
  G4double GetBeampipeThickness(); 

  G4double GetDeltaIntersection();
  G4double GetDeltaChord();

  G4double GetChordStepMinimum();

  G4double GetThresholdCutCharged();
  G4double GetThresholdCutPhotons();
  G4double GetTrackWeightFactor();

  G4bool GetSynchRadOn();
  G4bool GetSynchRescale();
  G4bool GetSynchTrackPhotons();
  G4double GetSynchLowX();
  G4double GetSynchLowGamE();

  G4bool GetPlanckOn();
  G4bool GetBDSeBremOn();

  G4bool GetVerboseStep();

  G4double GetLaserwireWavelength();
  G4ThreeVector GetLaserwireDir();
  G4bool GetLaserwireTrackPhotons();
  G4bool GetLaserwireTrackElectrons();

  G4bool GetReadBunchFile();
  G4bool GetWriteBunchFile();
  G4bool GetExtractBunchFile();

  G4bool GetTurnOnInteractions();
  G4bool GetUseTimer();
  G4bool GetUseEMHadronic();
  G4bool GetUseMuonPairProduction();
  G4bool GetStoreMuonTrajectories();
  G4bool GetUseMuonShowers();

  G4bool GetIncludeIronMagFields();

  G4bool GetUseHaloRadius();
  G4double GetHaloInnerRadius();
  G4double GetHaloOuterRadius();

  G4double GetInnerHaloX();
  G4double GetOuterHaloX();
  G4double GetInnerHaloY();
  G4double GetOuterHaloY();
  /*
  ifstream* GetInputBunchFile();
  ifstream* GetInputExtractBunchFile();
  ofstream* GetOutputBunchFile();
  */

  G4double GetLengthSafety();
  G4double GetEnergyOffset();
  G4double GetMuonProductionScaleFactor();
  G4double GetHadronInelasticScaleFactor();
  G4int GetVerboseEventNumber();

  G4Timer* GetTimer();

  ofstream GetEventOutput();

  G4long GetRandomSeed();
  G4bool GetUseBatch();
  G4int GetNumberToGenerate();

  G4int GetNumberOfEventsPerNtuple();

  // Internally transmitted variables:
  void SetGlobalBeamlineRotationY(G4double aYRotation);
  G4double GetGlobalBeamlineRotationY();

  G4double GetLWCalWidth();
  G4double GetLWCalOffset();
  G4String GetLWCalMaterial();


  G4FieldManager* GetZeroFieldManager();

protected:
private:
  // Data Members for Class Attributes
  ifstream ifs;
  ostream* log;

  BDSAcceleratorType* itsAccelerator;
  G4double itsBackgroundScaleFactor;

  G4double itsComponentBoxSize;
  G4double itsMagnetPoleSize;
  G4double itsMagnetPoleRadius;

  G4double itsTunnelRadius;
  G4double itsHorizontalBeamlineOffset;
  G4double itsVerticalBeamlineOffset;
 
  G4double itsBeampipeRadius; 
  G4double itsBeampipeThickness; 
  G4double itsWorldSizeZ; 
  G4double itsTotalS;

  G4double itsDeltaIntersection;
  G4double itsDeltaChord;
 
  G4double itsChordStepMinimum;

  G4double itsThresholdCutCharged;
  G4double itsThresholdCutPhotons;

  G4bool itsSynchRadOn;
  G4bool itsSynchRescale;
  G4bool itsSynchTrackPhotons;
  G4double itsSynchLowX;
  G4double itsSynchLowGamE;

  G4bool itsPlanckOn;
  G4bool itsBDSeBremOn;

  G4double itsLaserwireWavelength;
  G4ThreeVector itsLaserwireDir;
  G4bool itsLaserwireTrackPhotons;
  G4bool itsLaserwireTrackElectrons;
  G4bool itsTurnOnInteractions;

  G4bool itsReadBunchFile;
  G4bool itsWriteBunchFile;
  G4bool itsExtractBunchFile;
  G4bool itsVerboseStep;
  G4bool itsUseTimer;
  G4bool itsUseEMHadronic;
  G4bool itsUseMuonPairProduction;
  G4bool itsStoreMuonTrajectories;
  G4bool itsUseMuonShowers;

  G4bool itsIncludeIronMagFields;

  G4double itsLengthSafety;
  G4double itsTrackWeightFactor;
  G4double itsMuonProductionScaleFactor;
  G4double itsHadronInelasticScaleFactor;

  G4Timer* itsTimer;

  G4double itsEnergyOffset;
  G4int itsVerboseEventNumber;

  G4bool itsUseHaloRadius;
  G4double itsHaloInnerRadius;
  G4double itsHaloOuterRadius;

  G4double itsInnerHaloX;
  G4double itsOuterHaloX;
  G4double itsInnerHaloY;
  G4double itsOuterHaloY;

  G4bool itsUseBatch;
  G4long itsRandomSeed;
  G4int itsNumberToGenerate;

  G4String itsOutputNtupleFileName;
  G4int itsNumberOfEventsPerNtuple;

  // Internally used variables:
  G4double GlobalBeamlineRotationY;

  G4FieldManager* itsZeroFieldManager;

  G4double itsLWCalWidth;
  G4double itsLWCalOffset;
  G4String itsLWCalMaterial;

};

inline void BDSGlobalConstants::SetLogFile(ofstream & os)
{
  log=&os;
}

inline BDSAcceleratorType* BDSGlobalConstants::GetAcceleratorType()
{return itsAccelerator;}

inline G4double BDSGlobalConstants::GetBackgroundScaleFactor()
{return itsBackgroundScaleFactor;}

inline G4double BDSGlobalConstants::GetComponentBoxSize()
{return itsComponentBoxSize;}

inline G4double BDSGlobalConstants::GetMagnetPoleSize()
{return itsMagnetPoleSize;}

inline G4double BDSGlobalConstants::GetMagnetPoleRadius()
{return itsMagnetPoleRadius;}

inline G4double BDSGlobalConstants::GetTunnelRadius()
{return itsTunnelRadius;}

inline  G4double BDSGlobalConstants::GetHorizontalBeamlineOffset()
{return itsHorizontalBeamlineOffset;}
inline  G4double BDSGlobalConstants::GetVerticalBeamlineOffset()
{return itsVerticalBeamlineOffset;}



inline G4double BDSGlobalConstants::GetBeampipeRadius() 
{return itsBeampipeRadius;}

inline G4double BDSGlobalConstants::GetBeampipeThickness() 
{return itsBeampipeThickness;}

inline G4double BDSGlobalConstants::GetDeltaChord() 
{return itsDeltaChord;}

inline G4double BDSGlobalConstants::GetDeltaIntersection() 
{return itsDeltaIntersection;}

inline G4double BDSGlobalConstants::GetChordStepMinimum() 
{return itsChordStepMinimum;}

inline G4double BDSGlobalConstants::GetThresholdCutCharged() 
{return itsThresholdCutCharged;}
inline G4double BDSGlobalConstants::GetThresholdCutPhotons() 
{return itsThresholdCutPhotons;}

inline G4double BDSGlobalConstants::GetWorldSizeZ() 
{return itsWorldSizeZ;}
inline void BDSGlobalConstants::SetWorldSizeZ(G4double WorldSizeZ) 
{itsWorldSizeZ=WorldSizeZ;}

inline G4double BDSGlobalConstants::GetTotalS() 
{return itsTotalS;}
inline void BDSGlobalConstants::SetTotalS(G4double TotalS) 
{itsTotalS=TotalS;}

inline G4bool BDSGlobalConstants::GetSynchRadOn()
{return itsSynchRadOn;}

inline G4bool BDSGlobalConstants::GetSynchRescale()
{return itsSynchRescale;}

inline G4bool BDSGlobalConstants::GetSynchTrackPhotons()
{return itsSynchTrackPhotons ;}

inline G4double BDSGlobalConstants::GetSynchLowX()
{return itsSynchLowX ;}

inline G4double BDSGlobalConstants::GetSynchLowGamE()
{return itsSynchLowGamE ;}

inline G4bool BDSGlobalConstants::GetPlanckOn()
{return itsPlanckOn;}

inline G4bool BDSGlobalConstants::GetBDSeBremOn()
{return itsBDSeBremOn;}

inline G4bool BDSGlobalConstants::GetVerboseStep()
{return itsVerboseStep;}

inline G4double BDSGlobalConstants::GetLaserwireWavelength()
{return itsLaserwireWavelength ;}

inline G4ThreeVector BDSGlobalConstants::GetLaserwireDir()
{return itsLaserwireDir ;}

inline G4bool BDSGlobalConstants::GetLaserwireTrackPhotons()
{return itsLaserwireTrackPhotons ;}

inline G4bool BDSGlobalConstants::GetLaserwireTrackElectrons()
{return itsLaserwireTrackElectrons ;}

inline G4bool BDSGlobalConstants::GetReadBunchFile()
{return itsReadBunchFile;}

inline G4bool BDSGlobalConstants::GetExtractBunchFile()
{return itsExtractBunchFile;}

inline G4bool BDSGlobalConstants::GetWriteBunchFile()
{return itsWriteBunchFile ;}

inline G4double BDSGlobalConstants::GetLengthSafety()
{return itsLengthSafety;}

inline G4bool BDSGlobalConstants::GetTurnOnInteractions()
{return itsTurnOnInteractions;}

inline G4bool BDSGlobalConstants::GetUseTimer()
{return itsUseTimer;}

inline G4bool BDSGlobalConstants::GetUseEMHadronic()
{return itsUseEMHadronic;}

inline G4bool BDSGlobalConstants::GetIncludeIronMagFields()
{return itsIncludeIronMagFields;}

inline G4bool BDSGlobalConstants::GetUseMuonPairProduction()
{return itsUseMuonPairProduction;}

inline G4bool BDSGlobalConstants::GetStoreMuonTrajectories()
{return itsStoreMuonTrajectories;}

inline G4bool BDSGlobalConstants::GetUseMuonShowers()
{return itsUseMuonShowers;}

inline G4Timer* BDSGlobalConstants::GetTimer()
{return itsTimer;}

inline G4double BDSGlobalConstants::GetEnergyOffset()
{return itsEnergyOffset;}

inline  G4int BDSGlobalConstants::GetVerboseEventNumber()
{return itsVerboseEventNumber;}

inline  G4double BDSGlobalConstants::GetTrackWeightFactor()
{return itsTrackWeightFactor;}

inline  G4double BDSGlobalConstants::GetMuonProductionScaleFactor()
{return itsMuonProductionScaleFactor;}

inline  G4double BDSGlobalConstants::GetHadronInelasticScaleFactor()
{return itsHadronInelasticScaleFactor;}

inline  G4bool BDSGlobalConstants::GetUseHaloRadius()
{return itsUseHaloRadius;}
inline  G4double BDSGlobalConstants::GetHaloInnerRadius()
{return itsHaloInnerRadius;}
inline  G4double BDSGlobalConstants::GetHaloOuterRadius()
{return itsHaloOuterRadius;}


inline G4double BDSGlobalConstants::GetInnerHaloX()
{return itsInnerHaloX;}
inline G4double BDSGlobalConstants::GetOuterHaloX()
{return itsOuterHaloX;}
inline G4double BDSGlobalConstants::GetInnerHaloY()
{return itsInnerHaloY;}
inline G4double BDSGlobalConstants::GetOuterHaloY()
{return itsOuterHaloY;}

inline G4long BDSGlobalConstants::GetRandomSeed()
{return itsRandomSeed;}
inline  G4bool BDSGlobalConstants::GetUseBatch()
{return itsUseBatch;}
inline G4int BDSGlobalConstants::GetNumberToGenerate()
{return itsNumberToGenerate;}

inline G4String BDSGlobalConstants::GetOutputNtupleFileName()
{return itsOutputNtupleFileName;}
inline G4int BDSGlobalConstants::GetNumberOfEventsPerNtuple()
{return itsNumberOfEventsPerNtuple;}

// Internally transmitted variables
  // Internally transmitted variables:
inline void BDSGlobalConstants::SetGlobalBeamlineRotationY(G4double aYRotation)
{GlobalBeamlineRotationY=aYRotation;}
inline G4double BDSGlobalConstants::GetGlobalBeamlineRotationY()
{return GlobalBeamlineRotationY;}

inline G4FieldManager* BDSGlobalConstants::GetZeroFieldManager()
{return itsZeroFieldManager;}

inline  G4double BDSGlobalConstants::GetLWCalWidth()
{return itsLWCalWidth;}
inline  G4double BDSGlobalConstants::GetLWCalOffset()
{return itsLWCalOffset;}
inline  G4String BDSGlobalConstants::GetLWCalMaterial()
{return itsLWCalMaterial;}

extern BDSGlobalConstants* BDSGlobals;
#endif

