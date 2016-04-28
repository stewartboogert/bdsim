#ifndef BDSEXECOPTIONS_H
#define BDSEXECOPTIONS_H

#include <getopt.h>

#include "globals.hh" // geant4 types / globals

#include "parser/options.h"

/**
 * @brief Executable option processing for BDSIM
 */

class BDSExecOptions
{
public: 
  BDSExecOptions(int argc, char **argv);

  ~BDSExecOptions(){;}

  /// Print out the available executable commands.
  void Usage() const;

  /// Print out the commands and their set values.
  void Print() const;

  /// Accessor for options generated by command line parsing.
  const GMAD::Options& Options() const {return options;}

  /// Acessor for convenience for the one thing that's needed before the parser options.
  inline G4String InputFileName() const {return options.inputFileName;}
  
private:
  /// Private default constructor to force use of provided one.
  BDSExecOptions();
  
  /** Parse the command line options
   * @param[in] argc command line input number of variables
   * @param[in] argv array of char* of the input parameters
   */
  void Parse(int arcg, char **argv);
  
  /// Helper method to set the BDSIMPath correctly
  G4String GetPath(G4String filename);

  /// The options instance that is populated by parsing the command line options.
  GMAD::Options options;
};

#endif
