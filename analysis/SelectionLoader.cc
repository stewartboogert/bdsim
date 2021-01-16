/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2021.

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
#include "SelectionLoader.hh"

#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>

std::string RBDS::LoadSelection(const std::string& selectionFile)
{
  std::ifstream f(selectionFile.c_str());
  if(!f)
    {throw std::string("LoadSelection> could not open file " + selectionFile);}

  int lineCounter = 0;
  std::string line;

  // unique patterns to match
  // match a line starting with #
  std::regex comment("^\\#.*");

  std::string selection;  
  while(std::getline(f, line))
    {
      lineCounter++;
      if (std::all_of(line.begin(), line.end(), isspace))
	{continue;} // skip empty lines
      else if (std::regex_search(line, comment))
	{continue;} // skip lines starting with '#'
      else
	{
	  std::vector<std::string> results;
	  std::regex wspace("\\s+"); // any whitepsace
	  // -1 here makes it point to the suffix, ie the word rather than the wspace
	  std::sregex_token_iterator iter(line.begin(), line.end(), wspace, -1);
	  std::sregex_token_iterator end;
	  for (; iter != end; ++iter)
	    {
	      std::string res = (*iter).str();
	      results.push_back(res);
	    }
	  if (results.size() != 1)
	    {throw std::string("More than one word on line -> no white space allowed");}
	  else
	    {selection = results[0];}
	}
    }
  
  f.close();

  return selection;
}
