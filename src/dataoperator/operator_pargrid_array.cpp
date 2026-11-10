/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2025 Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include <vlsv_common.h>

#include "operator_pargrid_array.h"

using namespace std;

static const string regionName = "OperatorPargridData";

OperatorPargridArray::OperatorPargridArray(): DataOperator() {
   arrayInfoRead = false;
   #if PROFILE_LEVEL > 0
      profTotalTime = -1;
   #endif
}

OperatorPargridArray::~OperatorPargridArray() { }

bool OperatorPargridArray::addConfigFileItems(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
   cr.addComposed(regionName+".array_pargrid_name","Name of ParGrid array that is written to output file (string)");
   cr.addComposed(regionName+".array_output_name","Name of array in output file, must be compatible with VisIt (string)");
   cr.addComposed(regionName+".array_units","Units for the data in array (string)");
   return true;
}

bool OperatorPargridArray::finalize() { return true; }

std::string OperatorPargridArray::getName() const { return "OperatorPargridArrays"; }

bool OperatorPargridArray::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   // Init base class:
   bool success = DataOperator::initialize(cr,sim,simClasses);
   if (success == false) { return success; }
   // Read config file:
   if (addConfigFileItems(sim,simClasses,cr) == false) { return false; }
   cr.parse();
   cr.get(regionName+".array_pargrid_name",arrayPargridNames);
   cr.get(regionName+".array_output_name",arrayOutputNames);
   cr.get(regionName+".array_units",arrayUnits);
   if (arrayPargridNames.size() != arrayOutputNames.size()) { success = false; }
   if (arrayPargridNames.size() != arrayUnits.size()) { success = false; }
   return success;
}

/** Fetch information of ParGrid arrays and copy information relevant
 * to writing data to vector arrayInfo.
 * Return value: if true, array info was read successfully.*/
bool OperatorPargridArray::readArrayInfo() {
   // Get static ParGrid array info:
   vector<pargrid::DataID> dataIDs; /**< ParGrid ID of array.*/
   vector<string> names;            /**< Name of ParGrid array.*/
   vector<unsigned int> elements;   /**< Number of elements per cell, i.e. size of data vector.*/
   vector<string> datatypes;        /**< Array datatype as "int", "uint", or "float".*/
   vector<unsigned int> byteSizes;  /**< Byte size of each element.*/
   vector<const char*> pointers;    /**< Pointer to each array.*/
   simClasses->pargrid.getStaticUserDataInfo(dataIDs,names,elements,datatypes,byteSizes,pointers);
   // Search static ParGrid arrays for names defined in arrayNames:
   for (size_t j=0; j<names.size(); ++j) {
      for (size_t i=0; i<arrayPargridNames.size(); ++i) {
	 // Skip if array is not written to file:
	 if (arrayPargridNames[i] != names[j]) { continue; }
	 // Store array info to vector arrayInfo:
	 ArrayInfo info;
	 info.byteSize   = byteSizes[j];
	 info.dataID     = dataIDs[j];
	 info.datatype   = datatypes[j];
	 info.name       = arrayOutputNames[i];
	 info.units      = arrayUnits[i];
	 info.vectorSize = elements[j];
	 arrayInfo.push_back(info);
      }
   }
   // Deallocate vector arrayNames:
   vector<string> dummy;
   arrayPargridNames.swap(dummy);
   vector<string> dummy1;
   arrayOutputNames.swap(dummy1);
   vector<string> dummy2;
   arrayUnits.swap(dummy2);
   return true;
}

bool OperatorPargridArray::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
   bool success = true;
   // Check that class initialized successfully:
   if (getInitialized() == false) { return false; }
   #if PROFILE_LEVEL > 0
      profile::start("Pargrid Array(s)",profTotalTime);
   #endif
   // Ensure that array info has been read:
   if (arrayInfoRead == false) {
      if (readArrayInfo() == false) {
	 simClasses->logger << "(OPERATOR PARGRID ARRAY) ERROR: Failed to get array info" << endl << write;
         #if PROFILE_LEVEL > 0
	    profile::stop();
	 #endif
	 return false;
      }
   }
   // Write all arrays to output file:
   for (size_t i=0; i<arrayInfo.size(); ++i) {
      // Create map containing XML attributes:
      const uint64_t arraySize = simClasses->pargrid.getNumberOfLocalCells()*block::SIZE;
      map<string,string> attributes;
      attributes["name"] = arrayInfo[i].name;
      attributes["mesh"] = spatMeshName;
      attributes["centering"] = "zone";
      attributes["geometry"] = vlsv::getMeshGeometry(sim->meshGeometry);
      attributes["units"] = arrayInfo[i].units;
      // Get pointer to user data array:
      const char* ptr = simClasses->pargrid.getUserDataStatic<const char>(arrayInfo[i].dataID);
      if (ptr == NULL) {
	 simClasses->logger << "(OPERATOR PARGRID ARRAY) ERROR: Failed to write array '" << arrayInfo[i].name << "'" << endl << write;
	 success = false;
	 continue;
      }
      int blockSize = block::SIZE;
      const int vectorSize = arrayInfo[i].vectorSize / blockSize;
      // Write array to output file:
      if (simClasses->vlsv.writeArray("VARIABLE",attributes,arrayInfo[i].datatype,arraySize,vectorSize,arrayInfo[i].byteSize,ptr) == false) {
	 simClasses->logger << "(OPERATOR PARGRID ARRAY) ERROR: Failed to write array '" << arrayInfo[i].name << "'" << endl << write;
	 success = false;
      }
   }
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}
