/** This file is part of Corsair simulation.
 *
 *  Copyright 2011, 2012 Finnish Meteorological Institute
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

#include "operator_cellid.h"

using namespace std;

CellIDOP::CellIDOP(): DataOperator() { 
   cellIDprofileID = -1;
}

CellIDOP::~CellIDOP() { }

bool CellIDOP::finalize() {return true;}

std::string CellIDOP::getName() const {return "CellID";}

bool CellIDOP::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   return DataOperator::initialize(cr,sim,simClasses);
}

bool CellIDOP::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
   if (getInitialized() == false) return false;
   profile::start("CellID",cellIDprofileID);

   // Allocate output buffer:
   const vector<pargrid::CellID>& globalIDs = simClasses->pargrid.getGlobalIDs();
   const uint64_t N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();
   vector<pargrid::CellID> blockIDs(N_localBlocks*block::SIZE);
   for (size_t i=0; i<N_localBlocks; ++i) {
      for (int j=0; j<block::SIZE; ++j) blockIDs[i*block::SIZE+j] = globalIDs[i];
   }

   // Write buffer to output file:
   const string arrayType = "celldata";
   const uint64_t arraySize = N_localBlocks*block::SIZE;
   const uint64_t vectorSize = 1;

   map<string,string> attributes;
   attributes["name"] = getName();
   attributes["mesh"] = spatMeshName;
   attributes["type"] = arrayType;

   bool success = true;
   if (simClasses->vlsv.writeArray("VARIABLE",attributes,arraySize,vectorSize,&(blockIDs[0])) == false) success = false;
   profile::stop();
   return success;
}


