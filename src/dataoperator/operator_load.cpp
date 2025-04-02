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

#include "operator_load.h"

using namespace std;

LoadOP::LoadOP(): DataOperator() { 
   profileID = -1;
}

LoadOP::~LoadOP() { }

bool LoadOP::finalize() {return true;}

std::string LoadOP::getName() const {return "Load";}

bool LoadOP::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   return DataOperator::initialize(cr,sim,simClasses);
}

bool LoadOP::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
   if (getInitialized() == false) return false;
   profile::start("Load",profileID);

   const string arrayType = "celldata";
   const uint64_t arraySize = simClasses->pargrid.getNumberOfLocalCells()*block::SIZE;
   const uint64_t vectorSize = 1;

   // Allocate output buffer:
   vector<Real> weights(arraySize);

   // Calculate computational load for each block and store values to buffer:
   for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
      const Real weight = simClasses->pargrid.getCellWeights()[block];
      for (int i=0; i<block::SIZE; ++i) weights[block*block::SIZE+i] = weight;
   }

   // Write buffer to output file:
   map<string,string> attributes;
   attributes["name"] = getName();
   attributes["mesh"] = spatMeshName;
   attributes["type"] = arrayType;

   bool success = true;
   if (simClasses->vlsv.writeArray("VARIABLE",attributes,arraySize,vectorSize,&(weights[0])) == false) success = false;
   profile::stop();
   return success;
}


