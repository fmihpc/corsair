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

#include "operator_mpirank.h"

using namespace std;

MPIRank::MPIRank(): DataOperator() { 
   mpiRankID = -1;
}

MPIRank::~MPIRank() { }

bool MPIRank::finalize() {return true;}

std::string MPIRank::getName() const {return "MPI_rank";}

bool MPIRank::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   return DataOperator::initialize(cr,sim,simClasses);
}

bool MPIRank::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
   if (getInitialized() == false) return false;
   profile::start("MPIRank",mpiRankID);

   const string arrayType = "celldata";
   const uint64_t arraySize = simClasses->pargrid.getNumberOfLocalCells()*block::SIZE;
   const uint64_t vectorSize = 1;

   // Allocate output buffer:
   int* buffer = new int[arraySize*vectorSize];

   // Copy data to buffer:
   const int rank = sim->mpiRank;
   for (size_t cell=0; cell<simClasses->pargrid.getNumberOfLocalCells(); ++cell) {
      for (int i=0; i<block::SIZE; ++i) buffer[cell*block::SIZE+i] = rank;
   }

   // Write buffer to output file:
   map<string,string> attributes;
   attributes["name"] = getName();
   attributes["mesh"] = spatMeshName;
   attributes["type"] = arrayType;

   bool success = true;
   if (simClasses->vlsv.writeArray("VARIABLE",attributes,arraySize,vectorSize,buffer) == false) success = false;
   delete [] buffer;
   profile::stop();
   return success;
}

