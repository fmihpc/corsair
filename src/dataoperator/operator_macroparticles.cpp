/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2014 Finnish Meteorological Institute
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

#include <main.h>
#include "operator_macroparticles.h"

using namespace std;

MacroParticlesOperator::MacroParticlesOperator(): DataOperator() {
   #ifdef PROFILE
      totalTimeID = -1;
   #endif
}

MacroParticlesOperator::~MacroParticlesOperator() {finalize();}

bool MacroParticlesOperator::finalize() {return true;}

std::string MacroParticlesOperator::getName() const {
   return "N_macroparticles";
}

bool MacroParticlesOperator::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   return DataOperator::initialize(cr,sim,simClasses);
}

bool MacroParticlesOperator::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
   bool success = true;
   #ifdef PROFILE
      profile::start("MacroParticles",totalTimeID);
   #endif

   const pargrid::CellID N_cells = corsair::getObjectWrapper().simClasses.pargrid.getNumberOfLocalCells();
   uint32_t* buffer = new uint32_t[N_cells*block::SIZE];
   for (pargrid::CellID blockLID=0; blockLID<N_cells*block::SIZE; ++blockLID) {
      buffer[blockLID] = 0;
   }

   // Sum up total number of macroparticles in each block:
   for (size_t species=0; species<particleLists.size(); ++species) {
      const uint32_t* N_particles = particleLists[species]->getParticleNumberArray();

      for (pargrid::CellID blockLID=0; blockLID<N_cells; ++blockLID) {
	 for (int32_t c=0; c<block::SIZE; ++c) {
	    buffer[blockLID*block::SIZE+c] += N_particles[blockLID];
	 }
      }      
   }

   // Write data to output file:
   map<string,string> attribs;
   attribs["name"] = getName();
   attribs["mesh"] = spatMeshName;
   attribs["type"] = "celldata";
   if (simClasses->vlsv.writeArray("VARIABLE",attribs,N_cells,1,buffer) == false) {
      simClasses->logger << "\t ERROR failed to write number of macroparticles" << endl;
      success = false;
   }

   delete [] buffer; buffer = NULL;

   #ifdef PROFILE
      profile::stop();
   #endif
   return success;
}

