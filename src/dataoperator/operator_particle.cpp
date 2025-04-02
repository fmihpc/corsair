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

#include "operator_particle.h"

using namespace std;

ParticleOperator::ParticleOperator(): DataOperator() {
   #ifdef PROFILE
      totalTimeID = -1;
   #endif
}

ParticleOperator::~ParticleOperator() {finalize();}

bool ParticleOperator::finalize() {return true;}

std::string ParticleOperator::getName() const {
   return "particles";
}

bool ParticleOperator::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   return DataOperator::initialize(cr,sim,simClasses);
}

bool ParticleOperator::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
   bool success = true;
   #ifdef PROFILE
      profile::start("Particles",totalTimeID);
   #endif

   // TODO: Write only the species defined in config file:
   //size_t counter = 0;
   for (size_t species=0; species<particleLists.size(); ++species) {
      if (particleLists[species]->writeParticles(spatMeshName) == false) success = false;
      /* const size_t N_particles = particles[species].size();
      Real* buffer = new Real[N_particles*3];
      const double* coordinates = getBlockCoordinateArray(*sim,*simClasses);

      counter = 0;
      const unsigned int* sizes = NULL;
      const Particle** particleLists = NULL;
      particles[species].getParticles(sizes,particleLists);
      for (size_t cell=0; cell<simClasses->pargrid.getNumberOfLocalCells(); ++cell) {
	 for (unsigned int p=0; p<sizes[cell]; ++p) {
	    buffer[counter*3+0] = particleLists[cell][p].state[0] + coordinates[3*cell+0];
	    buffer[counter*3+1] = particleLists[cell][p].state[1] + coordinates[3*cell+1];
	    buffer[counter*3+2] = particleLists[cell][p].state[2] + coordinates[3*cell+2];
	    ++counter;
	 }
      }

      map<string,string> attribs;
      attribs["name"] = particles[species].name();
      attribs["type"] = VLSV::MESH_POINT;
      if (simClasses->vlsv.writeArray("MESH",attribs,N_particles,3,buffer) == false) {
	 simClasses->logger << "\t ERROR failed to write particle species!" << endl;
	 success = false;
      }
      delete [] buffer; buffer = NULL;*/
   }

   #ifdef PROFILE
      profile::stop();
   #endif
   return success;
}

