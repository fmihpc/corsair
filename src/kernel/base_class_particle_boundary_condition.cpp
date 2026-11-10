/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2013 Finnish Meteorological Institute
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
#include <limits>
#include <cmath>
#include <mpi.h>

#include <base_class_particle_boundary_condition.h>

using namespace std;

ParticleBoundaryCondBase* ParticleBoundaryCondBaseMaker() {return new ParticleBoundaryCondBase();}

ParticleBoundaryCondBase::ParticleBoundaryCondBase() { }

ParticleBoundaryCondBase::~ParticleBoundaryCondBase() { }

bool ParticleBoundaryCondBase::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
   return true;
}

bool ParticleBoundaryCondBase::apply(pargrid::DataID particleDataID,unsigned int* N_particles,
				     const std::vector<pargrid::CellID>& boundaryBlocks) {
   // Get particles:
   pargrid::DataWrapper<char> wrapper = simClasses->pargrid.getUserDataDynamic<char>(particleDataID);

   // Remove particles from given blocks:
   Real t_propag = 0.0;
   for (size_t b=0; b<boundaryBlocks.size(); ++b) {
      // Measure computation time if we are testing for repartitioning:
      if (sim->countPropagTime == true) t_propag = MPI_Wtime();

      const pargrid::CellID blockLID = boundaryBlocks[b];
      wrapper.resize(blockLID,0);
      N_particles[blockLID] = 0;

      // Store block injection time:
      if (sim->countPropagTime == true) {
	 t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	 simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
      }
   }

   return true;
}

bool ParticleBoundaryCondBase::finalize() {
   return true;
}

bool ParticleBoundaryCondBase::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
					  const std::string& regionName,const ParticleListBase* plist) {
   this->sim = &sim;
   this->simClasses = &simClasses;
   return true;
}

