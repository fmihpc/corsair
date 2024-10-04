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

#include <particle_list_base.h>
#include <base_class_particle_accumulator.h>
#include <base_class_particle_boundary_condition.h>
#include <base_class_particle_injector.h>
#include <base_class_particle_propagator.h>

using namespace std;

/** Default constructor for ParticleListBase. Initializes some 
 * internal variables to dummy values. ParticleListBase::initialize must be 
 * called before ParticleListBase can be used.*/
ParticleListBase::ParticleListBase() {
   initialized = false;
   particleDataID = pargrid::INVALID_DATAID;
   N_particles = NULL;
   sim = NULL;
   simClasses = NULL;
   cr = NULL;
   
   accumulator = NULL;
   boundaryCond = NULL;
   injector = NULL;
   propagator = NULL;
   
   #ifdef PROFILE
      accumID = -1;
      accumArrayClears = -1;
      accumMPIWaits = -1;
      accumMPIOverheads = -1;
      bcondsID = -1;
      totalTimeID = -1;
      mpiBufferCopyID = -1;
      mpiOverheadID = -1;
      mpiWaitRecvsID = -1;
      mpiWaitSendsID = -1;
      particleInjectionID = -1;
      particlePropagationID = -1;
      particleWriteID = -1;
      particleTotalID = -1;
   #endif
}

/** Virtual destructor for base class ParticleListBase.
 * The destructor calls ParticleListBase::finalize().*/
ParticleListBase::~ParticleListBase() {
   finalize();
}

/** Finalizer for ParticleListBase. This function deallocates ParGrid data arrays
 * and frees committed MPI datatypes.
 * @return If true, ParticleListBase has finalized successfully.*/
bool ParticleListBase::finalize() {
   bool success = true;
   if (initialized == false) return success;
   
   // Deallocate ParGrid static data arrays:
   if (particleDataID != pargrid::INVALID_DATAID) {
      if (simClasses->pargrid.removeUserData(particleDataID) == false) success = false;
      particleDataID = pargrid::INVALID_DATAID;
   }
   delete [] N_particles; N_particles = NULL;

   // Deallocate particle injector etc:
   if (accumulator != NULL) accumulator->finalize();
   if (boundaryCond != NULL) boundaryCond->finalize();
   if (injector != NULL) injector->finalize();
   if (propagator != NULL) propagator->finalize();
   delete accumulator; accumulator = NULL;
   delete boundaryCond; boundaryCond = NULL;
   delete injector; injector = NULL;
   delete propagator; propagator = NULL;
   
   // Free MPI datatype that was used to transfer the particles:
   if (particleType != MPI_DATATYPE_NULL) MPI_Type_free(&particleType);
   
   initialized = false;
   return success;
}

ParticleAccumulatorBase* ParticleListBase::getAccumulator() const {return accumulator;}

void ParticleListBase::getAccumulator(std::string& name,std::string& configRegion) const {
   name = accumulatorName;
   configRegion = accumulatorParams;
}

ParticleBoundaryCondBase* ParticleListBase::getBoundaryCondition() const {return boundaryCond;}

void ParticleListBase::getBoundaryCondition(std::string& name,std::string& configRegion) const {
   name = bcondName;
   configRegion = bcondParams;
}

/** Query the initialization status of ParticleList.
 * @return If true, ParticleList has initialized successfully and is ready for use.*/
bool ParticleListBase::getInitialized() const {
   return initialized;
}

ParticleInjectorBase* ParticleListBase::getInjector() const {return injector;}

void ParticleListBase::getInjector(std::string& name,std::string& configRegion) const {
   name = injectorName;
   configRegion = injectorParams;
}

string ParticleListBase::getInjectorType() const {
   return injector->getType();
}

/** Get the name of particle species propagated by this particle list.
 * @return Name of propagated particle species.*/
const std::string& ParticleListBase::getName() const {return speciesName;}

/** Get pointer to array that contains number of particles in each block.
 * @return Pointer to array.*/
uint32_t* ParticleListBase::getParticleNumberArray() const {return N_particles;}

/** Get ParGrid data ID associated with the particles propagated
 * by this particle list.
 * @param particleDataID ParGrid DataID of particles is copied here.
 * @return If true, particleDataID contains a valid data ID upon exit.*/
bool ParticleListBase::getParticles(pargrid::DataID& particleDataID) const {
   particleDataID = this->particleDataID;
   return getInitialized();
}

ParticlePropagatorBase* ParticleListBase::getPropagator() const {return propagator;}

void ParticleListBase::getPropagator(std::string& name,std::string& configRegion) const {
   name = propagatorName;
   configRegion = propagatorParams;
}

/** Get the name of particle species.
 * @return Name of particle species propagated by this particle list.*/
const std::string& ParticleListBase::getSpeciesName() const {
   return speciesName;
}

/** Initialize ParticleListBase.
 * @param sim Generic simulation variables.
 * @param simClasses Generic simulation classes, including the parallel mesh.
 * @param cr Configuration file reader.
 * @return If true, ParticleListBase initialized successfully and is ready for use.*/
bool ParticleListBase::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
				  const ObjectFactories& objectFactories,
				  const std::string& speciesName) {
   this->sim = &sim;
   this->simClasses = &simClasses;
   this->cr = &cr;
   this->speciesName = speciesName;
   return true;
}
