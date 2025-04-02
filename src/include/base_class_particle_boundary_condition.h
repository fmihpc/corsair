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

#ifndef BASE_CLASS_PARTICLE_BOUNDARY_CONDITION_H
#define BASE_CLASS_PARTICLE_BOUNDARY_CONDITION_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <particle_list_base.h>

/** Base class for particle injectors. All user-defined injectors must 
 * implement the interface defined here. In addition to calling base class 
 * constructor, all classes inheriting ParticleBoundaryCondBase should also call 
 * ParticleBoundaryCondBase::initialize in their initialize function to 
 * set member variables sim and simClasses point to correct places.*/
class ParticleBoundaryCondBase {
 public:
   /** Default constructor (empty).*/
   ParticleBoundaryCondBase();

   /** Default virtual destructor (empty).*/
   virtual ~ParticleBoundaryCondBase();

   /** Add all configuration file items read by this class to given config file reader.
    * Base class implementation simply returns value 'true'.
    * @param cr Configuration file reader.
    * @param regionName Name of configuration file region that contains 
    * parameters for this particle injector.
    * @return If true, configuration file items were added successfully.*/
   virtual bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);

   /** Apply boundary conditions to all particles in given cells. Typically one should remove 
    * all particles on simulation boundary cells, i.e., set N_particles[blockLID]=0 and remove 
    * the particles from ParGrid dynamic data array containing particles.
    * Base class apply function simply removes all particles from simulation boundary blocks.
    * @param particleDataID Data ID of ParGrid dynamic array Containing particles.
    * @param N_particles Current number of particles in each cell.
    * @param cells List of simulation boundary cells on this process, identified by their local IDs.
    * @return If true, boundary conditions were successfully applied to all particles.*/
   virtual bool apply(pargrid::DataID particleDataID,unsigned int* N_particles,const std::vector<pargrid::CellID>& cells);

   /** Finalize class. This function should deallocate all internal memory.
    * Base class version simply returns value 'true'.
    * @return If true, class was finalized successfully.*/
   virtual bool finalize();

   /** Initialize particle injector. All inheriting classes should call base class initialize
    * function to set member variables ParticleBoundaryCondBase::sim and ParticleBoundaryCondBase::simClasses.
    * @param sim Generic simulation control variables.
    * @param simClasses Generic simulation classes.
    * @param cr Configuration file reader.
    * @param regionName Name of configuration file region that contains parameters for this class.
    * @param plist Particle list class that uses this class.
    * @return If true, class was successfully initialized and is ready for use.*/
   virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			   const std::string& regionName,const ParticleListBase* plist);

 protected:
   Simulation* sim;               /**< Generic simulation control variables.*/
   SimulationClasses* simClasses; /**< Generic simulation classes.*/
};

/** Maker function that returns a new object of type ParticleBoundaryCondBase.
 * This function can be registered as a Maker to an ObjectFactoryGeneric object factory.
 * @return New ParticleBoundaryCondBase object.*/
ParticleBoundaryCondBase* ParticleBoundaryCondBaseMaker();

#endif
