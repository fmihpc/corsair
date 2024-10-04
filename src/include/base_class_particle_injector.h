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

#ifndef BASE_CLASS_PARTICLE_INJECTOR_H
#define BASE_CLASS_PARTICLE_INJECTOR_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <particle_list_base.h>

/** Base class for particle injectors. All user-defined injectors must 
 * implement the interface defined here. In addition to calling base class 
 * constructor, all classes inheriting ParticleInjectorBase should also call 
 * ParticleInjectorBase::initialize in their initialize function to 
 * set member variables sim and simClasses point to correct places.*/
class ParticleInjectorBase {
 public:
   /** Default constructor (empty).*/
   ParticleInjectorBase();
   
   /** Default virtual destructor (empty).*/
   virtual ~ParticleInjectorBase();
   
   /** Add all configuration file items read by this class to given config file reader.
    * @param cr Configuration file reader.
    * @param regionName Name of configuration file region that contains 
    * parameters for this particle injector.
    * @return If true, configuration file items were added successfully.*/
   virtual bool addConfigFileItems(ConfigReader& cr,const std::string& regionName) = 0;
   
   /** Finalize class. This function should deallocate all internal memory.
    * @return If true, class was finalized successfully.*/
   virtual bool finalize() = 0;

   virtual Real getMaximumSpeed();
   virtual void setType(std::string type);
   virtual std::string getType();
   
   /** Initialize particle injector. All inheriting classes should call base class initialize
    * function to set member variables ParticleInjectorBase::sim and ParticleInjectorBase::simClasses.
    * @param sim Generic simulation control variables.
    * @param simClasses Generic simulation classes.
    * @param cr Configuration file reader.
    * @param regionName Name of configuration file region that contains parameters for this class.
    * @param plist Particle list class that uses this class.
    * @return If true, class was successfully initialized and is ready for use.*/
   virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			   const std::string& regionName,const ParticleListBase* plist);
   
   /** Inject new particles to mesh.
    * @param particleDataID Data ID of ParGrid dynamic array where new particles
    * should be added.
    * @param N_particles Array containing number of existing particles in each cell.
    * Add the number of injected particles per cell to this array.
    * @return If true, no errors occurred in particle injection.*/
   virtual bool inject(pargrid::DataID particleDataID,unsigned int* N_particles) = 0;
   
 protected:
   Simulation* sim;               /**< Generic simulation control variables.*/
   SimulationClasses* simClasses; /**< Generic simulation classes.*/
   std::string type;              // Injector type
};

#endif
