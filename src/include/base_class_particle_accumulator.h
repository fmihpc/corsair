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

#ifndef BASE_CLASS_PARTICLE_ACCUMULATOR_H
#define BASE_CLASS_PARTICLE_ACCUMULATOR_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <particle_list_base.h>

/** Base class for particle injectors. All user-defined accumulators must 
 * implement the interface defined here. In addition to calling base class 
 * constructor, all classes inheriting ParticleAccumulatorBase should also call 
 * ParticleAccumulatorBase::initialize in their initialize function to 
 * set member variables sim and simClasses point to correct places.*/
class ParticleAccumulatorBase {
 public:
   /** Default constructor (empty).*/
   ParticleAccumulatorBase();

   /** Default virtual destructor (empty).*/
   virtual ~ParticleAccumulatorBase();

   /** Accumulate particles on process boundary cells.
    * @param particleDataID Data ID of a ParGrid dynamic data array that contains accumulated particles.
    * @param N_particles Array that contains the number of particles in each cell.
    * @return If true, particles were successfully accumulated.*/
   virtual bool accumulateBoundaryCells(pargrid::DataID particleDataID,const unsigned int* N_particles) = 0;

   /** Accumulate particles on cells that do not have any remote neighbors on other processes.
    * @param particleDataID Data ID of a ParGrid dynamic data array that contains accumulated particles.
    * @param N_particles Array that contains the number of particles in each cell.
    * @return If true, particles were successfully accumulated.*/
   virtual bool accumulateInnerCells(pargrid::DataID particleDataID,const unsigned int* N_particles) = 0;

   /** Add all configuration file items read by this class to given config file reader.
    * @param cr Configuration file reader.
    * @param regionName Name of configuration file region that contains 
    * parameters for this particle accumulator.
    * @return If true, configuration file items were added successfully.*/
   virtual bool addConfigFileItems(ConfigReader& cr,const std::string& regionName) = 0;

   /** Add updates, received from neighbor processes, to local accumulation array(s).
    * @return If true, all remote updates were successfully added to local array(s).*/
   virtual bool addRemoteUpdates() = 0;

   /** Set all array(s) used by this accumulator to zero values.*/
   virtual bool clearAccumulationArrays() = 0;

   /** Finalize class. This function should deallocate all internal memory.
    * @return If true, class was finalized successfully.*/
   virtual bool finalize() = 0;

   /** Initialize particle accumulator. All inheriting classes should call base class initialize
    * function to set member variables ParticleAccumulatorBase::sim and ParticleAccumulatorBase::simClasses.
    * @param sim Generic simulation control variables.
    * @param simClasses Generic simulation classes.
    * @param cr Configuration file reader.
    * @param regionName Name of configuration file region that contains parameters for this class.
    * @param plist Particle list class that uses this accumulator.
    * @return If true, class was successfully initialized and is ready for use.*/
   virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			   const std::string& regionName,const ParticleListBase* plist);

   /** Send updates to neighbor processes' accumulation array(s).
    * @return If true, all updates were successfully sent.*/
   virtual bool sendUpdates() = 0;

   /** Wait for remote updates to accumulation array(s) to arrive to this process.
    * @return If true, all updates were successfully received.*/
   virtual bool wait() = 0;

 protected:
   Simulation* sim;               /**< Generic simulation control variables.*/
   SimulationClasses* simClasses; /**< Generic simulation classes.*/
};

#endif
