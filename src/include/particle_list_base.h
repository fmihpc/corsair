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

#ifndef PARTICLE_LIST_BASE_H
#define PARTICLE_LIST_BASE_H

#include <mpi.h>
#include <vector>

#include <pargrid.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

// Forward declarations of classes stored in object factories:
struct ObjectFactories;
class ParticleAccumulatorBase;
class ParticleBoundaryCondBase;
class ParticleInjectorBase;
class ParticlePropagatorBase;

/** Base class for particle propagators.
 * 
 * Note that this base class cannot contain any information on the particles 
 * being propagated because all particle lists are stored as pointers on a 
 * vector and thus must have the same type. For example, ParticleListSkeleton 
 * pointers cannot be stored on a vector unless they all have the same template 
 * parameters (see particle_list_skeleton.h).
 * 
 * This class replaces class ParticleList.*/
class ParticleListBase {
 public: 
   ParticleListBase();
   virtual ~ParticleListBase();

   virtual bool finalize();
   bool getInitialized() const;
   virtual const std::string& getName() const;
   ParticleAccumulatorBase* getAccumulator() const;
   virtual void getAccumulator(std::string& name,std::string& configRegion) const;
   ParticleBoundaryCondBase* getBoundaryCondition() const;
   virtual void getBoundaryCondition(std::string& name,std::string& configRegion) const;
   ParticleInjectorBase* getInjector() const;
   virtual void getInjector(std::string& name,std::string& configRegion) const;
   std::string getInjectorType() const;
   uint32_t* getParticleNumberArray() const;
   virtual bool getParticles(pargrid::DataID& particleDataID) const;
   ParticlePropagatorBase* getPropagator() const;
   virtual void getPropagator(std::string& name,std::string& configRegion) const;
   virtual const std::string& getSpeciesName() const;
   virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			   const ObjectFactories& objectFactories,const std::string& speciesName);

   // ******************************************************************** //
   // ***** PURE VIRTUAL FUNCTIONS. THESE MUST BE DEFINED BY CLASSES ***** //
   // *****             THAT INHERIT PARTICLELISTBASE                ***** //
   // ******************************************************************** //
   
   /** Accumulate particles on boundary cells.
    * @return If true, all particles on boundary cells were accumulated successfully.*/
   virtual bool accumulateBoundaryCells() = 0;
   
   /** Accumulate particles on inner cells.
    * @return If true, all particles on inner cells were accumulated successfully.*/
   virtual bool accumulateInnerCells() = 0;
   
   /** Apply boundary conditions to all particles on exterior cells, i.e. cells 
    * that are on the boundary of the simulation domain. Typically particles in 
    * these cells are removed.
    * @return If true, boundary conditions were successfully applied.*/
   virtual bool applyBoundaryConditions() = 0;
   
   /** Clear arrays used in particle accumulation. This step can be done while 
    * we are waiting for particle transfers to complete.
    * @return If true, accumulation arrays were cleared successfully.*/
   virtual bool clearAccumulationArrays() = 0;

   /** Get class defining the propagated particle species.
    * @return Particle species.*/
   virtual const void* getSpecies() const = 0;
   
   /** Get a description of the type of species propagated by this particle list.
    * @return Description of species.*/
   virtual const std::string& getSpeciesType() const = 0;
   
   /** Inject new particles.
    * @param If true, particle injection was completed successfully.*/
   virtual bool injectParticles() = 0;
   
   /** Propagate particles on boundary cells.
    * @return If true, all particles were propagated successfully.*/
   virtual bool propagateBoundaryCellParticles() = 0;
   
   /** Propagate particles on inner cells.
    * @return If true, all particles were propagated successfully.*/
   virtual bool propagateInnerCellParticles() = 0;

   /** Get the size of ParticleListBase, i.e. the total number of particles 
    * stored on all local cells.
    * @return Total number of particles stored on this ParticleList on this process.*/
   virtual size_t size() const = 0;
   
   virtual bool updateAfterRepartition() = 0;
   
   /** Wait for all particle sends to complete.
    * @return If true, all particles were sent successfully.*/
   virtual bool waitParticleSends() = 0;
   
   /** Write particles to output file.
    * @param spatMeshName Name of the spatial mesh. This parameter can be ignored.
    * @return If true, particles were written successfully.*/
   virtual bool writeParticles(const std::string& spatMeshName) = 0;
   
 protected:
   
   // ***** VARIABLES RELATED TO MPI ***** //
   MPI_Datatype particleType;                           /**< MPI datatype used to transfer particles.*/
   
   // ***** VARIABLES RELATED TO PARTICLE DATA ***** //
   pargrid::DataID particleDataID;                      /**< ID of ParGrid dynamic data array that contains the 
							 * particles this ParticleList propagates.*/
   uint32_t* N_particles;
   
   // ***** MISCELLANEOUS CLASS VARIABLES ****** //
   bool initialized;                                    /**< If true, class implementing ParticleListBase has initialized successfully.*/
   Simulation* sim;                                     /**< Struct containing generic simulation variables.*/
   SimulationClasses* simClasses;                       /**< Struct containing generic simulation classes, incuding parallel grid.*/
   ConfigReader* cr;                                    /**< Pointer to configuration file reader.*/
   std::string accumulatorName;                         /**< Name of accumulator class used by this particle list.*/
   std::string accumulatorParams;                       /**< Name of configuration file region containing parameters for accumulator class.*/
   std::string bcondName;                               /**< Name of particle boundary condition class used by this particle list.*/
   std::string bcondParams;                             /**< Name of configuration file region containing parameters for boundary condition class.*/
   std::string injectorName;                            /**< Name of injector used by this particle list.*/
   std::string injectorParams;                          /**< Name of configuration file region containing parameters for injector class.*/
   std::string propagatorName;                          /**< Name of propagator used by this particle list.*/
   std::string propagatorParams;                        /**< Name of configuration file region containing parameters for propagator class.*/
   std::string speciesName;                             /**< Name of particle species.*/

   ParticleAccumulatorBase* accumulator;                /**< Particle Accumulator.*/
   ParticleBoundaryCondBase* boundaryCond;              /**< Particle boundary condition.*/
   ParticleInjectorBase* injector;                      /**< Particle injector.*/
   ParticlePropagatorBase* propagator;                  /**< Particle propagator.*/
   
   // ***** VARIABLES RELATED TO PROFILING ***** //
   #if PROFILE_LEVEL > 0
      std::string profileName;       /**< Profiler label for this ParticleList.*/
   
      int accumID;                   /**< Total time spent in particle accumulation.*/
      int accumArrayClears;          /**< Time spen in clearing accumulation array(s).*/
      int accumMPIWaits;             /**< Time spent in MPI waits during accumulation.*/
      int accumMPIOverheads;         /**< Overhead due to MPI during accumulation.*/
      int bcondsID;                  /**< Time spent in applying particle boundary conditions.*/
      int totalTimeID;               /**< Total time used by ParticleList.*/
      int mpiOverheadID;             /**< Time spent in miscellaneous stuff related to parallelization.*/
      int mpiBufferCopyID;           /**< Time spent in copying particles from receive buffer(s).*/
      int mpiWaitRecvsID;            /**< Time spent in waiting for particles to arrive from remote processes.*/
      int mpiWaitSendsID;            /**< Time spent in waiting for MPI sends to complete.*/
      int particleInjectionID;       /**< Time spent in injecting particles.*/
      int particlePropagationID;     /**< Time spent in propagating particles.*/
      int particleTotalID;           /**< Total time spent in propagating particles.*/
      int particleWriteID;           /**< Time spent in writing particles to file.*/
   #endif
};

#endif
