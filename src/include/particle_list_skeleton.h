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

#ifndef PARTICLE_LIST_SKELETON_H
#define PARTICLE_LIST_SKELETON_H

#include <cstdlib>
#include <list>
#include <map>
#include <set>

#include <pargrid_buffers.h>
#include <pargrid_copy_protocol_mpi.h>
#include <particle_list_base.h>
#include <base_class_particle_injector.h>

// ******************************************************** //
// Functions that class Species must implement.             //
// Species is a template parameter to ParticleListSkeleton. //
// It defined parameters common to all particles,           //
// such as mass, charge, etc.                               //
//                                                          //
// bool finalize();                                         //                                                         
// string getName() const;                                  //
// bool readParameters(Simulation& sim,                     //
//                     SimulationClasses& simClasses,       //
//                     ConfigReader& cr,                    //
//                     const std::string& speciesName);     //
//                                                          //
// ******************************************************** //

/** Class ParticleListSkeleton is a minimal implementation of the interface described 
 * in class ParticleListBase. This class is a general-use class in the sense that all 
 * physics is done by template parameter classes.
 * 
 * Template class PARTICLE must contain a public member variable (array) state whose 
 * first three elements are the x, y, and z-coordinates of the particle. Most important
 * reason for this requirement is that this seems to give the best performance with g++ 
 * compiler by a far margin.
 * 
 * Template parameters are:
 * SPECIES Class that defines parameters common to all particles, such as mass, charge, etc.
 * PARTICLE Class that defines propagated particles.
 * INJ Class that injects particles.
 * PROP Class that propagates particles forward in time.
 * ACCUM Class that accumulates particle-related quantities.
 */
template<class SPECIES,class PARTICLE>
class ParticleListSkeleton: public ParticleListBase {
 public:
   ParticleListSkeleton();
   virtual ~ParticleListSkeleton();

   virtual bool accumulateBoundaryCells();
   virtual bool accumulateInnerCells();
   virtual bool applyBoundaryConditions();
   virtual bool clearAccumulationArrays();
   virtual bool finalize();
   virtual const void* getSpecies() const;
   virtual const std::string& getSpeciesType() const;
   virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			   const ObjectFactories& objectFactories,const std::string& speciesName);
   virtual bool injectParticles();
   virtual bool propagateBoundaryCellParticles();
   virtual bool propagateInnerCellParticles();
   virtual size_t size() const;
   virtual bool updateAfterRepartition();
   virtual bool writeParticles(const std::string& spatMeshName);
   virtual bool waitParticleSends();

 protected:

   // ***** VARIABLES RELATED TO SENDING AND RECEIVING PARTICLES ***** //
   typedef pargrid::Buffer<PARTICLE> BUFFER;
   std::vector<BUFFER> inBuffers;
   std::vector<BUFFER> outBuffers;
   std::vector<pargrid::CopyProtocolMPI<BUFFER> > inBufferCopy;
   std::vector<pargrid::CopyProtocolMPI<BUFFER> > outBufferCopy;

   static unsigned int LIDsCalculatedTimestep;
   static std::vector<std::vector<pargrid::CellID> > recvLIDs;
   static std::vector<std::vector<pargrid::CellID> > sendLIDs;

   // ***** MISCELLANEOUS CLASS VARIABLES ***** //
   SPECIES species;                   /**< Particle species, i.e. parameters common to all particles.*/
   unsigned int meshVariablesUpdated; /**< Time step at which updatePartitioning was last called, this variable 
				       * is used to prevent multiple calls at same time step.*/

   void propagateCell(pargrid::CellID cellID,const double* coordinates,
		      pargrid::DataWrapper<PARTICLE>& wrapper,unsigned int* N_particles);
   bool updatePartitioning();
};

const int XPOS = 0;
const int YPOS = 1;
const int ZPOS = 2;

const Real bufferIncrementFactor = 1.2;
const int TAG_FIRST_MESSAGE = 0;
const int TAG_SECOND_MESSAGE = 1;
const int TAG_THIRD_MESSAGE = 2;

// Declare static member functions:
template<class SPECIES,class PARTICLE>
unsigned int ParticleListSkeleton<SPECIES,PARTICLE>::LIDsCalculatedTimestep = std::numeric_limits<unsigned int>::max();
template<class SPECIES,class PARTICLE>
std::vector<std::vector<pargrid::CellID> > ParticleListSkeleton<SPECIES,PARTICLE>::recvLIDs;
template<class SPECIES,class PARTICLE> 
std::vector<std::vector<pargrid::CellID> > ParticleListSkeleton<SPECIES,PARTICLE>::sendLIDs;

/** Default constructor for ParticleListSkeleton.*/
template<class SPECIES,class PARTICLE> inline
ParticleListSkeleton<SPECIES,PARTICLE>::ParticleListSkeleton(): ParticleListBase() { 
   meshVariablesUpdated = std::numeric_limits<unsigned int>::max();
}

/** Destructor for ParticleListSkeleton. Calls ParticleListSkeleton::finalize().*/
template<class SPECIES,class PARTICLE> inline
ParticleListSkeleton<SPECIES,PARTICLE>::~ParticleListSkeleton() {
   finalize();
   ParticleListBase::finalize();
}

/** Accumulate all particles on boundary cells. The actual accumulation is done 
 * by the template parameter class ACCUM.
 * @return If true, particles on boundary cells were accumulated successfully.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::accumulateBoundaryCells() {
   bool success = true;
   if (initialized == false) return false;

   #if PROFILE_LEVEL > 0
      profile::start(profileName,totalTimeID);
   #endif
   if (sim->meshRepartitioned == true) updatePartitioning();

   #if PROFILE_LEVEL > 0
      profile::start("Accumulator",accumID);
   #endif

   // Accumulate particle quantities:
   if (accumulator->accumulateBoundaryCells(particleDataID,N_particles) == false) success = false;

   // Send locally calculated update(s) to remote neighbour(s):
   #if PROFILE_LEVEL > 0
      profile::start("MPI overhead",accumMPIOverheads);
   #endif
   if (accumulator->sendUpdates() == false) success = false;
   #if PROFILE_LEVEL > 0
      profile::stop(); // MPI overhead
   #endif

   #if PROFILE_LEVEL > 0
      profile::stop(); // accumulation step
      profile::stop(); // total time
   #endif

   return success;
}

/** Accumulate all particles on inner cells. The actual accumulation 
 * is done by the template parameter class ACCUM.
 * @return If true, all particles on inner cells were accumulated successfully.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::accumulateInnerCells() {
   bool success = true;
   if (initialized == false) return false;
   #if PROFILE_LEVEL > 0
      profile::start(profileName,totalTimeID);
   #endif
   if (sim->meshRepartitioned == true) updatePartitioning();

   // Accumulate particle quantities:
   #if PROFILE_LEVEL > 0
      profile::start("Accumulator",accumID);
   #endif

   if (accumulator->accumulateInnerCells(particleDataID,N_particles) == false) success = false;

   // Wait for updates to arrive from remote neighbours:
   #if PROFILE_LEVEL > 0
      profile::start("MPI waits",accumMPIWaits);
   #endif
   if (accumulator->wait() == false) success = false;
   #if PROFILE_LEVEL > 0
      profile::stop(); // MPI Waits
   #endif

   // Add updates received from remote neighbour(s) to local mesh variable(s):
   #if PROFILE_LEVEL > 0
      profile::start("MPI overhead",accumMPIOverheads);
   #endif
   if (accumulator->addRemoteUpdates() == false) success = false;
   #if PROFILE_LEVEL > 0
      profile::stop(); // MPI overhead
   #endif

   #if PROFILE_LEVEL > 0
      profile::stop(); // Accumulator
      profile::stop(); // Total time
   #endif
   return success;
}

/** Apply boundary conditions on particles on cells that are on the boundary 
 * of the simulation domain, i.e. exterior cells. This function simply removes 
 * such particles.
 * @return If true, boundary conditions were applied successfully.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::applyBoundaryConditions() {
   bool success = true;
   if (initialized == false) return false;

   #if PROFILE_LEVEL > 0
      profile::start(profileName,totalTimeID);
      profile::start("Particle Step",particleTotalID);
      profile::start("boundary conds.",bcondsID);
   #endif

   if (sim->meshRepartitioned == true) updatePartitioning();
   if (boundaryCond->apply(particleDataID,N_particles,simClasses->pargrid.getExteriorCells()) == false) success = false;

   #if PROFILE_LEVEL > 0
      profile::stop(); // boundary conds.
      profile::stop(); // particle total
      profile::stop(); // total time
   #endif
   return success;
}

/** Clear ParGrid data arrays used in accumulation to zero value. This step can be completed 
 * while were are waiting for particles to arrive from remote processes. 
 * The actual clearing is done by the template parameter class ACCUM.
 * @return If true, accumulation arrays were cleared successfully.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::clearAccumulationArrays() {
   if (initialized == false) return false;

   #if PROFILE_LEVEL > 0
      profile::start(profileName,totalTimeID);
      profile::start("Accumulator",accumID);
      profile::start("array clearing",accumArrayClears);
   #endif

   if (sim->meshRepartitioned == true) updatePartitioning();
   const bool success = accumulator->clearAccumulationArrays();

   #if PROFILE_LEVEL > 0
      profile::stop();
      profile::stop();
      profile::stop();
   #endif
   return success;
}

/** Finalizer function. This function calls finalize()-functions of 
 * template parameter classes and ParticleListBase::finalize().
 * @return If true, ParticleListSkeleton has successfully deallocated all internal memory.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::finalize() {
   bool success = true;
   if (initialized == false) return true;

   // Cleanly deallocate buffer copy protocols:
     {
	std::vector<pargrid::CopyProtocolMPI<BUFFER> > dummy1;
	std::vector<pargrid::CopyProtocolMPI<BUFFER> > dummy2;
	inBufferCopy.swap(dummy1);
	outBufferCopy.swap(dummy2);
     }

   // Cleanly deallocate buffers:
     {
	std::vector<BUFFER> dummy1;
	std::vector<BUFFER> dummy2;
	inBuffers.swap(dummy1);
	outBuffers.swap(dummy2);
     }

   // Deallocate memory and call finalizers:
   if (species.finalize() == false) success = false;

   // ParticleListBase finalize() deallocates ParGrid data arrays 
   // and frees MPI datatype.
   if (ParticleListBase::finalize() == false) success = false;
   return success;
}

/** Get propagated particle species.
 * @return Pointer to particle species.*/
template<class SPECIES,class PARTICLE> inline
const void* ParticleListSkeleton<SPECIES,PARTICLE>::getSpecies() const {
   return &species;
}

/** Get a description of particle species propagated by this particle list.
 * @return Description of species.*/
template<class SPECIES,class PARTICLE> inline
const std::string& ParticleListSkeleton<SPECIES,PARTICLE>::getSpeciesType() const {
   return species.getSpeciesType();
}

/** Initializer for ParticleListBase. This function reads parameters from 
 * config file, allocates ParGrid data arrays, and calls initialize-functions of 
 * template parameter classes. Particle injector is also called and it is assumed that 
 * the injector will create the initial state.
 * @param sim Generic simulation variables.
 * @param simClasses Generic simulation classes, including the parallel mesh.
 * @param cr Configuration file reader.
 * @param speciesName Name of config file region that contains parameters read by ParticleListSkeleton.
 * @return If true, ParticleListSkeleton initialized successfully and is ready for use.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								   const ObjectFactories& objectFactories,const std::string& speciesName) {
   bool success = true;
   if (initialized == true) return true;
   if (ParticleListBase::initialize(sim,simClasses,cr,objectFactories,speciesName) == false) success = false;

   // Read species' parameters from config file. Since we do not know what 
   // the parameters are (species is given as a template parameter), we 
   // ask the species class to read them:
   if (species.readParameters(sim,simClasses,cr,speciesName) == false) success = false;

   cr.add(speciesName+".accumulator.name","Name of the particle accumulator.",std::string(""));
   cr.add(speciesName+".boundary_condition.name","Name of particle boundary condition.",std::string("DefaultBoundaryCond"));
   cr.add(speciesName+".injector.name","Name of the particle injector.",std::string(""));
   cr.add(speciesName+".propagator.name","Name of the particle propagator",std::string(""));
   cr.add(speciesName+".accumulator.parameters","Name of the config region that defines parameters for accumulator.",std::string(""));
   cr.add(speciesName+".boundary_condition.parameters","Name of config file region that defined parameters for boundary condition.",std::string(""));
   cr.add(speciesName+".injector.parameters","Name of the config region that defines parameters for injector.",std::string(""));
   cr.add(speciesName+".propagator.parameters","Name of the config region that defines parameters for propagator.",std::string(""));
   cr.parse();

   cr.get(speciesName+".injector.name",injectorName);
   cr.get(speciesName+".injector.parameters",injectorParams);
   cr.get(speciesName+".propagator.name",propagatorName);
   cr.get(speciesName+".propagator.parameters",propagatorParams);
   cr.get(speciesName+".boundary_condition.name",bcondName);
   cr.get(speciesName+".boundary_condition.parameters",bcondParams);
   cr.get(speciesName+".accumulator.name",accumulatorName);
   cr.get(speciesName+".accumulator.parameters",accumulatorParams);

   // Attempt to create particle injector:
   injector = objectFactories.particleInjectors.create(injectorName);
   if (injector == NULL) {
      simClasses.logger << "(PARTICLELIST) ERROR: Particle injector '" << injectorName << "' does not exist!" << std::endl << write;
      success = false;
   }

   // Read injector's parameters from config file:
   if (injector != NULL) {
      if (injector->addConfigFileItems(cr,injectorParams) == false) {
	 simClasses.logger << "(PARTICLELIST) ERROR: Injector failed to add config file items" << std::endl << write;
	 success = false;
      }
      if (injector->initialize(sim,simClasses,cr,injectorParams,this) == false) {
	 simClasses.logger << "(PARTICLELIST) ERROR: Failed to init injector!" << std::endl << write;
	 success = false;
      }
   }

   // Attempt to create particle boundary condition:
   boundaryCond = objectFactories.particleBoundaryConditions.create(bcondName);
   if (boundaryCond == NULL) {
      simClasses.logger << "(PARTICLELIST) ERROR: Particle boundary condition '" << bcondName << "' does not exist!" << std::endl << write;
      success = false;
   }

   // Read particle propagator's parameters from config file:
   if (boundaryCond != NULL) {
      if (boundaryCond->addConfigFileItems(cr,bcondParams) == false) {
	 simClasses.logger << "(PARTICLELIST) ERROR: Boundary condition failed to add config file items" << std::endl << write;
	 success = false;
      }
      if (boundaryCond->initialize(sim,simClasses,cr,bcondParams,this) == false) {
	 simClasses.logger << "(PARTICLELIST) ERROR: Failed to init boundary condition!" << std::endl << write;
	 success = false;
      }
   }

   // Attempt to create particle propagator:
   propagator = objectFactories.particlePropagators.create(propagatorName);
   if (propagator == NULL) {
      simClasses.logger << "(PARTICLELIST) ERROR: Particle propagator '" << propagatorName << "' does not exist!" << std::endl << write;
      success = false;
   }

   // Read particle propagator's parameters from config file:
   if (propagator != NULL) {
      if (propagator->addConfigFileItems(cr,propagatorParams) == false) {
	 simClasses.logger << "(PARTICLELIST) ERROR: Propagator failed to add config file items" << std::endl << write;
	 success = false;
      }      
      if (propagator->initialize(sim,simClasses,cr,propagatorParams,this) == false) {
	 simClasses.logger << "(PARTICLELIST) ERROR: Failed to init propagator!" << std::endl << write;
	 success = false;
      }
   }

   // Attempt to create particle accumulator:
   accumulator = objectFactories.particleAccumulators.create(accumulatorName);
   if (accumulator == NULL) {
      simClasses.logger << "(PARTICLELIST) ERROR: Particle accumulator '" << accumulatorName << "' does not exist!" << std::endl << write;
      success = false;
   }

   // Read accumulator's parameters from config file:
   if (accumulator != NULL) {
      if (accumulator->addConfigFileItems(cr,accumulatorParams) == false) {
	 simClasses.logger << "(PARTICLELIST) ERROR: Accumulator failed to add config file items!" << std::endl << write;
	 success = false;
      }
      if (accumulator->initialize(sim,simClasses,cr,accumulatorParams,this) == false) {
	 simClasses.logger << "(PARTICLELIST) ERROR: Failed to init accumulator!" << std::endl << write;
	 success = false;
      }
   }

   // Clear N_particles array to zero value unless simulation was restarted.
   // In the latter case contents of N_particles are correct:
   delete [] N_particles;
   N_particles = new unsigned int[simClasses.pargrid.getNumberOfLocalCells()];
   if (success == true && sim.restarted == false) {
      for (pargrid::CellID block=0; block<simClasses.pargrid.getNumberOfLocalCells(); ++block) N_particles[block] = 0;
   }

   // Create a ParGrid dynamic data array for the particles:
   particleDataID = simClasses.pargrid.addUserData<PARTICLE>(speciesName+"_particles",0,true);
   if (particleDataID == pargrid::INVALID_DATAID) {
      simClasses.logger << "(PARTICLELIST) ERROR: Failed to create dynamic data array for particle species '";
      simClasses.logger << speciesName << "' !" << std::endl << write;
      success = false;
   }

   // Inject particles unless simulation was restarted. In the 
   // latter case ParGrid dynamic array contains already the particles:
   if (success == true && sim.restarted == false) {
      if (injector->inject(particleDataID,N_particles) == false) success = false;
   }

   // Get an MPI datatype for transferring particles and commit it:
   if (success == true) {
      PARTICLE::getDatatype(particleType);
      MPI_Type_commit(&particleType);
   }

   #if PROFILE_LEVEL > 0
      // Create profiler name:
      profileName = "ParticleList "+species.getName();
   #endif

   if (success == true) initialized = true;
   return initialized;
}

/** Inject new particles. The actual injection is done by 
 * template parameter class INJ. 
 * @return If true, new particles were injected successfully.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::injectParticles() {
   bool success = true;
   if (initialized == false) return false;

   #if PROFILE_LEVEL > 0
      profile::start(profileName,totalTimeID);
      profile::start("Particle Step",particleTotalID);
      profile::start("injection",particleInjectionID);
   #endif

   if (sim->meshRepartitioned == true) updatePartitioning();

   // Inject new particles:
   if (injector->inject(particleDataID,N_particles) == false) success = false;

   #if PROFILE_LEVEL > 0
      profile::stop();
      profile::stop();
      profile::stop();
   #endif
   return success;
}

/** Propagate all particles on boundary cells. The actual propagation is done 
 * by template parameter class PROP, but ParticleListSkeleton does all the work 
 * related to moving particles between processes etc. Specifically, this function 
 * posts MPI receives for incoming particles at the start, and MPI sends for outgoing 
 * particles at the end.
 * @return If true, particles were propagated successfully.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::propagateBoundaryCellParticles() {
   bool success = true;
   if (initialized == false) return false;

   #if PROFILE_LEVEL > 0
      profile::start(profileName,totalTimeID);
      profile::start("Particle Step",particleTotalID);
   #endif

   #if PROFILE_LEVEL > 0
      profile::start("MPI overhead",mpiOverheadID);
   #endif
   // Check if mesh has been repartitioned:
   if (sim->meshRepartitioned == true) updatePartitioning();

   // Post receives for incoming particles:
   MPITypes::rank bufferCounter = 0;
   const std::set<pargrid::MPI_processID>& neighbourProcesses = simClasses->pargrid.getNeighbourProcesses();
   for (std::set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
      #ifndef NDEBUG
         if (bufferCounter >= static_cast<MPITypes::rank>(inBuffers.size())) {
	    std::stringstream ss;
	    ss << "STEP " << sim->timestep << " P#" << sim->mpiRank << " bufferCounter=" << bufferCounter;
	    ss << " is too large, inBuffers.size()=" << inBuffers.size() << std::endl;
	    std::cerr << ss.str();
	    exit(1);
	 }
      #endif

      if (inBufferCopy[bufferCounter].start() == false) {
	 #ifndef NDEBUG
	    simClasses->logger << "(PLIST SKELETON) ERROR: Failed to start data receive" << std::endl << write;
	 #endif
	 success = false;
      }
      ++bufferCounter;
   }

   // Clear send buffers:
   for (size_t i=0; i<outBuffers.size(); ++i) outBuffers[i].getBuffer().clear();

   #if PROFILE_LEVEL > 0
      profile::stop(); // MPI overhead
   #endif

   // Propagate particles on boundary cells:
   Real t_propag = 0.0;
   const double* coordinates = reinterpret_cast<double*>(simClasses->pargrid.getUserData(Simulation::crdsDataID));
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   const std::vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(pargrid::DEFAULT_STENCIL);
   #ifndef NDEBUG
      bool allValid = true;
      if (coordinates == NULL) {std::cerr << "'coordinates' is NULL!" << std::endl; allValid = false;}
      if (N_particles == NULL) {std::cerr << "'N_particles' is NULL!" << std::endl; allValid = false;}
      if (wrapper.valid() == false) {std::cerr << "'wrapper' is not valid!" << std::endl; allValid = false;}
      if (allValid == false) exit(1);
   #endif

   #if PROFILE_LEVEL > 0
      profile::start("propagation",particlePropagationID);
   #endif
   for (size_t c=0; c<boundaryBlocks.size(); ++c) {
      // Measure cell propagation time if we are preparing for repartitioning:
      if (sim->countPropagTime == true) {
	 t_propag = MPI_Wtime();
      }
      #ifndef NDEBUG
         if (boundaryBlocks[c] >= simClasses->pargrid.getNumberOfLocalCells()) {
	    std::stringstream ss;
	    ss << "STEP " << sim->timestep << " P#" << sim->mpiRank << " boundary block LID#";
	    ss << boundaryBlocks[c] << " is out of bounds!" << std::endl;
	    std::cerr << ss.str();
	    exit(1);
	 }
      #endif

      propagateCell(boundaryBlocks[c],coordinates,wrapper,N_particles);

      if (sim->countPropagTime == true) {
	 t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	 simClasses->pargrid.getCellWeights()[boundaryBlocks[c]] += t_propag;
      }
   }   
   #if PROFILE_LEVEL > 0
      profile::stop(); // propagation
      profile::start("MPI overhead",mpiOverheadID);
   #endif

   // Prepare buffers for sending to remote neighbours:
   MPITypes::rank processCounter = 0;
   for (typename std::map<pargrid::MPI_processID,std::set<pargrid::CellID> >::const_iterator 
	it  = simClasses->pargrid.getSends(sim->inverseStencilID).begin();
	it != simClasses->pargrid.getSends(sim->inverseStencilID).end(); 
	++it) {

      #ifndef NDEBUG
      if (processCounter >= static_cast<MPITypes::rank>(outBuffers.size())) {
	 std::stringstream ss;
	 ss << "STEP " << sim->timestep << " P#" << sim->mpiRank << " processCounter=" << processCounter;
	 ss << " is too large, outBuffers.size()=" << outBuffers.size() << std::endl;
	 std::cerr << ss.str();
	 exit(1);
      }
      #endif

      // Iterate over all blocks sent to neighbour process nbrHostID:
      int counter = 0;
      size_t blockCounter = 0;
      for (std::set<pargrid::CellID>::const_iterator block=it->second.begin(); block!=it->second.end(); ++block) {
	 const pargrid::CellID blockLID = sendLIDs[processCounter][blockCounter];

	 #ifndef NDEBUG
	    if (blockLID >= simClasses->pargrid.getNumberOfAllCells()) {
	       std::stringstream ss;
	       ss << "STEP " << sim->timestep << " P#" << sim->mpiRank << " block LID#" << blockLID;
	       ss << " is out of bounds!" << std::endl;
	       std::cerr << ss.str();
	       exit(1);
	    }
	 #endif

	 outBuffers[processCounter].getBlockSizes()[counter] = wrapper.size()[blockLID];
	 ++counter;
	 ++blockCounter;

	 // Copy particles from remote cell to output buffer:
	 std::vector<PARTICLE>& buffer = outBuffers[processCounter].getBuffer();
	 for (pargrid::ArraySizetype p=0; p<wrapper.size()[blockLID]; ++p) {
	    buffer.push_back( wrapper.data()[blockLID][p] );
	 }

	 // Clear remote cell:
	 wrapper.resize(blockLID,0);
      }

      // Post send for particles sent to process nbrHostID:
      if (outBufferCopy[processCounter].start() == false) {
	 #ifndef NDEBUG
	    simClasses->logger << "(PLIST SKELETON) ERROR: Failed to start data send" << std::endl << write;
	 #endif
	 success = false;
      }
      ++processCounter;
   }

   #if PROFILE_LEVEL > 0
      profile::stop(); // MPI overhead
      profile::stop(); // particle step
      profile::stop(); // total time
   #endif

   return success;
}

template<class SPECIES,class PARTICLE> inline
void ParticleListSkeleton<SPECIES,PARTICLE>::propagateCell(pargrid::CellID cellID,
								      const double* coordinates,
								      pargrid::DataWrapper<PARTICLE>& wrapper,
								      unsigned int* N_particles) {
   PARTICLE* __restrict__ particleList = wrapper.data()[cellID];
   #ifndef NDEBUG
      if (particleList == NULL && N_particles[cellID] > 0) {
	 std::cerr << "(PARTICLE LIST SKELETON) ERROR: Received NULL pointer for cell LID#" << cellID << std::endl;
	 exit(1);      
      }
   #endif

   // Call external propagator function:
   #if PROFILE_LEVEL > 1
      static int externalProp = -1;
      profile::start("external propagator",externalProp);
   #endif

   propagator->propagateCell(cellID,particleDataID,coordinates+3*cellID,N_particles[cellID]);

   #if PROFILE_LEVEL > 1
      profile::stop();
   #endif

   // Loop over the particles again and check if particles changed cells. If a 
   // particle changed cell, it is either copied to the new cell or copied to 
   // a transfer buffer. The last condition occurs when the new cell is on another process.
   // The algorithm below does not leave any gaps in the particle list:
   Real blockSize[3];
   const pargrid::CellID* const __restrict__ nbrIDs = simClasses->pargrid.getCellNeighbourIDs(cellID);
   getBlockSize(*simClasses,*sim,cellID,blockSize);
   #ifndef NDEBUG
      if (nbrIDs == NULL) {
	 std::cerr << "(PARTICLE LIST SKELETON) ERROR: Received NULL nbrIDs pointer!" << std::endl;
	 exit(0);
      }
   #endif

   int current = 0;
   int end     = N_particles[cellID]-1;
   while (current <= end) {
      #ifdef USE_INDEX_OPERATOR
         const int i_off = static_cast<int>(floor(particleList[current][XPOS] / blockSize[0]));
         const int j_off = static_cast<int>(floor(particleList[current][YPOS] / blockSize[1]));
         const int k_off = static_cast<int>(floor(particleList[current][ZPOS] / blockSize[2]));
         particleList[current][XPOS] -= i_off * blockSize[0];
         particleList[current][YPOS] -= j_off * blockSize[1];
         particleList[current][ZPOS] -= k_off * blockSize[2];
      #else
         const int i_off = static_cast<int>(floor(particleList[current].state[XPOS] / blockSize[0]));
         const int j_off = static_cast<int>(floor(particleList[current].state[YPOS] / blockSize[1]));
         const int k_off = static_cast<int>(floor(particleList[current].state[ZPOS] / blockSize[2]));
         particleList[current].state[XPOS] -= i_off * blockSize[0];
         particleList[current].state[YPOS] -= j_off * blockSize[1];
         particleList[current].state[ZPOS] -= k_off * blockSize[2];
      #endif

      const pargrid::NeighbourID nbrOffset = simClasses->pargrid.calcNeighbourTypeID(i_off,j_off,k_off);
      const pargrid::CellID newCell = nbrIDs[nbrOffset];
      if (newCell != cellID) {
	 // Particle was propagated to non-existing cell, remove it
	 // and replace the hole with an unpropagated particle from the end of list:
	 if (newCell == pargrid::INVALID_CELLID) {
	    particleList[current] = particleList[end];
	    --end;
	    continue;
	 }

	 // Particle needs to be moved to another process. Copy it to remote cell:
	 wrapper.push_back(newCell,particleList[current]);
	 particleList[current] = particleList[end];
	 --end;
	 continue;
      }
      ++current;
   }

   // Particle list may have a gap in it between particles that still remain on this 
   // cell and that have moved into this cell. The gap is removed here:
   int N_holes     = N_particles[cellID] - current;
   int N_copied    = wrapper.size(cellID)-N_particles[cellID];
   const int N_remaining = N_particles[cellID];
   N_particles[cellID]   = end+1;

   if (N_copied > N_holes) {
      end = wrapper.size(cellID)-1;
      while (N_holes > 0) {
	 particleList[current] = particleList[end];
	 ++current;
	 --end;
	 --N_holes;
      }
      wrapper.resize(cellID,end+1);
   } else {
      end = N_remaining;
      while (N_copied > 0) {
	 particleList[current] = particleList[end];
	 ++current;
	 ++end;
	 --N_copied;
      }
      wrapper.resize(cellID,current);
   }
}

/* Propagate all particles on boundary cells. The actual propagation is done 
 * by template parameter class PROP, but ParticleListSkeleton does all the work 
 * related to moving particles between processes etc. Specifically, this function 
 * waits for incoming particle transfers to complete after propagation, and then 
 * copies received particles from transfer buffer(s) to local cells.
 * @return If true, particles were propagated successfully.*/
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::propagateInnerCellParticles() {
   bool success = true;
   if (initialized == false) return false;
   #if PROFILE_LEVEL > 0
      profile::start(profileName,totalTimeID);
      profile::start("Particle Step",particleTotalID);
   #endif

   if (sim->meshRepartitioned == true) updatePartitioning();

   // Create aliases to particle list and its size stored in cell it->first:
   const double* coordinates = reinterpret_cast<double*>(simClasses->pargrid.getUserData(Simulation::crdsDataID));
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   const std::vector<pargrid::CellID>& innerBlocks = simClasses->pargrid.getInnerCells(pargrid::DEFAULT_STENCIL);

   #ifndef NDEBUG
      if (coordinates == NULL) {
	 std::cerr << "'coordinates' array is NULL!" << std::endl;
	 exit(1);
      }
      if (wrapper.valid() == false) {
	 std::cerr << "wrapper received from ParGrid is invalid!" << std::endl;
	 exit(1);
      }
   #endif

   #if PROFILE_LEVEL > 0
      profile::start("propagation",particlePropagationID);
   #endif

   // Propagate particles on inner cells:
   Real t_propag = 0.0;
   for (size_t c=0; c<innerBlocks.size(); ++c) {
      // Measure cell propagation time if we are preparing for repartitioning:
      if (sim->countPropagTime == true) {
	 t_propag = MPI_Wtime();
      }

      propagateCell(innerBlocks[c],coordinates,wrapper,N_particles);

      if (sim->countPropagTime == true) {
	 t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	 simClasses->pargrid.getCellWeights()[innerBlocks[c]] += t_propag;
      }
   }
   #if PROFILE_LEVEL > 0
      profile::stop(); // propagation
   #endif

   // Set N_particles array to have the correct number of particles per block:
   for (pargrid::CellID blockID=0; blockID<simClasses->pargrid.getNumberOfLocalCells(); ++blockID) {
      N_particles[blockID] = wrapper.size(blockID);
   }

   #if PROFILE_LEVEL > 0
      profile::start("MPI wait recvs",mpiWaitRecvsID);
   #endif

   // Wait for incoming particles:
   size_t bufferCounter = 0;
   const std::set<pargrid::MPI_processID>& neighbourProcesses = simClasses->pargrid.getNeighbourProcesses();
   for (std::set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
      #ifndef NDEBUG
         if (bufferCounter >= inBufferCopy.size()) {
	    std::cerr << "buffercounter too large!" << std::endl;
	    exit(1);
	 }
      #endif

      inBufferCopy[bufferCounter].wait();
      ++bufferCounter;
   }

   #if PROFILE_LEVEL > 0
      profile::stop(); // MPI wait recvs
      profile::start("MPI buffer copy",mpiBufferCopyID);
   #endif

   // Copy particles from inbound buffers to cells' particle lists:
   bufferCounter = 0;
   for (std::set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
      unsigned int offset = 0;
      size_t blockCounter = 0;
      for (size_t b=0; b<inBuffers[bufferCounter].getNumberOfBlocks(); ++b) {
	 const pargrid::CellID blockLID = recvLIDs[bufferCounter][blockCounter];

	 const std::vector<PARTICLE>& buffer = inBuffers[bufferCounter].getBuffer();
	 const uint32_t* const blockSizes = inBuffers[bufferCounter].getBlockSizes();
	 for (uint32_t p=0; p<blockSizes[b]; ++p) {
	    wrapper.push_back(blockLID,buffer[offset+p]);
	 }
	 N_particles[blockLID] += blockSizes[b];
	 offset += blockSizes[b];
	 ++blockCounter;
      }
      ++bufferCounter;
   }

   #if PROFILE_LEVEL > 0
      profile::stop(); // MPI buffer copy
      profile::stop(); // particle step
      profile::stop(); // total time
   #endif
   return success;
}

template<class SPECIES,class PARTICLE> inline
size_t ParticleListSkeleton<SPECIES,PARTICLE>::size() const {
   size_t sum = 0;
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
      sum += wrapper.size()[block];
   }
   return sum;
}

template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::updateAfterRepartition() {
   return updatePartitioning();
}

/**
 * This function is included in "MPI overhead" profiler section.
 * */
template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::updatePartitioning() {
   bool success = true;

   // Prevent multiple calls to updatePartitioning at same time step:
   if (meshVariablesUpdated == sim->timestep) return success;

   // Reallocate MPI transfer buffers:
   const std::set<pargrid::MPI_processID>& neighbourProcesses = simClasses->pargrid.getNeighbourProcesses();

   // Deallocate old in/out buffers and copy protocols:
     {
	std::vector<BUFFER> dummy1;
	std::vector<BUFFER> dummy2;
	inBuffers.swap(dummy1);
	outBuffers.swap(dummy2);

	std::vector<pargrid::CopyProtocolMPI<BUFFER> > dummyCopy1;
	std::vector<pargrid::CopyProtocolMPI<BUFFER> > dummyCopy2;
	inBufferCopy.swap(dummyCopy1);
	inBufferCopy.swap(dummyCopy2);
     }

   // Allocate new buffers:
   inBuffers.resize(neighbourProcesses.size());
   outBuffers.resize(neighbourProcesses.size());
   inBufferCopy.resize(neighbourProcesses.size());
   outBufferCopy.resize(neighbourProcesses.size());

   // Reallocate recvLIDs, sendLIDs vectors:
   if (LIDsCalculatedTimestep != sim->timestep) {
      std::vector<std::vector<pargrid::CellID> > dummy1;
      std::vector<std::vector<pargrid::CellID> > dummy2;
      recvLIDs.swap(dummy1);
      sendLIDs.swap(dummy2);

      recvLIDs.resize(neighbourProcesses.size());
      sendLIDs.resize(neighbourProcesses.size());
   }

   // Reallocate N_particles array:
   size_t sum = 0;
   delete [] N_particles;
   N_particles = new unsigned int[simClasses->pargrid.getNumberOfLocalCells()];
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
      N_particles[block] = wrapper.size()[block];
      sum += wrapper.size()[block];
   }

   size_t counter = 0;
   for (std::set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
      #ifndef NDEBUG
         if (counter >= inBuffers.size()) {
	    std::stringstream ss;
	    ss << "STEP " << sim->timestep << " P#" << sim->mpiRank << " counter=" << counter << " is too large, ";
	    ss << "inBuffers.size()=" << inBuffers.size() << " in updatePartitioning." << std::endl;
	    std::cerr << ss.str();
	    exit(1);
	 }
         if (counter >= outBuffers.size()) {
	    std::stringstream ss;
	    ss << "STEP " << sim->timestep << " P#" << sim->mpiRank << " counter=" << counter << " is too large, ";
	    ss << "outBuffers.size()=" << outBuffers.size() << " in updatePartitioning." << std::endl;
	    std::cerr << ss.str();
	    exit(1);
	 }
      #endif

      const std::size_t N_recvs = simClasses->pargrid.getNumberOfReceives(sim->inverseStencilID,*it);
      inBuffers[counter].resize(N_recvs);

      if (inBufferCopy[counter].set(&(inBuffers[counter]),false,sim->comm,*it) == false) {
	 #ifndef NDEBUG
	    simClasses->logger << "(PLIST SKELETON) ERROR: Failed to set source process rank" << std::endl << write;
	 #endif
	 success = false;
      }

      const std::size_t N_sends = simClasses->pargrid.getNumberOfSends(sim->inverseStencilID,*it);
      outBuffers[counter].resize(N_sends);

      if (outBufferCopy[counter].set(&(outBuffers[counter]),true,sim->comm,*it) == false) {
	 #ifndef NDEBUG
	    simClasses->logger << "(PLIST SKELETON) ERROR: Failed to set destination process rank" << std::endl << write;
	 #endif
	 success = false;
      }

      // Recalculate vector recvLIDs, sendLIDs contents:
      if (LIDsCalculatedTimestep != sim->timestep) {
	 const std::map<pargrid::MPI_processID,std::set<pargrid::CellID> >::const_iterator recvs
	   = simClasses->pargrid.getReceives(sim->inverseStencilID).find(*it);
	 recvLIDs[counter].resize(recvs->second.size());
	 size_t c=0;
	 for (std::set<pargrid::CellID>::const_iterator cellGID=recvs->second.begin(); cellGID!=recvs->second.end(); ++cellGID) {
	    const pargrid::CellID cellLID = simClasses->pargrid.getLocalID(*cellGID);
	    recvLIDs[counter][c] = cellLID;
	    ++c;
	 }

	 const std::map<pargrid::MPI_processID,std::set<pargrid::CellID> >::const_iterator sends 
	   = simClasses->pargrid.getSends(sim->inverseStencilID).find(*it);
	 sendLIDs[counter].resize(sends->second.size());
	 c=0;
	 for (std::set<pargrid::CellID>::const_iterator cellGID=sends->second.begin(); cellGID!=sends->second.end(); ++cellGID) {
	    const pargrid::CellID cellLID = simClasses->pargrid.getLocalID(*cellGID);
	    sendLIDs[counter][c] = cellLID;
	    ++c;
	 }
      }

      ++counter;
   }
   LIDsCalculatedTimestep = sim->timestep;
   meshVariablesUpdated = sim->timestep;
   return success;
}

template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::waitParticleSends() {
   bool success = true;
   if (getInitialized() == false) return false;

   // Wait for all particle sends to complete:
   #ifdef PROFILE
      profile::start(profileName,totalTimeID);
      profile::start("Particle Step",particleTotalID);
      profile::start("MPI wait sends",mpiWaitSendsID);
   #endif

   // Wait for sends:
   for (size_t i=0; i<outBufferCopy.size(); ++i) outBufferCopy[i].wait();

   #ifdef PROFILE
      profile::stop(); // wait sends
      profile::stop(); // particle step
      profile::stop(); // total time
   #endif
   return success;
}

template<class SPECIES,class PARTICLE> inline
bool ParticleListSkeleton<SPECIES,PARTICLE>::writeParticles(const std::string& spatMeshName) {
   bool success = true;
   if (initialized == false) return false;

   #if PROFILE_LEVEL > 0
      profile::start(speciesName+" writing",particleWriteID);
   #endif

   const size_t sum_particles = size();
   Real* buffer = new Real[sum_particles*3];
   const double* coordinates = getBlockCoordinateArray(*sim,*simClasses);
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);   

   size_t counter = 0;
   PARTICLE** particleLists = wrapper.data();
   for (size_t block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
      for (unsigned int p=0; p<wrapper.size()[block]; ++p) {
	 #ifdef USE_INDEX_OPERATOR
	    buffer[counter*3+XPOS] = particleLists[block][p][XPOS] + coordinates[3*block+XPOS];
	    buffer[counter*3+YPOS] = particleLists[block][p][YPOS] + coordinates[3*block+YPOS];
	    buffer[counter*3+ZPOS] = particleLists[block][p][ZPOS] + coordinates[3*block+ZPOS];
	 #else
	    buffer[counter*3+XPOS] = particleLists[block][p].state[XPOS] + coordinates[3*block+XPOS];
	    buffer[counter*3+YPOS] = particleLists[block][p].state[YPOS] + coordinates[3*block+YPOS];
	    buffer[counter*3+ZPOS] = particleLists[block][p].state[ZPOS] + coordinates[3*block+ZPOS];
	 #endif
	 ++counter;
      }
   }

   std::map<std::string,std::string> attribs;
   attribs["name"] = speciesName;
   attribs["type"] = vlsv::mesh::STRING_POINT;
   if (simClasses->vlsv.writeArray("MESH",attribs,sum_particles,3,buffer) == false) {
      simClasses->logger << "\t ERROR failed to write particle species!" << std::endl;
      success = false;
   }
   delete [] buffer; buffer = NULL;

   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}

#endif
