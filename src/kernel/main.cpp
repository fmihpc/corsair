/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2015 Finnish Meteorological Institute
 *  Copyright 2016 Arto Sandroos
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
#include <sstream>
#include <cmath>

#ifdef _OPENMP
   #include <omp.h>
#endif

#include <main.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <datawriter.h>
#include <dataoperatorcontainer.h>
#include <dataoperator.h>
#include <user.h>
#include <restart_writer.h>
#include <restart_builder.h>
#include <object_factories.h>

using namespace std;

const int MASTER_RANK = 0;
static corsair::ObjectWrapper objectWrapper;
const string logFileName = "logfile.txt";

bool continueSimulation(Simulation& sim) {
   bool cont = true;
   if (sim.maximumTime >= 0.0) if (sim.t >= sim.maximumTime) cont = false;
   if (sim.maximumTimesteps >= 0) if (sim.timestep >= sim.maximumTimesteps) cont = false;
   return cont;
}

corsair::ObjectWrapper& corsair::getObjectWrapper() {
   return objectWrapper;
}

void corsair::evaluateComputationalLoad(corsair::ObjectWrapper& owrapper,corsair::ComputationalLoad& global) {
   
   Real localLoad[compuload::SIZE];
   for (int i=0; i<compuload::SIZE; ++i) localLoad[i] = 0.0;

   // Calculate computational load of this process:
   for (pargrid::CellID block=0; block<owrapper.simClasses.pargrid.getNumberOfLocalCells(); ++block) {
      localLoad[compuload::SECONDS_PER_CELL] += owrapper.simClasses.pargrid.getCellWeights()[block];
   }

   // Sum up number of macroparticles on this process:
   uint64_t N_particles = 0;
   for (size_t l=0; l<owrapper.particleLists.size(); ++l) {
      N_particles += owrapper.particleLists[l]->size();
   }
   localLoad[compuload::N_PARTICLES] = N_particles;

   // Reduce avg,min,max values to master process:
   MPI_Reduce(localLoad,global.loadMin,compuload::SIZE,MPI_Type<Real>(),MPI_MIN,owrapper.sim.MASTER_RANK,owrapper.sim.comm);
   MPI_Reduce(localLoad,global.loadMax,compuload::SIZE,MPI_Type<Real>(),MPI_MAX,owrapper.sim.MASTER_RANK,owrapper.sim.comm);
   MPI_Reduce(localLoad,global.loadAvg,compuload::SIZE,MPI_Type<Real>(),MPI_SUM,owrapper.sim.MASTER_RANK,owrapper.sim.comm);

   if (owrapper.sim.mpiRank == owrapper.sim.MASTER_RANK) {
      // Average load on processes:
      for (int i=0; i<compuload::SIZE; ++i) global.loadAvg[i] /= owrapper.sim.mpiProcesses;
      
      // Calculate imbalance ratio(s):
      for (int i=0; i<compuload::SIZE; ++i) {
	 Real minRatio = 0.0, maxRatio = 0.0;
	 if (fabs(global.loadMin[i]) > 1e-6) minRatio = global.loadAvg[i] / global.loadMin[i];
	 if (fabs(global.loadMax[i]) > 1e-6) maxRatio = global.loadMax[i] / global.loadAvg[i];
	 global.ratios[i] = max(minRatio,maxRatio);
      }
   }

   // Broadcast imbalance ratios to all processes:
   MPI_Bcast(global.ratios,compuload::SIZE,MPI_Type<Real>(),owrapper.sim.MASTER_RANK,owrapper.sim.comm);

   // Write logfile message:
   owrapper.simClasses.logger
     << "\t Current average load is:  " << global.loadAvg[compuload::SECONDS_PER_CELL] << " seconds." << endl
     << "\t Current minimum load is:  " << global.loadMin[compuload::SECONDS_PER_CELL] << " seconds." << endl
     << "\t Current maximum load is:  " << global.loadMax[compuload::SECONDS_PER_CELL] << " seconds." << endl
     << "\t Current load imbalance is: " << global.ratios[compuload::SECONDS_PER_CELL] << endl
     << "\t Threshold for repartitioning is: " << owrapper.sim.maximumLoadImbalance << endl;
}

bool checkFileExists(const string fileName) {
   bool rvalue = false;
   ifstream ifs(fileName.c_str(),ifstream::in);
   if (ifs.good() == true) { rvalue = true; }
   return rvalue;
}

int main(int argn,char* args[]) {
   bool initialized = true;
   MPI_Init(&argn,&args);
   if (objectWrapper.sim.initialize(argn,args) == false) initialized = false;
   MPI_Comm_dup(MPI_COMM_WORLD,&objectWrapper.sim.comm);
   MPI_Comm_size(objectWrapper.sim.comm,&objectWrapper.sim.mpiProcesses);
   MPI_Comm_rank(objectWrapper.sim.comm,&objectWrapper.sim.mpiRank);

   int profileInitID = -1;
   int profilePropagID = -1;
   int profileSaveID = -1;
   int profileMainLoopID = -1;
   int profileRepartitionID = -1;
   int profileTestsID = -1;

   // Do not start if a log file exists already:
   int logFileExists = 0;
   if (objectWrapper.sim.mpiRank == objectWrapper.sim.MASTER_RANK) {
      if (checkFileExists(logFileName) == true) {
	 logFileExists = 1;
	 cerr << "(MAIN) Log file (" << logFileName << ") exist already, exiting.." << endl << flush;
      }
   }
   MPI_Bcast(&logFileExists,1,MPI_Type<int>(),objectWrapper.sim.MASTER_RANK,objectWrapper.sim.comm);
   if (logFileExists == 1) {
      MPI_Finalize();
      return 0;
   }

   #if PROFILE_LEVEL > 0
      profile::start("corsair/initialization",profileInitID);
   #endif

   // Open a log file:
   objectWrapper.simClasses.logger.open(objectWrapper.sim.comm,objectWrapper.sim.MASTER_RANK,logFileName,true);
   objectWrapper.simClasses.logger << "(MAIN) Corsair starting initialisation." << endl;
   objectWrapper.simClasses.logger << "Using " << objectWrapper.sim.mpiProcesses << " MPI processes and ";
   #ifndef _OPENMP
      objectWrapper.simClasses.logger << 1;
   #else
      objectWrapper.simClasses.logger << omp_get_max_threads();
   #endif
   objectWrapper.simClasses.logger << " OpenMP threads per process" << endl << write;

   // Init configuration file reader and parse all parameteres needed in main:
   if (objectWrapper.configReader.initialize("CORSAIR_",argn,args,objectWrapper.sim.comm,objectWrapper.sim.MASTER_RANK) == false) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: ConfigReader initialization failed!" << endl << write;
      initialized = false;
   }

   const Real NaN = NAN;
   string doRestart,runTestsString;
   objectWrapper.configReader.add("Simulation.time_initial","Value of time at the start of simulation (float).",(Real)0.0);
   objectWrapper.configReader.add("Simulation.maximum_time","Simulation exits after current simulation time exceeds this value. A negative value indicates that maximum time is not used as an end condition (float).",(Real)-1.0);
   objectWrapper.configReader.add("Simulation.maximum_timesteps","Simulation exits after this many timesteps have been computed. A negative value indicates that maximum timesteps is not used as an end condition (int).",(int)-1);
   objectWrapper.configReader.add("Simulation.dt","Time step (float).",(Real)1.0);
   objectWrapper.configReader.add("Simulation.data_save_interval_unit","Unit in which 'data_save_interval' is given (time/timestep) (string).",string(""));
   objectWrapper.configReader.add("Simulation.data_save_interval","Interval between data saves (int or float).",NaN);
   objectWrapper.configReader.add("Simulation.data_save_initial_state","Save initial state (yes or no).",string("yes"));
   objectWrapper.configReader.add("Simulation.data_save_during_first_timesteps","Save state during first N timesteps (int).",0);
   objectWrapper.configReader.add("Simulation.restart_write_interval","Interval, in time steps, between restart file writes. Defaults to zero value, i.e. restart files are not written (int).",0);
   objectWrapper.configReader.add("Simulation.random_number_generator.seed","Seed value for random number generator, with value zero it is calculated from system clock (int).",(int)0);
   objectWrapper.configReader.addComposed("LoadBalance.methods","Load balancing method for each hierarchical level (NONE / BLOCK / RANDOM / RCB / RIB / GRAPH / HYPERGRAPH).");
   objectWrapper.configReader.addComposed("LoadBalance.tolerances","Imbalance tolerance for each hierarchical level (float).");
   objectWrapper.configReader.addComposed("LoadBalance.processes_per_partition","Processes per hierarchical partition (int).");
   objectWrapper.configReader.add("gridbuilder","Name of builder that is used to create the initial simulation grid (string).","");
   objectWrapper.configReader.add("Simulation.repartition_check_interval","Time interval, in number of time steps, for checking if mesh needs to be repartitioned (int).",-1);
   objectWrapper.configReader.add("Simulation.maximum_load_imbalance","Maximum computational load imbalance, defaults to 1.3 (real)",(Real)1.30);
   objectWrapper.configReader.add("Simulation.restart","Is simulation restarted from a file (yes/no), defaults to 'no' (string)",string("no"));
   objectWrapper.configReader.add("Simulation.restart_filename_prefix","Restart filename prefix, defaults to 'restart' (string).",string("restart"));
   objectWrapper.configReader.add("Simulation.restart_major_store_interval","Restart files written at these intervals are never deleted, defaults to zero (int)",0);
   objectWrapper.configReader.add("Simulation.restart_minor_store_amount","Keep only N newest restart files where N is the value of this parameter, unless the restart file timestep corresponds to major store interval. Defaults to zero value, i.e. all files are kept (int).",0);
   objectWrapper.configReader.add("Simulation.mesh_always_written","If 'yes' mesh is written to every VLSV file, default value is 'no' (string).",string("no"));
   objectWrapper.configReader.add("Simulation.run_tests","If 'yes' corsair will run user-specified test suite, default value is 'no' (string).",string("no"));
   objectWrapper.configReader.parse();
   objectWrapper.sim.timestep = 0;
   objectWrapper.configReader.get("Simulation.time_initial",objectWrapper.sim.t);
   objectWrapper.configReader.get("Simulation.maximum_time",objectWrapper.sim.maximumTime);
   objectWrapper.configReader.get("Simulation.maximum_timesteps",objectWrapper.sim.maximumTimesteps);
   objectWrapper.configReader.get("Simulation.repartition_check_interval",objectWrapper.sim.repartitionCheckInterval);
   objectWrapper.configReader.get("Simulation.maximum_load_imbalance",objectWrapper.sim.maximumLoadImbalance);
   objectWrapper.configReader.get("Simulation.dt",objectWrapper.sim.dt);
   objectWrapper.configReader.get("Simulation.restart_write_interval",objectWrapper.sim.restartWriteInterval);
   objectWrapper.configReader.get("Simulation.restart",doRestart);
   objectWrapper.configReader.get("Simulation.restart_filename_prefix",objectWrapper.sim.restartFilenamePrefix);
   objectWrapper.configReader.get("Simulation.restart_major_store_interval",objectWrapper.sim.restartMajorInterval);
   objectWrapper.configReader.get("Simulation.restart_minor_store_amount",objectWrapper.sim.restartMinorFileAmount);
   objectWrapper.configReader.get("Simulation.run_tests",runTestsString);
   objectWrapper.sim.restarted = false;
   if (doRestart == "yes") objectWrapper.sim.restarted = true;
   if (runTestsString == "yes") objectWrapper.sim.runTests = true;
   
   if (objectWrapper.sim.dt <= 0.0) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Timestep given in parameter 'Simulation.dt' must be > 0." << endl << write;
      initialized = false;
   }
   if (objectWrapper.sim.maximumTime < 0.0 && objectWrapper.sim.maximumTimesteps < 0) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Maximum number of timesteps, 'Simulation.maximum_timesteps', and maximum simulation time, 'Simulation.maximum_time', cannot both have negative values." << endl << write;
      initialized = false;
   }
   
   // Call user-defined "early" initialization function:
   if (initialized == true) if (userEarlyInitialization(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.configReader,objectWrapper.particleLists) == false) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: User-defined early initialization function failed!" << endl << write;
      initialized = false;
   }

   // Init random number generator
   int randomSeed;
   objectWrapper.configReader.get("Simulation.random_number_generator.seed",randomSeed);
   if (randomSeed == 0) randomSeed = time(NULL);
   if (initialized == true) {
      // Master process calculates a random number seed for each process:
      int* seeds = NULL;
      if (objectWrapper.sim.mpiRank == objectWrapper.sim.MASTER_RANK) {
         objectWrapper.simClasses.logger << "(MAIN) Master process generating random number generator seeds ";
         objectWrapper.simClasses.logger << " using initial seed " << randomSeed << endl;
         seeds = new int[objectWrapper.sim.mpiProcesses];
         srand(randomSeed);
         for (int p=0; p<objectWrapper.sim.mpiProcesses; ++p) {
            seeds[p] = rand();
            objectWrapper.simClasses.logger << "\t P#" << p << "\t seed: " << seeds[p] << endl;
         }
      }
      MPI_Scatter(seeds,1,MPI_Type<int>(),&randomSeed,1,MPI_Type<int>(),objectWrapper.sim.MASTER_RANK,objectWrapper.sim.comm);
      
      objectWrapper.simClasses.logger << "(MAIN) Received random seed " << randomSeed << " from master process" << endl;
      if (objectWrapper.simClasses.random.initialize(randomSeed) == false) {
         objectWrapper.simClasses.logger << "(MAIN) ERROR: Failed to initialize random number generator!" << endl << write;
         initialized = false;
      }      
      delete [] seeds; seeds = NULL;
   }
   
   // Init parallel grid
   vector<string> loadBalancingMethods;
   vector<string> imbalanceTolerances;
   vector<string> procsPerPartition;
   vector<map<pargrid::InputParameter,string> > lb_parameters;
   objectWrapper.configReader.get("LoadBalance.methods",loadBalancingMethods);
   objectWrapper.configReader.get("LoadBalance.tolerances",imbalanceTolerances);
   objectWrapper.configReader.get("LoadBalance.processes_per_partition",procsPerPartition);

   if (loadBalancingMethods.size() != imbalanceTolerances.size()) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Different number of load balancing methods and tolerances in config file!" << endl << write;
      initialized = false;
   }

   if (procsPerPartition.size() != loadBalancingMethods.size()) {
      if (loadBalancingMethods.size() > 1) {
         objectWrapper.simClasses.logger << "(MAIN) ERROR: Parameter 'processes_per_partition' has incorrect size, hierarchical partitioning requested!" << endl << write;
         initialized = false;
      } else {
         procsPerPartition.push_back("1");
      }
   }
   
   if (initialized == true) for (size_t i=0; i<loadBalancingMethods.size(); ++i) {
      map<pargrid::InputParameter,string> parameters;
      parameters[pargrid::cellWeightScale]       = "1.0";
//      parameters[pargrid::edgeWeightScale]       = "0.25";
      parameters[pargrid::imbalanceTolerance]    = imbalanceTolerances[i];
      parameters[pargrid::loadBalancingMethod]   = loadBalancingMethods[i];
      parameters[pargrid::processesPerPartition] = procsPerPartition[i];
      lb_parameters.push_back(parameters);
   }
   
   if (initialized == true) if (objectWrapper.simClasses.pargrid.initialize(objectWrapper.sim.comm,lb_parameters) == false) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Failed to initialize ParGrid!" << endl << write;
      initialized = false;
   }

   // Add all known GridBuilders to GridBuilderFactory:
   if (objectWrapper.objectFactories.gridBuilders.registerMaker("Restart",RestartCreator) == false) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Failed to register RestartBuilder!" << endl << write;
      initialized = false;
   }
   if (initialized == true) if (registerObjectMakers(objectWrapper.objectFactories) == false) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Failed to register GridBuilder(s)!" << endl << write;
      initialized = false;
   }

   // Create the requested GridBuilder, and use it to create initial simulation grid:
   #if PROFILE_LEVEL > 0
      static int profGridBuilding = -1;
      profile::start("Grid Builder",profGridBuilding);
   #endif
   string gridBuilderName;
   objectWrapper.configReader.get("gridbuilder",gridBuilderName);
   GridBuilder* builder = NULL;
   objectWrapper.sim.restartTimestep = numeric_limits<unsigned int>::max();
   if (initialized == true) {
      // If simulation is restarted the mesh is read from a restart file
      // instead of being generated by the regular GridBuilder:
      if (objectWrapper.sim.restarted == true) {
         GridBuilder* builder = objectWrapper.objectFactories.gridBuilders.create("Restart");
         if (builder == NULL) {
            objectWrapper.simClasses.logger << "(MAIN) ERROR: Failed to create RestartBuilder!" << endl << write;
            initialized = false;
         }
         if (initialized == true) {
            if (builder->initialize(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.configReader) == false) {
               objectWrapper.simClasses.logger << "(MAIN) ERROR: RestartBuilder initialization failed!" << endl << write;
               initialized = false;
            }
         }
         if (initialized == true) {
            if (builder->build(objectWrapper.sim,objectWrapper.simClasses) == false) {
               objectWrapper.simClasses.logger << "(MAIN) ERROR: Failed to rebuild mesh!" << endl << write;
               initialized = false;
            }
         }
         delete builder; builder = NULL;
         objectWrapper.sim.restartTimestep = objectWrapper.sim.timestep;
      }
      
      // Regular GridBuilder is always created and initialized:
      builder = objectWrapper.objectFactories.gridBuilders.create(gridBuilderName);
      if (builder == NULL) {
         objectWrapper.simClasses.logger << "(MAIN) ERROR: Failed to create GridBuilder with name '" << gridBuilderName << "'" << endl << write;
         initialized = false;
      }
      if (initialized == true) if (builder->initialize(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.configReader) == false) {
         objectWrapper.simClasses.logger << "(MAIN) ERROR: GridBuilder initialization failed!" << endl << write;
         initialized = false;
      }
      // If simulation was not restarted the regular GridBuilder is used
      // to create the simulation mesh:
      if (initialized == true && objectWrapper.sim.restarted == false) if (builder->build(objectWrapper.sim,objectWrapper.simClasses) == false) {
         objectWrapper.simClasses.logger << "(MAIN) ERROR: GridBuilder failed!" << endl << write;
         initialized = false;
      }
   }
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif

   // Create inverse of default stencil:
   objectWrapper.sim.defaultStencilID = pargrid::DEFAULT_STENCIL;
   vector<pargrid::NeighbourID> nbrTypeIDs;
   for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
      if (i==0 && (j== 0 && k == 0)) continue;
      nbrTypeIDs.push_back(objectWrapper.simClasses.pargrid.calcNeighbourTypeID(i,j,k));
   }
   objectWrapper.sim.inverseStencilID = objectWrapper.simClasses.pargrid.addStencil(pargrid::remoteToLocalUpdates,nbrTypeIDs);
   nbrTypeIDs.clear();
   if (objectWrapper.sim.inverseStencilID == pargrid::INVALID_STENCILID) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Failed to create inverse Stencil!" << endl << write;
      initialized = false;
   }
   
   // Initialize data operator container:
   if (initialized == true) if (objectWrapper.dataOperatorContainer.initialize(objectWrapper.configReader,objectWrapper.sim,objectWrapper.simClasses) == false) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: DataOperatorContainer initialization failed!" << endl << write;
      initialized = false;
   }
   // Call user-defined function which inserts data operators into the container:
   if (initialized == true) if (registerDataOperators() == false) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Registration of one or more DataOperators failed!" << endl << write;
      initialized = false;
   }

   // Read data save interval from config file:
   string dataSaveIntervalUnit;
   objectWrapper.sim.t_previousDataSave = objectWrapper.sim.t;   
   objectWrapper.configReader.get("Simulation.data_save_interval_unit",dataSaveIntervalUnit);
   if (dataSaveIntervalUnit == "time") objectWrapper.sim.dataIntervalIsTime = true;
   else if (dataSaveIntervalUnit == "timestep") objectWrapper.sim.dataIntervalIsTime = false;
   else {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Data save interval unit was not given with parameter 'Simulation.data_save_interval_unit'" << endl;
      objectWrapper.simClasses.logger << "\t or it has an unknown value. Accepted values are 'time' and 'timestep'." << endl << write;
      initialized = false;
   }
   
   if (objectWrapper.sim.dataIntervalIsTime == true) {
      objectWrapper.configReader.get("Simulation.data_save_interval",objectWrapper.sim.dataIntervalFloat);
      if (objectWrapper.sim.dataIntervalFloat != objectWrapper.sim.dataIntervalFloat) {
         objectWrapper.simClasses.logger << "(MAIN) ERROR: Data save interval was not given with parameter 'Simulation.data_save_interval' !" << endl;
         initialized = false;
      }
   } else {
      objectWrapper.configReader.get("Simulation.data_save_interval",objectWrapper.sim.dataIntervalInteger);
      if (objectWrapper.sim.dataIntervalInteger > 100000) {
         objectWrapper.simClasses.logger << "(MAIN) ERROR: Data save interval was not given with parameter 'Simulation.data_save_interval'" << endl;
         objectWrapper.simClasses.logger << "\t or it has a very large value!" << endl << write;
         initialized = false;
      }
   }

   // Read if output files should be saved after initialization
   string dataSaveInitial = "yes";
   objectWrapper.configReader.get("Simulation.data_save_initial_state",dataSaveInitial);

   // Read if output files should be saved during first N timesteps
   int dataSaveFirstTimesteps = 0;
   objectWrapper.configReader.get("Simulation.data_save_during_first_timesteps",dataSaveFirstTimesteps);

   // Read if mesh should be written to every VLSV file:
   string meshAlwaysWritten;
   objectWrapper.sim.meshAlwaysWritten = false;
   objectWrapper.configReader.get("Simulation.mesh_always_written",meshAlwaysWritten);
   if (meshAlwaysWritten == "yes") objectWrapper.sim.meshAlwaysWritten = true;
   
   // Call user-defined "late" initialization function:
   if (initialized == true) {
      bool status = userLateInitialization(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.configReader,
                                           objectWrapper.objectFactories,objectWrapper.particleLists);
      if (status == false) {
         objectWrapper.simClasses.logger << "(MAIN) ERROR: User-defined late initialization function failed " << endl;
         objectWrapper.simClasses.logger << "              on one or more MPI process" << endl << write;
         initialized = false;
      }
   }
   
   // If simulation was restarted, ask particle lists update their internal state(s):
   if (initialized == true && objectWrapper.sim.restarted == true) {
      objectWrapper.simClasses.logger << "(MAIN) Updating particle lists after restart" << endl;
      for (size_t plist=0; plist<objectWrapper.particleLists.size(); ++plist) {
         if (objectWrapper.particleLists[plist]->updateAfterRepartition() == false) {
            objectWrapper.simClasses.logger << "(MAIN) ERROR: Particle list failed to update after restart" << endl << write;
            initialized = false;
         }
      }
   }

   // Check that all processes have succeeded in initialization:
   if (objectWrapper.simClasses.pargrid.checkSuccess(initialized) == false) {
      objectWrapper.simClasses.logger << "(MAIN) ERROR: Simulation did not initialize correctly" << endl << write;
      initialized = false;
   }
   
   // Distribute initialization status to all processes. If any process 
   // failed to initialize, all processes will know about it:
   int initSum = 0, globalInitSum;
   if (initialized == false) ++initSum;
   MPI_Allreduce(&initSum,&globalInitSum,1,MPI_Type<int>(),MPI_SUM,objectWrapper.sim.comm);
   if (globalInitSum > 0) initialized = false;
   
   // Check that initialisation succeeded. If not, exit program (well we just skip 
   // the propagation loop for clean exit):
   if (initialized == true) {
      objectWrapper.simClasses.logger << "(MAIN) Initialization completed." << endl << write;
   } else {
      objectWrapper.simClasses.logger << "(MAIN) Initialization failed, exiting!" << endl << write;
      if (objectWrapper.sim.mpiRank == objectWrapper.sim.MASTER_RANK) {
         cerr
	   << endl << "One or more ERRORS occurred during initialization, aborting." << endl << endl
	   << "Here are configuration file parameters read by this compilation of CORSAIR." << endl
	   << "Check logfile for specific error(s) that were encountered durint init." << endl << endl;
         objectWrapper.configReader.printHelpMessage();
      }
   }
   // Stop profiling of corsair initialization:
   #if PROFILE_LEVEL > 0
      profile::stop();
      profile::start("corsair/tests",profileTestsID);
   #endif

   // Run test suite:
   if (initialized == true && objectWrapper.sim.runTests == true) {
      objectWrapper.simClasses.logger << "(MAIN) Starting to run test suite" << endl << write;
      if (userRunTests(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.particleLists) == false) {
         objectWrapper.simClasses.logger << "(MAIN) WARNING: User test suite failed" << endl << write;
      } else {
         objectWrapper.simClasses.logger << "(MAIN) Test suite run successfully" << endl << write;
      }
   }      
   #if PROFILE_LEVEL > 0
      profile::stop();
      profile::start("corsair/main loop",profileMainLoopID);
   #endif

   if (initialized == true) {
      objectWrapper.sim.initializing = false;
      objectWrapper.sim.countPropagTime = false;
      objectWrapper.simClasses.pargrid.setPartitioningMode(pargrid::repartition);

      if(dataSaveInitial == "yes") {
      #if PROFILE_LEVEL > 0
	  profile::start("corsair/save",profileSaveID);
      #endif
	  objectWrapper.simClasses.logger << "(MAIN) Saving initial state" << endl;
	  saveState(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.dataOperatorContainer,objectWrapper.particleLists,builder);
      
      #if PROFILE_LEVEL > 0
	  profile::stop();
      #endif
      }
      
      fstream loads("timeseries_loads.txt", fstream::out);
      if (objectWrapper.sim.mpiRank == objectWrapper.sim.MASTER_RANK) 
        loads << "#step load(min) load(max) load(avg) imbalance_ratio" << endl;
      
      objectWrapper.simClasses.logger << "(MAIN) Starting to run simulation" << endl << write;
      while (continueSimulation(objectWrapper.sim) == true) {
         // If partitioning is checked on this time step, set objectWrapper.sim.countPropagTime to true:
         if (objectWrapper.sim.repartitionCheckInterval > 0) {
            if ((objectWrapper.sim.timestep+1) % objectWrapper.sim.repartitionCheckInterval == 0) {
               objectWrapper.sim.countPropagTime = true;
               objectWrapper.simClasses.pargrid.clearCellWeights();
            }
         }

         // If data is saved on this time step, set objectWrapper.sim.atDataSaveStep to true:
         objectWrapper.sim.atDataSaveStep = false;
         if (objectWrapper.sim.dataIntervalIsTime == true) {
            if ((objectWrapper.sim.t + objectWrapper.sim.dt) > objectWrapper.sim.t_previousDataSave + objectWrapper.sim.dataIntervalFloat) {
               objectWrapper.sim.atDataSaveStep = true;
               objectWrapper.sim.t_previousDataSave = objectWrapper.sim.t + objectWrapper.sim.dt;
            }
         } else {
            if ((objectWrapper.sim.timestep+1) % objectWrapper.sim.dataIntervalInteger == 0) {
               objectWrapper.sim.atDataSaveStep = true;
            }
         }
	 if (objectWrapper.sim.timestep < dataSaveFirstTimesteps) {
	    objectWrapper.sim.atDataSaveStep = true;
	 }
	 if(objectWrapper.sim.atDataSaveStep == true) {
	    objectWrapper.simClasses.logger << "(MAIN) At data save step (time step = " << objectWrapper.sim.timestep << ", time = " << objectWrapper.sim.t << ")" << endl << write;
	 }
         
         // Check if restart file should be written:
         if (objectWrapper.sim.restartWriteInterval > 0 && objectWrapper.sim.timestep != objectWrapper.sim.restartTimestep) {
            if (objectWrapper.sim.timestep % objectWrapper.sim.restartWriteInterval == 0 || objectWrapper.sim.timestep == 0) {
               objectWrapper.simClasses.logger << "(MAIN) Writing restart file (time step = " << objectWrapper.sim.timestep << ", time = " << objectWrapper.sim.t << ")" << endl;
               if (writeRestart(objectWrapper.sim,objectWrapper.simClasses) == false) {
                  objectWrapper.simClasses.logger << "(MAIN) Failed to write restart file!" << endl << write;
               }
            }
         }

         // Call user-defined propagation function:
         #if PROFILE_LEVEL > 0
            profile::start("corsair/propagation",profilePropagID);
         #endif

         if (propagate(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.particleLists) == false) {
            objectWrapper.simClasses.logger << "(MAIN) ERROR: User-defined propagate-function returned value 'false', exiting." << endl << write;
            initialized = false;
            break;
         }

         #if PROFILE_LEVEL > 0
            profile::stop();
         #endif
	 
         ++objectWrapper.sim.timestep;
         objectWrapper.sim.t += objectWrapper.sim.dt;

	 if(continueSimulation(objectWrapper.sim) == false) {
	    objectWrapper.simClasses.logger << "(MAIN) At final time step (time step = " << objectWrapper.sim.timestep << ", time = " << objectWrapper.sim.t << ")" << endl << write;
	 }
	 
         if (objectWrapper.sim.atDataSaveStep == true || continueSimulation(objectWrapper.sim) == false) {
            #if PROFILE_LEVEL > 0
               profile::start("corsair/save",profileSaveID);
            #endif

            saveState(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.dataOperatorContainer,objectWrapper.particleLists,builder);

            #if PROFILE_LEVEL > 0
               profile::stop();
            #endif
         }
	 
         // Check if new partitioning should be calculated:
         objectWrapper.sim.meshRepartitioned = false;
         if (objectWrapper.sim.repartitionCheckInterval > 0 && objectWrapper.sim.mpiProcesses > 1) {
            if (objectWrapper.sim.timestep % objectWrapper.sim.repartitionCheckInterval == 0) {
	           #if PROFILE_LEVEL > 0
                  profile::start("corsair/repartition check",profileRepartitionID);
               #endif
	       
               objectWrapper.simClasses.logger << endl << "(MAIN) Testing for repartitioning (time step = " << objectWrapper.sim.timestep << ", time = " << objectWrapper.sim.t << ")" << endl;

               // Evaluate computational load:
               corsair::ComputationalLoad compuLoad;
               corsair::evaluateComputationalLoad(objectWrapper,compuLoad);
               
               // Write statistics of computational load to time series file (master process only):
               if (objectWrapper.sim.mpiRank == objectWrapper.sim.MASTER_RANK) {
                  loads << objectWrapper.sim.timestep << '\t'; 
                  
                  // Write min,max,avg,imbalance ratio:
                  loads << compuLoad.loadMin[corsair::compuload::SECONDS_PER_CELL] << '\t';
                  loads << compuLoad.loadMax[corsair::compuload::SECONDS_PER_CELL] << '\t';
                  loads << compuLoad.loadAvg[corsair::compuload::SECONDS_PER_CELL] << '\t';
                  loads << compuLoad.ratios[corsair::compuload::SECONDS_PER_CELL] << '\t';
                  
                  // Write min,max,avg number of particles and imbalance ratio:
                  loads << compuLoad.loadMin[corsair::compuload::N_PARTICLES] << '\t';
                  loads << compuLoad.loadMax[corsair::compuload::N_PARTICLES] << '\t';
                  loads << compuLoad.loadAvg[corsair::compuload::N_PARTICLES] << '\t';
                  loads << compuLoad.ratios[corsair::compuload::N_PARTICLES] << '\t';
                  loads << endl;
               }

               #if PROFILE_LEVEL > 0
                  profile::stop();
               #endif
	       
               // If current imbalance ratio exceeds maximum tolerance, repartition mesh:
               if (compuLoad.ratios[corsair::compuload::SECONDS_PER_CELL] >= objectWrapper.sim.maximumLoadImbalance) {
                  objectWrapper.simClasses.logger << "\t Repartitioning mesh" << endl;
                  if (objectWrapper.simClasses.pargrid.balanceLoad() == false) {
                     objectWrapper.simClasses.logger << "(MAIN) Load Balancing failed!" << endl << write;
                  }

                  // Set variable values to reflect the state of newly repartitioned simulation mesh:
                  objectWrapper.sim.meshChangedStep   = objectWrapper.sim.timestep;
                  objectWrapper.sim.meshRepartitioned = true;
                  
                  // Ask particle lists to update their internal state:
                  for (size_t plist=0; plist<objectWrapper.particleLists.size(); ++plist) {
                     if (objectWrapper.particleLists[plist]->updateAfterRepartition() == false) {
                        objectWrapper.simClasses.logger << "(MAIN) Particle list failed to update after repartition" << endl << write;
                        initialized = false;
                     }
                  }
               }
	       else {
		  objectWrapper.simClasses.logger << "\t No repartitioning needed" << endl;
	       }
               objectWrapper.sim.countPropagTime = false;
	       objectWrapper.simClasses.logger << endl;
            }
         }
      }
      loads.close();
   }
   #if PROFILE_LEVEL > 0
      profile::stop();
      if (initialized == true) profile::print(objectWrapper.sim.comm,objectWrapper.sim.MASTER_RANK);
   #endif
   
   // Deallocate memory for clean exit:
   if (initialized == true) objectWrapper.simClasses.logger << write << "(MAIN) Exiting simulation after successful run." << endl << write;
   else objectWrapper.simClasses.logger << "(MAIN) Exiting simulation due to error(s)." << endl << write;
   
   if (userFinalization(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.particleLists) == false) {
      objectWrapper.simClasses.logger << "(MAIN) WARNING: User finalization function returned 'false'" << endl << write;
   }
   
   for (size_t i=0; i<objectWrapper.particleLists.size(); ++i) {
      objectWrapper.particleLists[i]->finalize();
      delete objectWrapper.particleLists[i]; objectWrapper.particleLists[i] = NULL;
   }
   if (objectWrapper.dataOperatorContainer.finalize() == false) {
      objectWrapper.simClasses.logger << "(MAIN) WARNING: Data operator container finalization function returned 'false'" << endl;
      objectWrapper.simClasses.logger << "\t This most likely means that one or more user-defined data operators" << endl;
      objectWrapper.simClasses.logger << "\t failed to finalize properly." << endl << write;
   }
   objectWrapper.configReader.finalize();
   objectWrapper.simClasses.logger.close();
   objectWrapper.simClasses.pargrid.finalize();
   objectWrapper.simClasses.random.finalize();
   objectWrapper.sim.finalize();
   delete builder; builder = NULL;
   
   // Exit program:
   MPI_Finalize();
   return 0;
}
