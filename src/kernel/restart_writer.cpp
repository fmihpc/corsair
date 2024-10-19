#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <stdint.h>
#include <iomanip>
#include <dirent.h>
#include <stdio.h>
#include <mpi.h>

#include <restart_writer.h>

using namespace std;

#if PROFILE_LEVEL > 0
   static int writeRestartID = -1;
#endif

/** Stop profiling and return given exit status. This is a small helper 
 * function used to clean up code in writeRestart function below.
 * @param status Exit status of writeRestart function.
 * @return Value given in parameter status.*/
bool exitStatus(const bool& status) {
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return status;
}

void removeOldRestartFiles(Simulation& sim) {
   if (sim.mpiRank != sim.MASTER_RANK) return;
   if (sim.restartMinorFileAmount == 0) return;
   
   map<size_t,string> restartFiles;
   
   const string directory = ".";
   const string prefix = sim.restartFilenamePrefix;
   const string suffix = ".vlsv";
   
   DIR* dir = opendir(directory.c_str());
   if (dir == NULL) return;
   
   struct dirent* entry = readdir(dir);
   while (entry != NULL) {
      // Check that the directory entry is a valid restart file:
      const string entryName = entry->d_name;
      const size_t prefixStart = entryName.find(prefix);
      if (prefixStart == string::npos) {
	 entry = readdir(dir);
	 continue;
      }
      if (entryName.find(suffix) == string::npos) {
	 entry = readdir(dir);
	 continue;
      }

      // Parse the timestep at which restart file was written:
      const size_t size = entryName.size() - prefix.size() - suffix.size();
      const size_t timestep = atoi(entryName.substr(prefix.size(),size).c_str());
      
      // Add the timestep and filename to map 'restartFiles':
      restartFiles[timestep] = entryName;
      entry = readdir(dir);
   }
   closedir(dir);
   
   // Skip filenames that are kept:
   if (restartFiles.size() <= (sim.restartMinorFileAmount-1)) return;
   map<size_t,string>::reverse_iterator it=restartFiles.rbegin();
   for (size_t i=0; i<(sim.restartMinorFileAmount-1); ++i) ++it;

   // Iterate over the rest of files and remove them unless 
   // the timestep is at 'major interval' in which case the file is always kept:
   while (it != restartFiles.rend()) {
      if (sim.restartMajorInterval != 0 && it->first % sim.restartMajorInterval == 0) {
	 ++it;
	 continue;
      }
      
      remove(it->second.c_str());
      ++it;
   }
}

bool writeRestart(Simulation& sim,SimulationClasses& simClasses) {
   bool success = true;
   #if PROFILE_LEVEL > 0
      profile::start("Write restart",writeRestartID);
   #endif

   simClasses.logger << endl << "(RESTART WRITER) Starting to write restart data" << endl;
   removeOldRestartFiles(sim);
   
   // Name of the simulation mesh:
   const string meshName = "SpatialGrid";

   // Create output file name:
   string fileName;
   stringstream ss;
   ss.fill('0');
   ss << sim.restartFilenamePrefix << setw(8) << sim.timestep << ".vlsv";
   ss >> fileName;
   
   // Open output file for parallel writing:
   if (simClasses.vlsv.open(fileName,sim.comm,sim.MASTER_RANK) == false) {
      simClasses.logger << "\t Failed to open output file '" << fileName << "'!" << endl << write;
      return exitStatus(false);
   }

   // Simulation mesh is written as an old-style VLSV mesh:
   map<string,string> attribs;
   
   const uint64_t N_blocks = simClasses.pargrid.getNumberOfLocalCells();
   
   // Write necessary simulation parameters (only master process ends up writing these):
   if (simClasses.vlsv.writeParameter("time",&sim.t) == false) success = false;
   if (simClasses.vlsv.writeParameter("time step",&sim.timestep) == false) success = false;
   if (simClasses.vlsv.writeParameter("dt",&sim.dt) == false) success = false;
   
   // Write name of mesh (master process only):
   attribs["name"] = meshName;
   attribs["type"] = "mesh name";
   char dummyPtr = 'a';
   if (sim.mpiRank == sim.MASTER_RANK) {
      if (simClasses.vlsv.writeArray("MESH_NAME",attribs,1,1,&dummyPtr) == false) {
	 simClasses.logger << "\t Failed to write mesh name!" << endl << write;
	 success = false;
      }
   } else {
      if (simClasses.vlsv.writeArray("MESH_NAME",attribs,0,0,&dummyPtr) == false) {
	 simClasses.logger << "\t Failed to write mesh name!" << endl << write;
	 success = false;
      }
   }
   
   // Write global IDs of all local blocks as a variable:
   attribs["name"] = "CellID";
   attribs["mesh"] = meshName;
   attribs["type"] = "celldata";
   if (simClasses.vlsv.writeArray("VARIABLE",attribs,N_blocks,1,&(simClasses.pargrid.getGlobalIDs()[0])) == false) {
      simClasses.logger << "\t Failed to write block global IDs!" << endl << write;
      success = false;
   }
   
   // Write global IDs of blocks' neighbours. Note that ParGrid stores neighbour local IDs
   // so we need to convert them to global IDs before writing:
   pargrid::NeighbourID* nbrTypeIDs = new pargrid::NeighbourID[N_blocks*(pargrid::N_neighbours-1)];
   pargrid::CellID* nbrGlobalIDs = new pargrid::CellID[N_blocks*(pargrid::N_neighbours-1)];
   for (pargrid::CellID block=0; block<N_blocks; ++block) {
      // Get pointer to array containing cell neighbour local IDs:
      pargrid::CellID* nbrIDs = simClasses.pargrid.getCellNeighbourIDs(block);
      
      int counter = 0;
      for (pargrid::NeighbourID n=0; n<pargrid::N_neighbours; ++n) {
	 // Block is not its own neighbour here:
	 if (n == simClasses.pargrid.calcNeighbourTypeID(0,0,0)) continue;
	 
	 nbrTypeIDs[block*(pargrid::N_neighbours-1)+counter]   = n;
	 if (nbrIDs[n] == pargrid::INVALID_CELLID) {
	    nbrGlobalIDs[block*(pargrid::N_neighbours-1)+counter] = pargrid::INVALID_CELLID;
	 } else {
	    nbrGlobalIDs[block*(pargrid::N_neighbours-1)+counter] = simClasses.pargrid.getGlobalIDs()[nbrIDs[n]];
	 }
	 ++counter;
      }
   }
   
   attribs["name"] = "NbrTypeID";
   attribs["mesh"] = meshName;
   attribs["type"] = "celldata";
   if (simClasses.vlsv.writeArray("VARIABLE",attribs,N_blocks,pargrid::N_neighbours-1,nbrTypeIDs) == false) {
      simClasses.logger << "\t Failed to write cell neighbour type IDs!" << endl << write;
      success = false;
   }
   delete [] nbrTypeIDs; nbrTypeIDs = NULL;
   
   attribs["name"] = "NbrGlobalID";
   attribs["mesh"] = meshName;
   attribs["type"] = "celldata";
   if (simClasses.vlsv.writeArray("VARIABLE",attribs,N_blocks,pargrid::N_neighbours-1,nbrGlobalIDs) == false) {
      simClasses.logger << "\t Failed to write cell neighbour type IDs!" << endl << write;
      success = false;
   }
   delete [] nbrGlobalIDs; nbrGlobalIDs = NULL;

   // Write all static user data arrays to restart file:
   vector<pargrid::DataID> arrayDataIDs;
   vector<string> arrayNames;
   vector<unsigned int> arrayElements;
   vector<string> arrayDatatypes;
   vector<unsigned int> arrayByteSizes;
   vector<const char*> arrayPointers;
   simClasses.pargrid.getStaticUserDataInfo(arrayDataIDs,arrayNames,arrayElements,arrayDatatypes,arrayByteSizes,arrayPointers);
   for (size_t i=0; i<arrayDataIDs.size(); ++i) {
      attribs["name"] = arrayNames[i];
      attribs["mesh"] = meshName;
      attribs["type"] = "celldata";      
      const uint64_t vectorSize = arrayElements[i];
      
      if (simClasses.vlsv.writeArray("STATIC",attribs,arrayDatatypes[i],N_blocks,vectorSize,arrayByteSizes[i],arrayPointers[i]) == false) {
	 simClasses.logger << "\t Failed to write static user data array '" << arrayNames[i] << "'!" << endl << write;
	 success = false;
      }
   }

   // Write all dynamic user data arrays to restart file. This is a bit more 
   // complicated process than with static user data because array elements are 
   // not stored in a single continuous array.
   vector<pargrid::ArraySizetype*> arraySizes;
   vector<const char**> arrayDynamicPointers;
   simClasses.pargrid.getDynamicUserDataInfo(arrayDataIDs,arrayNames,arraySizes,arrayDatatypes,arrayByteSizes,arrayDynamicPointers);
   for (size_t i=0; i<arrayDataIDs.size(); ++i) {
      // Write array that contains number of stored elements per block:
      attribs["name"] = arrayNames[i];
      attribs["mesh"] = meshName;
      attribs["type"] = "celldata";
      if (simClasses.vlsv.writeArray("SIZE(DYNAMIC)",attribs,N_blocks,1,&(arraySizes[i][0])) == false) {
	 simClasses.logger << "\t Failed to write dynamic user data '" << arrayNames[i] << "' size array!" << endl << write;
	 success = false;
      }

      // Count the number of array elements written by this process, currently 
      // this needs to be known before starting VLSV multiwrite:
      uint64_t myArraySize = 0;
      for (pargrid::CellID c=0; c<simClasses.pargrid.getNumberOfLocalCells(); ++c) {
	 myArraySize += arraySizes[i][c];
      }

      // Init multiwrite mode:
      if (simClasses.vlsv.startMultiwrite(arrayDatatypes[i],myArraySize,1,arrayByteSizes[i]) == false) {
	 simClasses.logger << "\t Failed to start multiwrite!" << endl << write;
	 success = false;
      }

      // Write elements stored in each block in multiwrite mode:
      for (pargrid::CellID c=0; c<simClasses.pargrid.getNumberOfLocalCells(); ++c) {
	 // NOTE: addMultiwriteUnit function call here binds to the template wrapper function
	 // which creates an MPI datatype for char pointer, when it should be a continuous
	 // array of chars with arrayByteSizes[i] elements. This is corrected by writing
	 // arraySizes[i][c]*arrayByteSizes[i] elements instead of arraySizes[i][c].
	 if (simClasses.vlsv.addMultiwriteUnit(arrayDynamicPointers[i][c],arraySizes[i][c]*arrayByteSizes[i]) == false) {
	    simClasses.logger << "\t Failed to add a multiwrite unit!" << endl << write;
	    success = false;
	 }
      }

      // End multiwrite, the actual file writing occurs here:
      if (simClasses.vlsv.endMultiwrite("DYNAMIC",attribs) == false) {
	 simClasses.logger << "\t Failed to end multiwrite!" << endl << write;
	 success = false;
      }
   }

   // Close VLSV file:
   if (simClasses.vlsv.close() == false) success = false;

   // Check that all processes succeeded:
   if (simClasses.pargrid.checkSuccess(success) == false) {
      simClasses.logger << "ERROR: One or more processes failed to write restart data" << endl;
      success = false;
   } else {
      simClasses.logger << "\t Restart file written successfully" << endl;
   }
   // Flush log file message and exit:
   simClasses.logger << write;
   return exitStatus(success);
}
