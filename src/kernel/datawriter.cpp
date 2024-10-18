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
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include "datawriter.h"

using namespace std;

#if PROFILE_LEVEL > 0
   static int fileCloseID = -1;
   static int fileOpenID = -1;
   static int writeMeshID = -1;
   static int dataOperatorsID = -1;
#endif

#if PROFILE_LEVEL > 1
   static int allReduceID = -1;
#endif

bool saveState(Simulation& sim,SimulationClasses& simClasses,DataOperatorContainer& dataOperatorContainer,
	       vector<ParticleListBase*>& particles,GridBuilder* builder) {
   bool success = true;
   #if PROFILE_LEVEL > 0
      profile::start("File Open",fileOpenID);
   #endif
   simClasses.logger << "(DATAWRITER) Starting to write simulation data" << endl;

   // Create name for output file:
   stringstream ss;
   ss.fill('0');
   ss << "state" << setw(8) << sim.timestep << ".vlsv";
   string fileName;
   ss >> fileName;
   
   const string meshName = "SpatialGrid";
   
   // Open output file for parallel writing:
   if (simClasses.vlsv.open(fileName,sim.comm,sim.MASTER_RANK) == false) {
      simClasses.logger << "\t Failed to open output file '" << fileName << "'" << endl;
      success = false;
   }
   #if PROFILE_LEVEL > 0
      profile::stop(); // file open
   #endif
   
   // Write time and timestep:
   const Real t_start = MPI_Wtime();
   map<string,string> attributes;
   if (simClasses.pargrid.getRank() == sim.MASTER_RANK) {
      attributes["name"] = "time";
      if (simClasses.vlsv.writeArray("PARAMETER",attributes,1,1,&sim.t) == false) success = false;
      attributes["name"] = "timestep";
      if (simClasses.vlsv.writeArray("PARAMETER",attributes,1,1,&sim.timestep) == false) success = false;
   } else {
      if (simClasses.vlsv.writeArray("PARAMETER",attributes,0,0,&sim.t) == false) success = false;
      if (simClasses.vlsv.writeArray("PARAMETER",attributes,0,0,&sim.timestep) == false) success = false;
   }

   // Determine if mesh should be (re)written to output file:
   bool writeMesh = false;
   if (sim.meshChangedStep > sim.meshWrittenStep) writeMesh = true;
   if (sim.meshWrittenStep == numeric_limits<unsigned int>::max()) writeMesh = true;
   if (sim.meshAlwaysWritten == true) writeMesh = true;

   // Write mesh:
   #if PROFILE_LEVEL > 0
      profile::start("Mesh Write",writeMeshID);
   #endif
   if (builder->writeMesh(sim,simClasses,meshName,writeMesh) == false) {
      simClasses.logger << "(DATAWRITER) ERROR: Failed to write simulation mesh!" << endl << write;
      success = false;
   }
   if (success == true && writeMesh == true) {
      sim.meshFileName = fileName;
      sim.meshWrittenStep = sim.timestep;
   }
   #if PROFILE_LEVEL > 0
      profile::stop(); // mesh writing
   #endif
   
   // Write all DataOperators:
   #if PROFILE_LEVEL > 0
     profile::start("Data Operators",dataOperatorsID);
   #endif
   if (dataOperatorContainer.writeData(meshName,particles) == false) {
      simClasses.logger << "\t ERROR Failed to write data operators to output file!" << endl;
      success = false;
   }
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif

   const Real t_total = MPI_Wtime() - t_start;
   const uint64_t bytesWritten = simClasses.vlsv.getBytesWritten();
   
   // Close output file:
   #if PROFILE_LEVEL > 0
      profile::start("File Close",fileCloseID);
   #endif
   
   if (simClasses.vlsv.close() == false) {
      simClasses.logger << "\t Error occurred while closing output file '" << fileName << "'" << endl;
      success = false;
   }
   
   // Check that all processes succeeded in data writing:
   unsigned int successSum = 0;
   unsigned int mySuccess = 0;
   if (success == false) ++mySuccess;
   
   #if PROFILE_LEVEL > 1
      profile::start("allreduce",allReduceID);
   #endif
   
   MPI_Allreduce(&mySuccess,&successSum,1,MPI_Type<unsigned int>(),MPI_SUM,sim.comm);
   
   #if PROFILE_LEVEL > 1
      profile::stop(); // allreduce
   #endif
   
   if (successSum > 0) {
      simClasses.logger << "\t " << successSum << " processes failed to write data!" << endl;
      success = false;
   }

   // Flush log message and exit:
   Real divider = 1.0e3;
   string units = "kB";
   if (bytesWritten >= 1.0e9) {
      divider = 1.0e9;
      units = "GB";
   } else if (bytesWritten >= 1.0e6) {
      divider = 1.0e6;
      units = "MB";
   }
   
   Real datarate = bytesWritten/t_total;
   string datarateUnits = "kB/s";
   if (datarate >= 1.0e9) {
      datarateUnits = "GB/s";
      datarate /= 1.0e9;
   } else if (datarate >= 1.0e6) {
      datarateUnits = "MB/s";
      datarate /= 1.0e6;
   } else {
      datarate /= 1.0e3;
   }
   
   simClasses.logger
     << "\t wrote " << bytesWritten/divider << ' ' << units << " in " << t_total << " seconds, " << "throughput " << datarate << ' ' << datarateUnits << "." << endl;

   if(successSum < 1) { simClasses.logger << "\t Writing successful." << endl; }
   
   #if PROFILE_LEVEL > 0
      profile::stop(); // datawriter finalize
   #endif
   return success;
}
