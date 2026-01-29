/** This file is part of Corsair simulation.
 *
 *  Copyright 2011,2012,2015 Finnish Meteorological Institute
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
#include <stdint.h>
#include <map>
#include <list>

#include <vlsv_reader.h>
#include <vlsv_reader_parallel.h>
#include "restart_builder.h"

using namespace std;

static vlsv::ParallelReader* vlsvReader = NULL;

RestartBuilder::RestartBuilder(): GridBuilder() { 
   initialized = false;
   vlsvReader = NULL;
}

RestartBuilder::~RestartBuilder() { 
   finalize();
}

bool RestartBuilder::finalize() {
   if (vlsvReader != NULL) vlsvReader->close();
   delete vlsvReader; vlsvReader = NULL;
   initialized = false;
   return true;
}

/** This function is called by datawriter. RestartBuilder version always returns with a 
 * failed status as this function should never get called.
 * @return This function always returns false.*/
bool RestartBuilder::getCellIndices(unsigned int N_cells,const builder::ID* const cellIDs,builder::ID*& indices) {
   if (initialized == false) return false;
   return false;
}

/** Each process reads its cells' global IDs and number of neighbours from restart file.
 * @return If true cell global IDs were read successfully.*/
bool RestartBuilder::getCellInfo(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,unsigned char* N_neighbours) {
   if (initialized == false) return false;
   bool success = true;

   uint64_t arraySize;
   uint64_t vectorSize;
   uint64_t dataSize;
   vlsv::datatype::type datatype;

   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","CellID"));
   attribs.push_back(make_pair("mesh",meshName));
   if (vlsvReader->getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,datatype,dataSize) == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read cell ID array info!" << endl << write;
      success = false;
   }
   if (datatype != vlsv::datatype::UINT) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Unsupported cell ID datatype!" << endl << write;
      success = false;
   }
   if (vectorSize != 1) {
      simClasses->logger << "(RESTART BUILDER) ERROR: cell ID vectorsize is not unity!" << endl << write;
      success = false;
   }
   if (success == false) return success;

   char* buffer = new char[(cellOffsetEnd-cellOffsetStart)*vectorSize*dataSize];   
   if (vlsvReader->readArray("VARIABLE",attribs,cellOffsetStart,cellOffsetEnd-cellOffsetStart,buffer) == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read cell IDs!" << endl << write;
      success = false;
   }
   for (builder::ID i=0; i<cellOffsetEnd-cellOffsetStart; ++i) {
      char* ptr = buffer+i*dataSize;
      switch (dataSize) {
       case sizeof(uint16_t):
	 cellIDs[i] = *reinterpret_cast<uint16_t*>(ptr);
	 break;
       case sizeof(uint32_t):
	 cellIDs[i] = *reinterpret_cast<uint32_t*>(ptr);
	 break;
       case sizeof(uint64_t):
	 cellIDs[i] = *reinterpret_cast<uint64_t*>(ptr);
	 break;
       default:
	 success = false;
	 break;
      }
   }
   delete [] buffer; buffer = NULL;

   for (builder::ID i=0; i<cellOffsetEnd-cellOffsetStart; ++i) N_neighbours[i] = pargrid::N_neighbours-1;
   return success;
}

bool RestartBuilder::getCellNeighbours(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,
				       unsigned char* N_neighbours,builder::ID* neighbourIDs,unsigned char* neighbourTypeIDs) {
   if (initialized == false) return false;
   bool success = true;

   uint64_t arraySize;
   uint64_t vectorSize;
   uint64_t dataSize;
   vlsv::datatype::type datatype;

   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("name","NbrGlobalID"));
   attribs.push_back(make_pair("mesh",meshName));
   if (vlsvReader->getArrayInfo("VARIABLE",attribs,arraySize,vectorSize,datatype,dataSize) == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read cell neighbour ID array info!" << endl << write;
      success = false;
   }
   if (datatype != vlsv::datatype::UINT) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Unsupported neighbour global ID datatype!" << endl << write;
      success = false;
   }
   if (vectorSize != (pargrid::N_neighbours-1)) {
      simClasses->logger << "(RESTART BUILDER) ERROR: NbrGlobalID vectorsize is not 26!" << endl << write;
      success = false;
   }
   if (success == false) return success;

   char* buffer = new char[(cellOffsetEnd-cellOffsetStart)*vectorSize*dataSize];
   if (vlsvReader->readArray("VARIABLE",attribs,cellOffsetStart,cellOffsetEnd-cellOffsetStart,buffer) == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read cell neighbour IDs!" << endl << write;
      success = false;
   }
   for (builder::ID i=0; i<cellOffsetEnd-cellOffsetStart; ++i) {
      for (uint64_t j=0; j<vectorSize; ++j) {
	 char* ptr = buffer+i*dataSize*vectorSize + j*dataSize;
	 switch (dataSize) {
	  case sizeof(uint16_t):
	    neighbourIDs[i*vectorSize+j] = *reinterpret_cast<uint16_t*>(ptr);
	    break;
	  case sizeof(uint32_t):
	    neighbourIDs[i*vectorSize+j] = *reinterpret_cast<uint32_t*>(ptr);
	    break;
	  case sizeof(uint64_t):
	    neighbourIDs[i*vectorSize+j] = *reinterpret_cast<uint64_t*>(ptr);
	    break;
	  default:
	    success = false;
	    break;
	 }
      }
   }
   delete [] buffer; buffer = NULL;

   for (builder::ID i=0; i<cellOffsetEnd-cellOffsetStart; ++i) {
      int counter = 0;
      for (uint64_t j=0; j<pargrid::N_neighbours; ++j) {
	 if (j == simClasses->pargrid.calcNeighbourTypeID(0,0,0)) continue;
	 neighbourTypeIDs[i*(pargrid::N_neighbours-1)+counter] = j;
	 ++counter;
      }
   }
   return success;
}

/** This function is called by datawriter. RestartBuilder version always returns with a 
 * failed status as this function should never get called.
 * @return This function always returns false.*/
bool RestartBuilder::getMeshBoundingBox(Real& x_min,Real& y_min,Real& z_min,Real& dx0,Real& dy0,Real& dz0) {
   if (initialized == false) return false;
   return false;
}

bool RestartBuilder::getInitialState(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,double* coordinates) {
   if (initialized == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Class not initialized!" << endl << write;
      return false;
   }
   bool success = true;

   uint64_t arraySize;
   uint64_t vectorSize;
   uint64_t dataSize;
   vlsv::datatype::type datatype;

   list<pair<string,string> > attribs;
   attribs.push_back(make_pair("mesh",meshName));
   attribs.push_back(make_pair("name","cell_coordinates"));
   if (vlsvReader->getArrayInfo("STATIC",attribs,arraySize,vectorSize,datatype,dataSize) == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read cell coordinate array info!" << endl << write;
      success = false;
   }   
   if (datatype != vlsv::datatype::FLOAT) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Unsupported cell coordinate datatype!" << endl << write;
      success = false;
   }
   if (vectorSize != 3) {
      simClasses->logger << "(RESTART BUILDER) ERROR: cell coordinate vectorsize is not 3!" << endl << write;
      success = false;
   }
   if (success == false) return success;

   char* buffer = new char[(cellOffsetEnd-cellOffsetStart)*vectorSize*dataSize];
   if (vlsvReader->readArray("STATIC",attribs,cellOffsetStart,cellOffsetEnd-cellOffsetStart,buffer) == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read cell coordinate array!" << endl << write;
      success = false;
   }

   for (builder::ID i=0; i<cellOffsetEnd-cellOffsetStart; ++i) {
      for (uint64_t j=0; j<vectorSize; ++j) {
	 char* ptr = buffer+i*dataSize*vectorSize + j*dataSize;
	 switch (dataSize) {
	  case sizeof(float):
	    coordinates[i*vectorSize+j] = *reinterpret_cast<float*>(ptr);
	    break;
	  case sizeof(double):
	    coordinates[i*vectorSize+j] = *reinterpret_cast<double*>(ptr);
	    break;
	  case sizeof(long double):
	    coordinates[i*vectorSize+j] = *reinterpret_cast<long double*>(ptr);
	    break;
	  default:
	    success = false;
	    break;
	 }
      }
   }
   delete [] buffer; buffer = NULL;

   // Values of generic simulation parameters are read here
   // and copied to struct Simulation:
   if (vlsvReader->readParameter("time",sim->t) == false) success = false;
   if (vlsvReader->readParameter("timestep",sim->timestep) == false) success = false;
   if (vlsvReader->readParameter("dt",sim->dt) == false) success = false;
   if (success == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read t, timestep, or dt from restart file" << endl << write;
   }

   // ParGrid static and dynamic data arrays are read from 
   // restart file, allocated, and their data is copied to memory:
   set<string> arrayNames;
   if (vlsvReader->getUniqueAttributeValues("STATIC","name",arrayNames) == false) success = false;
   if (success == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read static array names from restart file" << endl << write;
      return success;
   }
   for (set<string>::const_iterator it=arrayNames.begin(); it!=arrayNames.end(); ++it) {
      // Get info of inserted static array:
      attribs.clear();
      attribs.push_back(make_pair("name",*it));
      attribs.push_back(make_pair("mesh",meshName));
      if (vlsvReader->getArrayInfo("STATIC",attribs,arraySize,vectorSize,datatype,dataSize) == false) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read static array info" << endl << write;
	 success = false;
	 break;
      }

      string datatypeString;
      switch (datatype) {
       case vlsv::datatype::INT:
	 datatypeString = "int";
	 break;
       case vlsv::datatype::UINT:
	 datatypeString = "uint";
	 break;
       case vlsv::datatype::FLOAT:
	 datatypeString = "float";
	 break;
       default:
	 success = false;
	 simClasses->logger << "(RESTART BUILDER) ERROR: Unsupported datatype in static array" << endl << write;
	 break;
      }

      // Create ParGrid array:
      pargrid::DataID dataID = simClasses->pargrid.addUserData(*it,vectorSize,datatypeString,dataSize,false);
      if (dataID == pargrid::INVALID_DATAID) {
	 success = false;
	 simClasses->logger << "(RESTART BUILDER) ERROR: Failed to add static data array to ParGrid" << endl << write;
	 break;
      }

      // Copy contents from file to array. Here we can just request 
      // a char pointer to the array since we are doing a byte copy of values:
      char* arrayPtr = simClasses->pargrid.getUserDataStatic<char>(dataID);
      if (vlsvReader->readArray("STATIC",attribs,cellOffsetStart,cellOffsetEnd-cellOffsetStart,arrayPtr) == false) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read array data from file" << endl << write;
	 success = false;
	 break;
      }
   }

   // Recreate ParGrid dynamic data arrays:
   arrayNames.clear();
   if (vlsvReader->getUniqueAttributeValues("DYNAMIC","name",arrayNames) == false) success = false;
   if (success == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read dynamic array names from restart file" << endl << write;
      return success;
   }
   for (set<string>::const_iterator it=arrayNames.begin(); it!=arrayNames.end(); ++it) {
      // Get array datatype info:
      attribs.clear();
      attribs.push_back(make_pair("name",*it));
      attribs.push_back(make_pair("mesh",meshName));
      if (vlsvReader->getArrayInfo("DYNAMIC",attribs,arraySize,vectorSize,datatype,dataSize) == false) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read dynamic array '" << *it << "' info" << endl << write;
	 success = false;
	 break;
      }

      string datatypeString;
      switch (datatype) {
       case vlsv::datatype::UNKNOWN:
	 datatypeString = "unknown";
	 break;
       case vlsv::datatype::INT:
	 datatypeString = "int";
	 break;
       case vlsv::datatype::UINT:
	 datatypeString = "uint";
	 break;
       case vlsv::datatype::FLOAT:
	 datatypeString = "float";
	 break;
       default:
	 success = false;
	 simClasses->logger << "(RESTART BUILDER) ERROR: Unsupported datatype in dynamic array" << endl << write;
	 break;
      }

      // Create ParGrid array:
      pargrid::DataID dataID = simClasses->pargrid.addUserData(*it,vectorSize,datatypeString,dataSize,true);
      if (dataID == pargrid::INVALID_DATAID) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Failed to create dynamic array '" << *it << "' !" << endl << write;
	 success = false; break;
      }
      pargrid::DataWrapper<char> wrapper = simClasses->pargrid.getUserDataDynamic<char>(dataID);

      // Get sizes of dynamic array elements:
      attribs.clear();
      attribs.push_back(make_pair("name",*it));
      attribs.push_back(make_pair("mesh",meshName));
      if (vlsvReader->getArrayInfo("SIZE_DYNAMIC",attribs,arraySize,vectorSize,datatype,dataSize) == false) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Failed to obtain dynamic data '" << *it << "' size array info!" << endl << write;
	 success = false; break;
      }
      if (vectorSize != 1) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Dynamic data '" << *it << "' size array has incorrect vectorsize!" << endl << write;
	 success = false; break;
      }
      if (datatype != vlsv::datatype::UINT) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Dynamic data '" << *it << "' size array has unsupported datatype!" << endl << write;
	 success = false; break;
      }

      // Read dynamic data size array:
      char* buffer = new char[arraySize*vectorSize*dataSize];
      if (vlsvReader->readArray("SIZE_DYNAMIC",attribs,cellOffsetStart,cellOffsetEnd-cellOffsetStart,buffer) == false) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read dynamic array data from file" << endl << write;
	 success = false;
      }

      // Actual dynamic data is written as a consecutive array in VLSV file. Each process should
      // reads its portion of that array but we do not know the correct file offsets yet.
      // Here we compute them:
      const uint8_t*  ptr8 = NULL;
      const uint16_t* ptr16 = NULL;
      const uint32_t* ptr32 = NULL;
      const uint64_t* ptr64 = NULL;
      uint64_t elementSum = 0;
      switch (dataSize) {
       case sizeof(uint8_t):
	 ptr8  = reinterpret_cast<uint8_t*>(buffer);
	 for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block)
	   elementSum += ptr8[block];
	 break;
       case sizeof(uint16_t):
	 ptr16 = reinterpret_cast<uint16_t*>(buffer);
	 for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) 
	   elementSum += ptr16[block];
	 break;
       case sizeof(uint32_t):
	 ptr32 = reinterpret_cast<uint32_t*>(buffer);
	 for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) 
	   elementSum += ptr32[block];
	 break;
       case sizeof(uint64_t):
	 ptr64 = reinterpret_cast<uint64_t*>(buffer);
	 for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) 
	   elementSum += ptr64[block];
	 break;
      }

      // Collect element sums to master process who then computes correct offsets for everyone:
      uint64_t* globalOffsets = NULL;
      if (sim->mpiRank == sim->MASTER_RANK) globalOffsets = new uint64_t[sim->mpiProcesses];
      MPI_Gather(&elementSum,1,MPI_Type<uint64_t>(),globalOffsets,1,MPI_Type<uint64_t>(),sim->MASTER_RANK,sim->comm);
      if (sim->mpiRank == sim->MASTER_RANK) {
	 uint64_t offset = 0;
	 for (int p=0; p<sim->mpiProcesses; ++p) {
	    const uint64_t tmp = globalOffsets[p];
	    globalOffsets[p] = offset;
	    offset += tmp;
	 }
      }
      MPI_Scatter(globalOffsets,1,MPI_Type<uint64_t>(),&elementSum,1,MPI_Type<uint64_t>(),sim->MASTER_RANK,sim->comm);
      delete [] globalOffsets; globalOffsets = NULL;

      // Resize dynamic array elements to correct sizes and add multiread units to 
      // vlsv::ParallelReader for reading the data from VLSV file to memory:
      char* arrayPtr = NULL;
      attribs.clear();
      attribs.push_back(make_pair("name",*it));
      attribs.push_back(make_pair("mesh",meshName));
      if (vlsvReader->startMultiread("DYNAMIC",attribs) == false) {
         simClasses->logger << "(RESTART BUILDER) ERROR: Failed to start dynamic array '" << *it << "' multiread mode!" << endl << write;
         success = false;
      }
      switch (dataSize) {
       case sizeof(uint8_t):
         ptr8 = reinterpret_cast<uint8_t*>(buffer);
         for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
            wrapper.resize(block,ptr8[block]);
            arrayPtr = wrapper.data()[block];
            if (vlsvReader->addMultireadUnit(arrayPtr,ptr8[block]) == false) success = false;
         }
         break;
       case sizeof(uint16_t):
         ptr16 = reinterpret_cast<uint16_t*>(buffer);
         for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
            wrapper.resize(block,ptr16[block]);
            arrayPtr = wrapper.data()[block];
            if (vlsvReader->addMultireadUnit(arrayPtr,ptr16[block]) == false) success = false;
         }
         break;
       case sizeof(uint32_t):
         ptr32 = reinterpret_cast<uint32_t*>(buffer);
         for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
            wrapper.resize(block,ptr32[block]);
            arrayPtr = wrapper.data()[block];
            if (vlsvReader->addMultireadUnit(arrayPtr,ptr32[block]) == false) success = false;
         }
         break;
       case sizeof(uint64_t):
         ptr64 = reinterpret_cast<uint64_t*>(buffer);
         for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
            wrapper.resize(block,ptr64[block]);
            arrayPtr = wrapper.data()[block];
            if (vlsvReader->addMultireadUnit(arrayPtr,ptr64[block]) == false) success = false;
         }
         break;
      }
      delete [] buffer; buffer = NULL;

      if (vlsvReader->endMultiread(elementSum) == false) {
         simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read dynamic array '" << *it << "' !" << endl << write;
         success = false; break;
      }
   }
   return success;
}

bool RestartBuilder::getNumberOfCells(builder::ID& N_cells) {
   if (initialized == false) return initialized;
   bool success = true;

   N_cells = 0;
   map<string,string> attribsOut;
   list<pair<string,string> > attribsIn;
   attribsIn.push_back(make_pair("name","CellID"));
   attribsIn.push_back(make_pair("mesh",meshName));
   if (vlsvReader->getArrayAttributes("VARIABLE",attribsIn,attribsOut) == false) {
      simClasses->logger << "(RESTART BUILDER) ERROR: Failed to read number of cells in mesh!" << endl << write;
      success = false;
   } else {
      map<string,string>::const_iterator it=attribsOut.find("arraysize");
      if (it == attribsOut.end()) {
	 simClasses->logger << "(RESTART BUILDER) ERROR: Variable 'CellID' did not contain array size!" << endl << write;
	 success = false;
      } else {
	 N_cells = atoi(it->second.c_str());
      }
   }
   return success;
}

bool RestartBuilder::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
   initialized = true;
   this->sim = &sim;
   this->simClasses = &simClasses;
   simClasses.logger << "(RESTART BUILDER) Starting initialization." << endl;

   string restartFileName;
   const string regionName = "Restart";
   cr.add(regionName+".filename","Name of the restart file (string).",string(""));
   cr.parse();
   cr.get(regionName+".filename",restartFileName);

   vlsvReader = new vlsv::ParallelReader();
   if (vlsvReader->open(restartFileName,sim.comm,sim.MASTER_RANK,MPI_INFO_NULL) == false) {
      simClasses.logger << "\t Failed to open file '" << restartFileName << "' for restarting!" << endl;
      initialized = false;
   }

   // Attempt to read mesh name:
   map<string,string> attribsOut;
   list<pair<string,string> > attribsIn;
   attribsIn.push_back(make_pair("type","mesh name"));
   if (initialized == true) if (vlsvReader->getArrayAttributes("MESH_NAME",attribsIn,attribsOut) == false) {
      simClasses.logger << "(RESTART BUILDER) ERROR: Failed to get mesh name from restart file!" << endl;
      initialized = false;
   }
   if (initialized == true) {
      if (attribsOut.find("name") == attribsOut.end()) {
	 simClasses.logger << "(RESTART BUILDER) ERROR: Array 'MESH_NAME' did not contain mesh name!" << endl;
	 initialized = false;
      } else {
	 meshName = attribsOut["name"];
      }
   }
   if (initialized == true) {simClasses.logger << "\t Initialization successful." << endl; }
   simClasses.logger << write;
   return initialized;
}

// Registered to ObjectFactoryGeneric:
GridBuilder* RestartCreator() {return new RestartBuilder();}

