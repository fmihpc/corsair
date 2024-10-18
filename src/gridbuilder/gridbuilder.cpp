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

#include "gridbuilder.h"

using namespace std;

GridBuilder::GridBuilder() { }

GridBuilder::~GridBuilder() { }

bool GridBuilder::build(Simulation& sim,SimulationClasses& simClasses) {
   bool success = true;
   
   simClasses.logger << "(GRIDBUILDER) Starting to build grid." << endl;
   
   // Get the total number of cells in the grid:
   builder::ID N_totalCells = 0;
   if (getNumberOfCells(N_totalCells) == false) {
      simClasses.logger << "\t Failed to get number of cells to create!" << endl;
      success = false;
   }
   
   // Do an initial partitioning of cells amongs the processes. 
   // Every process does the partitioning and arrives at the same result:
   builder::ID* cellOffsets = new builder::ID[sim.mpiProcesses+1];
   builder::ID* cellsPerProcess = new builder::ID[sim.mpiProcesses];
   cellOffsets[0] = 0;
   for (int p=0; p<sim.mpiProcesses; ++p) {
      cellsPerProcess[p] = N_totalCells / sim.mpiProcesses;
      if (p < static_cast<int>(N_totalCells % sim.mpiProcesses)) ++cellsPerProcess[p];
      cellOffsets[p+1] = cellOffsets[p] + cellsPerProcess[p];
   }
   const builder::ID N_cells = cellsPerProcess[sim.mpiRank];
   const builder::ID myCellsStart = cellOffsets[sim.mpiRank];
   const builder::ID myCellsEnd   = cellOffsets[sim.mpiRank+1];
   delete [] cellOffsets; cellOffsets = NULL;
   delete [] cellsPerProcess; cellsPerProcess = NULL;
   
   // Get cell global IDs and their number of neighbours:
   builder::ID* cellIDs = new builder::ID[N_cells];
   unsigned char* N_neighbours = new unsigned char[N_cells];
   if (getCellInfo(myCellsStart,myCellsEnd,cellIDs,N_neighbours) == false) {
      simClasses.logger << "\t Failed to get cell info!" << endl;
      success = false;
   }

   // Prepare to add cells to parallel grid:
   builder::ID neighbourSum = 0;
   builder::ID* neighbourIDs = NULL;
   unsigned char* neighbourTypeIDs = NULL;
   if (success == true) {
      // Count the total number of neighbours (needed for mem alloc):
      for (builder::ID i=0; i<N_cells; ++i) neighbourSum += N_neighbours[i];
      neighbourIDs = new builder::ID[neighbourSum];
      neighbourTypeIDs = new unsigned char[neighbourSum];
      // Get neighbour global IDs and type IDs:
      if (getCellNeighbours(myCellsStart,myCellsEnd,cellIDs,N_neighbours,neighbourIDs,neighbourTypeIDs) == false) {
         simClasses.logger << "\t Failed to get cell neighbour IDs!" << endl;
         success = false;
      }
   }
   
   // Add cells to parallel grid:
   #if PROFILE_LEVEL > 0
      static int profCreateCells = -1;
      profile::start("Mesh Cell Creation",profCreateCells);
   #endif
   if (success == true) {
      vector<pargrid::CellID> nbrIDs;
      vector<unsigned char> nbrTypeIDs;

      builder::ID counter = 0;
      for (builder::ID i=0; i<N_cells; ++i) {
         nbrIDs.clear();
         nbrTypeIDs.clear();

         const pargrid::CellID cellID = cellIDs[i];
         for (int j=0; j<N_neighbours[i]; ++j) {
            nbrIDs.push_back(neighbourIDs[counter+j]);
            nbrTypeIDs.push_back(neighbourTypeIDs[counter+j]);
         }
         
         if (simClasses.pargrid.addCell(cellID,nbrIDs,nbrTypeIDs) == false) success = false;
         counter += N_neighbours[i];
      }
   }
   if (success == true) {
      if (simClasses.pargrid.addCellFinished() == false) {
         simClasses.logger << "\t ParGrid addCellFinished failed!" << endl;
         success = false;
      }
   } else {
      simClasses.logger << "\t One ore more cell addition failed!" << endl;
   }   
   delete [] N_neighbours; N_neighbours = NULL;
   delete [] neighbourIDs; neighbourIDs = NULL;
   delete [] neighbourTypeIDs; neighbourTypeIDs = NULL;
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   
   // ParGrid uses local IDs, but we gave global IDs when add cells. 
   // Request local IDs from ParGrid:
   vector<pargrid::CellID> localIDs(N_cells);
   for (pargrid::CellID i=0; i<N_cells; ++i) {
      localIDs[i] = simClasses.pargrid.getLocalID(cellIDs[i]);
   }

   // Create ParGrid data array that contains cell coordinates:
   Simulation::crdsDataID = simClasses.pargrid.addUserData<Real>("cell_coordinates",3);
   if (sim.crdsDataID == simClasses.pargrid.invalidDataID()) {
      simClasses.logger << "\t Failed to create a ParGrid static data array for cell coordinates!" << endl;
      success = false;
   }
   Simulation::crdsStencilID = pargrid::DEFAULT_STENCIL;
   if (simClasses.pargrid.addDataTransfer(Simulation::crdsDataID,Simulation::crdsStencilID) == false) {
      simClasses.logger << "\t Failed to add a ParGrid transfer for cell coordinates array!" << endl;
      success = false;
   }
   if (success == true) {
      double* coordinates = reinterpret_cast<double*>(simClasses.pargrid.getUserData(Simulation::crdsDataID));      
      if (getInitialState(myCellsStart,myCellsEnd,cellIDs,coordinates) == false) {
         simClasses.logger << "\t Initial state calculation failed!" << endl;
         success = false;
      }
   }
   delete [] cellIDs; cellIDs = NULL;
   
   // Balance load:
   if (success == true) for (size_t cell=0; cell<N_cells; ++cell) {
      simClasses.pargrid.getCellWeights()[cell] = 1.0;
   }
   if (success == true) {
      if (simClasses.pargrid.balanceLoad() == false) {
         simClasses.logger << "\t Load balancing failed!" << endl;
         //success = false;
      }
      sim.meshChangedStep   = sim.timestep;
      sim.meshRepartitioned = true;
   }
   
   // Check that all processes succeeded in grid building:
   unsigned int successSum = 0;
   unsigned int mySuccess = 0;
   if (success == false) ++mySuccess;
   MPI_Allreduce(&mySuccess,&successSum,1,MPI_Type<unsigned int>(),MPI_SUM,sim.comm);
   if (successSum > 0) {
      simClasses.logger << "\t " << successSum << " processes failed in grid building!" << endl;
      success = false;
   }

   if(success == true) { simClasses.logger << "\t Build successful." << endl; }
   
   // Flush log message and exit:
   simClasses.logger << write;
   return success;
}

/** Write simulation mesh to output file.
 * @param sim Struct containing generic simulation variables.
 * @param simClasses Struct containing generic simulation classes.
 * @param meshName Name of the simulation mesh.
 * @param writeMesh If true, the mesh data is written out. If false, mesh exists in another VLSV file. 
 * Only the information that is necessary to find the mesh is written out.
 * @return If true, all processes succeeded in mesh writing.*/
bool GridBuilder::writeMesh(Simulation& sim,SimulationClasses& simClasses,const std::string& meshName,bool writeMesh) {
   bool success = true;
   
   // XML tag called MESH must appear on each VLSV file even if the mesh is not written out:
   map<string,string> attributes;
   if (writeMesh == false) {
      int* array = NULL;
      attributes["name"] = meshName;
      attributes["file"] = sim.meshFileName;
      attributes["type"] = vlsv::mesh::STRING_QUAD_MULTI;
      if (simClasses.vlsv.writeArray("MESH",attributes,0,1,array) == false) success = false;
      
      // If mesh is read from another file we can exit now:
      return success;
   }
   
   // Write mesh bounding box (master process only):
   Real array[6];
   if (sim.mpiRank == sim.MASTER_RANK) {
      if (getMeshBoundingBox(array[0],array[1],array[2],array[3],array[4],array[5]) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to obtain mesh bounding box!" << endl << write;
	 success = false;
      }
   }
   array[3] /= block::WIDTH_X;
   array[4] /= block::WIDTH_Y;
   array[5] /= block::WIDTH_Z;
   
   attributes.clear();
   attributes["mesh"] = meshName;
   if (sim.mpiRank == sim.MASTER_RANK) {
      if (simClasses.vlsv.writeArray("MESH_BBOX",attributes,6,1,array) == false) success = false;
   } else {
      if (simClasses.vlsv.writeArray("MESH_BBOX",attributes,0,1,array) == false) success = false;
   }
   
   // Check that everything is ok:
   if (success == false) return success;

   // Create ParGrid array that contains boundary cell local IDs and exchange data with neighbours:
   pargrid::DataID lidArrayID = simClasses.pargrid.addUserData<pargrid::CellID>("localIDs",block::SIZE);
   if (lidArrayID == simClasses.pargrid.invalidDataID()) {
      simClasses.logger << "(GRIDBUILDER) ERROR: Failed to create array for local IDs!" << endl << write;
      success = false;
   }   
   if (success == true) {
      bool ok = simClasses.pargrid.addDataTransfer(lidArrayID,pargrid::DEFAULT_STENCIL);
      if (ok == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to add transfer for local ID array!" << endl << write;
	 success = false;
      }
   }
   
   // Write block local IDs to all boundary blocks:
   pargrid::CellID* lidArray = NULL;
   if (success == true) {
      lidArray = simClasses.pargrid.getUserDataStatic<pargrid::CellID>(lidArrayID);
      const vector<pargrid::CellID>& boundaryBlocks = simClasses.pargrid.getBoundaryCells(pargrid::DEFAULT_STENCIL);
      for (size_t block=0; block<boundaryBlocks.size(); ++block) {
	 const pargrid::CellID blockLID = boundaryBlocks[block]*block::SIZE;
	 for (int i=0; i<block::SIZE; ++i) lidArray[blockLID+i] = blockLID + i;
      }
   }
   
   // Start transfer:
   if (success == true) {
      if (simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,lidArrayID) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to exchange local IDs!" << endl << write;
	 success = false;
      }
   }
   
   // Calculate mesh block (i,j,k) indices:
   builder::ID* indices = NULL;
   const pargrid::CellID N_blocks = simClasses.pargrid.getNumberOfAllCells();
   const std::vector<pargrid::CellID>& globalIDs = simClasses.pargrid.getGlobalIDs();
   if (getCellIndices(N_blocks,&(globalIDs[0]),indices) == false) {
      delete [] indices; indices = NULL; 
      success = false;
      simClasses.logger << "(GRIDBUILDER) ERROR: Failed to calculate mesh block (i,j,k) indices!" << endl << write;
   }
   
   // Calculate cell indices for all cells in the blocks:
   pargrid::CellID* indicesOut = NULL;
   if (success == true) {
      size_t counter = 0;
      indicesOut = new pargrid::CellID[N_blocks*3*block::SIZE];
      for (pargrid::CellID block=0; block<N_blocks; ++block) {
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    indicesOut[counter+0] = indices[block*3+0]*block::WIDTH_X + i;
	    indicesOut[counter+1] = indices[block*3+1]*block::WIDTH_Y + j;
	    indicesOut[counter+2] = indices[block*3+2]*block::WIDTH_Z + k;
	    counter += 3;
	 }
      }
   }
   delete [] indices; indices = NULL;
   
   // Write cell indices to file:
   if (success == true) {
      attributes.clear();
      attributes["name"] = meshName;
      attributes["type"] = vlsv::mesh::STRING_QUAD_MULTI;
      if (simClasses.vlsv.writeArray("MESH",attributes,N_blocks*block::SIZE,3,indicesOut) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write cell indices!" << endl << write;
	 success = false;
      }
   }
   delete [] indicesOut; indicesOut = NULL;
   
   // Write number of cells in this multimesh:
   attributes.clear();
   attributes["mesh"] = meshName;
   uint32_t N_zones[2] = {N_blocks*block::SIZE,(N_blocks-simClasses.pargrid.getNumberOfLocalCells())*block::SIZE};
   if (simClasses.vlsv.writeArray("MESH_ZONES",attributes,1,2,N_zones) == false) {
      simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write mesh domain size!" << endl << write;
      success = false;
   }
   
   // Wait for block local IDs:
   if (simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,lidArrayID) == false) {
      simClasses.logger << "(GRIDBUILDER) ERROR: Failed to exchange block local IDs!" << endl << write;
      success = false;
   }
   
   if (success == true) {
      // Write N_ghosts local IDs:
      const uint64_t N_ghosts = N_blocks-simClasses.pargrid.getNumberOfLocalCells();
      pargrid::CellID* ptr_localIDs = lidArray + simClasses.pargrid.getNumberOfLocalCells()*block::SIZE;
      if (simClasses.vlsv.writeArray("MESH_GHOST_LOCALIDS",attributes,N_ghosts*block::SIZE,1,ptr_localIDs) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write block local IDs!" << endl << write;
	 success = false;
      }
      
      // Write N_ghosts host IDs:  
      pargrid::CellID counter = 0;
      pargrid::MPI_processID* hosts = new pargrid::MPI_processID[N_ghosts*block::SIZE];
      for (pargrid::CellID block=simClasses.pargrid.getNumberOfLocalCells(); block<simClasses.pargrid.getNumberOfAllCells(); ++block) {
	 for (int i=0; i<block::SIZE; ++i) {
	    hosts[counter+i] = simClasses.pargrid.getHosts()[block];
	 }
	 counter += block::SIZE;
      }
      if (simClasses.vlsv.writeArray("MESH_GHOST_DOMAINS",attributes,N_ghosts*block::SIZE,1,hosts) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write block ghost domain IDs!" << endl << write;
	 success = false;
      }
      delete [] hosts; hosts = NULL;
   }
   
   // Remove ParGrid localID array created above:
   simClasses.pargrid.removeDataTransfer(pargrid::DEFAULT_STENCIL,lidArrayID);
   simClasses.pargrid.removeUserData(lidArrayID);
   return success;
}

