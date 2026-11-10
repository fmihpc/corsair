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

#ifndef SIMULATIONCLASSES_H
#define SIMULATIONCLASSES_H

#include <vector>
#include <mpi.h>

#include <pargrid.h>
#include <vlsv_writer.h>

#include <constants.h>
#include <mpilogger.h>
#include <simulation.h>
#include <randomnumber.h>
#include <cell.h>

struct SimulationClasses {
   Constants constants;
   MPILogger logger;                             /**< Log file writer.*/
   pargrid::ParGrid<Cell<Real> > pargrid;        /**< Parallel grid.*/
   RandomNumber random;                          /**< Random number generator.*/
   vlsv::Writer vlsv;                            /**< Parallel file writer.*/
};

namespace block {
   template<typename INT>
   INT arrayIndex(INT i,INT j,INT k);

   template<typename INT>
   INT arrayIndex_XY(INT i,INT j);

   template<typename INT>
   INT index(INT i,INT j,INT k);

   template<typename INT>
   INT index(INT i,INT j);

   template<typename T>
   void addValues2D_XZ(SimulationClasses& simClasses,pargrid::CellID block,const T* const array,T* const data);
   template<typename T>
   void addValues3D(SimulationClasses& simClasses,pargrid::CellID block,const T* const array,T* const data,const int vectorDim=1);

   template<typename INT> inline
   INT calculateGlobalIndex(Simulation& sim,INT& i,INT& j,INT& k);

   template<typename T>
   void fetchValues2D_XZ(SimulationClasses& simClasses,pargrid::CellID block,T* const dest,const T* const source);
   template<typename T>
   void fetchValues3D(SimulationClasses& simClasses,pargrid::CellID block,T* const dest,const T* const source);
   template<typename T>
   void fetchValues3D_face(SimulationClasses& simClasses,pargrid::CellID block,T* const array,const T* const data);
}

/** Get size of mesh block. blockSize[0] is the sum of individual cell's x-sizes, 
 * i.e. blockSize[0]=sum(dx_i), where dx_i is the size of cell i (within the block) 
 * in x-direction. yz-sizes are calculated similarly. This definition works for 
 * uniform and non-uniform cell sizes.
 * @param simClasses Struct containing generic simulation classes.
 * @param sim
 * @param blockLID Local ID of the block.
 * @param blockSize Array in which block (x,y,z) sizes are written to.
 */
template<typename REAL> inline
void getBlockSize(const SimulationClasses& simClasses,const Simulation& sim,pargrid::CellID blockLID,REAL* blockSize) {
   blockSize[0] = sim.dx_block[0];
   blockSize[1] = sim.dy_block[0];
   blockSize[2] = sim.dz_block[0];
}

template<typename REAL> inline
void getBlockCellSize(const SimulationClasses& simClasses,const Simulation& sim,pargrid::CellID blockLID,REAL* cellSize) {
   cellSize[0] = sim.dx_cell[0];
   cellSize[1] = sim.dy_cell[0];
   cellSize[2] = sim.dz_cell[0];
}

inline const double* getBlockCoordinateArray(const Simulation& sim,SimulationClasses& simClasses) {
   return reinterpret_cast<const double*>(simClasses.pargrid.getUserData(Simulation::crdsDataID));
}

namespace block {

   /** Two-dimensional accumulator. This function accumulates values from given array 
    * to x- and z-neighbours of the given block.
    * Note that this function works for scalar fields only.
    * @param simClasses Struct containing generic simulation classes.
    * @param block Local ID of the block.
    * @param array Array containing the accumulated values.
    * @param data Array in which values from 'array' are added to.*/
   template<typename T> inline
   void addValues2D_XZ(SimulationClasses& simClasses,pargrid::CellID block,const T* const array,T* const data) {
      pargrid::CellID nbrLID = simClasses.pargrid.invalid();
      const pargrid::CellID* const nbrs = simClasses.pargrid.getCellNeighbourIDs(block);

      // ***** ACCUMULATE TO FACE NEIGHBOURS ***** //
      // -x neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k)
	  data[nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,k)] += array[block::arrayIndex(0,1,k+1)];
      // +x neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k)
	  data[nbrLID*block::SIZE+block::index(0,0,k)] += array[block::arrayIndex(block::WIDTH_X+1,1,k+1)];
      // -z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int i=0; i<block::WIDTH_X; ++i)
	  data[nbrLID*block::SIZE+block::index(i,0,block::WIDTH_Z-1)] += array[block::arrayIndex(i+1,1,0)];
      // +z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int i=0; i<block::WIDTH_X; ++i)
	  data[nbrLID*block::SIZE+block::index(i,0,0)] += array[block::arrayIndex(i+1,1,block::WIDTH_Z+1)];

      // ***** ACCUMULATE TO EDGE NEIGHBOURS ***** //
      // -x,-z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j)
	  data[nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,block::WIDTH_Z-1)] += array[block::arrayIndex(0,j+1,0)];
      // +x,-z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j)
	  data[nbrLID*block::SIZE+block::index(0,j,block::WIDTH_Z-1)] += array[block::arrayIndex(block::WIDTH_X+1,j+1,0)];
      // -x,+z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j)
	  data[nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,0)] += array[block::arrayIndex(0,j+1,block::WIDTH_Z+1)];
      // +x,+z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j)
	  data[nbrLID*block::SIZE+block::index(0,j,0)] += array[block::arrayIndex(block::WIDTH_X+1,j+1,block::WIDTH_Z+1)];

      // ***** ACCUMULATE TO THIS BLOCK ***** //
      for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) {
	 data[block*block::SIZE+block::index(i,0,k)] += array[block::arrayIndex(i+1,+1,k+1)];
      }
   }

   /** Three-dimensional accumulator. This function accumulates values from the given 
    * array to all existing neighbours of the given block.
    * Note that this function works for scalar fields only.
    * @param simClasses Struct containing generic simulation classes.
    * @param block Local ID of the block.
    * @param array Array containing the values to be accumulated.
    * @param data Array in which values from 'array' are added to.*/
   template<typename T> inline
   void addValues3D(SimulationClasses& simClasses,pargrid::CellID block,const T* const array,T* const data,const int vectorDim) {
      pargrid::CellID nbrLID = simClasses.pargrid.invalid();
      const pargrid::CellID* const nbrs = simClasses.pargrid.getCellNeighbourIDs(block);

      // ***** ACCUMULATE TO FACE NEIGHBOURS ***** //
      // Add values to -x neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,k))*vectorDim+l] += array[block::arrayIndex(0,j+1,k+1)*vectorDim+l];
      // Add values to -y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,k))*vectorDim+l] += array[block::arrayIndex(i+1,0,k+1)*vectorDim+l];
      // Add values to -z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(i,j,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(i+1,j+1,0)*vectorDim+l];
      // Add values to +x neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,j,k))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,j+1,k+1)*vectorDim+l];
      // Add values to +y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(i,0,k))*vectorDim+l] += array[block::arrayIndex(i+1,block::WIDTH_Y+1,k+1)*vectorDim+l];
      // Add values to +z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(i,j,0))*vectorDim+l] += array[block::arrayIndex(i+1,j+1,block::WIDTH_Z+1)*vectorDim+l];

      // ***** ACCUMULATE TO EDGE NEIGHBOURS ***** //
      // Add values to -x,-y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,k))*vectorDim+l] += array[block::arrayIndex(0,0,k+1)*vectorDim+l];
      // Add values to +x,-y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,k))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,0,k+1)*vectorDim+l];
      // Add values to -x,+y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,k))*vectorDim+l] += array[block::arrayIndex(0,block::WIDTH_Y+1,k+1)*vectorDim+l];
      // Add values to +x,+y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+0)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int k=0; k<block::WIDTH_Z; ++k) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,0,k))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,k+1)*vectorDim+l];

      // Add values to -y,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(i+1,0,0)*vectorDim+l];
      // Add values to +y,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(i,0,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(i+1,block::WIDTH_Y+1,0)*vectorDim+l];
      // Add values to -y,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,0))*vectorDim+l] += array[block::arrayIndex(i+1,0,block::WIDTH_Z+1)*vectorDim+l];
      // Add values to +y,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(i,0,0))*vectorDim+l] += array[block::arrayIndex(i+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorDim+l];

      // Add values to -x,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(0,j+1,0)*vectorDim+l];
      // Add values to +x,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,j,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,j+1,0)*vectorDim+l];
      // Add values to -x,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,0))*vectorDim+l] += array[block::arrayIndex(0,j+1,block::WIDTH_Z+1)*vectorDim+l];
      // Add values to +x,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int j=0; j<block::WIDTH_Y; ++j) for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,j,0))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,j+1,block::WIDTH_Z+1)*vectorDim+l];

      // ***** ACCUMULATE TO CORNER NEIGHBOURS ***** //
      // Add value to -x,-y,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(0,0,0)*vectorDim+l];
      // Add value to +x,-y,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,0,0)*vectorDim+l];
      // Add value to -x,+y,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(0,block::WIDTH_Y+1,0)*vectorDim+l];
      // Add value to +x,+y,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,-1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,0,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,0)*vectorDim+l];
      // Add value to -x,-y,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,0))*vectorDim+l] += array[block::arrayIndex(0,0,block::WIDTH_Z+1)*vectorDim+l];
      // Add value to +x,-y,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,0))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,0,block::WIDTH_Z+1)*vectorDim+l];
      // Add value to -x,+y,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,0))*vectorDim+l] += array[block::arrayIndex(0,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorDim+l];
      // Add value to +x,+y,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+1)];
      if (nbrLID != simClasses.pargrid.invalid())
	for (int l=0;l<vectorDim;++l)
	  data[(nbrLID*block::SIZE+block::index(0,0,0))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorDim+l];

      // ***** ACCUMULATE TO THIS BLOCK ***** //
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (int l=0;l<vectorDim;++l)
	 data[(block*block::SIZE+block::index(i,j,k))*vectorDim+l] += array[block::arrayIndex(i+1,j+1,k+1)*vectorDim+l];
   }

   /** Calculate index to temporary computation array. It is assumed here that the temporary 
    * array is two cells wider per coordinate direction than the block.
    * Note that this function works for scalar fields only.
    * @param i i-index of the cell in temporary array.
    * @param j j-index of the cell in temporary array.
    * @param k k-index of the cell in temporary array.*/
   template<typename INT> inline
   INT arrayIndex(INT i,INT j,INT k) {return k*(block::WIDTH_Y+2)*(block::WIDTH_X+2) + j*(block::WIDTH_X+2) + i;}

   template<typename INT> inline
   INT arrayIndex_XY(INT i,INT j) {return j*(block::WIDTH_X+2) + i;}

   template<typename INT> inline
   INT calculateGlobalIndex(Simulation& sim,INT& i,INT& j,INT& k) {
      return k*sim.y_blocks*sim.x_blocks + j*sim.x_blocks + i;
   }

   template<typename INT> inline
   void calculateBlockIndices(Simulation& sim,pargrid::CellID blockGID,INT& i,INT& j,INT& k) {
      k = blockGID / (sim.y_blocks*sim.x_blocks);
      blockGID -= k*(sim.y_blocks*sim.x_blocks);
      j = blockGID / sim.x_blocks;
      blockGID -= j*sim.x_blocks;
      i = blockGID;
   }

   template<typename INT> inline
   INT index_XY(INT i,INT j) {return j*block::WIDTH_X + i;}

   /** Calculate index into block. Note that this function works for scalar fields only.
    * @param i i-index of the cell in the block.
    * @param j j-index of the cell in the block.
    * @param k k-index of the cell in the block.*/
   template<typename INT> inline
   INT index(INT i,INT j,INT k) {return k*block::WIDTH_Y*block::WIDTH_X + j*block::WIDTH_X + i;}

   /** Fetch data from simulation mesh to a temporary computation array. This
    * function copies data from block's xz-face neighbours (two-dimensional version).
    * Note that this function works for scalar fields only.
    * @param simClasses Struct containing generic simulation classes.
    * @param block Local ID of the block.
    * @param dest Pointer to temporary computation array in which data is copied.
    * @param source Pointer to array containing the fetched data.*/
   template<typename T> inline
   void fetchValues2D_XZ(SimulationClasses& simClasses,pargrid::CellID block,T* const dest,const T* const source) {
      const T* srcTmp = NULL;

      // ***** FETCH DATA FROM THIS BLOCK ***** //
      srcTmp = source + block*block::SIZE;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i)
	dest[block::arrayIndex(i+1,1,k+1)] = srcTmp[block::index(i,0,k)];

      pargrid::CellID nbrLID = simClasses.pargrid.invalidCellID();
      const pargrid::CellID* const nbrs = simClasses.pargrid.getCellNeighbourIDs(block);

      // ***** FETCH VALUES FROM FACE NEIGHBOURS ***** //
      // -x neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k)
	   dest[block::arrayIndex(0,1,k+1)] = srcTmp[block::index(block::WIDTH_X-1,0,k)];
      }
      // +x neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k)
	   dest[block::arrayIndex(block::WIDTH_X+1,1,k+1)] = srcTmp[block::index(0,0,k)];
      }
      // -z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,1,0)] = srcTmp[block::index(i,0,block::WIDTH_Z-1)];
      }
      // +z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,1,block::WIDTH_Z+1)] = srcTmp[block::index(i,0,0)];
      }

      // ***** FETCH VALUES FROM CORNER NEIGHBOURS ***** //
      // -x,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(0,1,0)] = srcTmp[block::index(block::WIDTH_X-1,0,block::WIDTH_Z-1)];
      }
      // +x,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(block::WIDTH_X+1,1,0)] = srcTmp[block::index(0,0,block::WIDTH_Z-1)];
      }
      // -x,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(0,1,block::WIDTH_Z+1)] = srcTmp[block::index(block::WIDTH_X-1,0,0)];
      }
      // +x,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(block::WIDTH_X+1,1,block::WIDTH_Z+1)] = srcTmp[block::index(0,0,0)];
      }
   }

   template<typename T> inline
   void fetchValues3D(SimulationClasses& simClasses,pargrid::CellID blockLID,T* RESTRICT const dest,const T* RESTRICT const source) {
      pargrid::CellID nbrLID = simClasses.pargrid.invalidCellID();
      const pargrid::CellID* const nbrs = simClasses.pargrid.getCellNeighbourIDs(blockLID);
      const T* srcTmp = NULL;

      // ***** FETCH DATA FROM THIS BLOCK ***** //
      srcTmp = source + blockLID*block::SIZE;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
	dest[block::arrayIndex(i+1,j+1,k+1)] = srcTmp[block::index(i,j,k)];

      // ***** LOAD VALUES FROM 6 FACE NEIGHBOURS ***** //
      // Copy from -x neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j)
	   dest[block::arrayIndex(0,j+1,k+1)] = srcTmp[block::index(block::WIDTH_X-1,j,k)];
      }
      // Copy from +x neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j)
	   dest[block::arrayIndex(block::WIDTH_X+1,j+1,k+1)] = srcTmp[block::index(0,j,k)];
      }
      // Copy from -y neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,0,k+1)] = srcTmp[block::index(i,block::WIDTH_Y-1,k)];
      }
      // Copy from +y neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,block::WIDTH_Y+1,k+1)] = srcTmp[block::index(i,0,k)];
      }
      // Copy from -z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,j+1,0)] = srcTmp[block::index(i,j,block::WIDTH_Z-1)];
      }
      // Copy from +z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,j+1,block::WIDTH_Z+1)] = srcTmp[block::index(i,j,0)];
      }

      // ***** LOAD VALUES FROM 12 EDGE NEIGHBOURS ***** //
      // -x,-y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k)
	   dest[block::arrayIndex(0,0,k+1)] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,k)];
      }
      // +x,-y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k)
	   dest[block::arrayIndex(block::WIDTH_X+1,0,k+1)] = srcTmp[block::index(0,block::WIDTH_Y-1,k)];
      }
      // -x,+y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) 
	   dest[block::arrayIndex(0,block::WIDTH_Y+1,k+1)] = srcTmp[block::index(block::WIDTH_X-1,0,k)];
      }
      // +x,+y neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k)
	   dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,k+1)] = srcTmp[block::index(0,0,k)];
      }
      // -x,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j)
	   dest[block::arrayIndex(0,j+1,0)] = srcTmp[block::index(block::WIDTH_X-1,j,block::WIDTH_Z-1)];
      }
      // +x,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j)
	   dest[block::arrayIndex(block::WIDTH_X+1,j+1,0)] = srcTmp[block::index(0,j,block::WIDTH_Z-1)];
      }
      // -x,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j)
	   dest[block::arrayIndex(0,j+1,block::WIDTH_Z+1)] = srcTmp[block::index(block::WIDTH_X-1,j,0)];
      }
      // +x,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j)
	   dest[block::arrayIndex(block::WIDTH_X+1,j+1,block::WIDTH_Z+1)] = srcTmp[block::index(0,j,0)];
      }
      // -y,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,0,0)] = srcTmp[block::index(i,block::WIDTH_Y-1,block::WIDTH_Z-1)];
      }
      // +y,-z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,block::WIDTH_Y+1,0)] = srcTmp[block::index(i,0,block::WIDTH_Z-1)];
      }
      // -y,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,0,block::WIDTH_Z+1)] = srcTmp[block::index(i,block::WIDTH_Y-1,0)];
      }
      // +y,+z neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i)
	   dest[block::arrayIndex(i+1,block::WIDTH_Y+1,block::WIDTH_Z+1)] = srcTmp[block::index(i,0,0)];
      }

      // ***** LOAD VALUES FROM 8 CORNER NEIGHBOURS ***** //
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(0,0,0)] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,block::WIDTH_Z-1)];
      }
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(block::WIDTH_X+1,0,0)] = srcTmp[block::index(0,block::WIDTH_Y-1,block::WIDTH_Z-1)];
      }
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(0,block::WIDTH_Y+1,0)] = srcTmp[block::index(block::WIDTH_X-1,0,block::WIDTH_Z-1)];
      }
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,0)] = srcTmp[block::index(0,0,block::WIDTH_Z-1)];
      }
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(0,0,block::WIDTH_Z+1)] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,0)];
      }
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(block::WIDTH_X+1,0,block::WIDTH_Z+1)] = srcTmp[block::index(0,block::WIDTH_Y-1,0)];
      }
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(0,block::WIDTH_Y+1,block::WIDTH_Z+1)] = srcTmp[block::index(block::WIDTH_X-1,0,0)];
      }
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*block::SIZE;
	 dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,block::WIDTH_Z+1)] = srcTmp[block::index(0,0,0)];
      }
   }

   /** Fetch data from simulation mesh to a temporary computation array. This 
    * function copies data from block's face neighbours (three-dimensional version).
    * Note that this function works for scalar fields only.
    * @param simClasses Struct containing generic simulation classes.
    * @param block Local ID of the block.
    * @param array Pointer to temporary computation array in which data is copied.
    * @param data Pointer to array containing the fetched data.*/
   template<typename T> inline
   void fetchValues3D_face(SimulationClasses& simClasses,pargrid::CellID block,T* __restrict__ const array,const T* __restrict__ const data) {
      pargrid::CellID nbrLID = simClasses.pargrid.invalidCellID();
      const pargrid::CellID* const nbrs = simClasses.pargrid.getCellNeighbourIDs(block);
      const T* dataTmp = NULL;

      // ***** FETCH DATA FROM THIS BLOCK ***** //
      dataTmp = data + block*block::SIZE;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
	array[block::arrayIndex(i+1,j+1,k+1)] = dataTmp[block::index(i,j,k)];

      // ***** FETCH DATA FROM FACE NEIGHBOURS ***** //
      // Copy from -x neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 dataTmp = data + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j)
	   array[block::arrayIndex(0,j+1,k+1)] = dataTmp[block::index(block::WIDTH_X-1,j,k)];
      }
      // Copy from +x neighbour:
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 dataTmp = data + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j)
	   array[block::arrayIndex(block::WIDTH_X+1,j+1,k+1)] = dataTmp[block::index(0,j,k)];
      }
      // Copy from -y neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 dataTmp = data + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i)
	   array[block::arrayIndex(i+1,0,k+1)] = dataTmp[block::index(i,block::WIDTH_Y-1,k)];
      }
      // Copy from +y neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 dataTmp = data + nbrLID*block::SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i)
	   array[block::arrayIndex(i+1,block::WIDTH_Y+1,k+1)] = dataTmp[block::index(i,0,k)];
      }
      // Copy from -z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 dataTmp = data + nbrLID*block::SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
	   array[block::arrayIndex(i+1,j+1,0)] = dataTmp[block::index(i,j,block::WIDTH_Z-1)];
      }
      // Copy from +z neighbour
      nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1)];
      if (nbrLID != simClasses.pargrid.invalidCellID()) {
	 dataTmp = data + nbrLID*block::SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
	   array[block::arrayIndex(i+1,j+1,block::WIDTH_Z+1)] = dataTmp[block::index(i,j,0)];
      }
   }

   template<typename REAL> inline
   void getBlockSizes(Simulation& sim,pargrid::CellID blockGID,REAL* RESTRICT blockSize) {
      // Calculate block i,j,k indices:
      unsigned int i_block,j_block,k_block;
      calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);

      // Copy values to output array:
      blockSize[0] = sim.dx_block[i_block];
      blockSize[1] = sim.dx_block[j_block];
      blockSize[2] = sim.dx_block[k_block];
   }

   template<typename REAL> inline
   void getBlockSizes(Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockLID,REAL* RESTRICT blockSize) {
      // Get block global ID:
      const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];

      // Calculate block i,j,k indices:
      unsigned int i_block,j_block,k_block;
      calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);

      // Copy values to output array:
      blockSize[0] = sim.dx_block[i_block];
      blockSize[1] = sim.dx_block[j_block];
      blockSize[2] = sim.dx_block[k_block];
   }

   /** Get cell (x,y,z) sizes in a mesh block. Output pointers may be indexed as
    * dx[i], where 0 < i_cell < block::WIDTH_X, and similarly for dy,dz arrays.
    * @param sim 
    * @param blockGID Global ID of the mesh block.
    * @param dx Pointer that can be indexed to get mesh block's cells' x-sizes.
    * @param dy Pointer that can be indexed to get mesh block's cells' y-sizes.
    * @param dz Pointer that can be indexed to get mesh block's cells' z-sizes.*/
   template<typename REAL> inline
   void getCellSizes(Simulation& sim,pargrid::CellID blockGID,REAL*& RESTRICT dx,REAL*& RESTRICT dy,REAL*& RESTRICT dz) {
      // Calculate block i,j,k indices;
      unsigned int i_block,j_block,k_block;
      calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);

      // Set output pointers:
      dx = sim.dx_cell + i_block*block::WIDTH_X;
      dy = sim.dy_cell + j_block*block::WIDTH_Y;
      dz = sim.dz_cell + k_block*block::WIDTH_Z;
   }

   /** Get cell (x,y,z) sizes in a mesh block. Output pointers may be indexed as
    * dx[i], where 0 < i_cell < block::WIDTH_X, and similarly for dy,dz arrays.
    * @param sim 
    * @param simClasses Struct containing generic simulation classes.
    * @param blockLID Local ID of the mesh block.
    * @param dx Pointer that can be indexed to get mesh block's cells' x-sizes.
    * @param dy Pointer that can be indexed to get mesh block's cells' y-sizes.
    * @param dz Pointer that can be indexed to get mesh block's cells' z-sizes.*/
   template<typename REAL> inline
   void getCellSizes(Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockLID,
		     REAL*& RESTRICT dx,REAL*& RESTRICT dy,REAL*& RESTRICT dz) {
      // Get block global ID:
      const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];

      // Calculate block i,j,k indices;
      unsigned int i_block,j_block,k_block;
      calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);

      // Set output pointers:
      dx = sim.dx_cell + i_block*block::WIDTH_X;
      dy = sim.dy_cell + j_block*block::WIDTH_Y;
      dz = sim.dz_cell + k_block*block::WIDTH_Z;
   }

} // namespace block

#endif
