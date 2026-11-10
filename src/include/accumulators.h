/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2014 Finnish Meteorological Institute
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

// This file contains some general-use particle accumulation methods, 
// such as nearest grid point (NGP), cloud-in-cell (CIC), and 
// triangular-shaped cloud (TSC).
// 
// User is assumed to accumulate all particles on a given block into a 
// temporary array (parameter "array in all the methods below"), whose 
// contents are copied to a ParGrid data array afterwards. One can use 
// one of the data copying methods in simulationclasses.h here.

#ifndef ACCUMULATORS_H
#define ACCUMULATORS_H

#include <cmath>

#include <definitions.h>
#include <simulationclasses.h>

namespace sfactor {
   enum {W_X_N, W_X_C, W_X_P, W_Y_N, W_Y_C, W_Y_P, W_Z_N, W_Z_C, W_Z_P};
}

template<typename T,class PARTICLE>
void accumulateScalarCIC_XYZ(const Real* cellSizes,T* array,T value,const PARTICLE& particle);
template<typename T,class PARTICLE>
void accumulateScalarNGP_XYZ(const Real* cellSizes,T* array,T value,const PARTICLE& particle);
template<typename T,class PARTICLE>
void accumulateScalarTSC_XYZ(const Real* cellSizes,T* array,T value,const PARTICLE& particle);

template<typename T>
void addValues3D(SimulationClasses& simClasses,pargrid::CellID block,const T* const array,T* const data,const int vectorDim,const int stride);

template<typename REAL>
void evaluateShapeFactorsNGP_3D(const REAL* const sizes,const REAL* const pos,REAL* shapeFactors,int i,int j,int k);
template<typename REAL>
void evaluateShapeFactorsTSC_3D(const REAL* const sizes,const REAL* const pos,REAL* shapeFactors,int i,int j,int k);

template<typename T>
void fetchValues3D(SimulationClasses& simClasses,pargrid::CellID blockLID,T* RESTRICT const dest,const T* RESTRICT const source,int vectorsize);

/** Three-dimensional accumulator. This function accumulates values from the given 
 * array to all existing neighbours of the given block.
 * Note that this function works for scalar fields only.
 * @param simClasses Struct containing generic simulation classes.
 * @param block Local ID of the block.
 * @param array Array containing the values to be accumulated.
 * @param data Array in which values from 'array' are added to.
 * @param vectorDim Size of the data vector.
 * @param stride Element of data vector which is accumulated.*/
template<typename T> inline
void addValues3D(SimulationClasses& simClasses,pargrid::CellID block,const T* const array,T* const data,const int vectorDim,const int stride) {
  pargrid::CellID nbrLID = simClasses.pargrid.invalid();
  const pargrid::CellID* const nbrs = simClasses.pargrid.getCellNeighbourIDs(block);
  const int l = stride;

  // ***** ACCUMULATE TO FACE NEIGHBOURS ***** //
  // Add values to -x neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j)
      data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,k))*vectorDim+l] += array[block::arrayIndex(0,j+1,k+1)];
  // Add values to -y neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i)
      data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,k))*vectorDim+l] += array[block::arrayIndex(i+1,0,k+1)];
  // Add values to -z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
      data[(nbrLID*block::SIZE+block::index(i,j,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(i+1,j+1,0)];
  // Add values to +x neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j)
      data[(nbrLID*block::SIZE+block::index(0,j,k))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,j+1,k+1)];
  // Add values to +y neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i)
      data[(nbrLID*block::SIZE+block::index(i,0,k))*vectorDim+l] += array[block::arrayIndex(i+1,block::WIDTH_Y+1,k+1)];
  // Add values to +z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
      data[(nbrLID*block::SIZE+block::index(i,j,0))*vectorDim+l] += array[block::arrayIndex(i+1,j+1,block::WIDTH_Z+1)];

  // ***** ACCUMULATE TO EDGE NEIGHBOURS ***** //
  // Add values to -x,-y neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+0)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int k=0; k<block::WIDTH_Z; ++k)
      data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,k))*vectorDim+l] += array[block::arrayIndex(0,0,k+1)];
  // Add values to +x,-y neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+0)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int k=0; k<block::WIDTH_Z; ++k)
      data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,k))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,0,k+1)];
  // Add values to -x,+y neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+0)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int k=0; k<block::WIDTH_Z; ++k)
      data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,k))*vectorDim+l] += array[block::arrayIndex(0,block::WIDTH_Y+1,k+1)];
  // Add values to +x,+y neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+0)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int k=0; k<block::WIDTH_Z; ++k)
      data[(nbrLID*block::SIZE+block::index(0,0,k))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,k+1)];

  // Add values to -y,-z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int i=0; i<block::WIDTH_X; ++i)
      data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(i+1,0,0)];
  // Add values to +y,-z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int i=0; i<block::WIDTH_X; ++i)
      data[(nbrLID*block::SIZE+block::index(i,0,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(i+1,block::WIDTH_Y+1,0)];
  // Add values to -y,+z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int i=0; i<block::WIDTH_X; ++i)
      data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,0))*vectorDim+l] += array[block::arrayIndex(i+1,0,block::WIDTH_Z+1)];
  // Add values to +y,+z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int i=0; i<block::WIDTH_X; ++i)
      data[(nbrLID*block::SIZE+block::index(i,0,0))*vectorDim+l] += array[block::arrayIndex(i+1,block::WIDTH_Y+1,block::WIDTH_Z+1)];

  // Add values to -x,-z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int j=0; j<block::WIDTH_Y; ++j)
      data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(0,j+1,0)];
  // Add values to +x,-z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int j=0; j<block::WIDTH_Y; ++j)
      data[(nbrLID*block::SIZE+block::index(0,j,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,j+1,0)];
  // Add values to -x,+z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int j=0; j<block::WIDTH_Y; ++j) 
      data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,0))*vectorDim+l] += array[block::arrayIndex(0,j+1,block::WIDTH_Z+1)];
  // Add values to +x,+z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    for (int j=0; j<block::WIDTH_Y; ++j)
      data[(nbrLID*block::SIZE+block::index(0,j,0))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,j+1,block::WIDTH_Z+1)];

  // ***** ACCUMULATE TO CORNER NEIGHBOURS ***** //
  // Add value to -x,-y,-z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(0,0,0)];
  // Add value to +x,-y,-z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,0,0)];
  // Add value to -x,+y,-z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(0,block::WIDTH_Y+1,0)];
  // Add value to +x,+y,-z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,-1)];
  if (nbrLID != simClasses.pargrid.invalid())
    data[(nbrLID*block::SIZE+block::index(0,0,block::WIDTH_Z-1))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,0)];
  // Add value to -x,-y,+z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,0))*vectorDim+l] += array[block::arrayIndex(0,0,block::WIDTH_Z+1)];
  // Add value to +x,-y,+z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,0))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,0,block::WIDTH_Z+1)];
  // Add value to -x,+y,+z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,0))*vectorDim+l] += array[block::arrayIndex(0,block::WIDTH_Y+1,block::WIDTH_Z+1)];
  // Add value to +x,+y,+z neighbour:
  nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+1)];
  if (nbrLID != simClasses.pargrid.invalid())
    data[(nbrLID*block::SIZE+block::index(0,0,0))*vectorDim+l] += array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,block::WIDTH_Z+1)];

  // ***** ACCUMULATE TO THIS BLOCK ***** //
  for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
    data[(block*block::SIZE+block::index(i,j,k))*vectorDim+l] += array[block::arrayIndex(i+1,j+1,k+1)];
}

template<typename T,class PARTICLE> inline
void accumulateScalarCIC_XYZ(const Real* cellSizes,T* array,T value,const PARTICLE& particle) {
   const Real x0 = particle.state[0] + 0.5*cellSizes[0];
   const Real y0 = particle.state[1] + 0.5*cellSizes[1];
   const Real z0 = particle.state[2] + 0.5*cellSizes[2];
   const int i = static_cast<int>(floor(x0 / cellSizes[0]));
   const int j = static_cast<int>(floor(y0 / cellSizes[1]));
   const int k = static_cast<int>(floor(z0 / cellSizes[2]));

   const Real w_x = (x0 - i*cellSizes[0]) / cellSizes[0];
   const Real w_y = (y0 - j*cellSizes[1]) / cellSizes[1];
   const Real w_z = (z0 - k*cellSizes[2]) / cellSizes[2];

   array[block::arrayIndex(i+0,j+0,k+0)] += (1-w_x)*(1-w_y)*(1-w_z) * value;
   array[block::arrayIndex(i+1,j+0,k+0)] +=   w_x  *(1-w_y)*(1-w_z) * value;
   array[block::arrayIndex(i+0,j+1,k+0)] += (1-w_x)*  w_y  *(1-w_z) * value;
   array[block::arrayIndex(i+1,j+1,k+0)] +=   w_x  *  w_y  *(1-w_z) * value;
   array[block::arrayIndex(i+0,j+0,k+1)] += (1-w_x)*(1-w_y)*  w_z  * value;
   array[block::arrayIndex(i+1,j+0,k+1)] +=   w_x  *(1-w_y)*  w_z  * value;
   array[block::arrayIndex(i+0,j+1,k+1)] += (1-w_x)*  w_y  *  w_z  * value;
   array[block::arrayIndex(i+1,j+1,k+1)] +=   w_x  *  w_y  *  w_z  * value;
}

template<typename T,class PARTICLE> inline
void accumulateScalarNGP_XYZ(const Real* cellSizes,T* array,T value,const PARTICLE& particle) {
   const int i = 1 + static_cast<int>(floor(particle.state[0] / cellSizes[0]));
   const int j = 1 + static_cast<int>(floor(particle.state[1] / cellSizes[1]));
   const int k = 1 + static_cast<int>(floor(particle.state[2] / cellSizes[2]));
   array[block::arrayIndex(i,j,k)] += value;
}

template<typename T,class PARTICLE> inline
void accumulateScalarTSC_XYZ(const Real* cellSizes,T* array,T value,const PARTICLE& particle) {
   const Real x0 = particle.state[0];
   const Real y0 = particle.state[1];
   const Real z0 = particle.state[2];
   const int i = 1 + static_cast<int>(floor(x0 / cellSizes[0]));
   const int j = 1 + static_cast<int>(floor(y0 / cellSizes[1]));
   const int k = 1 + static_cast<int>(floor(z0 / cellSizes[2]));

   // Calculate shape factors for TSC method. The method is from 
   // Hockney & Eastwood, p. 131. Note that a constant 0.25*dx^2 is
   // used here to make shape factors positive definitive:
   Real x;
   const Real dx2 = cellSizes[0]*cellSizes[0];
   Real CONST = 0.25*dx2;
   x = x0 - (i-1.5)*cellSizes[0];
   const Real W_x_n = 0.5*(x*x - 3.0*cellSizes[0]*x + 2.0*dx2 + CONST)/dx2;
   x = x0 - (i-0.5)*cellSizes[0];
   const Real W_x_c = 1.0 - (x*x + CONST)/dx2;
   x = x0 - (i+0.5)*cellSizes[0];
   const Real W_x_p = 0.5*(x*x + 3.0*cellSizes[0]*x + 2.0*dx2 + CONST)/dx2;

   Real y;
   const Real dy2 = cellSizes[1]*cellSizes[1];
   CONST = 0.25*dy2;
   y = y0 - (j-1.5)*cellSizes[1];
   const Real W_y_n = 0.5*(y*y - 3.0*cellSizes[1]*y + 2.0*dy2 + CONST)/dy2;
   y = y0 - (j-0.5)*cellSizes[1];
   const Real W_y_c = 1.0 - (y*y + CONST)/dy2;
   y = y0 - (j+0.5)*cellSizes[1];
   const Real W_y_p = 0.5*(y*y + 3.0*cellSizes[1]*y + 2.0*dy2 + CONST)/dy2;

   Real z;
   const Real dz2 = cellSizes[2]*cellSizes[2];
   CONST = 0.25*dz2;
   z = z0 - (k-1.5)*cellSizes[2];
   const Real W_z_n = 0.5*(z*z - 3.0*cellSizes[2]*z + 2.0*dz2 + CONST)/dz2;
   z = z0 - (k-0.5)*cellSizes[2];
   const Real W_z_c = 1.0 - (z*z + CONST)/dz2;
   z = z0 - (k+0.5)*cellSizes[2];
   const Real W_z_p = 0.5*(z*z + 3.0*cellSizes[2]*z + 2.0*dz2 + CONST)/dz2;

   #ifndef NDEBUG
   bool factorsOK = true;
   const Real EPS = 1.0e-14;
   if (W_x_n < -EPS || W_x_n > 1.0) factorsOK = false;
   if (W_x_c < -EPS || W_x_c > 1.0) factorsOK = false;
   if (W_x_p < -EPS || W_x_p > 1.0) factorsOK = false;
   if (fabs(W_x_n+W_x_c+W_x_p-1.0) > EPS) factorsOK = false;
   if (W_y_n < -EPS || W_y_n > 1.0) factorsOK = false;
   if (W_y_c < -EPS || W_y_c > 1.0) factorsOK = false;
   if (W_y_p < -EPS || W_y_p > 1.0) factorsOK = false;
   if (fabs(W_y_n+W_y_c+W_y_p-1.0) > EPS) factorsOK = false;
   if (W_z_n < -EPS || W_z_n > 1.0) factorsOK = false;
   if (W_z_c < -EPS || W_z_c > 1.0) factorsOK = false;
   if (W_z_p < -EPS || W_z_p > 1.0) factorsOK = false;
   if (fabs(W_z_n+W_z_c+W_z_p-1.0) > EPS) factorsOK = false;
   if (factorsOK == false) {
      std::cerr << "(ACCUMULATOR.H) ERROR in TSC shape factors!" << std::endl;
      std::cerr << "x,y,z: " << x0 << '\t' << y0 << '\t' << z0 << std::endl;
      std::cerr << "x,y,z (relative): " << x0/cellSizes[0] << '\t' << y0/cellSizes[1] << '\t' << z0/cellSizes[2] << std::endl;
      std::cerr << "X: " << W_x_n << '\t' << W_x_c << '\t' << W_x_p << std::endl;
      std::cerr << "Y: " << W_y_n << '\t' << W_y_c << '\t' << W_y_p << std::endl;
      std::cerr << "Z: " << W_z_n << '\t' << W_z_c << '\t' << W_z_p << std::endl;
      std::cerr.precision(14); std::cerr << std::scientific;
      std::cerr << "SUMS: " << W_x_n+W_x_c+W_x_p-1.0 << '\t' << W_y_n+W_y_c+W_y_p-1.0 << '\t' << W_z_n+W_z_c+W_z_p-1.0 << std::endl;
      exit(1);
   }
   #endif

   array[block::arrayIndex(i-1,j-1,k-1)] += W_x_n*W_y_n*W_z_n * value;
   array[block::arrayIndex(i  ,j-1,k-1)] += W_x_c*W_y_n*W_z_n * value;
   array[block::arrayIndex(i+1,j-1,k-1)] += W_x_p*W_y_n*W_z_n * value;
   array[block::arrayIndex(i-1,j  ,k-1)] += W_x_n*W_y_c*W_z_n * value;
   array[block::arrayIndex(i  ,j  ,k-1)] += W_x_c*W_y_c*W_z_n * value;
   array[block::arrayIndex(i+1,j  ,k-1)] += W_x_p*W_y_c*W_z_n * value;
   array[block::arrayIndex(i-1,j+1,k-1)] += W_x_n*W_y_p*W_z_n * value;
   array[block::arrayIndex(i  ,j+1,k-1)] += W_x_c*W_y_p*W_z_n * value;
   array[block::arrayIndex(i+1,j+1,k-1)] += W_x_p*W_y_p*W_z_n * value;

   array[block::arrayIndex(i-1,j-1,k  )] += W_x_n*W_y_n*W_z_c * value;
   array[block::arrayIndex(i  ,j-1,k  )] += W_x_c*W_y_n*W_z_c * value;
   array[block::arrayIndex(i+1,j-1,k  )] += W_x_p*W_y_n*W_z_c * value;
   array[block::arrayIndex(i-1,j  ,k  )] += W_x_n*W_y_c*W_z_c * value;
   array[block::arrayIndex(i  ,j  ,k  )] += W_x_c*W_y_c*W_z_c * value;
   array[block::arrayIndex(i+1,j  ,k  )] += W_x_p*W_y_c*W_z_c * value;
   array[block::arrayIndex(i-1,j+1,k  )] += W_x_n*W_y_p*W_z_c * value;
   array[block::arrayIndex(i  ,j+1,k  )] += W_x_c*W_y_p*W_z_c * value;
   array[block::arrayIndex(i+1,j+1,k  )] += W_x_p*W_y_p*W_z_c * value;

   array[block::arrayIndex(i-1,j-1,k+1)] += W_x_n*W_y_n*W_z_p * value;
   array[block::arrayIndex(i  ,j-1,k+1)] += W_x_c*W_y_n*W_z_p * value;
   array[block::arrayIndex(i+1,j-1,k+1)] += W_x_p*W_y_n*W_z_p * value;
   array[block::arrayIndex(i-1,j  ,k+1)] += W_x_n*W_y_c*W_z_p * value;
   array[block::arrayIndex(i  ,j  ,k+1)] += W_x_c*W_y_c*W_z_p * value;
   array[block::arrayIndex(i+1,j  ,k+1)] += W_x_p*W_y_c*W_z_p * value;
   array[block::arrayIndex(i-1,j+1,k+1)] += W_x_n*W_y_p*W_z_p * value;
   array[block::arrayIndex(i  ,j+1,k+1)] += W_x_c*W_y_p*W_z_p * value;
   array[block::arrayIndex(i+1,j+1,k+1)] += W_x_p*W_y_p*W_z_p * value;
}

template<typename REAL>
void evaluateShapeFactorsNGP_3D(const REAL* const sizes,const REAL* const pos,REAL* shapeFactors,int i,int j,int k) {
   shapeFactors[sfactor::W_X_N] = 0.0;
   shapeFactors[sfactor::W_X_C] = 1.0;
   shapeFactors[sfactor::W_X_P] = 0.0;
   shapeFactors[sfactor::W_Y_N] = 0.0;
   shapeFactors[sfactor::W_Y_C] = 1.0;
   shapeFactors[sfactor::W_Y_P] = 0.0;
   shapeFactors[sfactor::W_Z_N] = 0.0;
   shapeFactors[sfactor::W_Z_C] = 1.0;
   shapeFactors[sfactor::W_Z_P] = 0.0;
}

template<typename REAL>
void evaluateShapeFactorsTSC_3D(const REAL* const sizes,const REAL* const pos,REAL* shapeFactors,int i,int j,int k) {
   Real x;
   const Real dx2 = sizes[0]*sizes[0];
   Real CONST = 0.25*dx2;
   x = pos[0] - (i-1.5)*sizes[0];   
   shapeFactors[sfactor::W_X_N] = 0.5*(x*x - 3.0*sizes[0]*x + 2.0*dx2 + CONST)/dx2;
   x = pos[0] - (i-0.5)*sizes[0];   
   shapeFactors[sfactor::W_X_C] = 1.0 - (x*x + CONST)/dx2;
   x = pos[0] - (i+0.5)*sizes[0];   
   shapeFactors[sfactor::W_X_P] = 0.5*(x*x + 3.0*sizes[0]*x + 2.0*dx2 + CONST)/dx2;

   Real y;
   const Real dy2 = sizes[1]*sizes[1];
   CONST = 0.25*dy2;
   y = pos[1] - (j-1.5)*sizes[1];
   shapeFactors[sfactor::W_Y_N] = 0.5*(y*y - 3.0*sizes[1]*y + 2.0*dy2 + CONST)/dy2;
   y = pos[1] - (j-0.5)*sizes[1];
   shapeFactors[sfactor::W_Y_C] = 1.0 - (y*y + CONST)/dy2;
   y = pos[1] - (j+0.5)*sizes[1];
   shapeFactors[sfactor::W_Y_P] = 0.5*(y*y + 3.0*sizes[1]*y + 2.0*dy2 + CONST)/dy2;

   Real z;
   const Real dz2 = sizes[2]*sizes[2];
   CONST = 0.25*dz2;
   z = pos[2] - (k-1.5)*sizes[2];
   shapeFactors[sfactor::W_Z_N] = 0.5*(z*z - 3.0*sizes[2]*z + 2.0*dz2 + CONST)/dz2;
   z = pos[2] - (k-0.5)*sizes[2];
   shapeFactors[sfactor::W_Z_C] = 1.0 - (z*z + CONST)/dz2;
   z = pos[2] - (k+0.5)*sizes[2];
   shapeFactors[sfactor::W_Z_P] = 0.5*(z*z + 3.0*sizes[2]*z + 2.0*dz2 + CONST)/dz2;

   #ifndef NDEBUG
   bool factorsOK = true;
   const Real EPS = -1e-14;
   if (shapeFactors[sfactor::W_X_N] < EPS || shapeFactors[sfactor::W_X_N] > 1.0) factorsOK = false;
   if (shapeFactors[sfactor::W_X_C] < EPS || shapeFactors[sfactor::W_X_C] > 1.0) factorsOK = false;
   if (shapeFactors[sfactor::W_X_P] < EPS || shapeFactors[sfactor::W_X_P] > 1.0) factorsOK = false;
   if (shapeFactors[sfactor::W_Y_N] < EPS || shapeFactors[sfactor::W_Y_N] > 1.0) factorsOK = false;
   if (shapeFactors[sfactor::W_Y_C] < EPS || shapeFactors[sfactor::W_Y_C] > 1.0) factorsOK = false;
   if (shapeFactors[sfactor::W_Y_P] < EPS || shapeFactors[sfactor::W_Y_P] > 1.0) factorsOK = false;
   if (shapeFactors[sfactor::W_Z_N] < EPS || shapeFactors[sfactor::W_Z_N] > 1.0) factorsOK = false;
   if (shapeFactors[sfactor::W_Z_C] < EPS || shapeFactors[sfactor::W_Z_C] > 1.0) factorsOK = false;
   if (shapeFactors[sfactor::W_Z_P] < EPS || shapeFactors[sfactor::W_Z_P] > 1.0) factorsOK = false;
   if (factorsOK == false) {
      std::cerr << "(ACCUMULATORS.H) ERROR in TSC evaluateShapeFactors!" << std::endl;
      std::cerr << "\t pos= " << pos[0] << '\t' << pos[1] << '\t' << pos[2] << std::endl;
      std::cerr << "\t W_X: " << shapeFactors[sfactor::W_X_N] << '\t' << shapeFactors[sfactor::W_X_C] << '\t' << shapeFactors[sfactor::W_X_P] << std::endl;
      std::cerr << "\t W_X: " << shapeFactors[sfactor::W_Y_N] << '\t' << shapeFactors[sfactor::W_Y_C] << '\t' << shapeFactors[sfactor::W_Y_P] << std::endl;
      std::cerr << "\t W_X: " << shapeFactors[sfactor::W_Z_N] << '\t' << shapeFactors[sfactor::W_Z_C] << '\t' << shapeFactors[sfactor::W_Z_P] << std::endl;
      exit(1);
   }
   #endif
}

template<typename T> inline
void fetchValues3D(SimulationClasses& simClasses,pargrid::CellID blockLID,T* RESTRICT const dest,const T* RESTRICT const source,int vectorsize) {
   pargrid::CellID nbrLID = simClasses.pargrid.invalidCellID();
   const pargrid::CellID* const nbrs = simClasses.pargrid.getCellNeighbourIDs(blockLID);
   const T* srcTmp = NULL;

   // ***** FETCH DATA FROM THIS BLOCK ***** //
   srcTmp = source + blockLID*block::SIZE*vectorsize;
   for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
     dest[block::arrayIndex(i+1,j+1,k+1)*vectorsize+v] = srcTmp[block::index(i,j,k)*vectorsize+v];

   // ***** LOAD VALUES FROM 6 FACE NEIGHBOURS ***** //
   // Copy from -x neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(0,j+1,k+1)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,j,k)*vectorsize+v];
   }
   // Copy from +x neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(block::WIDTH_X+1,j+1,k+1)*vectorsize+v] = srcTmp[block::index(0,j,k)*vectorsize+v];
   }
   // Copy from -y neighbour
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(i+1,0,k+1)*vectorsize+v] = srcTmp[block::index(i,block::WIDTH_Y-1,k)*vectorsize+v];
   }
   // Copy from +y neighbour
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(i+1,block::WIDTH_Y+1,k+1)*vectorsize+v] = srcTmp[block::index(i,0,k)*vectorsize+v];
   }
   // Copy from -z neighbour
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(i+1,j+1,0)*vectorsize+v] = srcTmp[block::index(i,j,block::WIDTH_Z-1)*vectorsize+v];
   }
   // Copy from +z neighbour
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(i+1,j+1,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(i,j,0)*vectorsize+v];
   }

   // ***** LOAD VALUES FROM 12 EDGE NEIGHBOURS ***** //
   // -x,-y neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+0)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(0,0,k+1)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,k)*vectorsize+v];
   }
   // +x,-y neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+0)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(block::WIDTH_X+1,0,k+1)*vectorsize+v] = srcTmp[block::index(0,block::WIDTH_Y-1,k)*vectorsize+v];
   }
   // -x,+y neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+0)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(0,block::WIDTH_Y+1,k+1)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,0,k)*vectorsize+v];
   }
   // +x,+y neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+0)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,k+1)*vectorsize+v] = srcTmp[block::index(0,0,k)*vectorsize+v];
   }
   // -x,-z neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int j=0; j<block::WIDTH_Y; ++j) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(0,j+1,0)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,j,block::WIDTH_Z-1)*vectorsize+v];
   }
   // +x,-z neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int j=0; j<block::WIDTH_Y; ++j) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(block::WIDTH_X+1,j+1,0)*vectorsize+v] = srcTmp[block::index(0,j,block::WIDTH_Z-1)*vectorsize+v];
   }
   // -x,+z neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int j=0; j<block::WIDTH_Y; ++j) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(0,j+1,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,j,0)*vectorsize+v];
   }
   // +x,+z neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int j=0; j<block::WIDTH_Y; ++j) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(block::WIDTH_X+1,j+1,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(0,j,0)*vectorsize+v];
   }
   // -y,-z neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(i+1,0,0)*vectorsize+v] = srcTmp[block::index(i,block::WIDTH_Y-1,block::WIDTH_Z-1)*vectorsize+v];
   }
   // +y,-z neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(i+1,block::WIDTH_Y+1,0)*vectorsize+v] = srcTmp[block::index(i,0,block::WIDTH_Z-1)*vectorsize+v];
   }
   // -y,+z neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(i+1,0,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(i,block::WIDTH_Y-1,0)*vectorsize+v];
   }
   // +y,+z neighbour:
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int i=0; i<block::WIDTH_X; ++i) for (int v=0; v<vectorsize; ++v)
	dest[block::arrayIndex(i+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(i,0,0)*vectorsize+v];
   }

   // ***** LOAD VALUES FROM 8 CORNER NEIGHBOURS ***** //
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int v=0; v<vectorsize; ++v)
      dest[block::arrayIndex(0,0,0)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,block::WIDTH_Z-1)*vectorsize+v];
   }
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int v=0; v<vectorsize; ++v)
      dest[block::arrayIndex(block::WIDTH_X+1,0,0)*vectorsize+v] = srcTmp[block::index(0,block::WIDTH_Y-1,block::WIDTH_Z-1)*vectorsize+v];
   }
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int v=0; v<vectorsize; ++v)
      dest[block::arrayIndex(0,block::WIDTH_Y+1,0)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,0,block::WIDTH_Z-1)*vectorsize+v];
   }
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,-1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int v=0; v<vectorsize; ++v)
      dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,0)*vectorsize+v] = srcTmp[block::index(0,0,block::WIDTH_Z-1)*vectorsize+v];
   }
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int v=0; v<vectorsize; ++v)
      dest[block::arrayIndex(0,0,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,0)*vectorsize+v];
   }
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int v=0; v<vectorsize; ++v)
      dest[block::arrayIndex(block::WIDTH_X+1,0,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(0,block::WIDTH_Y-1,0)*vectorsize+v];
   }
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int v=0; v<vectorsize; ++v)
      dest[block::arrayIndex(0,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(block::WIDTH_X-1,0,0)*vectorsize+v];
   }
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+1)];
   if (nbrLID != simClasses.pargrid.invalidCellID()) {
      srcTmp = source + nbrLID*block::SIZE*vectorsize;
      for (int v=0; v<vectorsize; ++v)
      dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorsize+v] = srcTmp[block::index(0,0,0)*vectorsize+v];
   }
}

#endif
