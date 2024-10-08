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
#include <limits>
#include <cmath>

#include "rect_cuboid_builder.h"
#include <mpilogger.h>

using namespace std;

static const builder::ID MAX_INDEX = numeric_limits<builder::ID>::max() - 3;

RectCuboidBuilder::RectCuboidBuilder(): GridBuilder() {
   initialized = false;
   periodic_x = false;
   periodic_y = false;
   periodic_z = false;
}

RectCuboidBuilder::~RectCuboidBuilder() {
   
}

void RectCuboidBuilder::calculateCellIndices(builder::ID cellID,builder::ID& i,builder::ID& j,builder::ID& k) {
   builder::ID tmp = cellID;
   k = tmp / (size_y*size_x);
   tmp -= k*size_y*size_x;
   j = tmp / size_x;
   tmp -= j*size_x;
   i = tmp;
}

builder::ID RectCuboidBuilder::calculateNeighbourID(builder::ID I,builder::ID J,builder::ID K,int i,int j,int k) {
   // Check that the given neighbour is within the simulation box:
   builder::ID i_out = I + i;
   builder::ID j_out = J + j;
   builder::ID k_out = K + k;
   // Check that neighbour i-index is within the simulation volume:
   if (i_out > size_x-1) {
      if (periodic_x == false) {
	 return numeric_limits<builder::ID>::max();
      } else {
	 if (size_x == 1) i_out = 0;
	 else if (i_out > MAX_INDEX) i_out = size_x-1 - (numeric_limits<builder::ID>::max()-i_out);
	 else i_out -= size_x;
      }
   }
   // Check that neighbour j-index is within the simulation volume:
   if (j_out > size_y-1) {
      if (periodic_y == false) {
	 return numeric_limits<builder::ID>::max();
      } else {
	 if (size_y == 1) j_out = 0;
	 else if (j_out > MAX_INDEX) j_out = size_y-1 - (numeric_limits<builder::ID>::max()-j_out);
	 else j_out -= size_y;
      }
   }
   // Check that neighbour k-index is within the simulation volume:
   if (k_out > size_z-1) {
      if (periodic_z == false) {
	 return numeric_limits<builder::ID>::max();
      } else {
	 if (size_z == 1) k_out = 0;
	 else if (k_out > MAX_INDEX) k_out = size_z-1 - (numeric_limits<builder::ID>::max()-k_out);
	 else k_out -= size_z;
      }
   }
   return k_out*size_y*size_x + j_out*size_x + i_out;
}

unsigned char RectCuboidBuilder::calculateNeighbourTypeID(int i,int j,int k) {
   ++i; ++j; ++k;
   return k*9+j*3+i;
}

unsigned char RectCuboidBuilder::countNeighbours(builder::ID i,builder::ID j,builder::ID k) {
   unsigned char N_neighbours = 0;
   for (int ii=-1; ii<2; ++ii) for (int jj=-1; jj<2; ++jj) for (int kk=-1; kk<2; ++kk) {
      if (ii == 0 && (jj == 0 && kk == 0)) continue; // cell is not its own neighbour
      if (i + ii > size_x-1 && periodic_x == false) continue;
      if (j + jj > size_y-1 && periodic_y == false) continue;
      if (k + kk > size_z-1 && periodic_z == false) continue;
      ++N_neighbours;
   }
   return N_neighbours;
}

bool RectCuboidBuilder::finalize() {
   bool success = true;
 
   return success;
}

bool RectCuboidBuilder::getCellIndices(unsigned int N_cells,const builder::ID* const cellIDs,builder::ID*& indices) {
   if (initialized == false) return false;
   bool success = true;

   builder::ID i,j,k;
   indices = new builder::ID[N_cells*3];

   // Calculate (i,j,k) indices for all given cell global IDs and copy them to indices:
   for (builder::ID c=0; c<N_cells; ++c) {
      const builder::ID cellID = cellIDs[c];
      calculateCellIndices(cellID,i,j,k);
      
      indices[c*3+0] = i;
      indices[c*3+1] = j;
      indices[c*3+2] = k;
   }
   return success;
}

bool RectCuboidBuilder::getCellInfo(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,
				    builder::ID* cellIDs,unsigned char* N_neighbours) {
   if (initialized == false) return false;
   bool success = true;
   
   builder::ID i,j,k;
   const builder::ID N_cells = cellOffsetEnd-cellOffsetStart;
   for (builder::ID c=0; c<N_cells; ++c) {
      builder::ID cellID = cellOffsetStart + c;
      calculateCellIndices(cellID,i,j,k);
      cellIDs[c]      = cellID;
      N_neighbours[c] = countNeighbours(i,j,k);
   }
   return success;
}

bool RectCuboidBuilder::getCellNeighbours(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,
					  unsigned char* N_neighbours,builder::ID* neighbourIDs,unsigned char* neighbourTypeIDs) {
   if (initialized == false) return false;
   bool success = true;
   
   builder::ID counter = 0;
   builder::ID I,J,K;
   const builder::ID N_cells = cellOffsetEnd-cellOffsetStart;
   for (builder::ID c=0; c<N_cells; ++c) {
      const builder::ID cellID = cellIDs[c];
      calculateCellIndices(cellID,I,J,K);
      
      // Add neighbour IDs that exist in a 3x3x3 cube centered at this cell:
      for (int i=-1; i<2; ++i) for (int j=-1; j<2; ++j) for (int k=-1; k<2; ++k) {
	 if (i == 0 && (j == 0 && k == 0)) continue;
	 builder::ID nbrID     = calculateNeighbourID(I,J,K,i,j,k);
	 if (nbrID == numeric_limits<builder::ID>::max()) continue;
	 
	 const unsigned char nbrTypeID = calculateNeighbourTypeID(i,j,k);	 
	 neighbourIDs[counter]     = nbrID;
	 neighbourTypeIDs[counter] = nbrTypeID;
	 ++counter;
      }
   }
   //cerr << "RectCuboid counter: " << counter << endl;
   return success;
}

bool RectCuboidBuilder::getInitialState(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,double* coordinates) {
   if (initialized == false) return false;
   bool success = true;
   
   builder::ID i,j,k;
   const builder::ID N_cells = cellOffsetEnd-cellOffsetStart;
   for (builder::ID cell=0; cell<N_cells; ++cell) {
      builder::ID cellID = cellIDs[cell];
      calculateCellIndices(cellID,i,j,k);
      
      coordinates[3*cell+0] = x_min + i*dx;
      coordinates[3*cell+1] = y_min + j*dy;
      coordinates[3*cell+2] = z_min + k*dz;
   }
   return success;
}

bool RectCuboidBuilder::getMeshBoundingBox(Real& x_min,Real& y_min,Real& z_min,Real& dx0,Real& dy0,Real& dz0) {
   if (initialized == false) return false;
   bool success = true;
   x_min = this->x_min;
   y_min = this->y_min;
   z_min = this->z_min;
   dx0   = this->dx;
   dy0   = this->dy;
   dz0   = this->dz;
   return success;
}

bool RectCuboidBuilder::getNumberOfCells(builder::ID& N_cells) {
   if (initialized == false) return false;
   bool success = true;
   N_cells = size_x*size_y*size_z;
   return success;
}

bool RectCuboidBuilder::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
   bool success = true;
   simClasses.logger << "(RECTCUBOID) Starting to initialize" << endl;

   Real NaN = NAN;
   string xPeriodic,yPeriodic,zPeriodic;
   builder::ID MAX = numeric_limits<builder::ID>::max();
   cr.add("RectCuboid.x_min","Minimum value for x-coordinate (float).",NaN);
   cr.add("RectCuboid.y_min","Minimum value for y-coordinate (float).",NaN);
   cr.add("RectCuboid.z_min","Minimum value for z-coordinate (float).",NaN);
   cr.add("RectCuboid.x_max","Minimum value for x-coordinate (float).",NaN);
   cr.add("RectCuboid.y_max","Minimum value for y-coordinate (float).",NaN);
   cr.add("RectCuboid.z_max","Minimum value for z-coordinate (float).",NaN);
   cr.add("RectCuboid.x_size","Number of blocks in x-direction (int).",MAX);
   cr.add("RectCuboid.y_size","Number of blocks in y-direction (int).",MAX);
   cr.add("RectCuboid.z_size","Number of blocks in z-direction (int).",MAX);
   cr.add("RectCuboid.x_periodic","Is grid periodic in x-direction (yes/no) ?","no");
   cr.add("RectCuboid.y_periodic","Is grid periodic in y-direction (yes/no) ?","no");
   cr.add("RectCuboid.z_periodic","Is grid periodic in z-direction (yes/no) ?","no");
   cr.parse();
   cr.get("RectCuboid.x_min",x_min);
   cr.get("RectCuboid.y_min",y_min);
   cr.get("RectCuboid.z_min",z_min);
   cr.get("RectCuboid.x_max",x_max);
   cr.get("RectCuboid.y_max",y_max);
   cr.get("RectCuboid.z_max",z_max);
   cr.get("RectCuboid.x_size",size_x);
   cr.get("RectCuboid.y_size",size_y);
   cr.get("RectCuboid.z_size",size_z);
   cr.get("RectCuboid.x_periodic",xPeriodic);
   cr.get("RectCuboid.y_periodic",yPeriodic);
   cr.get("RectCuboid.z_periodic",zPeriodic);
   
   if (x_min != x_min) success = false;
   if (y_min != y_min) success = false;
   if (z_min != z_min) success = false;
   if (x_max != x_max) success = false;
   if (y_max != y_max) success = false;
   if (z_max != z_max) success = false;
   if (size_x == MAX) success = false;
   if (size_y == MAX) success = false;
   if (size_z == MAX) success = false;
   if (xPeriodic == "yes") periodic_x = true;
   if (yPeriodic == "yes") periodic_y = true;
   if (zPeriodic == "yes") periodic_z = true;

   dx = (x_max - x_min) / size_x;
   dy = (y_max - y_min) / size_y;
   dz = (z_max - z_min) / size_z;

   // Calculate mesh block sizes:
   sim.dx_block = new Real[1];
   sim.dy_block = new Real[1];
   sim.dz_block = new Real[1];
   sim.dx_block[0] = dx;
   sim.dy_block[0] = dy;
   sim.dz_block[0] = dz;
   
   // Calculate cell sizes:
   sim.dx_cell = new Real[1];
   sim.dy_cell = new Real[1];
   sim.dz_cell = new Real[1];
   sim.dx_cell[0] = sim.dx_block[0] / block::WIDTH_X;
   sim.dy_cell[0] = sim.dy_block[0] / block::WIDTH_Y;
   sim.dz_cell[0] = sim.dz_block[0] / block::WIDTH_Z;

   sim.meshGeometry = vlsv::geometry::CARTESIAN;
   
   sim.x_periodic = periodic_x;
   sim.y_periodic = periodic_y;
   sim.z_periodic = periodic_z;

   sim.x_min = x_min;
   sim.x_max = x_max;
   sim.y_min = y_min;
   sim.y_max = y_max;
   sim.z_min = z_min;
   sim.z_max = z_max;
   
   if (success == true) simClasses.logger << "\t Initialization successful." << endl << write;
   else simClasses.logger << "\t Initialization failed." << endl << write;
   initialized = success;
   return success;
}


// Register to ObjectFactoryGeneric:
GridBuilder* RCBCreator() {return new RectCuboidBuilder();}

