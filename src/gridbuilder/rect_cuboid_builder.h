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

#ifndef RECT_CUBOID_BUILDER_H
#define RECT_CUBOID_BUILDER_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <gridbuilder.h>

class RectCuboidBuilder: public GridBuilder {
 public:
   RectCuboidBuilder();
   ~RectCuboidBuilder();

   bool finalize();
   bool getCellIndices(unsigned int N_cells,const builder::ID* const cellIDs,builder::ID*& indices);
   bool getCellInfo(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,
		    builder::ID* cellIDs,unsigned char* N_neighbours);
   bool getCellNeighbours(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,
			  unsigned char* N_neighbours,builder::ID* neighbourIDs,unsigned char* neighbourTypeIDs);
   bool getInitialState(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,double* coordinates);
   bool getMeshBoundingBox(Real& x_min,Real& y_min,Real& z_min,Real& dx0,Real& dy0,Real& dz0);
   bool getNumberOfCells(builder::ID& N_cells);
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);

 private:
   bool initialized;
   bool periodic_x;                            /**< If true, x-coordinate is periodic.*/
   bool periodic_y;                            /**< If true, y-coordinate is periodic.*/
   bool periodic_z;                            /**< If true, z-coordinate is periodic.*/
   builder::ID size_x;                         /**< Number of cells in x-direction in initial grid.*/
   builder::ID size_y;                         /**< Number of cells in y-direction in initial grid.*/
   builder::ID size_z;                         /**< Number of cells in z-direction in initial grid.*/
   Real x_min;                                 /**< Minimum value for x-coordinate in initial grid.*/
   Real y_min;                                 /**< Minimum value for y-coordinate in initial grid.*/
   Real z_min;                                 /**< Minimum value for z-coordinate in initial grid.*/
   Real x_max;                                 /**< Maximum value for x-coordinate in initial grid.*/
   Real y_max;                                 /**< Maximum value for y-coordinate in initial grid.*/
   Real z_max;                                 /**< Maximum value for z-coordinate in initial grid.*/
   Real dx;                                    /**< Cell size in x-direction.*/
   Real dy;                                    /**< Cell size in y-direction.*/
   Real dz;                                    /**< Cell size in z-direction.*/

   void calculateCellIndices(builder::ID cellID,builder::ID& i,builder::ID& j,builder::ID& k);
   builder::ID calculateNeighbourID(builder::ID I,builder::ID J,builder::ID K,int i,int j,int k);
   unsigned char calculateNeighbourTypeID(int i,int j,int k);
   unsigned char countNeighbours(builder::ID i,builder::ID j,builder::ID k);
};

GridBuilder* RCBCreator();

#endif
