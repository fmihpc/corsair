/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2013 Finnish Meteorological Institute
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

#ifndef LOGICALLY_CARTESIAN_BUILDER_H
#define LOGICALLY_CARTESIAN_BUILDER_H

#include <simulation.h>
#include <simulationclasses.h>
#include <gridbuilder.h>

class LogicallyCartesianBuilder: public GridBuilder {
 public:
   LogicallyCartesianBuilder();
   virtual ~LogicallyCartesianBuilder();

   virtual bool defineParameters(ConfigReader& cr);

   virtual bool finalize();
   virtual bool getCellIndices(unsigned int N_cells,const builder::ID* const globalIDs,builder::ID*& indices);
   virtual bool getCellInfo(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,unsigned char* N_neighbours);
   virtual bool getCellNeighbours(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,
				  unsigned char* N_neighbours,builder::ID* neighbourIDs,unsigned char* neighbourTypeIDs);
   virtual bool getMeshBoundingBox(Real& x_min,Real& y_min,Real& z_min,Real& dx0,Real& dy0,Real& dz0);
   virtual bool getInitialState(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,double* coordinates);
   virtual bool getNumberOfCells(builder::ID& N_cells);
   virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);

   virtual bool writeMesh(Simulation& sim,SimulationClasses& simClasses,const std::string& meshName,bool meshWrite);

 protected:
   bool initialized;       /**< If true, LogicallyCartesianBuilder initialized correctly.*/
   const std::string prefix;

   bool calcCellVolumes;   /**< If true, calculate cell volumes and insert them to ParGrid.*/

   enum Geometry {
      Cartesian,
      Cylindrical,
      Spherical
   };

   Geometry geometry;

   std::string name;
   builder::ID xCells;     /**< Number of cells in x-direction.*/
   builder::ID yCells;     /**< Number of cells in y-direction.*/
   builder::ID zCells;     /**< Number of cells in z-direction.*/
   std::string xLabel;     /**< Name of x-coordinate, default is "x-coordinate"..*/
   std::string yLabel;     /**< Name of y-coordinate, default is "y-coordinate".*/
   std::string zLabel;     /**< Name of z-coordinate, default is "z-coordinate"..*/
   Real x_max;             /**< Maximum coordinate value of mesh in x-direction.*/
   Real y_max;             /**< Maximum coordinate value of mesh in y-direction.*/
   Real z_max;             /**< Maximum coordinate value of mesh in z-direction.*/
   Real x_min;             /**< Minimum coordinate value of mesh in x-direction.*/
   Real y_min;             /**< Minimum coordinate value of mesh in y-direction.*/
   Real z_min;             /**< Minimum coordinate value of mesh in z-direction.*/
   bool xPeriodic;         /**< If true, mesh is periodic in x-direction.*/
   bool yPeriodic;         /**< If true, mesh is periodic in y-direction.*/
   bool zPeriodic;         /**< If true, mesh is periodic in z-direction.*/
   bool xUniform;          /**< If true, cell size in x-direction is uniform, i.e. dx_cell is constant.*/
   bool yUniform;          /**< If true, cell size in y-direction is uniform, i.e. dy_cell is constant.*/
   bool zUniform;          /**< If true, cell size in z-direction is uniform, i.e. dz_cell is constant.*/
   std::string xUnits;     /**< Unit of measure for x-coordinate, defaults to "".*/
   std::string yUnits;     /**< Unit of measure for y-coordinate, defaults to "".*/
   std::string zUnits;     /**< Unit of measure for z-coordinate, defaults to "".*/

   Simulation* sim;

   void calculateCellIndices(builder::ID cellID,builder::ID& i,builder::ID& j,builder::ID& k);
   virtual bool calculateBlockSizes();
   virtual bool calculateCellSizes();
   builder::ID calculateNeighbourID(builder::ID I,builder::ID J,builder::ID K,int i,int j,int k);
   unsigned char calculateNeighbourTypeID(int i,int j,int k);
   unsigned char countNeighbours(builder::ID i,builder::ID j,builder::ID k);
   virtual bool writeInitStatus(bool success,SimulationClasses& simClasses,const std::string& reason);
};

GridBuilder* LCCreator();

#endif
