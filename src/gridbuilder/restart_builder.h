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

#ifndef RESTART_BUILDER_H
#define RESTART_BUILDER_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <gridbuilder.h>

class RestartBuilder: public GridBuilder {
 public:
   RestartBuilder();
   ~RestartBuilder();

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
   bool initialized;                /**< If true, RestartBuilder has initialized correctly.*/
   std::string meshName;            /**< Name of the mesh RestartBuilder should re-create.*/
   Simulation* sim;                 /**< Pointer to struct containing generic simulation variables.*/
   SimulationClasses* simClasses;   /**< Pointer to struct containing generic simulation classes.*/
};

GridBuilder* RestartCreator();

#endif
