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

#ifndef GRIDBUILDER_H
#define GRIDBUILDER_H

#include <mpi.h>
#include <vector>
#include <stdint.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

namespace builder {
   typedef uint32_t ID;
}

class GridBuilder {
 public:
   GridBuilder();
   virtual ~GridBuilder();

   bool build(Simulation& sim,SimulationClasses& simClasses);
   virtual bool writeMesh(Simulation& sim,SimulationClasses& simClasses,const std::string& name,bool write);

   // ********************************************************* //
   // ***** DECRALATIONS OF PURE VIRTUAL MEMBER FUNCTIONS ***** //
   // *****   IMPLEMENTING CLASSES MUST DEFINE ALL THESE  ***** //
   // ********************************************************* //

   virtual bool finalize() = 0;
   virtual bool getCellIndices(unsigned int N_cells,const builder::ID* const cellIDs,builder::ID*& indices) = 0;
   virtual bool getCellInfo(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,unsigned char* N_neighbours) = 0;
   virtual bool getCellNeighbours(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,
				  unsigned char* N_neighbours,builder::ID* neighbourIDs,unsigned char* neighbourTypeIDs) = 0;
   virtual bool getMeshBoundingBox(Real& x_min,Real& y_min,Real& z_min,Real& dx0,Real& dy0,Real& dz0) = 0;
   virtual bool getInitialState(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,double* coordinates) = 0;
   virtual bool getNumberOfCells(builder::ID& N_cells) = 0;
   virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) = 0;

 protected:

};

#endif
