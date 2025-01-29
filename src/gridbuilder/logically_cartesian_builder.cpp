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

#include <cstdlib>
#include <iostream>

#include "logically_cartesian_builder.h"

using namespace std;

static const builder::ID MAX_INDEX = numeric_limits<builder::ID>::max() - 3;

// Function that creates new LogicallyCartesianBuilder:
GridBuilder* LCCreator() {return new LogicallyCartesianBuilder();}

LogicallyCartesianBuilder::LogicallyCartesianBuilder(): prefix("LogicallyCartesian") {
   xUniform = true;
   yUniform = true;
   zUniform = true;
   xPeriodic = false;
   yPeriodic = false;
   zPeriodic = false;
   initialized = false;
}

LogicallyCartesianBuilder::~LogicallyCartesianBuilder() { }

void LogicallyCartesianBuilder::calculateCellIndices(builder::ID cellID,builder::ID& i,builder::ID& j,builder::ID& k) {
   builder::ID tmp = cellID;
   k = tmp / (sim->y_blocks*sim->x_blocks);
   tmp -= k*sim->y_blocks*sim->x_blocks;
   j = tmp / sim->x_blocks;
   tmp -= j*sim->x_blocks;
   i = tmp;
}

bool LogicallyCartesianBuilder::calculateBlockSizes() {
   if (sim->dx_cell == NULL) return false;
   if (sim->dy_cell == NULL) return false;
   if (sim->dz_cell == NULL) return false;

   sim->dx_block = new Real[sim->x_blocks];
   for (unsigned int i_block=0; i_block<sim->x_blocks; ++i_block) {
      sim->dx_block[i_block] = 0.0;
      for (int i_cell=0; i_cell<block::WIDTH_X; ++i_cell) {
	 sim->dx_block[i_block] += sim->dx_cell[i_block*block::WIDTH_X+i_cell];
      }
   }

   sim->dy_block = new Real[sim->y_blocks];
   for (unsigned int j_block=0; j_block<sim->y_blocks; ++j_block) {
      sim->dy_block[j_block] = 0.0;
      for (int j_cell=0; j_cell<block::WIDTH_Y; ++j_cell) {
	 sim->dy_block[j_block] += sim->dy_cell[j_block*block::WIDTH_Y+j_cell];
      }
   }

   sim->dz_block = new Real[sim->z_blocks];
   for (unsigned int k_block=0; k_block<sim->z_blocks; ++k_block) {
      sim->dz_block[k_block] = 0.0;
      for (int k_cell=0; k_cell<block::WIDTH_Z; ++k_cell) {
	 sim->dz_block[k_block] += sim->dz_cell[k_block*block::WIDTH_Z+k_cell];
      }
   }

   return true;
}

bool LogicallyCartesianBuilder::calculateCellSizes() {
   // Allocate memory for cell size arrays:
   sim->dx_cell = new Real[xCells];
   if (xUniform == true) {
      for (builder::ID i=0; i<xCells; ++i) sim->dx_cell[i] = (x_max-x_min) / xCells;
   } else {
      if (xCells == 1) {
	 sim->dx_cell[0] = x_max-x_min;
      } else {
	 const Real dx0 = (x_max-x_min) / xCells / 3.0;
	 const Real dx = 2*((x_max-x_min) - xCells*dx0) / (xCells*(xCells-1));
	 for (builder::ID i=0; i<xCells; ++i) sim->dx_cell[i] = dx0 + i*dx;
      }
   }

   sim->dy_cell = new Real[yCells];
   if (yUniform == true) {
      for (builder::ID i=0; i<yCells; ++i) sim->dy_cell[i] = (y_max-y_min) / yCells;
   } else {
      if (yCells == 1) {
	 sim->dy_cell[0] = y_max-y_min;
      } else {
	 const Real dy0 = (y_max-y_min) / yCells / 3.0;
	 const Real dy = 2*((y_max-y_min) - yCells*dy0) / (yCells*(yCells-1));
	 for (builder::ID i=0; i<yCells; ++i) sim->dy_cell[i] = dy0 + i*dy;
      }
   }

   sim->dz_cell = new Real[zCells];
   if (zUniform == true) {
      for (builder::ID i=0; i<zCells; ++i) sim->dz_cell[i] = (z_max-z_min) / zCells;
   } else {
      if (zCells == 1) {
	 sim->dz_cell[0] = z_max-z_min;
      } else {
	 const Real dz0 = (z_max-z_min) / zCells / 3.0;
	 const Real dz = 2*((z_max-z_min) - zCells*dz0) / (zCells*(zCells-1));
	 for (builder::ID i=0; i<zCells; ++i) sim->dz_cell[i] = dz0 + i*dz;
      }
   }

   return true;
}

builder::ID LogicallyCartesianBuilder::calculateNeighbourID(builder::ID I,builder::ID J,builder::ID K,int i,int j,int k) {
   // Check that the given neighbour is within the simulation box:
   builder::ID i_out = I + i;
   builder::ID j_out = J + j;
   builder::ID k_out = K + k;

   // Check that neighbour i-index is within the simulation volume:
   if (i_out > sim->x_blocks-1) {
      if (xPeriodic == false) {
	 return numeric_limits<builder::ID>::max();
      } else {
	 if (sim->x_blocks == 1) i_out = 0;
	 else if (i_out > MAX_INDEX) i_out = sim->x_blocks-1 - (numeric_limits<builder::ID>::max()-i_out);
	 else i_out -= sim->x_blocks;
      }
   }

   // Check that neighbour j-index is within the simulation volume:
   if (j_out > sim->y_blocks-1) {
      if (yPeriodic == false) {
	 return numeric_limits<builder::ID>::max();
      } else {
	 if (sim->y_blocks == 1) j_out = 0;
	 else if (j_out > MAX_INDEX) j_out = sim->y_blocks-1 - (numeric_limits<builder::ID>::max()-j_out);
	 else j_out -= sim->y_blocks;
      }
   }

   // Check that neighbour k-index is within the simulation volume:
   if (k_out > sim->z_blocks-1) {
      if (zPeriodic == false) {
	 return numeric_limits<builder::ID>::max();
      } else {
	 if (sim->z_blocks == 1) k_out = 0;
	 else if (k_out > MAX_INDEX) k_out = sim->z_blocks-1 - (numeric_limits<builder::ID>::max()-k_out);
	 else k_out -= sim->z_blocks;
      }
   }
   return k_out*sim->y_blocks*sim->x_blocks + j_out*sim->x_blocks + i_out;
}

unsigned char LogicallyCartesianBuilder::calculateNeighbourTypeID(int i,int j,int k) {
   ++i; ++j; ++k;
   return k*9+j*3+i;
}

unsigned char LogicallyCartesianBuilder::countNeighbours(builder::ID i,builder::ID j,builder::ID k) {
   unsigned char N_neighbours = 0;
   for (int ii=-1; ii<2; ++ii) for (int jj=-1; jj<2; ++jj) for (int kk=-1; kk<2; ++kk) {
      if (ii == 0 && (jj == 0 && kk == 0)) continue; // cell is not its own neighbour
      if (i + ii > sim->x_blocks-1 && xPeriodic == false) continue;
      if (j + jj > sim->y_blocks-1 && yPeriodic == false) continue;
      if (k + kk > sim->z_blocks-1 && zPeriodic == false) continue;
      ++N_neighbours;
   }
   return N_neighbours;
}

bool LogicallyCartesianBuilder::defineParameters(ConfigReader& cr) {
   const Real DEFAULT = numeric_limits<Real>::infinity();
   const builder::ID ONE = 1;

   cr.add(prefix+".dx_uniform","Is cell size uniform in x-direction (yes/no) ?",string("yes"));
   cr.add(prefix+".dy_uniform","Is cell size uniform in y-direction (yes/no) ?",string("yes"));
   cr.add(prefix+".dz_uniform","Is cell size uniform in z-direction (yes/no) ?",string("yes"));
   cr.add(prefix+".x_min","Minimum value for x-coordinate (float).",DEFAULT);
   cr.add(prefix+".y_min","Minimum value for y-coordinate (float).",DEFAULT);
   cr.add(prefix+".z_min","Minimum value for z-coordinate (float).",DEFAULT);
   cr.add(prefix+".x_max","Maximum value for x-coordinate (float).",DEFAULT);
   cr.add(prefix+".y_max","Maximum value for y-coordinate (float).",DEFAULT);
   cr.add(prefix+".z_max","Maximum value for z-coordinate (float).",DEFAULT);
   cr.add(prefix+".x_size","Number of cells in x-direction",ONE);
   cr.add(prefix+".y_size","Number of cells in y-direction",ONE);
   cr.add(prefix+".z_size","Number of cells in z-direction",ONE);
   cr.add(prefix+".x_periodic","Is mesh periodic in x-direction (yes/no) ?",string("no"));
   cr.add(prefix+".y_periodic","Is mesh periodic in y-direction (yes/no) ?",string("no"));
   cr.add(prefix+".z_periodic","Is mesh periodic in z-direction (yes/no) ?",string("no"));

   cr.add(prefix+".x_units","Units of measure for x-coordinate (string).",string(""));
   cr.add(prefix+".y_units","Units of measure for y-coordinate (string).",string(""));
   cr.add(prefix+".z_units","Units of measure for z-coordinate (string).",string(""));
   cr.add(prefix+".x_label","Label for x-coordinate (string).","x-coordinate");
   cr.add(prefix+".y_label","Label for y-coordinate (string).","y-coordinate");
   cr.add(prefix+".z_label","Label for z-coordinate (string).","z-coordinate");

   cr.add(prefix+".input_units","Units in which x_min etc. are given, defaults to meters. Accepted m/km/AU/RE/RS/RV (string).",string("m"));
   cr.add(prefix+".geometry","Geometry of mesh (cartesian/cylindrical/spherical)",string("cartesian"));
   cr.add(prefix+".force_cubic_cells","Force Cartesian cells to be cubic by automatically setting coordinate limits of axes with only one cell (yes/no)",string("no"));
   return true;
}

bool LogicallyCartesianBuilder::finalize() {return false;}

bool LogicallyCartesianBuilder::getCellIndices(unsigned int N_cells,const builder::ID* const globalIDs,builder::ID*& indices) {
   indices = NULL;
   if (initialized == false) return initialized;

   builder::ID i,j,k;
   indices = new builder::ID[3*N_cells];

   // Calculate (i,j,k) indices for all given cell global IDs and copy them to indices:
   for (builder::ID c=0; c<N_cells; ++c) {
      const builder::ID cellID = globalIDs[c];
      calculateCellIndices(cellID,i,j,k);

      indices[c*3+0] = i;
      indices[c*3+1] = j;
      indices[c*3+2] = k;
   }

   return initialized;
}

bool LogicallyCartesianBuilder::getCellInfo(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,unsigned char* N_neighbours) {
   if (initialized == false) return initialized;

   builder::ID i,j,k;
   const builder::ID N_cells = cellOffsetEnd-cellOffsetStart;

   for (builder::ID cell=0; cell<N_cells; ++cell) {
      builder::ID cellID = cellOffsetStart + cell;
      calculateCellIndices(cellID,i,j,k);
      cellIDs[cell]      = cellID;
      N_neighbours[cell] = countNeighbours(i,j,k);
   }

   return true;
}

bool LogicallyCartesianBuilder::getCellNeighbours(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* cellIDs,
						  unsigned char* N_neighbours,builder::ID* neighbourIDs,unsigned char* neighbourTypeIDs) {
   if (initialized == false) return initialized;

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

   return true;
}

bool LogicallyCartesianBuilder::getInitialState(builder::ID cellOffsetStart,builder::ID cellOffsetEnd,builder::ID* globalIDs,double* coordinates) {
   if (initialized == false) return initialized;

   builder::ID i_block,j_block,k_block;
   const builder::ID N_blocks = cellOffsetEnd-cellOffsetStart;

   for (builder::ID block=0; block<N_blocks; ++block) {
      builder::ID blockGID = globalIDs[block];
      calculateCellIndices(blockGID,i_block,j_block,k_block);

      coordinates[3*block+0] = sim->x_crds_node[i_block*block::WIDTH_X + 0];
      coordinates[3*block+1] = sim->y_crds_node[j_block*block::WIDTH_Y + 0];
      coordinates[3*block+2] = sim->z_crds_node[k_block*block::WIDTH_Z + 0];
   }

   return true;
}

bool LogicallyCartesianBuilder::getMeshBoundingBox(Real& x_min,Real& y_min,Real& z_min,Real& dx0,Real& dy0,Real& dz0) {
   if (initialized == false) return initialized;

   // TEST
   x_min = this->x_min;
   y_min = this->y_min;
   z_min = this->z_min;
   dx0 = sim->dx_cell[0];
   dy0 = sim->dy_cell[0];
   dz0 = sim->dz_cell[0];
   // END TEST

   return initialized;
}

bool LogicallyCartesianBuilder::getNumberOfCells(builder::ID& N_blocks) { 
   N_blocks = sim->x_blocks*sim->y_blocks*sim->z_blocks;
   return initialized;
}

bool LogicallyCartesianBuilder::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) { 
   initialized = true;
   simClasses.logger << "(LC BUILDER) Starting to initialize" << endl;
   this->sim = &sim;

   // Insert config file parameters to ConfigReader:
   if (defineParameters(cr) == false) {
      initialized = false;
      return writeInitStatus(initialized,simClasses,"Failed to add parameters to config file reader");
   }

   // Parse parameter values:
   cr.parse();
   string xPeriodicString,yPeriodicString,zPeriodicString;
   string dxUniformString,dyUniformString,dzUniformString;
   string geometryString,inputUnitsString,forceCubicCellsString;
   cr.get(prefix+".dx_uniform",dxUniformString);
   cr.get(prefix+".dy_uniform",dyUniformString);
   cr.get(prefix+".dz_uniform",dzUniformString);
   cr.get(prefix+".x_min",x_min);
   cr.get(prefix+".y_min",y_min);
   cr.get(prefix+".z_min",z_min);
   cr.get(prefix+".x_max",x_max);
   cr.get(prefix+".y_max",y_max);
   cr.get(prefix+".z_max",z_max);
   cr.get(prefix+".x_size",xCells);
   cr.get(prefix+".y_size",yCells);
   cr.get(prefix+".z_size",zCells);
   cr.get(prefix+".x_periodic",xPeriodicString);
   cr.get(prefix+".y_periodic",yPeriodicString);
   cr.get(prefix+".z_periodic",zPeriodicString);
   cr.get(prefix+".x_label",xLabel);
   cr.get(prefix+".y_label",yLabel);
   cr.get(prefix+".z_label",zLabel);
   cr.get(prefix+".x_units",xUnits);
   cr.get(prefix+".y_units",yUnits);
   cr.get(prefix+".z_units",zUnits);
   cr.get(prefix+".geometry",geometryString);
   cr.get(prefix+".input_units",inputUnitsString);
   cr.get(prefix+".force_cubic_cells",forceCubicCellsString);

   // Process input values:
   if (dxUniformString == "yes") { xUniform = true; }
   else if (dxUniformString == "no") { xUniform = false; }
   else { simClasses.logger << "\t WARNING: Unknown value in parameter '" << prefix+".dx_uniform', should be 'yes' or 'no'." << endl; }
   if (dyUniformString == "yes") { yUniform = true; }
   else if (dyUniformString == "no") { yUniform = false; }
   else { simClasses.logger << "\t WARNING: Unknown value in parameter '" << prefix+".dy_uniform', should be 'yes' or 'no'." << endl; }
   if (dzUniformString == "yes") { zUniform = true; }
   else if (dzUniformString == "no") { zUniform = false; }
   else { simClasses.logger << "\t WARNING: Unknown value in parameter '" << prefix+".dz_uniform', should be 'yes' or 'no'." << endl; }

   if (xPeriodicString == "yes") { xPeriodic = true; }
   else if (xPeriodicString == "no") { xPeriodic = false; }
   else { simClasses.logger << "\t WARNING: Unknown value in parameter '" << prefix+".x_periodic', should be 'yes' or 'no'." << endl; }
   if (yPeriodicString == "yes") { yPeriodic = true; }
   else if (yPeriodicString == "no") { yPeriodic = false; }
   else { simClasses.logger << "\t WARNING: Unknown value in parameter '" << prefix+".y_periodic', should be 'yes' or 'no'." << endl; }
   if (zPeriodicString == "yes") { zPeriodic = true; }
   else if (zPeriodicString == "no") { zPeriodic = false; }
   else { simClasses.logger << "\t WARNING: Unknown value in parameter '" << prefix+".z_periodic', should be 'yes' or 'no'." << endl; }

   double inputUnitsScale = simClasses.constants.getDistanceInSI(inputUnitsString);
   if (inputUnitsScale == numeric_limits<double>::infinity()) {
      simClasses.logger << "\t ERROR: Invalid input units '" << inputUnitsString << "' given in confile file." << endl << write;
      initialized = false;
   }

   if (geometryString == "cartesian") {
      geometry = Cartesian;
      sim.meshGeometry = vlsv::geometry::CARTESIAN;
      x_min *= inputUnitsScale;
      x_max *= inputUnitsScale;
      y_min *= inputUnitsScale;
      y_max *= inputUnitsScale;
      z_min *= inputUnitsScale;
      z_max *= inputUnitsScale;
   } else if (geometryString == "cylindrical") {
      geometry = Cylindrical;
      sim.meshGeometry = vlsv::geometry::CYLINDRICAL;
      x_min *= inputUnitsScale;
      x_max *= inputUnitsScale;
      z_min *= inputUnitsScale;
      z_max *= inputUnitsScale;
   } else if (geometryString == "spherical") {
      geometry = Spherical;
      sim.meshGeometry = vlsv::geometry::SPHERICAL;
      x_min *= inputUnitsScale;
      x_max *= inputUnitsScale;
   } else {
      simClasses.logger << "\t ERROR: Unknown mesh geometry '" << geometryString << "' should be 'cartesian', 'cylindrical', or 'spherical'." << endl << write;
      initialized = false;
      sim.meshGeometry = vlsv::geometry::UNKNOWN;
   }

   // try to force cubic cells in Cartesian geometry
   if(forceCubicCellsString == "yes") {
      if(geometryString != "cartesian") {
	 simClasses.logger << "\t ERROR: force_cubic_cells = yes allowed only with geometry = cartesian" << endl << write;
	 initialized = false;
	 goto forceCubicCellsFailed;
      }
      if(dxUniformString != "yes" || dyUniformString != "yes" || dzUniformString != "yes") {
	 simClasses.logger << "\t ERROR: force_cubic_cells = yes allowed only with dx_uniform = dy_uniform = dz_uniform = yes" << endl << write;
	 initialized = false;
	 goto forceCubicCellsFailed;
      }
      // find out cell size of the first axis with more than one cell or otherwise use x axis
      Real dc = -1.0;
      if(xCells > 1)      { dc = (x_max - x_min)/xCells; }
      else if(yCells > 1) { dc = (y_max - y_min)/yCells; }
      else if(zCells > 1) { dc = (z_max - z_min)/zCells; }
      else                { dc = (x_max - x_min)/xCells; }
      // reset coordinate limits of single cell axes to [-dc/2,dc/2]
      const Real dcHalf = dc/2.0;
      if(dc > 0) {
	 if(xCells == 1) {
	    x_min = -dcHalf;
	    x_max = dcHalf;
	    simClasses.logger << "\t forcing Cartesian cubic cells: resetting coordinate limits of a single cell axis: x (x_min = " << x_min << ", x_max = " << x_max << ")" << endl;
	 }
	 if(yCells == 1) {
	    y_min = -dcHalf;
	    y_max = dcHalf;
	    simClasses.logger << "\t forcing Cartesian cubic cells: resetting coordinate limits of a single cell axis: y (y_min = " << y_min << ", y_max = " << y_max << ")" << endl;
	 }
	 if(zCells == 1) {
	    z_min = -dcHalf;
	    z_max = dcHalf;
	    simClasses.logger << "\t forcing Cartesian cubic cells: restting coordinate limits of a single cell axis: z (z_min = " << z_min << ", z_max = " << z_max << ")" << endl;
	 }
      }
      // check resulting cell sizes
      const Real dx = (x_max - x_min)/xCells;
      const Real dy = (y_max - y_min)/yCells;
      const Real dz = (z_max - z_min)/zCells;
      if(dx != dy || dx != dz || dy != dz) {
	 simClasses.logger << "\t ERROR: failed forcing Cartesian cells to be cubic (force_cubic_cells = yes), (dx,dy,dz) = (" << dx << "," << dy << "," << dz << ")" << endl << write;
	 initialized = false;
	 goto forceCubicCellsFailed;
      }
   }
forceCubicCellsFailed:

   // Check input parameters for sanity:
   const Real DEFAULT = numeric_limits<Real>::infinity();
   if (xCells % block::WIDTH_X != 0) {
      xCells += (block::WIDTH_X - xCells % block::WIDTH_X);
      simClasses.logger
	<< "\t WARNING: Given mesh x-size was not exactly divisible by block::WIDTH_X" << endl
	<< "\t          New x_size is " << xCells << endl << write;
   }
   if (yCells % block::WIDTH_Y != 0) {
      yCells += (block::WIDTH_Y - yCells % block::WIDTH_Y);
      simClasses.logger
	<< "\t WARNING: Given mesh y-size was not exactly divisible by block::WIDTH_Y" << endl
	<< "\t          New y_size is " << yCells << endl << write;
   }
   if (zCells % block::WIDTH_Z != 0) {
      zCells += (block::WIDTH_Z - zCells % block::WIDTH_Z);
      simClasses.logger
	<< "\t WARNING: Given mesh z-size was not exactly divisible by block::WIDTH_Z" << endl
	<< "\t          New z_size is " << zCells << endl << write;
   }
   if (x_min == DEFAULT) { initialized = false; }
   if (y_min == DEFAULT) { initialized = false; }
   if (z_min == DEFAULT) { initialized = false; }
   if (x_max == DEFAULT) { initialized = false; }
   if (y_max == DEFAULT) { initialized = false; }
   if (z_max == DEFAULT) { initialized = false; }
   if (x_max < x_min) { initialized = false; }
   if (y_max < y_min) { initialized = false; }
   if (z_max < z_min) { initialized = false; }

   if (initialized == false) {
      simClasses.logger
	<< "\t ERROR: One or more min/max coordinate limits were not found in config file" << endl
	<< "\t        or possibly min value is larger than max value" << endl
	<< "\t x min/max: " << x_min << '\t' << x_max << endl
	<< "\t y min/max: " << y_min << '\t' << y_max << endl
	<< "\t z min/max: " << z_min << '\t' << z_max << endl
	<< write;
   }

   sim.dx_uniform = xUniform;
   sim.dy_uniform = yUniform;
   sim.dz_uniform = zUniform;
   sim.x_periodic = xPeriodic;
   sim.y_periodic = yPeriodic;
   sim.z_periodic = zPeriodic;
   sim.x_min = x_min;
   sim.x_max = x_max;
   sim.y_min = y_min;
   sim.y_max = y_max;
   sim.z_min = z_min;
   sim.z_max = z_max;
   sim.x_blocks = xCells / block::WIDTH_X;
   sim.y_blocks = yCells / block::WIDTH_Y;
   sim.z_blocks = zCells / block::WIDTH_Z;

   // Check mesh periodicity in cylindrical and spherical geometries:
   switch (geometry) {
    case Cartesian:
      break;
    case Cylindrical:
      if (xPeriodic == true) {
	 simClasses.logger << "\t WARNING: Ignoring periodic radial coordinate in cylindrical geometry." << endl;
	 xPeriodic = false;
      }
      if (yPeriodic == true) {
	 y_min = 0.0;
	 y_max = 2.0*M_PI;
      }
      break;
    case Spherical:
      if (xPeriodic == true) {
	 simClasses.logger << "\t WARNING: Ignoring periodic radial coordinate in spherical geometry." << endl;
	 xPeriodic = false;
      }
      if (yPeriodic == true) {
	 simClasses.logger << "\t WARNING: Theta-coordinate set to periodic in spherical geometry" << endl;
	 //simClasses.logger << "\t WARNING: Ignoring periodic theta coordinate in spherical geometry." << endl << write;
	 //yPeriodic = false;
      }
      if (zPeriodic == true) {
	 simClasses.logger << "\t WARNING: Phi-coordinate set to periodic in spherical geometry" << endl;
	 //z_min = 0.0;
	 //z_max = 2.0*M_PI;
      }
      break;
    default:
      break;
   }

   if (calculateCellSizes() == false) {
      simClasses.logger << "\t WARNING: Failed to calculate cell sizes!" << endl;
      initialized = false;
   }

   if (calculateBlockSizes() == false) {
      simClasses.logger << "\t WARNING: Failed to calculate block sizes!" << endl;
      initialized = false;
   }

   sim.x_crds_node = new Real[xCells+1];
   sim.y_crds_node = new Real[yCells+1];
   sim.z_crds_node = new Real[zCells+1];
   sim.x_crds_node[0] = x_min;
   for (builder::ID i=1; i<xCells; ++i) sim.x_crds_node[i] = sim.x_crds_node[i-1] + sim.dx_cell[i-1];
   sim.x_crds_node[xCells] = x_max;
   sim.y_crds_node[0] = y_min;
   for (builder::ID i=1; i<yCells; ++i) sim.y_crds_node[i] = sim.y_crds_node[i-1] + sim.dy_cell[i-1];
   sim.y_crds_node[yCells] = y_max;
   sim.z_crds_node[0] = z_min;
   for (builder::ID i=1; i<zCells; ++i) sim.z_crds_node[i] = sim.z_crds_node[i-1] + sim.dz_cell[i-1];
   sim.z_crds_node[zCells] = z_max;

   return writeInitStatus(initialized,simClasses,"");
}

bool LogicallyCartesianBuilder::writeInitStatus(bool success,SimulationClasses& simClasses,const std::string& reason) {
   if (success == true) {
      simClasses.logger << "\t Initialization successful." << endl;
   } else {
      simClasses.logger
	<< "\t Initialization failed." << endl
	<< "\t " << reason << endl<< write;
   }
   return success;
}

bool LogicallyCartesianBuilder::writeMesh(Simulation& sim,SimulationClasses& simClasses,const std::string& meshName,bool writeMesh) {
   bool success = true;

   // XML tag called MESH must appear on each VLSV file even if the mesh is not written out:
   map<string,string> attributes;
   if (writeMesh == false) {
      int* array = NULL;
      attributes["name"] = meshName;
      attributes["file"] = sim.meshFileName;
      attributes["type"] = vlsv::mesh::STRING_UCD_MULTI;
      stringstream ss;
      ss << (pargrid::MPI_processID)sim.mpiProcesses;
      attributes["domains"] = ss.str();
      if (simClasses.vlsv.writeArray("MESH",attributes,0,1,array) == false) success = false;

      // If mesh is read from another file we can exit now:
      return success;
   }

   // Create ParGrid array that contains boundary cell local IDs and exchange data with neighbours:
   pargrid::DataID lidArrayID = simClasses.pargrid.addUserData<pargrid::CellID>("localIDs",1);
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
	 const pargrid::CellID blockLID = boundaryBlocks[block];
	 lidArray[blockLID] = blockLID;
      }
   }

   // Exchange block local IDs with neighbor processes:
   if (success == true) {
      if (simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,lidArrayID) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to exchange local IDs!" << endl << write;
	 success = false;
      }
   }

   // Write mesh bounding box and node coordinates (master process only):
   attributes.clear();
   attributes["mesh"] = meshName;
   if (sim.mpiRank == sim.MASTER_RANK) {
      unsigned int bbox[6];
      bbox[0] = sim.x_blocks;   // Number of blocks in x-direction in mesh bounding box.
      bbox[1] = sim.y_blocks;   // Number of blocks in y-direction in mesh bounding box.
      bbox[2] = sim.z_blocks;   // Number of blocks in z-direction in mesh bounding box.
      bbox[3] = block::WIDTH_X; // Number of cells in each block in x-direction.
      bbox[4] = block::WIDTH_Y; // Number of cells in each block in y-direction.
      bbox[5] = block::WIDTH_Z; // Number of cells in each block in z-direction.
      if (simClasses.vlsv.writeArray("MESH_BBOX",attributes,6,1,bbox) == false) success = false;

      if (simClasses.vlsv.writeArray("MESH_NODE_CRDS_X",attributes,xCells+1,1,sim.x_crds_node) == false) success = false;
      if (simClasses.vlsv.writeArray("MESH_NODE_CRDS_Y",attributes,yCells+1,1,sim.y_crds_node) == false) success = false;
      if (simClasses.vlsv.writeArray("MESH_NODE_CRDS_Z",attributes,zCells+1,1,sim.z_crds_node) == false) success = false;
   } else {
      unsigned int* ptr = NULL;
      if (simClasses.vlsv.writeArray("MESH_BBOX",attributes,0,1,ptr) == false) success = false;

      Real* crd = NULL;
      if (simClasses.vlsv.writeArray("MESH_NODE_CRDS_X",attributes,0,1,crd) == false) success = false;
      if (simClasses.vlsv.writeArray("MESH_NODE_CRDS_Y",attributes,0,1,crd) == false) success = false;
      if (simClasses.vlsv.writeArray("MESH_NODE_CRDS_Z",attributes,0,1,crd) == false) success = false;
   }
   if (success == false) {
      simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write mesh bounding box information" << endl << write;
   }

   // Count total number of blocks (local+remote) in the domain:
   const vector<pargrid::CellID>& globalIDs = simClasses.pargrid.getGlobalIDs();

   if (success == true) {
      // Flag periodic neighbors so that plugin knows to exclude them:
      const vector<pargrid::MPI_processID>& hosts = simClasses.pargrid.getHosts();
      const vector<pargrid::CellID>& boundaryBlocks = simClasses.pargrid.getBoundaryCells(sim.defaultStencilID);
      vector<bool> validity;
      validity.resize(simClasses.pargrid.getNumberOfAllCells());
      for (pargrid::CellID block=simClasses.pargrid.getNumberOfLocalCells(); block<simClasses.pargrid.getNumberOfAllCells(); ++block)
	validity[block] = false;

      for (size_t block=0; block<boundaryBlocks.size(); ++block) {
	 // Calculate local block's (i,j,k) indices:
	 const pargrid::CellID localBlockLID = boundaryBlocks[block];
	 const pargrid::CellID localBlockGID = simClasses.pargrid.getGlobalIDs()[localBlockLID];

	 builder::ID i_block,j_block,k_block;
	 calculateCellIndices(localBlockGID,i_block,j_block,k_block);

	 // Get local block's neighbors:
	 const pargrid::CellID* const neighborLIDs = simClasses.pargrid.getCellNeighbourIDs(localBlockLID);
	 for (size_t n=0; n<pargrid::N_neighbours; ++n) {
	    // Get neighbor local ID:
	    const pargrid::CellID neighborLID = neighborLIDs[n];

	    // Skip non-existing neighbors:
	    if (neighborLID == pargrid::INVALID_CELLID) continue;

	    // Skip neighbors that are local to this process:
	    if (hosts[neighborLID] == sim.mpiRank) continue;

	    // Calculate neighbor block's (i,j,k) indices:
	    const pargrid::CellID neighborGID = simClasses.pargrid.getGlobalIDs()[neighborLID];
	    builder::ID i_nbr,j_nbr,k_nbr;
	    calculateCellIndices(neighborGID,i_nbr,j_nbr,k_nbr);

	    bool periodic = false;
	    switch (geometry) {
	     case Cartesian:
	       // In Cartesian geometry periodic ghost zones are removed:
	       if (abs((int64_t)i_nbr-(int64_t)i_block) > 1) periodic = true;
	       if (abs((int64_t)j_nbr-(int64_t)j_block) > 1) periodic = true;
	       if (abs((int64_t)k_nbr-(int64_t)k_block) > 1) periodic = true;
	       break;
	     case Cylindrical:
	       // In cylindrical geometry phi-periodic (y-coordinate) ghost zones are accepted:
	       if (abs((int64_t)i_nbr-(int64_t)i_block) > 1) periodic = true;
	       if (abs((int64_t)k_nbr-(int64_t)k_block) > 1) periodic = true;
	       break;
	     case Spherical:
	       // In spherical geometry phi-periodic (z-coordinate) ghost zones are accepted:
	       if (abs((int64_t)i_nbr-(int64_t)i_block) > 1) periodic = true;
	       if (abs((int64_t)j_nbr-(int64_t)j_block) > 1) periodic = true;
	       if (abs((int64_t)k_nbr-(int64_t)k_block) > 1) periodic = true;
	       break;
	     default:
	       // In Cartesian geometry periodic ghost zones are removed:
	       if (abs((int64_t)i_nbr-(int64_t)i_block) > 1) periodic = true;
	       if (abs((int64_t)j_nbr-(int64_t)j_block) > 1) periodic = true;
	       if (abs((int64_t)k_nbr-(int64_t)k_block) > 1) periodic = true;
	       break;
	    }

	    // Flag periodic neighbor as invalid:
	    if (periodic == false) validity[neighborLID] = true;
	 }
      }

      const pargrid::CellID N_locals = simClasses.pargrid.getNumberOfLocalCells();

      // Push all valid block global IDs (local+ghost) to vector validGlobalIDs:
      vector<pargrid::CellID> validGlobalIDs;
      for (pargrid::CellID block=0; block<N_locals; ++block) {
	 validGlobalIDs.push_back( globalIDs[block] );
      }
      for (size_t i=N_locals; i<simClasses.pargrid.getNumberOfAllCells(); ++i) {
	 if (validity[i] == false) continue;
	 validGlobalIDs.push_back( globalIDs[i] );
      }

      const pargrid::CellID N_blocks = validGlobalIDs.size();
      const pargrid::CellID N_ghosts = N_blocks - N_locals;

      // Write mesh block global IDs and metadata:
      attributes.clear();
      attributes["name"] = meshName;
      attributes["type"] = vlsv::mesh::STRING_UCD_MULTI;
      stringstream ss;
      ss << (pargrid::MPI_processID)sim.mpiProcesses;
      attributes["domains"] = ss.str();
      attributes["xperiodic"] = "no";
      attributes["yperiodic"] = "no";
      attributes["zperiodic"] = "no";
      if (xPeriodic == true) attributes["xperiodic"] = "yes";
      if (yPeriodic == true) attributes["yperiodic"] = "yes";
      if (zPeriodic == true) attributes["zperiodic"] = "yes";
      if (xLabel.size() > 0) attributes["xlabel"] = xLabel;
      if (yLabel.size() > 0) attributes["ylabel"] = yLabel;
      if (zLabel.size() > 0) attributes["zlabel"] = zLabel;
      if (xUnits.size() > 0) attributes["xunits"] = xUnits;
      if (yUnits.size() > 0) attributes["yunits"] = yUnits;
      if (zUnits.size() > 0) attributes["zunits"] = zUnits;

      switch (geometry) {
       case Cartesian:
	 attributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	 break;
       case Cylindrical:
	 attributes["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
	 break;
       case Spherical:
	 attributes["geometry"] = vlsv::geometry::STRING_SPHERICAL;
	 break;
       default:
	 attributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	 break;
      }
      if (simClasses.vlsv.writeArray("MESH",attributes,N_blocks,1,&(validGlobalIDs[0])) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write block global IDs" << endl << write;
	 success = false;
      }

      // Write number of total and ghost blocks in the domain:
      attributes.clear();
      attributes["mesh"] = meshName;
      int domainSize[2];
      domainSize[0] = N_blocks;
      domainSize[1] = N_ghosts;
      if (simClasses.vlsv.writeArray("MESH_DOMAIN_SIZES",attributes,1,2,domainSize) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write domain size" << endl << write;
	 success = false;
      }

      // Wait for block local IDs:
      if (simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,lidArrayID) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to exchange block local IDs!" << endl << write;
	 success = false;
      }

      // Make lists of valid ghost block local IDs and domains:
      vector<pargrid::CellID> validNeighbors;
      vector<pargrid::MPI_processID> validHosts;
      for (size_t i=N_locals; i<simClasses.pargrid.getNumberOfAllCells(); ++i) {
	 if (validity[i] == false) continue;
	 validNeighbors.push_back(lidArray[i]);
	 validHosts.push_back(hosts[i]);
      }

      // Write ghost block local IDs:
      if (simClasses.vlsv.writeArray("MESH_GHOST_LOCALIDS",attributes,N_ghosts,1,&(validNeighbors[0])) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write ghost block local IDs!" << endl << write;
	 success = false;
      }

      // Write ghost block domain IDs:
      if (simClasses.vlsv.writeArray("MESH_GHOST_DOMAINS",attributes,N_ghosts,1,&(validHosts[0])) == false) {
	 simClasses.logger << "(GRIDBUILDER) ERROR: Failed to write block ghost domain IDs!" << endl << write;
	 success = false;
      }
   }

   // Remove ParGrid localID array created above:
   simClasses.pargrid.removeDataTransfer(pargrid::DEFAULT_STENCIL,lidArrayID);
   simClasses.pargrid.removeUserData(lidArrayID);

   return success;
}


