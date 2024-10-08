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

#include <limits>
#include <simulation.h>

using namespace std;

//unsigned int Simulation::N_particleSpecies = 0;
pargrid::DataID Simulation::crdsDataID       = numeric_limits<pargrid::DataID>::max();
pargrid::StencilID Simulation::crdsStencilID = numeric_limits<pargrid::StencilID>::max();

Simulation::Simulation(): MASTER_RANK(0) {
   argn = 0;
   args = NULL;
   t = 0.0;
   timestep = 0;
   meshWrittenStep = numeric_limits<unsigned int>::max();
   runTests = false;

   dx_block = NULL;
   dy_block = NULL;
   dz_block = NULL;
   dx_cell = NULL;
   dy_cell = NULL;
   dz_cell = NULL;
   dx_uniform = true;
   dy_uniform = true;
   dz_uniform = true;
   x_periodic = false;
   y_periodic = false;
   z_periodic = false;
   x_min = 0;
   x_max = 0;
   y_min = 0;
   y_max = 0;
   z_min = 0;
   z_max = 0;
   x_crds_node = NULL;
   y_crds_node = NULL;
   z_crds_node = NULL;
   meshGeometry = vlsv::geometry::UNKNOWN;
   initializing = true;
}

Simulation::~Simulation() { 
   delete [] dx_block; dx_block = NULL;
   delete [] dy_block; dy_block = NULL;
   delete [] dz_block; dz_block = NULL;
   delete [] dx_cell; dx_cell = NULL;
   delete [] dy_cell; dy_cell = NULL;
   delete [] dz_cell; dz_cell = NULL;
   delete [] x_crds_node; x_crds_node = NULL;
   delete [] y_crds_node; y_crds_node = NULL;
   delete [] z_crds_node; z_crds_node = NULL;
}

bool Simulation::finalize() {
   argn = 0;
   args = NULL;
   return true;
}

bool Simulation::initialize(int argn,char** args) {
   this->argn = argn;
   this->args = args;
   return true;
}
