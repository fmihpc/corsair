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
#include <map>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "operator_timeseries_load.h"

using namespace std;

LoadTSeriesOP::LoadTSeriesOP(): DataOperator() { 
   #if PROFILE_LEVEL > 0
      profileID = -1;
   #endif
}

LoadTSeriesOP::~LoadTSeriesOP() { }

bool LoadTSeriesOP::finalize() {return true;}

std::string LoadTSeriesOP::getName() const {return "TimeSeriesLoad";}

bool LoadTSeriesOP::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   return DataOperator::initialize(cr,sim,simClasses);
}

bool LoadTSeriesOP::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
   if (getInitialized() == false) return false;
   #if PROFILE_LEVEL > 0
      profile::start("Load Time Series",profileID);
   #endif

   // Write buffer to output file:
   map<string,string> attributes;
   //attributes["mesh"] = spatMeshName;
   attributes["xlabel"] = "Time";
   attributes["ylabel"] = "Computational load";
   attributes["xunit"] = "s";
   attributes["yunit"] = "s";

   Real load = 0.0;
   for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
      load += simClasses->pargrid.getCellWeights()[block];
   }

   bool success = true;
   attributes["name"] = "load_max";
   if (simClasses->vlsv.writeWithReduction("TIMESERIES",attributes,1,&load,MPI_MAX) == false) success = false;
   attributes["name"] = "load_min";
   if (simClasses->vlsv.writeWithReduction("TIMESERIES",attributes,1,&load,MPI_MIN) == false) success = false;
   attributes["name"] = "load_avg";   
   load /= sim->mpiProcesses;
   if (simClasses->vlsv.writeWithReduction("TIMESERIES",attributes,1,&load,MPI_SUM) == false) success = false;

   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}


