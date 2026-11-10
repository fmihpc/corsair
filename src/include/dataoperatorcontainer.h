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

#ifndef DATAOPERATORCONTAINER_H
#define DATAOPERATORCONTAINER_H

#include <stdint.h>
#include <set>
#include <vector>

#include <dataoperator.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <particle_list_base.h>

class DataOperatorContainer {
 public:

   DataOperatorContainer();
   ~DataOperatorContainer();

   bool finalize();
   bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   bool registerOperator(DataOperator* op);
   unsigned int size() const;
   bool writeData(const std::string& meshName,std::vector<ParticleListBase*>& particleLists);

 private:
   bool initialized;
   std::set<std::string> excludedOperators; /**< Names of DataOperators that are registered but 
					     * should not be called during simulation.*/
   std::vector<DataOperator*> operators;    /**< All registered DataOperators.*/
   ConfigReader* configReader;              /**< Pointer to configuration file reader.*/
   Simulation* sim;                         /**< Pointer to generic simulation variables.*/
   SimulationClasses* simClasses;           /**< Pointer to generic simulation classes.*/
};

#endif
