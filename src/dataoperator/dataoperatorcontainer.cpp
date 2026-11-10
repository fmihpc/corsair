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

#include "dataoperatorcontainer.h"

using namespace std;

DataOperatorContainer::DataOperatorContainer() {
   initialized = false;
}

DataOperatorContainer::~DataOperatorContainer() {
   finalize();
}

bool DataOperatorContainer::finalize() {
   // Do not finelize if DRC has not been initialized:
   if (initialized == false) return false;

   // Call finalize for all DataOperators and delete them:
   for (size_t i=0; i<operators.size(); ++i) {
      operators[i]->finalize();
      delete operators[i]; 
      operators[i] = NULL;
   }
   operators.resize(0);
   initialized = false;
   return true;
}

bool DataOperatorContainer::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   // Do not init twice:
   if (initialized == true) return true;

   bool success = true;
   this->sim = &sim;
   this->simClasses = &simClasses;
   this->configReader = &cr;

   // Read DataOperator exclude list from config file:
   const string prefix = "DataOperatorExcludes";
   vector<string> excludeList;
   cr.addComposed(prefix+".exclude_list","Names of DataOperators that will not be called during simulation (string).");
   cr.parse();
   cr.get(prefix+".exclude_list",excludeList);

   for (size_t i=0; i<excludeList.size(); ++i) {
      // Skip empty lines:
      if (excludeList[i].size() == 0) continue;

      // Add DataOperator's name to excludedOperators:
      excludedOperators.insert(excludeList[i]);
   }

   initialized = true;   
   return success;
}

bool DataOperatorContainer::registerOperator(DataOperator* op) {
   if (initialized == false) return false;
   bool success = true;

   // Check if an operator with the same name already exists:
   string name = op->getName();
   for (size_t i=0; i<operators.size(); ++i) {
      if (operators[i]->getName() == name) {
	 success = false;
	 simClasses->logger << "(DATAOPERATORCONTAINER) ERROR: DataOperator with name '" << op->getName() << " already exists!" << endl << write;
	 break;
      }
   }

   // Check if the operator is on exlude list. If yes, do not add it to loaded operators:
   set<string>::const_iterator it = excludedOperators.find(name);
   if (it != excludedOperators.end()) {
      simClasses->logger << "(DATAOPERATORCONTAINER) Added DataOperator '" << name << "' to exclude list" << endl;
   }

   // Operator does not exist yet, add it 
   // to vector operators and initialze:
   if (success == true && it == excludedOperators.end()) {
      operators.push_back(op);
      if (op->initialize(*configReader,*sim,*simClasses) == false) {
	 simClasses->logger << "(DATAOPERATORCONTAINER) ERROR: DataOperator '" << op->getName() << "' failed to initialize!" << endl << write;
	 success = false;
      }
   }
   return success;
}

unsigned int DataOperatorContainer::size() const {
   return operators.size();
}

bool DataOperatorContainer::writeData(const string& meshName,vector<ParticleListBase*>& particleLists) {
   bool success = true;
   for (size_t i=0; i<operators.size(); ++i) {
      if (operators[i]->writeData(meshName,particleLists) == false) {
	 simClasses->logger << "(DATAOPERATORCONTAINER) ERROR: DataOperator '" << operators[i]->getName();
	 simClasses->logger << "' failed to write its data!" << endl << write;
	 success = false;
      }
   }
   return success;
}
