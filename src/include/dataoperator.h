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

#ifndef DATAOPERATOR_H
#define DATAOPERATOR_H

#include <simulation.h>
#include <simulationclasses.h>
#include <particle_list_base.h>
#include <configreader.h>

// *********************************************** //
// ***** DATAOPERATOR BASE CLASS DECLARATION ***** //
// *********************************************** //

class DataOperator {
 public:
   DataOperator();
   virtual ~DataOperator();

   virtual bool finalize() = 0;
   bool getInitialized() const {return initialized;}
   virtual std::string getName() const = 0;
   virtual bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   virtual bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) = 0;

 protected:
   bool initialized;                     /**< If true, derived DataOperator initialized successfully.*/
   Simulation* sim;                      /**< Pointer to generic simulation variables.*/
   SimulationClasses* simClasses;        /**< Pointer to generic simulation classes.*/
};

#endif

