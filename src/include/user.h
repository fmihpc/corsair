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

#ifndef USER_H
#define USER_H

#include <vector>

#include <configreader.h>
#include <dataoperatorcontainer.h>
#include <gridbuilder.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <object_factories.h>
#include <particle_list_base.h>

bool propagate(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists);
bool userEarlyInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,std::vector<ParticleListBase*>& particleLists);
bool userLateInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			    const ObjectFactories& objectFactories,std::vector<ParticleListBase*>& particleLists);
bool userFinalization(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists);
bool userRunTests(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists);

/** User-defined function which registers all DataOperators to DataOperatorContainer. Only
 * registered DataOperators are visible in the simulation.
 * @param doc Container for DataOperators.
 * @return If true, all DataOperators were registered successfully.*/
bool registerDataOperators();

/** User-defined function which registers one or more GridBuilders 
 * to GridBuilderFactory.
 * @param gbf GridBuilderFactory in which GridBuilders are to be registered.
 * @return If true, all GridBuilders were registered successfully.*/
bool registerObjectMakers(ObjectFactories& objectFactories);

#endif
