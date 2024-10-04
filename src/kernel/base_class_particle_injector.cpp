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
#include <limits>
#include <cmath>

#include <base_class_particle_injector.h>

using namespace std;

ParticleInjectorBase::ParticleInjectorBase() { }

ParticleInjectorBase::~ParticleInjectorBase() { }

Real ParticleInjectorBase::getMaximumSpeed() {return 0.0;}
void ParticleInjectorBase::setType(string type) { this->type = type; }
string ParticleInjectorBase::getType() { return this->type; }

bool ParticleInjectorBase::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
				      const std::string& regionName,const ParticleListBase* plist) {
   this->sim = &sim;
   this->simClasses = &simClasses;
   this->type = "";
   return true;
}

