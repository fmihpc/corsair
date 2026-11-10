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

#ifndef OBJECT_FACTORIES_H
#define OBJECT_FACTORIES_H

#include <object_factory_generic.h>

#include <gridbuilder.h>
#include <base_class_particle_accumulator.h>
#include <base_class_particle_boundary_condition.h>
#include <base_class_particle_injector.h>
#include <base_class_particle_propagator.h>

/** Container class for all Object Factories used by Corsair kernel.*/
struct ObjectFactories {
   ObjectFactoryGeneric<GridBuilder> gridBuilders;                            /**< Factory that creates grid builders.*/
   ObjectFactoryGeneric<ParticleAccumulatorBase> particleAccumulators;        /**< Factory that creates particle accumulators.*/
   ObjectFactoryGeneric<ParticleBoundaryCondBase> particleBoundaryConditions; /**< Factory that creates particle boundary conditions.*/
   ObjectFactoryGeneric<ParticleInjectorBase> particleInjectors;              /**< Factory that creates particle injectors.*/
   ObjectFactoryGeneric<ParticlePropagatorBase> particlePropagators;          /**< Factory that creates particle propagators.*/
};

#endif
