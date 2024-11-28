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
#include <cmath>
#include <random>

#include "randomnumber.h"

using namespace std;

namespace corsair {
   random_device rd;
   mt19937 generator(rd());                          /**< Mersenne Twister 19937 random number generator.*/
   normal_distribution<Real> normal(0.0,1.0);        /**< Gaussian random number distribution with zero mean and unit variance.*/
   uniform_real_distribution<Real> uniform(0.0,1.0); /**< Uniform random number distribution in interval [0.0,1.0)*/
   uniform_real_distribution<Real> uniformPM(-1.0,1.0); /**< Uniform random number distribution in interval [-1.0,1.0)*/
}

/** Class default constructor. Sets RandomNumber::initialized to false.*/
RandomNumber::RandomNumber(): initialized(false) { }

/** Destructor. Calls RandomNumber::finalize().
 @see RandomNumber::finalize().*/
RandomNumber::~RandomNumber() {
   finalize();
}

/** Class finalizer function. Sets RandomNumber::initialized to false.
 * @return Return value true if finalization succeeded.*/
bool RandomNumber::finalize() {
   initialized = false;
   return true;
}

/** Return a normally distributed random number with zero mean and unit variance.
 * @return Generated random number.*/
Real RandomNumber::gaussian() {
   return corsair::normal(corsair::generator);
}

/** Initialize random number generator with given seed value.
 * @param seed Seed value.
 * @return Returns value true if initialization succeeded.*/
bool RandomNumber::initialize(int seed) {
   if (initialized == true) { return true; }
   corsair::generator.seed(abs(seed));
   initialized = true;
   return initialized;
}

/** Return a uniform random number in interval [0,1)
 * @return Generated random number.*/
Real RandomNumber::uniform() {
   return corsair::uniform(corsair::generator);
}

/** Return a uniform random number in interval [-1,1)
 * @return Generated random number.*/
Real RandomNumber::uniformPM() {
   return corsair::uniformPM(corsair::generator);
}
