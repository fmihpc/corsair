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

#include <root_finding.h>

using namespace std;

namespace rootf {

   Real bisection(Real (*function)(Real x),Real xmin,Real xmax,
		  Real tol_dif,Real tol_min,int max_iter) {

      int iterations = 0;
      Real a = xmin;
      Real b = xmax;
      Real root = 0.0, func_a, func_m;

      bool cont = true;
      while (cont == true) {
	 func_a = (*function)(a);
	 // Check if we're close enough to the root of equation
	 if (b - a <= tol_dif) {
	    root = 0.5 * (a + b);
	    cont = false;
	 }

	 // Check if function value is close enough to zero
	 Real m = 0.5 * (a + b);
	 func_m = (*function)(m);
	 if (fabs(func_m) <= tol_min) {
	    root = m;
	    cont = false;
	 }

	 // Split the interval between a and b. Choose new 
	 // interval so that a solution exists.
	 if (func_a * func_m < 0.0)
	   b = m;
	 else
	   a = m;

	 // Check if we've done too many iterations
	 ++iterations;
	 if (iterations > max_iter) {
	    cont = false;
	    root = xmax = 2.0;
	 }
      }
      return root;
   }
}
