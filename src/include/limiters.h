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

#ifndef LIMITERS_H
#define LIMITERS_H

#include <cmath>
#include <climits>

// ***** TWO-POINT LIMITERS ***** //

template<typename REAL> REAL minmod(REAL u_left,REAL u_rght);

// ***** THREE-POINT LIMITERS ***** //

template<typename REAL> REAL nolimiter(REAL u_left,REAL u_cent,REAL u_rght);
template<typename REAL> REAL MClimiter(REAL u_left,REAL u_cent,REAL u_rght);
template<typename REAL> REAL superbee(REAL u_left,REAL u_cent,REAL u_rght);
template<typename REAL> REAL vanLeer(REAL u_left,REAL u_cent,REAL u_rght);

template<typename REAL> inline
REAL minmod(REAL u_left,REAL u_rght) {
   const REAL ZERO = 0.0;
   if (u_left*u_rght < ZERO) return ZERO;
   if (fabs(u_left) < fabs(u_rght)) return u_left;
   else return u_rght;
}

template<typename REAL> inline
REAL nolimiter(REAL u_left,REAL u_cent,REAL u_rght) {
   const REAL HALF = 0.5;
   return HALF*(u_rght-u_left);
}

template<typename REAL> inline
REAL MClimiter(REAL u_left,REAL u_cent,REAL u_rght) {
   const REAL HALF = 0.5;
   const REAL TWO  = 2.0;   
   const REAL forw   = u_rght-u_cent;
   const REAL back   = u_cent-u_left;
   const REAL cntr   = HALF*(u_rght-u_left);
   const REAL slope1 = minmod(TWO*forw,cntr);
   const REAL slope2 = minmod(TWO*back,cntr);
   return minmod(slope1,slope2);
}

template<typename REAL> inline
REAL superbee(REAL u_left,REAL u_cent,REAL u_rght) {
   const REAL HALF = 0.5;
   const REAL ONE  = 1.0;
   const REAL forw = u_rght-u_cent;
   const REAL back = u_cent-u_left;

   REAL tmp = std::min(fabs(back),fabs(forw));
   tmp = std::min(tmp,HALF*fabs(back));
   tmp = std::min(tmp,fabs(forw));

   REAL sign_forw = ONE;
   if (forw < 0.0) sign_forw = -ONE;
   REAL sign_back = ONE;
   if (back < 0.0) sign_back = -ONE;

   return (sign_forw+sign_back)*tmp;
}

template<typename REAL> inline
REAL vanLeer(REAL u_left,REAL u_cent,REAL u_rght) {
   const REAL ZERO = 0.0;
   const REAL EPS  = std::numeric_limits<REAL>::min();
   return std::max((u_rght-u_cent)*(u_cent-u_left),ZERO)/(EPS+u_rght-u_left);
}

#endif
