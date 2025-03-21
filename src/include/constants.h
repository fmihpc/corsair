/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2025 Finnish Meteorological Institute
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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>
#include <map>

namespace constants {
   const double BOLTZMANN = 1.3806503e-23;            // Boltzmann's constant, J/K
   const double CHARGE_ELEMENTARY = 1.60217646e-19;   // Elementary charge, C
   const double DIST_ASTRONOMICAL_UNIT = 1.49598e+11; // Astronomical unit, m
   const double DIST_METER = 1.0;                     // Meter in meters (for completeness)
   const double DIST_SOLAR_RADIUS = 6.955e+08;        // Sun equatorial radius in meters
   const double DIST_MERCURY_RADIUS = 2439.7e3;       // Radius of Mercury in meters
   const double DIST_VENUS_RADIUS = 6051.8e3;         // Radius of Venus in meters
   const double DIST_EARTH_RADIUS = 6378.1e3;         // Radius of Earth in meters
   const double DIST_MARS_RADIUS = 3389.5e3;          // Radius of Mars in meters
   const double GRAVITY = 6.67300e-11;                // Gravitational constant, m^3/(kg s^2)
   const double MASS_ELECTRON = 9.10938188e-31;       // Electron mass, kg
   const double MASS_PROTON = 1.67262158e-27;         // Proton mass, kg
   const double MASS_ALPHA = 6.6446573450e-27;        // Alpha particle mass, kg
   const double MASS_OXYGEN = 2.6567625437e-26;       // Oxygen atom mass, kg
   const double MASS_SOLAR = 1.9891e+30;              // Solar mass, kg
   const double MASS_MERCURY = 0.3301e24;             // Mercury mass, kg
   const double MASS_VENUS = 4.8685e24;               // Venus mass, kg
   const double MASS_EARTH = 5.9722e24;               // Earth mass, kg
   const double MASS_MARS = 6.4185e23;                // Mars mass, kg
   const double OMEGA_SOLAR_EQUATOR = 2.972e-6;       // Solar angular velocity at equator, 14.713 degrees/day, 1/s
   const double PERMEABILITY = 1.25663706e-6;         // Permeability of free space, kg m/C^2
   const double PERMITTIVITY = 8.85418782e-12;        // Permittivity of free space, s^4 A^2/(kg m^3) 
   const double SPEED_LIGHT = 299792458.0;            // Speed of light, m/s
   const double EV_TO_KELVIN = 11604.5250061657;      // Electron volt to kelvin conversion factor, K/eV
}

class Constants {
 public:
   Constants();
   double get(const std::string& name) const;
   double getDistanceInSI(const std::string& name) const;
   double getEnergyInSI(const std::string& unit) const;
   double notFound() const;
 private:
   std::map<std::string,double> values;
};

#endif
