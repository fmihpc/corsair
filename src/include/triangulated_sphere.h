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

#ifndef TRIANGULATED_SPHERE_H
#define TRIANGULATED_SPHERE_H

#include <vector>
#include <stdint.h>

#include <definitions.h>
#include <ucd_mesh.h>

/** Using createSphere function declared in this file a user can create triangulated spherical 
 * surfaces. The generation method is fairly standard, i.e., we start with a triangulated polyhedron 
 * such as a tetrahedron, refine the triangles for chosen number of times, and renormalize the 
 * node coordinates to lie on a spherical surface.
 */

namespace ucdmesh {

   namespace spheregenerator {
      /** Supported polyhedra that can be used in generation of triangulated spherical surface.*/
      enum Method {
	 icosahedron,      // Generate spherical surface starting from an icosahedron
	 octahedron,       // Generate spherical surface starting from an octahedron
	 tetrahedron       // Generate spherical surface starting from a tetrahedron
      };

   } // namespace spheregenerator

   bool createSphere(spheregenerator::Method method,int refinements,Real centroid[3],Real radius,
		     std::vector<Node>& nodes,std::vector<TriangularFace>& faces,std::size_t& N_connections);

   bool getFaceData(const std::vector<Node>& nodes,const std::vector<TriangularFace>& faces,std::vector<Real>& faceData);

} // namespace ucdmesh

#endif
