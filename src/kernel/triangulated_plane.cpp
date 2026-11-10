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

#include <cmath>
#include <limits>

#include <linear_algebra.h>
#include <triangulated_plane.h>
#include <vlsv_common.h>

using namespace std;

namespace ucdmesh {

   bool getPlanarFaceData(const std::vector<Node>& nodes,const std::vector<TriangularFace>& faces,std::vector<Real>& faceData) {
      bool success = true;
      faceData.clear();
      faceData.resize(faces.size()*ucdmesh::facedataelement::SIZE);

      Real* ptrArea = NULL;
      Real* ptrNormal = NULL;
      Real vector1[3];
      Real vector2[3];
      ptrArea   = &(faceData[ucdmesh::facedataelement::AREA]);
      ptrNormal = &(faceData[ucdmesh::facedataelement::NORMAL_X]);
      for (size_t face=0; face<faces.size(); ++face) {
	 switch (faces[face].cellType) {
	  case vlsv::celltype::TRIANGLE:
	    // Create two triangle edge vectors:
	    vector1[0] = nodes[faces[face].node2].x - nodes[faces[face].node1].x;
	    vector1[1] = nodes[faces[face].node2].y - nodes[faces[face].node1].y;
	    vector1[2] = nodes[faces[face].node2].z - nodes[faces[face].node1].z;
	    vector2[0] = nodes[faces[face].node3].x - nodes[faces[face].node2].x;
	    vector2[1] = nodes[faces[face].node3].y - nodes[faces[face].node2].y;
	    vector2[2] = nodes[faces[face].node3].z - nodes[faces[face].node2].z;

	    // Calculate surface area and normal vector:
	    crossProduct(vector1,vector2,ptrNormal);
	    *ptrArea = vectorMagnitude<3>(ptrNormal);
	    for (int i=0; i<3; ++i) ptrNormal[i] /= *ptrArea;
	    *ptrArea *= 0.5;
	    break;
	  default:
	    *ptrArea = 0.0;
	    for (int i=0; i<3; ++i) ptrNormal[i] = numeric_limits<Real>::infinity();
	    success = false;
	    break;
	 }

	 ptrArea   += ucdmesh::facedataelement::SIZE;
	 ptrNormal += ucdmesh::facedataelement::SIZE;
      }      
      return success;
   }

   void initializePlanarSurface(const Real* corners,std::vector<Node>& nodes,
				std::vector<TriangularFace>& faces) {
      // Store current number of nodes, needed to add correct faces below
      // in case there are multiple planar surfaces:
      const size_t N = nodes.size();

      // Calculate additional node at centroid position:
      Real center[3];
      for (int i=0; i<3; ++i) center[i] = 0.0;
      for (int i=0; i<4; ++i) center[0] += corners[3*i+0];
      for (int i=0; i<4; ++i) center[1] += corners[3*i+1];
      for (int i=0; i<4; ++i) center[2] += corners[3*i+2];
      for (int i=0; i<3; ++i) center[i] /= 4;

      // Generate nodes:
      for (int i=0; i<4; ++i) nodes.push_back(Node(corners[3*i+0],corners[3*i+1],corners[3*i+2]));
      nodes.push_back(Node(center[0],center[1],center[2]));

      // Generate faces:
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,N+0,N+1,N+4));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,N+1,N+2,N+4));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,N+2,N+3,N+4));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,N+0,N+4,N+3));
   }

   /** Create a triangulated spherical surface, starting from a chosen polyhedron.
    * @param corners Vector containing 12 coordinates of the four corners of the plane. Vector can 
    * contain corner coordinates of multiple planes.
    * @param refinements Number of times triangles on chosen polyhedron are refined.
    * @param nodes List of nodes on generated spherical surface.
    * @param faces List of faces on generated spherical surface.
    * @param N_connections Total number of connectivity entries.
    * @return If true, surface was generated successfully.*/
   bool createPlane(const std::vector<Real>& corners,int refinements,std::vector<Node>& nodes,
		    std::vector<TriangularFace>& faces,std::size_t& N_connections) {
      bool result = true;

      // Sanity check on input parameters:
      if (refinements < 0) result = false;

      // Exit if error(s) have occurred:
      if (result == false) return result;

      nodes.clear();
      faces.clear();

      // Generate a triangulated planar surface(s):
      const size_t surfaces = corners.size() / 12;
      for (size_t s=0; s<surfaces; ++s) {
	 initializePlanarSurface(&(corners[12*s]),nodes,faces);
      }

      // Exit if error(s) have occurred:
      if (result == false) {
	 nodes.clear();
	 faces.clear();
	 N_connections = 0;
	 return result;
      }

      // Refine:
      for (int i=0; i<refinements; ++i) refineFaces(nodes,faces);

      N_connections = 5*faces.size();
      return result;
   }

} // namespace ucdmesh
