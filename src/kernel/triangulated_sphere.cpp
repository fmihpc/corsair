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
#include <triangulated_sphere.h>
#include <vlsv_common.h>

using namespace std;

namespace ucdmesh {

   /** Create a triangulated icosahedron.
    * @param nodes Vector where generated nodes are written to.
    * @param faces Vector where generated faces are written to.
    * @return If true, icosahedron was created successfully.*/
   bool createIcosahedron(std::vector<Node>& nodes,std::vector<TriangularFace>& faces) {
      // Generate nodes that make up an icosahedron:
      const Real t = (1.0 + sqrt(5))/2.0;
      const Real s = sqrt(1.0 + t*t);

      nodes.push_back(Node(t/s,1/s,0));
      nodes.push_back(Node(-t/s,1/s,0));
      nodes.push_back(Node(t/s,-1/s,0));

      nodes.push_back(Node(-t/s,-1/s,0));
      nodes.push_back(Node(1/s,0,t/s));
      nodes.push_back(Node(1/s,0,-t/s));

      nodes.push_back(Node(-1/s,0,t/s));
      nodes.push_back(Node(-1/s,0,-t/s));
      nodes.push_back(Node(0,t/s,1/s));

      nodes.push_back(Node(0,-t/s,1/s));
      nodes.push_back(Node(0,t/s,-1/s));
      nodes.push_back(Node(0,-t/s,-1/s));

      // Generate faces:
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,0,8,4));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,1,10,7));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,2,9,11));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,7,3,1));

      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,0,5,10));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,3,9,6));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,3,11,9));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,8,6,4));

      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,2,4,9));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,3,7,11));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,4,2,0));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,9,4,6));

      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,2,11,5));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,0,10,8));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,5,0,2));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,10,5,7));

      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,1,6,8));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,1,8,10));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,6,1,3));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,11,7,5));

      return true;
   }

   /** Create a triangulated octahedron.
    * @param nodes Vector where generated nodes are written to.
    * @param faces Vector where generated faces are written to.
    * @return If true, octahedron was created successfully.*/
   bool createOctahedron(std::vector<Node>& nodes,std::vector<TriangularFace>& faces) {
      // Generate nodes that make up an octahedron:
      nodes.push_back(Node(+1, 0, 0));
      nodes.push_back(Node(-1, 0, 0));
      nodes.push_back(Node( 0,+1, 0));
      nodes.push_back(Node( 0,-1, 0));
      nodes.push_back(Node( 0, 0,+1));
      nodes.push_back(Node( 0, 0,-1));

      // Generate faces:
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,4,0,2));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,5,2,0));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,4,2,1));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,5,1,2));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,4,1,3));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,5,3,1));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,4,3,0));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,5,0,3));

      return true;
   }

   /** Create a triangulated tetrahedron.
    * @param nodes Vector where generated nodes are written to.
    * @param faces Vector where generated faces are written to.
    * @return If true, tetrahedron was created successfully.*/
   bool createTetrahedron(std::vector<Node>& nodes,std::vector<TriangularFace>& faces) {
      // Generate nodes that make up a tetrahedron:
      nodes.push_back(Node(0,0,1));
      nodes.push_back(Node(2.0/3.0*sqrt(2),0,-1.0/3.0));
      nodes.push_back(Node(-sqrt(2.0)/3.0,sqrt(6.0)/3.0,-1.0/3.0));
      nodes.push_back(Node(-sqrt(2.0)/3.0,-sqrt(6.0)/3.0,-1.0/3.0));

      // Generate faces:
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,0,1,2));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,0,2,3));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,0,3,1));
      faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,1,3,2));

      return true;
   }

   bool getFaceData(const std::vector<Node>& nodes,const std::vector<TriangularFace>& faces,std::vector<Real>& faceData) {
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

   /** Normalize node coordinates so that they lie on a spherical surface of given 
    * radius, and shift (translate) the centroid to given position.
    * @param centroid Centroid of the sphere.
    * @param radius Radius of the spherical surface.
    * @param nodes List of nodes on triangulated surface.*/
   void normalizeRadiusAndShiftCentroid(Real centroid[3],Real radius,std::vector<Node>& nodes) {
      Real currentCentroid[3] = {0.0,0.0,0.0};
      for (size_t node=0; node<nodes.size(); ++node) {
	 // Calculate current distance to origin:
	 const Real r = sqrt(nodes[node].x*nodes[node].x + nodes[node].y*nodes[node].y + nodes[node].z*nodes[node].z);

	 // Normalize radius to radius0:
	 nodes[node].x *= (radius/r);
	 nodes[node].y *= (radius/r);
	 nodes[node].z *= (radius/r);

	 currentCentroid[0] += nodes[node].x;
	 currentCentroid[1] += nodes[node].y;
	 currentCentroid[2] += nodes[node].z;
      }

      for (int i=0; i<3; ++i) currentCentroid[i] /= nodes.size();

      // Shift centroid:
      for (size_t node=0; node<nodes.size(); ++node) {
	 nodes[node].x += (centroid[0]-currentCentroid[0]);
	 nodes[node].y += (centroid[1]-currentCentroid[1]);
	 nodes[node].z += (centroid[2]-currentCentroid[2]);
      }
   }

   /** Create a triangulated spherical surface, starting from a chosen polyhedron.
    * @param method Type of polyhedron that is used in mesh generation.
    * @param refinements Number of times triangles on chosen polyhedron are refined.
    * @param centroid Centroid of generated sphere.
    * @param radius Radius of generated sphere.
    * @param nodes List of nodes on generated spherical surface.
    * @param faces List of faces on generated spherical surface.
    * @param N_connections Total number of connectivity entries.
    * @return If true, spherical surface was generated successfully.*/
   bool createSphere(spheregenerator::Method method,int refinements,Real centroid[3],Real radius,
		     std::vector<Node>& nodes,std::vector<TriangularFace>& faces,std::size_t& N_connections) {
      bool result = true;

      // Sanity check on input parameters:
      if (refinements < 0) result = false;
      if (radius <= 0.0) result = false;

      // Exit if error(s) have occurred:
      if (result == false) return result;

      nodes.clear();
      faces.clear();

      // Generate a triangulated sphere using the chosen method:
      switch (method) {
       case spheregenerator::icosahedron:
	 result = createIcosahedron(nodes,faces);
	 break;
       case spheregenerator::octahedron:
	 result = createOctahedron(nodes,faces);
	 break;
       case spheregenerator::tetrahedron:
	 result = createTetrahedron(nodes,faces);
	 break;
       default:
	 // Unknown or unsupported method:
	 result = false;
	 break;
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

      // Shift centroid to given position and set sphere's radius to given one:
      normalizeRadiusAndShiftCentroid(centroid,radius,nodes);

      N_connections = 5*faces.size();
      return result;
   }

} // namespace ucdmesh
