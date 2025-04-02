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
#include <unordered_map>
#include <limits>
#include <vlsv_common.h>

#include <ucd_mesh.h>

using namespace std;

static const uint32_t maxHash = 65535;

namespace ucdmesh {

   /** Test if two edges are equal.
    * This function is used in triangular face refinement algorithm where an edge, 
    * defined by indices of the Nodes it connects to, is refined unless it has 
    * already been refined at an earlier stage in the algorithm.
    * @param first First edge.
    * @param second Second edge.
    * @return If true, edges first and second are equal. Otherwise they are different.*/
   bool EdgesEqual::operator()(const pair<uint32_t,uint32_t>& first,const pair<uint32_t,uint32_t>& second) const {
      if (first.first != second.first) return false;
      if (first.second != second.second) return false;
      return true;
   }

   /** Calculate a hash value for given edge. An edge is defined by 
    * indices of Nodes it connects to.
    * @param nodeIndices Indices of Nodes defining an edge.
    * @return 64bit hash value of the edge.*/
   uint64_t EdgeHash::operator()(const pair<uint32_t,uint32_t>& nodeIndices) const {
      uint64_t result = 0;
      uint64_t tmp = nodeIndices.first % maxHash;
      result = (result | tmp);

      tmp = nodeIndices.second % maxHash;
      return (result | (tmp << 16));
   }

   /** Constructor for struct Node. Given x,y,z coordinate values are 
    * copied into corresponding member variables.
    * @param x Node's x-coordinate.
    * @param y Node's y-coordinate.
    * @param z Node's z-coordinate.*/
   Node::Node(Real x,Real y,Real z): x(x),y(y),z(z) { }

   /** Contructor for a triangulare face element.
    * @param cellType Cell type, one of the values defined in vlsv_common.h header (vlsv::celltype::TRIANGLE for triangle).
    * @param N_nodes Number of nodes this face element connects to (three for triangular).
    * @param node1 Index of the first node.
    * @param node2 Index of the second node.
    * @param node3 Index of the third node.*/
   TriangularFace::TriangularFace(uint32_t cellType,uint32_t N_nodes,uint32_t node1,uint32_t node2,uint32_t node3):
     cellType(cellType),N_nodes(N_nodes),node1(node1),node2(node2),node3(node3) { }

   /** Refine triangular faces. After this function exits, each triangle has 
    * been replaced by four new triangles.
    * @param nodes List of nodes on the triangulated surface.
    * @param faces List of faces making up the triangulated surface.*/
   void refineFaces(std::vector<Node>& nodes,std::vector<TriangularFace>& faces) {
      // Hash table where each refined edge is stored. An edge is defined by the 
      // indices of Nodes it connects to, this is the pair<uint32_t,uint32_t> template parameter.
      // Second template parameter (int) is the index of the new Node, i.e., its position in vector nodes.
      unordered_map<pair<uint32_t,uint32_t>,int,EdgeHash,EdgesEqual> newNodeIndices;

      // Iterator to unordered_map newNodeIndices. Algorithm below attempts to refine 
      // all edges of each triangle, and sets this iterator point to the new Node. In 
      // case the new Node already exists, iterator will point to that.
      pair<unordered_map<pair<uint32_t,uint32_t>,int,EdgeHash,EdgesEqual>::iterator,bool> result;

      const uint32_t invalidValue = numeric_limits<uint32_t>::max();
      uint32_t newNode1,newNode2,newNode3;
      const size_t N_existingTriangles = faces.size();

      for (size_t face=0; face<N_existingTriangles; ++face) {
	 const uint32_t node1 = faces[face].node1;
	 const uint32_t node2 = faces[face].node2;
	 const uint32_t node3 = faces[face].node3;

	 // Attempt to create three new nodes:
	 result = newNodeIndices.insert(make_pair(make_pair(min(node1,node2),max(node1,node2)),invalidValue));
	 if (result.second == true) {
	    newNode1 = nodes.size();
	    result.first->second = newNode1;
	    nodes.push_back(Node(0.5*(nodes[node1].x+nodes[node2].x),
				 0.5*(nodes[node1].y+nodes[node2].y),
				 0.5*(nodes[node1].z+nodes[node2].z)));
	 } else {
	    newNode1 = result.first->second;
	 }

	 result = newNodeIndices.insert(make_pair(make_pair(min(node2,node3),max(node2,node3)),invalidValue));
	 if (result.second == true) {
	    newNode2 = nodes.size();
	    result.first->second = newNode2;
	    nodes.push_back(Node(0.5*(nodes[node2].x+nodes[node3].x),
				 0.5*(nodes[node2].y+nodes[node3].y),
				 0.5*(nodes[node2].z+nodes[node3].z)));
	 } else {
	    newNode2 = result.first->second;
	 }

	 result = newNodeIndices.insert(make_pair(make_pair(min(node3,node1),max(node3,node1)),invalidValue));
	 if (result.second == true) {
	    newNode3 = nodes.size();
	    result.first->second = newNode3;
	    nodes.push_back(Node(0.5*(nodes[node3].x+nodes[node1].x),
				 0.5*(nodes[node3].y+nodes[node1].y),
				 0.5*(nodes[node3].z+nodes[node1].z)));
	 } else {
	    newNode3 = result.first->second;
	 }

	 // Insert new faces:
	 faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,newNode1,node2,newNode2));
	 faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,newNode2,node3,newNode3));
	 faces.push_back(TriangularFace(vlsv::celltype::TRIANGLE,3,newNode1,newNode2,newNode3));
	 faces[face].node2 = newNode1;
	 faces[face].node3 = newNode3;
      }
   }

} // namespace ucdmesh
