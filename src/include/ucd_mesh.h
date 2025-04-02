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

#ifndef UCD_MESH_H
#define UCD_MESH_H

#include <vector>
#include <stdint.h>

#include <definitions.h>

/** Using createSphere function declared in this file a user can create triangulated spherical 
 * surfaces. The generation method is fairly standard, i.e., we start with a triangulated polyhedron 
 * such as a tetrahedron, refine the triangles for chosen number of times, and renormalize the 
 * node coordinates to lie on a spherical surface.
 */

namespace ucdmesh {

   namespace facedataelement {
      enum Element {
	 AREA,
	 NORMAL_X,
	 NORMAL_Y,
	 NORMAL_Z,
	 SIZE
      };
   } // namespace facedataelement

   /** Coordinates of a node on triangulated surface.*/
   struct Node {
      Node(Real x,Real y,Real z);

      Real x; /**< Node's x-coordinate.*/
      Real y; /**< Node's y-coordinate.*/
      Real z; /**< Node's z-coordinate.*/
   };

   /** Triangular face element. This struct contains all necessary information 
    * that is required to write the ucd mesh to VLSV file(s).*/
   struct TriangularFace {
      TriangularFace(uint32_t cellType,uint32_t N_nodes,uint32_t node1,uint32_t node2,uint32_t node3);

      uint32_t cellType; /**< Type of cell, one of the pre-defined types in vlsv_common.h file.*/
      uint32_t N_nodes;  /**< Number of nodes the cell is connected to.*/
      uint32_t node1;    /**< Index of first node.*/
      uint32_t node2;    /**< Index of second node.*/
      uint32_t node3;    /**< Index of third node.*/
   };

   /** Wrapper struct for testing if two edges are equal, i.e., if they connect to same Nodes.*/
   struct EdgesEqual {
      bool operator()(const std::pair<uint32_t,uint32_t>& first,const std::pair<uint32_t,uint32_t>& second) const;
   };

   /** Wrapper struct for calculating a hash value for an edge.
    * Calculated value is based on indices of nodes the edge connects to.*/
   struct EdgeHash {
      /** Calculate a hash value for given edge. An edge is defined by 
       * indices of Nodes it connects to.
       * @param nodeIndices Indices of Nodes defining an edge.
       * @return 64bit hash value of the edge.*/
      uint64_t operator()(const std::pair<uint32_t,uint32_t>& nodeIndices) const;
   };

   void refineFaces(std::vector<Node>& nodes,std::vector<TriangularFace>& faces);

} // namespace ucdmesh

#endif
