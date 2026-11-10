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

#ifndef PARTICLE_LIST_UTIL_H
#define PARTICLE_LIST_UTIL_H

#include <vector>
#include <pargrid_definitions.h>
#include <pargrid_datawrapper.h>
#include <limits>

namespace particlelist {

   template<class PARTICLE>
   void removeParticles(pargrid::CellID blockLID,const std::vector<pargrid::ArraySizetype>& indicesToRemove,
			pargrid::DataWrapper<PARTICLE>& wrapper,unsigned int* N_particles);


   /** Remove particles from given block. Particles to be removed are defined by 
    * their positions in block's particle list, i.e., by their indices in the 
    * particle array.
    * @param blockLID Local ID of the block.
    * @param indicesToRemove Indices of particles to be removed.
    * @param wrapper ParGrid DataWrapper containing particles.
    * @param N_particles Array passed by ParticleList class.*/
   template<class PARTICLE> inline
   void removeParticles(pargrid::CellID blockLID,const std::vector<pargrid::ArraySizetype>& indicesToRemove,
			pargrid::DataWrapper<PARTICLE>& wrapper,unsigned int* N_particles) {
      // Exit immediately if there are no particles to remove:
      if (indicesToRemove.size() == 0) return;

      // Check if all particles are removed:
      if (indicesToRemove.size() == wrapper.size(blockLID)) {
	 wrapper.resize(blockLID,0);
	 N_particles[blockLID] = 0;
	 return;
      }

      std::vector<pargrid::ArraySizetype>::const_iterator begin = indicesToRemove.begin();
      std::vector<pargrid::ArraySizetype>::const_iterator end = indicesToRemove.end(); --end;

      // First hole is given by the first element in indicesToRemove. Make lastExisting 
      // point to last particle in list that is not removed:
      pargrid::ArraySizetype hole = *begin;
      pargrid::ArraySizetype lastExisting = wrapper.size(blockLID)-1;
      while (lastExisting == *end && end != indicesToRemove.begin()) {
	 --end;
	 --lastExisting;
      }

      // Fill holes (=particles to be removed) with particles from the end of list:
      while (begin != indicesToRemove.end() && *begin < lastExisting) {
	 // Copy last remaining particle to hole:
	 wrapper.data()[blockLID][hole] = wrapper.data()[blockLID][lastExisting];

	 #ifndef NDEBUG
	    wrapper.data()[blockLID][lastExisting].state[0] = std::numeric_limits<Real>::infinity();
	    wrapper.data()[blockLID][lastExisting].state[1] = std::numeric_limits<Real>::infinity();
	    wrapper.data()[blockLID][lastExisting].state[2] = std::numeric_limits<Real>::infinity();
	 #endif

	 ++begin;
	 hole = *begin;

	 // Make lastExisting point to last particle in list that is not removed:
	 --lastExisting;
	 while (lastExisting == *end && lastExisting > hole) {
	    --end;
	    --lastExisting;
	 }
      }

      // Reduce array sizes:
      const pargrid::ArraySizetype oldSize = wrapper.size(blockLID);
      const pargrid::ArraySizetype newSize = oldSize-indicesToRemove.size();
      wrapper.resize(blockLID,newSize);
      N_particles[blockLID] = newSize;
   }

} // namespace particlelist

#endif
