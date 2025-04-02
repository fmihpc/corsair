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

#ifndef OPERATOR_PITCH_BINNER_H
#define OPERATOR_PITCH_BINNER_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include <dataoperator.h>

template<class FIELD,class SPECIES,class PARTICLE>
class PitchBinner: public DataOperator {
 public:
   PitchBinner();
   ~PitchBinner();

   bool finalize();
   std::string getName() const;
   bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   virtual bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);

 private:
   bool divideByBinWidth;           /**< If true, value in each bin is divided by bin width
				     * before data is written.*/
   int N_bins;                      /**< Number of bins used for [-1,+1] pitch interval.*/
   int pitchBinnerID;

   int blockIndex(int i,int j,int k);
};

template<class FIELD,class SPECIES,class PARTICLE>
PitchBinner<FIELD,SPECIES,PARTICLE>::PitchBinner(): DataOperator() { 
   pitchBinnerID = -1;
}

template<class FIELD,class SPECIES,class PARTICLE>
PitchBinner<FIELD,SPECIES,PARTICLE>::~PitchBinner() { }

template<class FIELD,class SPECIES,class PARTICLE>
int PitchBinner<FIELD,SPECIES,PARTICLE>::blockIndex(int i,int j,int k) {return k*block::WIDTH_Y*block::WIDTH_X+j*block::WIDTH_X+i;}

template<class FIELD,class SPECIES,class PARTICLE>
bool PitchBinner<FIELD,SPECIES,PARTICLE>::finalize() {return true;}

template<class FIELD,class SPECIES,class PARTICLE>
std::string PitchBinner<FIELD,SPECIES,PARTICLE>::getName() const {return "pitch";}

template<class FIELD,class SPECIES,class PARTICLE>
bool PitchBinner<FIELD,SPECIES,PARTICLE>::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   if (DataOperator::initialize(cr,sim,simClasses) == false) {
      simClasses.logger << "(PITCHBINNER) ERROR: DataOperator initialization failed!" << std::endl << write;
      return false;
   }
   bool success = true;

   // Parse configuration file:
   std::string divide;
   cr.add("PitchBinner.bins","Number of bins used for [-1,+1] pitch interval (integer).",(int)9);
   cr.add("PitchBinner.divide_by_bin_width","Should values in bins be divided by bin widths (yes/no) ?","");
   cr.parse();
   cr.get("PitchBinner.divide_by_bin_width",divide);
   cr.get("PitchBinner.bins",N_bins);

   // Sanity check on values:
   if (N_bins < 2) {
      simClasses.logger << "(PITCHBINNER) ERROR: Parameter N_bins >= 1" << std::endl << write;
      success = false;
   }
   if (divide == "yes") divideByBinWidth = true;
   else if (divide == "no") divideByBinWidth = false;
   else {
      simClasses.logger << "(PITCHBINNER) ERROR: Parameter 'divide_by_bin_width' must have value 'yes' or 'no'." << std::endl << write;
      success = false;
   }
   return success;
}

template<class FIELD,class SPECIES,class PARTICLE>
bool PitchBinner<FIELD,SPECIES,PARTICLE>::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
   if (getInitialized() == false) return false;
   profile::start("PitchBinner",pitchBinnerID);
   bool success = true;

   // Allocate arrays which are written to output file:
   const Real pitch_min = -1.0;
   const Real pitch_max = +1.0;
   const Real d_pitch = (pitch_max - pitch_min) / N_bins;

   // Get list of local cells. The list is ordered.
   pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();

   // Allocate arrays:
   Real** buffer = new Real*[N_bins];
   for (int i=0; i<N_bins; ++i) buffer[i] = new Real[N_localBlocks*block::SIZE];
   const double* blockCrds = getBlockCoordinateArray(*sim,*simClasses);

   for (size_t species=0; species<particleLists.size(); ++species) {
      // This operator can only bin guiding center "GC" particles:
      pargrid::DataID particleDataID = pargrid::INVALID_DATAID;
      if (particleLists[species]->getParticles(particleDataID) == false) continue;
      if (particleLists[species]->getSpeciesType() != "GC") continue;

      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE** particles = wrapper.data();
      const SPECIES* particleSpecies = reinterpret_cast<const SPECIES*>(particleLists[species]->getSpecies());

      for (size_t block=0; block<N_localBlocks; ++block) {
	 // Clear output arrays:
	 for (int i=0; i<N_bins; ++i) for (int j=0; j<block::SIZE; ++j) buffer[i][block*block::SIZE+j] = 0.0;

	 Real cellSize[3];
	 getBlockCellSize(*simClasses,*sim,block,cellSize);

	 for (size_t p=0; p<wrapper.size(block); ++p) {
	    // Get fields at particle position:
	    Real E[3];
	    Real B[3];
	    Real pos[3];
	    for (int j=0; j<3; ++j) pos[j] = blockCrds[3*block+j] + particles[block][p].state[j];
	    FIELD::getFields(block,sim->t,pos,E,B);
	    const Real B_mag = vectorMagnitude<3>(B);

	    // Calculate pitch and accumulate to output arrays:
	    //const Real V_gyro_2 = 2.0*particles[block][p].state[GC::MU]*B_mag/particleLists[species]->getSpecies().m;
	    const Real V_gyro_2 = 2.0*particles[block][p].state[GC::MU]*B_mag/particleSpecies->m;
	    const Real V_star   = sqrt(particles[block][p].state[GC::VPAR]*particles[block][p].state[GC::VPAR] + V_gyro_2);
	    const Real pitch    = particles[block][p].state[GC::VPAR] / V_star;

	    int bin = static_cast<int>(floor((pitch - pitch_min) / d_pitch));
	    if (bin < 0) bin = 0;
	    else if (bin > N_bins-1) bin = N_bins-1;

	    // Calculate index into cells in block:
	    const int i = static_cast<int>(floor(particles[block][p].state[GC::X]/cellSize[0]));
	    const int j = static_cast<int>(floor(particles[block][p].state[GC::Y]/cellSize[1]));
	    const int k = static_cast<int>(floor(particles[block][p].state[GC::Z]/cellSize[2]));

	    buffer[bin][block*block::SIZE+blockIndex(i,j,k)] += particles[block][p].state[GC::WEIGHT];
	 }

	 // Divide values by cell size to make units 1/m^3:
	 Real factor = 1.0 / (cellSize[0]*cellSize[1]*cellSize[2]);

	 // Divide by bin width if necessary:
	 if (divideByBinWidth == true) {
	    factor /= d_pitch;
	 }

	 for (int i=0; i<N_bins; ++i) for (int j=0; j<block::SIZE; ++j) buffer[i][block*block::SIZE+j] *= factor;
      }

      // Write pitch arrays to file:
      for (int i=0; i<N_bins; ++i) {
	 std::stringstream ss;
	 ss << "bin" << i;

	 std::map<std::string,std::string> attributes;
	 attributes["mesh"] = spatMeshName;
	 attributes["type"] = "celldata";
	 attributes["name"] = ss.str().c_str();
	 attributes["dir"] = "/PitchBinner/" + particleLists[species]->getName();
	 if (divideByBinWidth == false) attributes["unit"] = "1/m3";
	 else attributes["unit"] = "1/(m3 rad)";
	 const uint64_t arraySize = N_localBlocks*block::SIZE;
	 const uint64_t vectorSize = 1;
	 if (simClasses->vlsv.writeArray("VARIABLE",attributes,arraySize,vectorSize,buffer[i]) == false) success = false;
      }
   }

   // Deallocate memory and exit:
   for (int i=0; i<N_bins; ++i) {delete [] buffer[i]; buffer[i] = NULL;}
   delete [] buffer; buffer = NULL;
   profile::stop();
   return success;
}

#endif
