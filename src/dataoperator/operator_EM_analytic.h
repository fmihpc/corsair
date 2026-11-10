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

#ifndef OPERATOR_EM_ANALYTIC_H
#define OPERATOR_EM_ANALYTIC_H

#include <map>
#include <vector>

#include <dataoperator.h>

template<class C>
class EMAnalyticOperator: public DataOperator {
 public:
   EMAnalyticOperator();
   ~EMAnalyticOperator();

   bool finalize();
   std::string getName() const;
   bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists);

 private:
   #if PROFILE_LEVEL > 0
      int emAnalyticID;
   #endif
};

template<class C>
EMAnalyticOperator<C>::EMAnalyticOperator(): DataOperator() { 
   #if PROFILE_LEVEL > 0
     emAnalyticID = -1;
   #endif
}

template<class C>
EMAnalyticOperator<C>::~EMAnalyticOperator() {finalize();}

template<class C>
bool EMAnalyticOperator<C>::finalize() {
   if (initialized == false) return true;
   initialized=false;
   return true;
}  

template<class C>
std::string EMAnalyticOperator<C>::getName() const {return "EManalytic";}

template<class C>
bool EMAnalyticOperator<C>::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   initialized = true;

   if (DataOperator::initialize(cr,sim,simClasses) == false) {
      simClasses.logger << "(EManalytic) ERROR: Field function (given as a template parameter) failed to initialize!" << std::endl << write;
      initialized = false;
   }

   return initialized;
}

template<class C>
bool EMAnalyticOperator<C>::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
   if (getInitialized() == false) return false;
   bool success = true;

   #if PROFILE_LEVEL > 0
      profile::start("EMAnalytic",emAnalyticID);
   #endif

   // Get local cell IDs:
   const pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();  

   // Allocate arrays for E,B vectors:
   Real* B = new Real[3*N_localBlocks*block::SIZE];
   Real* E = new Real[3*N_localBlocks*block::SIZE];
   const double* blockCrds = getBlockCoordinateArray(*sim,*simClasses);

   // Calculate E & B vectors at local cell centroids:
   size_t index = 0;
   for (size_t block=0; block<N_localBlocks; ++block) {
      Real cellSize[3];
      getBlockCellSize(*simClasses,*sim,block,cellSize);

      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	 // Sample E,B fields in the cell with N_samples in each coordinate direction.
	 // Samples are divided in such a way that with one sample dx_sample=dx_cell,
	 // with two samples dx_sample=dx_cell/2, and so forth:
	 const int N_samples = 3;
	 const int N_samples3 = N_samples*N_samples*N_samples;
	 Real E_tmp[3];
	 Real B_tmp[3];
	 Real position[3];
	 const Real dx_sample = cellSize[0]/N_samples;
	 const Real dy_sample = cellSize[1]/N_samples;
	 const Real dz_sample = cellSize[2]/N_samples;

	 for (int n=0; n<3; ++n) E[3*index+n] = 0.0;
	 for (int n=0; n<3; ++n) B[3*index+n] = 0.0;
	 for (int kk=0; kk<N_samples; ++kk) for (int jj=0; jj<N_samples; ++jj) for (int ii=0; ii<N_samples; ++ii) {
	    position[0] = blockCrds[3*block+0] + i*cellSize[0] + (ii+0.5)*dx_sample;
	    position[1] = blockCrds[3*block+1] + j*cellSize[1] + (jj+0.5)*dy_sample;
	    position[2] = blockCrds[3*block+2] + k*cellSize[2] + (kk+0.5)*dz_sample;

	    C::getFields(block,sim->t,position,E_tmp,B_tmp);
	    for (int n=0; n<3; ++n) E[3*index+n] += E_tmp[n];
	    for (int n=0; n<3; ++n) B[3*index+n] += B_tmp[n];
	 }
	 for (int n=0; n<3; ++n) E[3*index+n] /= N_samples3;
	 for (int n=0; n<3; ++n) B[3*index+n] /= N_samples3;
	 ++index;
      }
   }

   // Write arrays to output file:
   std::map<std::string,std::string> attribs;
   attribs["name"] = "E";
   attribs["mesh"] = spatMeshName;
   attribs["type"] = "celldata";
   attribs["unit"] = "V/m";
   const uint64_t arraySize = N_localBlocks*block::SIZE;
   const uint64_t vectorSize = 3;
   if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,E) == false) success = false;

   attribs["name"] = "B";
   attribs["unit"] = "T";
   if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,B) == false) success = false;

   // Deallocate memory and exit:
   delete [] E; E = NULL;
   delete [] B; B = NULL;

   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}

#endif
