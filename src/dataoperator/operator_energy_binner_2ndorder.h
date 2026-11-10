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

#ifndef OPERATOR_ENERGY_BINNER_2NDORDER_H
#define OPERATOR_ENERGY_BINNER_2NDORDER_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include <simulationclasses.h>
#include <dataoperator.h>
#include <accumulators.h>

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
class EnergyBinner2nd: public DataOperator {
 public:
   EnergyBinner2nd();
   ~EnergyBinner2nd();

   bool finalize();
   std::string getName() const;
   bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   virtual bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists);

 private:

   struct Instrument {
      std::vector<std::string> channelNames;        /**< Name of each energy channel.*/
      bool divideByBinWidth;                        /**< If true, measured counts are divided by channel bin
						     *                                                      * width. Physical units will be 1/m^3/eV instead of 1/m^3.*/
      bool energyPerNucleon;                        /**< If true, energy limits in minValues, maxValues are
						     *                                                      * per nucleon (energy * proton mass / particle mass).*/
      std::vector<Real> binWidth;                   /**< Energy bin widths in eV.*/
      std::vector<Real> maxValues;                  /**< Maximum energy for each channel.*/
      std::vector<Real> minValues;                  /**< Minimum energy for each channel.*/
      std::string name;                             /**< Name of the Instrument.*/
      std::vector<size_t> speciesIndices;           /**< For each species the corresponding index into particle lists.*/
      std::vector<std::string> speciesNames;        /**< Names of accumulated particle species.*/
   };

   #if PROFILE_LEVEL > 0
      int accumulateID;
      int energyBinnerID;
      int dataWriteID;
      int mpiOverheadID;
      int mpiWaitID;
      std::string profileName;
   #endif

   std::vector<Instrument> instruments;

   int accCellIndex(int i,int j,int k) const;
   int cellIndex(int i,int j,int k) const;
   void accumulateCell(size_t instr,size_t channel,pargrid::CellID cell,const std::vector<ParticleListBase*>& particleLists,Real* const data);
   bool getScaleFactor(const std::string& energyUnit,Real& energyScaling) const;
};

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
EnergyBinner2nd<FIELD,SPECIES,PARTICLE,ORDER>::EnergyBinner2nd(): DataOperator() { 
   #if PROFILE_LEVEL > 0
      accumulateID = -1;
      energyBinnerID = -1;
      dataWriteID = -1;
      mpiOverheadID = -1;
      mpiWaitID = -1;
   #endif
}

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
EnergyBinner2nd<FIELD,SPECIES,PARTICLE,ORDER>::~EnergyBinner2nd() { }

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
void EnergyBinner2nd<FIELD,SPECIES,PARTICLE,ORDER>::accumulateCell(size_t instr,size_t channel,pargrid::CellID block,
								   const std::vector<ParticleListBase*>& particleLists,Real* const data) {
   Real E[3];
   Real B[3];
   Real V_drift[3];

   // This is a second-order accumulator, so we need two extra cells per coordinate direction:
   const int accBlockSize  = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real acc[accBlockSize];
   for (int i=0; i<accBlockSize; ++i) acc[i] = 0.0;

   const double* blockCrds = getBlockCoordinateArray(*sim,*simClasses);
   Real cellSize[3];
   getBlockCellSize(*simClasses,*sim,block,cellSize);

   // Accumulate all particle species that contribute to this instrument's measurements:
   for (size_t s=0; s<instruments[instr].speciesNames.size(); ++s) {
      // Find particle list propagating the species with name speciesNames[s]:
      ParticleListBase* particleList = NULL;
      for (size_t i=0; i<particleLists.size(); ++i) {
	 if (particleLists[i]->getName() == instruments[instr].speciesNames[s]) {
	    particleList = particleLists[i];
	    break;
	 }
      }
      if (particleList == NULL) continue;

      // Get particle species:
      const SPECIES* const species = reinterpret_cast<const SPECIES*>(particleList->getSpecies());

      // Either use the species' mass or proton mass when calculating energy,
      // depending on whether to count energies per nucleon or not:
      Real mass = species->m;
      if (instruments[instr].energyPerNucleon == true) mass = constants::MASS_PROTON;

      // Get pointer to particles:
      pargrid::DataID particleDataID = pargrid::INVALID_DATAID;
      if (particleList->getParticles(particleDataID) == false) {
	 simClasses->logger << "(ENERGYBINNER) ERROR: Failed to get particle data ID!" << std::endl << write;
	 continue;
      }
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE** const particles = wrapper.data();

      // Accumulate particles:
      for (unsigned int p=0; p<wrapper.size(block); ++p) {
	 Real pos[3];
	 pos[GC::X] = blockCrds[3*block+GC::X] + particles[block][p].state[GC::X];
	 pos[GC::Y] = blockCrds[3*block+GC::Y] + particles[block][p].state[GC::Y];
	 pos[GC::Z] = blockCrds[3*block+GC::Z] + particles[block][p].state[GC::Z];
	 FIELD::getFields(block,sim->t,pos,E,B);
	 const Real B_mag2 = vectorMagnitude2<3>(B);

	 // Calculate particle's non-relativistic energy:
	 crossProduct(E,B,V_drift);
	 for (int j=0; j<3; ++j) V_drift[j] /= B_mag2;
	 const Real v_mag2 = particles[block][p].state[GC::VPAR]*particles[block][p].state[GC::VPAR] + vectorMagnitude2<3>(V_drift);
	 const Real energy = particles[block][p].state[GC::MU]*sqrt(B_mag2) + 0.5*mass*v_mag2;
	 if (energy < instruments[instr].minValues[channel] || energy > instruments[instr].maxValues[channel]) continue;

	 const Real weight = particles[block][p].state[GC::WEIGHT];
	 switch (ORDER) {
	  case 0:
	    accumulateScalarNGP_XYZ(cellSize,acc,weight,particles[block][p]);
	    break;
	  case 1:
	    accumulateScalarCIC_XYZ(cellSize,acc,weight,particles[block][p]);
	    break;
	  case 2:
	    accumulateScalarTSC_XYZ(cellSize,acc,weight,particles[block][p]);
	    break;
	  default:
	    accumulateScalarCIC_XYZ(cellSize,acc,weight,particles[block][p]);
	    break;
	 }
      }
   }

   // All particles in this cell have been accumulated. Add values from
   // array acc to this cell and its existing neighbours:
   block::addValues3D(*simClasses,block,acc,data);
}

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
bool EnergyBinner2nd<FIELD,SPECIES,PARTICLE,ORDER>::finalize() {return true;}

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
std::string EnergyBinner2nd<FIELD,SPECIES,PARTICLE,ORDER>::getName() const {
   std::stringstream ss;
   ss << "EnergyBinner" << ORDER;
   return ss.str();
}

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
bool EnergyBinner2nd<FIELD,SPECIES,PARTICLE,ORDER>::getScaleFactor(const std::string& energyUnit,Real& energyScaling) const {
   bool success = true;
   if (energyUnit == "eV") {
      energyScaling = constants::CHARGE_ELEMENTARY;
   } else if (energyUnit == "keV") {
      energyScaling = 1.0e3 * constants::CHARGE_ELEMENTARY;
   } else if (energyUnit == "MeV") {
      energyScaling = 1.0e6 * constants::CHARGE_ELEMENTARY;
   } else if (energyUnit == "GeV") {
      energyScaling = 1.0e9 * constants::CHARGE_ELEMENTARY;
   } else {
      success = false;
   }
   return success;
}

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
bool EnergyBinner2nd<FIELD,SPECIES,PARTICLE,ORDER>::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   if (DataOperator::initialize(cr,sim,simClasses) == false) return false;

   bool success = true;

   // Read the number of energy binners from config file, and check the obtained value for sanity:
   int N_binners;
   std::string energyUnitIn = "";
   cr.add("EnergyBinner.binners","Number of energy binner instruments (int).",(int)0);
   cr.parse();
   cr.get("EnergyBinner.binners",N_binners);

   if (N_binners < 0) {
      simClasses.logger << "(ENERGYBINNER) ERROR: parameter 'binners' has negative value!" << std::endl << write;
      return false;
   }
   if (N_binners == 0) return true;

   // Read parameters of each instrument from config file:
   instruments.resize(N_binners);
   for (size_t i=0; i<instruments.size(); ++i) {
      std::stringstream ss;
      ss << '.' << i << '.';
      std::string middix;              // Since this will be in the middle, it can't be prefix or suffix :)
      ss >> middix;

      energyUnitIn = "";
      cr.add("EnergyBinner"+middix+"divide_by_bin_width","Should counts be divided by energy channel width (yes/no) ?","");
      cr.add("EnergyBinner"+middix+"name","Name of the instrument (string).","");
      cr.add("EnergyBinner"+middix+"energy_unit","Unit in which channel limits are given (eV/keV/MeV/GeV).",energyUnitIn);
      cr.add("EnergyBinner"+middix+"energy_per_amu","Should energy be calculated per nucleon, i.e. using proton energy (yes/no) ?","");
      cr.addComposed("EnergyBinner"+middix+"channels","Minimum and maximum energy for an energy channel, in units specified with 'energy_unit' (float,multiple).");
      cr.addComposed("EnergyBinner"+middix+"channel_names","Name of each energy channel (string,multiple).");
      cr.addComposed("EnergyBinner"+middix+"particle_species","Names of particle species this instrument bins (string,multiple).");
   }
   cr.parse();

   for (size_t i=0; i<instruments.size(); ++i) {
      std::stringstream ss;
      ss << '.' << i << '.';
      std::string middix;
      ss >> middix;

      std::string name,energyPerAmu,divideByBinWidth;
      std::vector<std::string> channels;
      std::vector<std::string> channelNames;
      std::vector<std::string> speciesNames;
      if (cr.get("EnergyBinner"+middix+"divide_by_bin_width",divideByBinWidth) == false) success = false;
      if (cr.get("EnergyBinner"+middix+"name",name) == false) success = false;
      if (cr.get("EnergyBinner"+middix+"energy_unit",energyUnitIn) == false) success = false;
      if (cr.get("EnergyBinner"+middix+"energy_per_amu",energyPerAmu) == false) success = false;
      if (cr.get("EnergyBinner"+middix+"channels",channels) == false) success = false;
      if (cr.get("EnergyBinner"+middix+"channel_names",channelNames) == false) success = false;
      if (cr.get("EnergyBinner"+middix+"particle_species",speciesNames) == false) success = false;
      if (success == false) break;

      // Remove empty entries from speciesNames:
      do {
	 bool finished = true;
	 for (std::vector<std::string>::iterator it=speciesNames.begin(); it!=speciesNames.end(); ++it) {
	    if (*it == "") {
	       speciesNames.erase(it);
	       finished = false;
	       break;
	    }
	 }
	 if (finished == true) break;
      } while (true);

      // Check obtained values for sanity:
      if (channels.size() != channelNames.size()) {
	 simClasses.logger << "(ENERGYBINNER) ERROR: Fields '" << "EnergyBinner"+middix+"channels" << "' and '";
	 simClasses.logger << "EnergyBinner"+middix+"channel_names" << "' have different number of entries!" << std::endl << write;
	 success = false;
      }

      instruments[i].name = name;
      instruments[i].channelNames = channelNames;
      instruments[i].minValues.resize(channels.size());
      instruments[i].maxValues.resize(channels.size());
      instruments[i].speciesIndices.resize(speciesNames.size());
      instruments[i].speciesNames = speciesNames;
      instruments[i].binWidth.resize(channels.size());

      // Determine input/output energy units:
      Real energyScalingInput;
      if (getScaleFactor(energyUnitIn,energyScalingInput) == false) {
	 simClasses.logger << "(ENERGYBINNER) ERROR: Input energy unit '" << energyUnitIn << "' given in parameter ";
	 simClasses.logger << "'EnergyBinner" << middix << "energy_unit'" << std::endl;
	 simClasses.logger << "does not correspond to any known unit. Valid values are eV,keV,MeV,GeV." << std::endl << write;
	 return false;
      }

      // Determine if energy should be binned per nucleon:
      if (energyPerAmu == "yes") instruments[i].energyPerNucleon = true;
      else if (energyPerAmu == "no") instruments[i].energyPerNucleon = false;
      else {
	 simClasses.logger << "(ENERGYBINNER) ERROR: Parameter 'EnergyBinner" << middix << "energy_per_amu' may only have value 'yes' or 'no'." << std::endl << write;
	 return false;
      }

      // Determine if counts should be divided by bin widths:
      if (divideByBinWidth == "yes") instruments[i].divideByBinWidth = true;
      else if (divideByBinWidth == "no") instruments[i].divideByBinWidth = false;
      else {
	 simClasses.logger << "(ENERGYBINNER) ERROR: Parameter 'EnergyBinner" << middix << "divide_by_bin_width' may only have value 'yes' or 'no'." << std::endl << write;
	 return false;
      }

      for (size_t j=0; j<channels.size(); ++j) {
	 // Check that each "channels" entry contains minimum and maximum energy, separated by a whitespace:
	 const size_t separator = channels[j].find(' ');
	 if (separator == std::string::npos) {
	    simClasses.logger << "(ENERGYBINNER) ERROR: channel entry '" << channels[j] << "' does not contain " << std::endl;
	    simClasses.logger << "minimum and maximum energies separated by whitespace" << std::endl << write;
	    success = false;
	 }
	 // Store min/max energy limits in simulation units:
	 instruments[i].minValues[j] = energyScalingInput * atof(channels[j].substr(0,separator).c_str());
	 instruments[i].maxValues[j] = energyScalingInput * atof(channels[j].substr(separator+1,channels[j].size()-separator).c_str());
	 instruments[i].binWidth[j]  = (instruments[i].maxValues[j] - instruments[i].minValues[j])/constants::CHARGE_ELEMENTARY;
      }

      // Exit if errors have come up:
      if (success == false) break;
   }

   #if PROFILE_LEVEL > 0
      std::stringstream ss;
      ss << "EnergyBinner" << ORDER;
      ss >> profileName;
   #endif   
   return success;
}

template<class FIELD,class SPECIES,class PARTICLE,int ORDER>
bool EnergyBinner2nd<FIELD,SPECIES,PARTICLE,ORDER>::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
   if (getInitialized() == false) return false;
   bool success = true;
   #if PROFILE_LEVEL > 0
      profile::start(profileName,energyBinnerID);
      profile::start("MPI overhead",mpiOverheadID);
   #endif

   // Allocate array for accumulation:
   const pargrid::CellID N_blocks = simClasses->pargrid.getNumberOfAllCells();
   const std::vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(0);

   // Create a ParGrid data array that is used in accumulation:
   const pargrid::DataID arrayID = simClasses->pargrid.addUserData<Real>("eb2nd",block::SIZE);
   if (arrayID == simClasses->pargrid.invalidDataID()) {
      simClasses->logger << "(ENERGYBINNER 2ND) ERROR: Failed to create a ParGrid data array!" << std::endl << write;
      return false;
   }

   Real* const data = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(arrayID));

   // Create a ParGrid Stencil for transferring array 'data':
   std::vector<pargrid::NeighbourID> recvNbrTypeIDs;
   recvNbrTypeIDs.reserve(26);
   for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
      const pargrid::NeighbourID nbrTypeID = simClasses->pargrid.calcNeighbourTypeID(i,j,k);
      if (nbrTypeID == 13) continue;
      recvNbrTypeIDs.push_back(nbrTypeID);
   }
   const pargrid::StencilID stencilID = simClasses->pargrid.addStencil(pargrid::remoteToLocalUpdates,recvNbrTypeIDs);

   // Tell ParGrid that 'array' is transferred with the stencil created above:
   simClasses->pargrid.addDataTransfer(arrayID,stencilID);
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif

   // Write all instrument data:
   for (size_t instr=0; instr<instruments.size(); ++instr) {
      // Write data from all channels:
      for (size_t channel=0; channel<instruments[instr].minValues.size(); ++channel) {
	 // Clear data array:
	 for (size_t i=0; i<N_blocks*block::SIZE; ++i) data[i] = 0.0;

	 // Accumulate particle data from boundary cells:
	 #if PROFILE_LEVEL > 0
	    profile::start("accumulation",accumulateID);
	 #endif
	 for (size_t block=0; block<boundaryBlocks.size(); ++block) {
	    accumulateCell(instr,channel,boundaryBlocks[block],particleLists,data);
	 }
	 #if PROFILE_LEVEL > 0
	    profile::stop();
	    profile::start("MPI overhead",mpiOverheadID);
	 #endif

	 // Post sends:
	 simClasses->pargrid.startNeighbourExchange(stencilID,arrayID);
	 #if PROFILE_LEVEL > 0
	    profile::stop();
	    profile::start("accumulation",accumulateID);
	 #endif

	 // Accumulate particle data from inner cells:
	 const std::vector<pargrid::CellID>& innerBlocks = simClasses->pargrid.getInnerCells(0);
	 for (size_t block=0; block<innerBlocks.size(); ++block) {
	    accumulateCell(instr,channel,innerBlocks[block],particleLists,data);
	 }
	 #if PROFILE_LEVEL > 0
	    profile::stop();
	    profile::start("MPI waits",mpiWaitID);
	 #endif

	 // Wait for incoming data:
	 simClasses->pargrid.wait(stencilID,arrayID);
	 #if PROFILE_LEVEL > 0
	    profile::stop();
	    profile::start("accumulation",accumulateID);
	 #endif

	 // Sum contributions from remote processes into boundary cell data:
	 unsigned int* offsets = NULL;
	 char* tmpBuffers = NULL;
	 Real* buffers = NULL;
	 simClasses->pargrid.getRemoteUpdates(stencilID,arrayID,offsets,tmpBuffers);
	 buffers = reinterpret_cast<Real*>(tmpBuffers);
	 for (size_t block=0; block<boundaryBlocks.size(); ++block) {
	    const pargrid::CellID blockLID = boundaryBlocks[block];
	    for (unsigned int b=offsets[block]; b<offsets[block+1]; ++b) {
	       for (int j=0; j<block::SIZE; ++j) data[blockLID*block::SIZE+j] += buffers[b*block::SIZE + j];
	    }
	 }

	 // Divide values by cell size to make units 1/m^3:
	 Real cellSize[3];
	 getBlockCellSize(*simClasses,*sim,0,cellSize);
	 Real factor = 1.0 / (cellSize[0]*cellSize[1]*cellSize[2]);

	 // Divide values in buffer by energy bin width (if necessary):
	 if (instruments[instr].divideByBinWidth == true) {
	    factor /= instruments[instr].binWidth[channel];
	 }
	 for (size_t i=0; i<simClasses->pargrid.getNumberOfLocalCells()*block::SIZE; ++i) data[i] *= factor;
	 #if PROFILE_LEVEL > 0
	    profile::stop();
	    profile::start("data writing",dataWriteID);
	 #endif

	 std::stringstream dirName;
	 dirName << "/EnergyBinner" << ORDER << '/' << instruments[instr].name;

	 std::stringstream channelName;
	 channelName << instruments[instr].channelNames[channel] << '_' << ORDER;

	 // Write data to file:
	 std::map<std::string,std::string> attribs;
	 attribs["name"] = channelName.str();
	 attribs["mesh"] = spatMeshName;
	 attribs["type"] = "celldata";
	 attribs["dir"] = dirName.str();
	 //attribs["dir"] = "/EnergyBinner2nd/" + instruments[instr].name;
	 if (instruments[instr].divideByBinWidth == true) attribs["unit"] = "1/(m3 eV)";
	 else attribs["unit"] = "1/m3";
	 const uint64_t arraySize = simClasses->pargrid.getNumberOfLocalCells()*block::SIZE;
	 const uint64_t vectorSize = 1;
	 if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,&(data[0])) == false) success = false;
	 #if PROFILE_LEVEL > 0
	    profile::stop();
	 #endif
      }
   }

   // Deallocate memory:
   simClasses->pargrid.removeStencil(stencilID);
   simClasses->pargrid.removeUserData(arrayID);
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}

#endif
