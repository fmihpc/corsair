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

#ifndef OPERATOR_PARGRID_ARRAY_H
#define OPERATOR_PARGRID_ARRAY_H

#include <vector>

#include <dataoperator.h>

class OperatorPargridArray: public DataOperator {
 public:
   OperatorPargridArray();
   ~OperatorPargridArray();

   bool finalize();
   std::string getName() const;
   bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   virtual bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists);

 private:

   /** Struct that contains all information that is required 
    * to correctly write a ParGrid array to output file.*/
   struct ArrayInfo {
      unsigned int byteSize;            /**< Byte size of single data vector element.*/
      pargrid::CellID dataID;           /**< ParGrid DataID of array.*/
      std::string datatype;             /**< String representation of the datatype, either "int", "uint", or "float".*/
      std::string name;                 /**< Output name of the array.*/
      std::string units;                /**< Physical units of variable in array.*/
      unsigned int vectorSize;          /**< Size of data vector.*/
   };

   bool arrayInfoRead;                         /**< If true, ParGrid array information for arrays defined in 
						* vector arrayNames have been read.*/
   std::vector<ArrayInfo> arrayInfo;           /**< Vector containing ArrayInfo for each ParGrid array that is written 
						* to output file by this DataOperator.*/
   std::vector<std::string> arrayOutputNames;  /**< Name of each array in output file, must be compatible with VisIt.
						May contain variable path.*/
   std::vector<std::string> arrayPargridNames; /**< Names of ParGrid arrays that are written to output files.*/   
   std::vector<std::string> arrayUnits;

   #if PROFILE_LEVEL > 0
      int profTotalTime;                       /**< Total wallclock time used by this DataOperator.*/
   #endif

   bool addConfigFileItems(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
   bool readArrayInfo();
};

#endif
