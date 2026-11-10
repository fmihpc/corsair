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

#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <mpi.h>
#include <stdint.h>
#include <vector>

namespace configreader {
   template<int SIZE,typename C>
   bool checkInputValues(C* values,const C comparedValue);
}

/** ConfigReader is a class meant for reading [option,value] pairs from configuration 
 * files, command line and environment variables. Typically these options are used 
 * as parameters for the program.
 * ConfigReader has to be initialized by calling ConfigReader::initialize, and 
 * similarly finalized by calling ConfigReader::finalize. It is erroneous to call 
 * other member functions before initialization or after finalization.
 * The intended way of using ConfigReader, after initialization, is to add requested 
 * values with ConfigReader::add member functions. After requested options have been added, 
 * ConfigReader::parse must be called for ConfigReader to actually (re)read configuration 
 * file(s) and fetch values for specified options. After parsing the obtained values 
 * may be queried with ConfigReader::get functions.
 * ConfigReader is a parallel reader. Option values are broadcasted to all processes 
 * after their values have been parsed.
 * ConfigReader uses two Boost libraries, boost/lexical_cast.hpp and boost/program_options.hpp. 
 * These libraries must be included when compiling ConfigReader. Additionally, library 
 * boost_program_options must be linked to the main program.
 */
class ConfigReader {
 public:
   ConfigReader();
   ~ConfigReader();

   bool add(const std::string& name,const std::string& descr,const bool& defValue);
   bool add(const std::string& name,const std::string& descr,const float& defValue);
   bool add(const std::string& name,const std::string& descr,const double& defValue);
   bool add(const std::string& name,const std::string& descr,const int32_t& defValue);
   bool add(const std::string& name,const std::string& descr,const uint32_t& defValue);
   bool add(const std::string& name,const std::string& descr,const int64_t& defValue);
   bool add(const std::string& name,const std::string& descr,const uint64_t& defValue);
   bool add(const std::string& name,const std::string& descr,const std::string& defValue);
   bool addComposed(const std::string& name,const std::string& descr);
   bool finalize();
   bool get(const std::string& name,bool& value);
   bool get(const std::string& name,float& value);
   bool get(const std::string& name,double& value);
   bool get(const std::string& name,int32_t& value);
   bool get(const std::string& name,uint32_t& value);
   bool get(const std::string& name,int64_t& value);
   bool get(const std::string& name,uint64_t& value);
   bool get(const std::string& name,std::string& value);
   bool get(const std::string& name,std::vector<bool>& value);
   bool get(const std::string& name,std::vector<float>& value);
   bool get(const std::string& name,std::vector<double>& value);
   bool get(const std::string& name,std::vector<int32_t>& value);
   bool get(const std::string& name,std::vector<uint32_t>& value);
   bool get(const std::string& name,std::vector<int64_t>& value);
   bool get(const std::string& name,std::vector<uint64_t>& value);
   bool get(const std::string& name,std::vector<std::string>& value);
   bool initialize(const std::string& envPrefix,int argn,char** args,MPI_Comm comm,int masterRank);
   bool isInitialized();
   bool parse();
   bool printHelpMessage();

 private:
   int argn;                                                 /**< Private copy of argn, as given in main.*/
   char** args;                                              /**< Private pointer to args, as given in main.*/
   MPI_Comm comm;                                            /**< MPI communicator used by ConfigReader.*/
   std::string envPrefix;                                    /**< Prefix for environment variables read by ConfigReader.*/
   std::string fileNameGlobalConfig;                         /**< Name of global configuration file, read from command line
							      * or environment variable.*/
   std::string fileNameRunConfig;                            /**< Name of runtime configuration file, read from command line 
							      * or environment variable.*/
   std::string fileNameUserConfig;                           /**< User-specific configuration file, read from command line 
							      * or environment variable.*/
   bool initialized;                                         /**< If true, ConfigReader has initialized successfully.*/
   int masterRank;                                           /**< MPI rank of master process in communicator comm.*/
   std::map<std::string,std::string> options;                /**< List of currently known [option,value] pairs, identified 
							      * by the option names.*/
   std::map<std::string,bool> optionBroadcasted;             /**< For each option, boolean status if the option's value 
							      * has been distributed to all processes in communicator comm.*/
   int rank;                                                 /**< MPI rank of this process in communicator comm.*/
   std::map<std::string,std::vector<std::string> > vectorOptions;
   std::map<std::string,bool> vectorOptionBroadcasted;
};

namespace configreader {
   /** Compare values in array against the given value. If any 
    * element in array is equal to given value, as compared by '==' 
    * operator, function returns 'false'. This is a helper function 
    * that is intended for checking if user gave required values in 
    * configuration file.
    * Template parameter SIZE is the number of elements in input array.
    * @param values Array containing input values.
    * @param comparedValue Each element in input array is compared against this value.
    * @return If true, one or more elements in input array are equal to comparedValue.*/
   template<int SIZE,typename C> inline
   bool checkInputValues(C* array,const C comparedValue) {
      for (int i=0; i<SIZE; ++i) {
	 if (array[i] == comparedValue) return false;
      }
      return true;
   }
}

#endif
