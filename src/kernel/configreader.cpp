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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "mpiconversion.h"
#include "configreader.h"

using namespace std;
using boost::lexical_cast;
namespace PO = boost::program_options;

// ***** INIT STATIC MEMBER VARIABLES *****
PO::options_description* descriptions = NULL;
PO::variables_map* variables = NULL;

/** Default constructor. ConfigReader::initialize must be 
 * called before ConfigReader can be used.
 */
ConfigReader::ConfigReader() {
   initialized = false;
}

/** Default destructor. ConfigReader::finalize must be called 
 * before destruction.*/
ConfigReader::~ConfigReader() {

}

/** Add an option whose value should be read from configuration file(s).
 * This function will fail if an option with the same name has already been added.
 * This function can be called on all processes, but it only does something useful
 * at master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @param defValue Default value for the option.
 * @return If true, option was successfully added.*/
bool ConfigReader::add(const string& name,const string& descr,const bool& defValue) {
   stringstream ss;
   ss << defValue;
   return add(name,descr,ss.str());
}

/** Add an option whose value should be read from configuration file(s). 
 * This function will fail if an option with the same name has already been added.
 * This function can be called on all processes, but it only does something useful 
 * at master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @param defValue Default value for the option.
 * @return If true, option was successfully added.*/
bool ConfigReader::add(const string& name,const string& descr,const float& defValue) {
   stringstream ss;
   ss << setprecision(numeric_limits<float>::digits10+1) << defValue;
   return add(name,descr,ss.str());
}

/** Add an option whose value should be read from configuration file(s).
 * This function will fail if an option with the same name has already been added.
 * This function can be called on all processes, but it only does something useful
 * at master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @param defValue Default value for the option.
 * @return If true, option was successfully added.*/
bool ConfigReader::add(const string& name,const string& descr,const double& defValue) {
   stringstream ss;
   ss << setprecision(numeric_limits<double>::digits10+1) << defValue;
   return add(name,descr,ss.str());
}

/** Add an option whose value should be read from configuration file(s).
 * This function will fail if an option with the same name has already been added.
 * This function can be called on all processes, but it only does something useful
 * at master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @param defValue Default value for the option.
 * @return If true, option was successfully added.*/
bool ConfigReader::add(const string& name,const string& descr,const int32_t& defValue) {
   stringstream ss;
   ss << defValue;
   return add(name,descr,ss.str());
}

/** Add an option whose value should be read from configuration file(s).
 * This function will fail if an option with the same name has already been added.
 * This function can be called on all processes, but it only does something useful
 * at master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @param defValue Default value for the option.
 * @return If true, option was successfully added.*/
bool ConfigReader::add(const string& name,const string& descr,const uint32_t& defValue) {
   stringstream ss;
   ss << defValue;
   return add(name,descr,ss.str());
}

/** Add an option whose value should be read from configuration file(s).
 * This function will fail if an option with the same name has already been added.
 * This function can be called on all processes, but it only does something useful
 * at master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @param defValue Default value for the option.
 * @return If true, option was successfully added.*/
bool ConfigReader::add(const string& name,const string& descr,const int64_t& defValue) {
   stringstream ss;
   ss << defValue;
   return add(name,descr,ss.str());
}

/** Add an option whose value should be read from configuration file(s).
 * This function will fail if an option with the same name has already been added.
 * This function can be called on all processes, but it only does something useful
 * at master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @param defValue Default value for the option.
 * @return If true, option was successfully added.*/
bool ConfigReader::add(const string& name,const string& descr,const uint64_t& defValue) {
   stringstream ss;
   ss << defValue;
   return add(name,descr,ss.str());
}

/** Add an option whose value should be read from configuration file(s).
 * This function will fail if an option with the same name has already been added.
 * This function can be called on all processes, but it only does something useful
 * at master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @param defValue Default value for the option.
 * @return If true, option was successfully added.*/
bool ConfigReader::add(const string& name,const string& descr,const string& defValue) {
   if (initialized == false) return false;
   if (rank != masterRank) return true;

   // Add new option only if it does not already exist:
   pair<map<string,string>::iterator,bool> result = options.insert(make_pair(name,""));
   if (result.second == false) return false;

   // Add it to Boost:
   optionBroadcasted.insert(make_pair(name,false));
   descriptions->add_options()(name.c_str(),PO::value<string>(&(result.first->second))->default_value(defValue),descr.c_str());   
   return true;
}

/** Add a composed option whose values should be read from configuration file(s).
 * This function will fail if an option with the same name has already been added.
 * This function can be called by all processes, but it only has an effect on 
 * master process.
 * @param name Name of the option.
 * @param descr Description for the option (optional).
 * @return If true, composed option was successfully added.*/
bool ConfigReader::addComposed(const string& name,const string& descr) {
   if (initialized == false) return false;
   if (rank != masterRank) return true;

   // Add new vector option only if it doesn't already exist:
   pair<map<string,vector<string> >::iterator,bool> result = vectorOptions.insert(make_pair(name,vector<string>()));
   if (result.second == false) return false;

   // Add vector option to Boost:
   vectorOptionBroadcasted.insert(make_pair(name,false));
   descriptions->add_options()(name.c_str(),PO::value<vector<string> >(&(result.first->second))->composing(),descr.c_str());
   return true;
}

/** Finalize ConfigReader. This function deallocates internal memory.
 * Calling other member functions except initialize after finalize 
 * will result in unspecified behaviour.
 * @return If true, ConfigReader was successfully finalized.*/
bool ConfigReader::finalize() {
   initialized = false;
   delete descriptions; descriptions = NULL;
   delete variables; variables = NULL;
   MPI_Comm_free(&comm);
   return true;
}

/** Get an options value. The option must have been previously added to ConfigReader by 
 * calling one of the ConfigReader::add functions. If successful, this function will 
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param value Variable in which option's value is copied.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,bool& value) {
   string s;
   if (get(name,s) == false) return false;
   value = lexical_cast<bool>(s);
   return true;
}

/** Get an options value. The option must have been previously added to ConfigReader by 
 * calling one of the ConfigReader::add functions. If successful, this function will 
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param value Variable in which option's value is copied.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,float& value) {
   string s;
   if (get(name,s) == false) return false;
   if (s == "nan") {value = NAN; return true;}
   else if (s == "inf") {value = numeric_limits<float>::infinity(); return true;}
   value = lexical_cast<float>(s);
   return true;
}

/** Get an options value. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param value Variable in which option's value is copied.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,double& value) {
   string s;
   if (get(name,s) == false) return false;
   if (s == "nan") {value = NAN; return true;}
   else if (s == "inf") {value = numeric_limits<double>::infinity(); return true;}
   value = lexical_cast<double>(s);
   return true;
}

/** Get an options value. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param value Variable in which option's value is copied.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,int32_t& value) {
   string s;
   if (get(name,s) == false) return false;
   value = lexical_cast<int32_t>(s);
   return true;
}

/** Get an options value. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param value Variable in which option's value is copied.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,uint32_t& value) {
   string s;
   if (get(name,s) == false) return false;
   value = lexical_cast<uint32_t>(s);
   return true;
}

/** Get an options value. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param value Variable in which option's value is copied.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,int64_t& value) {
   string s;
   if (get(name,s) == false) return false;
   value = lexical_cast<int64_t>(s);
   return true;
}

/** Get an options value. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param value Variable in which option's value is copied.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,uint64_t& value) {
   string s;
   if (get(name,s) == false) return false;
   value = lexical_cast<uint64_t>(s);
   return true;
}

/** Get an options value. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param value Variable in which option's value is copied.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,string& value) {
   // Attempt to find the requested option:
   map<string,string>::const_iterator it = options.find(name);
   if (it == options.end()) return false;
   // Copy its value to output variable:
   value = it->second;
   return true;
}

/** Get all values for a composed option. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param values Vector in which option's values are copied to.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,vector<bool>& values) {
   // Attempt to find the requested option:
   map<string,vector<string> >::const_iterator it = vectorOptions.find(name);
   if (it == vectorOptions.end()) return false;
   // Copy values to output vector:
   values.resize(it->second.size());
   for (size_t i=0; i<it->second.size(); ++i) values[i] = lexical_cast<bool>(it->second[i]);
   return true;
}

/** Get all values for a composed option. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param values Vector in which option's values are copied to.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,vector<float>& values) {
   // Attempt to find the requested option:
   map<string,vector<string> >::const_iterator it = vectorOptions.find(name);
   if (it == vectorOptions.end()) return false;
   // Copy values to output vector:
   values.resize(it->second.size());
   for (size_t i=0; i<it->second.size(); ++i) values[i] = lexical_cast<float>(it->second[i]);
   return true;
}

/** Get all values for a composed option. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param values Vector in which option's values are copied to.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,vector<double>& values) {
   // Attempt to find the requested option:
   map<string,vector<string> >::const_iterator it = vectorOptions.find(name);
   if (it == vectorOptions.end()) return false;
   // Copy values to output vector:
   values.resize(it->second.size());
   for (size_t i=0; i<it->second.size(); ++i) values[i] = lexical_cast<double>(it->second[i]);
   return true;
}

/** Get all values for a composed option. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param values Vector in which option's values are copied to.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,vector<int32_t>& values) {
   // Attempt to find the requested option:
   map<string,vector<string> >::const_iterator it = vectorOptions.find(name);
   if (it == vectorOptions.end()) return false;
   // Copy values to output vector:
   values.resize(it->second.size());
   for (size_t i=0; i<it->second.size(); ++i) values[i] = lexical_cast<int32_t>(it->second[i]);
   return true;
}

/** Get all values for a composed option. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param values Vector in which option's values are copied to.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,vector<uint32_t>& values) {
   // Attempt to find the requested option:
   map<string,vector<string> >::const_iterator it = vectorOptions.find(name);
   if (it == vectorOptions.end()) return false;
   // Copy values to output vector:
   values.resize(it->second.size());
   for (size_t i=0; i<it->second.size(); ++i) values[i] = lexical_cast<uint32_t>(it->second[i]);
   return true;
}

/** Get all values for a composed option. The option must have been previously added to ConfigReader by
 * calling one of the ConfigReader::add functions. If successful, this function will
 * return the same value on all processes.
 * @param name Name of the option whose value is requested.
 * @param values Vector in which option's values are copied to.
 * @return If true, the output value is valid.*/
bool ConfigReader::get(const string& name,vector<string>& values) {
   // Attempt to find the requested option:
   map<string,vector<string> >::const_iterator it = vectorOptions.find(name);
   if (it == vectorOptions.end()) return false;

   // Copy values to output vector:
   values = it->second;
   return true;
}

/** Initialize ConfigReader.
 * @param comm MPI communicator that ConfigReader should use to broadcast option values.
 * @param masterRank MPI master process rank in communicator comm.
 * @param argn argn as given to main function (number of command line parameters).
 * @param args args as given to main function (command line parameters).
 * @return If true, ConfigReader initialized successfully.
 */
bool ConfigReader::initialize(const std::string& envPrefix,int argn,char** args,MPI_Comm comm,int masterRank) {
   // Do not initialize, if initialization has already been done:
   if (initialized == true) return false;

   initialized = false;
   //this->sim = &sim;
   this->argn = argn;
   this->args = args;
   this->envPrefix = envPrefix;

   // Get the rank of this process in communicator comm
   MPI_Comm_dup(comm,&this->comm);
   this->masterRank = masterRank;
   MPI_Comm_rank(this->comm,&rank);
   //rank = sim.mpiRank;

   if (rank == masterRank) {
      descriptions = new PO::options_description("All known input file parameters (name,default value,description)");
      initialized = true;

      descriptions->add_options()
	("globalconfig",PO::value<string>(&fileNameGlobalConfig)->default_value(""),
	    "Name of file containing global configuration parameters, relative to current working directory (string).")
	("userconfig",PO::value<string>(&fileNameUserConfig)->default_value(""),
	    "Name of file containing per-user configuration parameters, relative to current working directory (string).")
	("runconfig",PO::value<string>(&fileNameRunConfig)->default_value(""),
	    "Name of file containing runtime configuration parameters, relative to current working directory (string).");
   }

   // Master broadcasts init status to all processes:
   MPI_Bcast(&initialized,sizeof(bool),MPI_BYTE,masterRank,comm);
   return initialized;
}

/** Query if ConfigReader has initialized successfully.
 * @return If true, ConfigReader is ready for use.
 */
bool ConfigReader::isInitialized() {return initialized;}

/** Parse values for options that have been added to ConfigReader by calling 
 * ConfigReader::add functions. After this function completes successfully, all 
 * processes have the same [option,value] pairs in their ConfigReaders.
 * @return If true, option values were successfully parsed.
 */
bool ConfigReader::parse() {
   if (initialized == false) return false;

   if (rank == masterRank) {
      const bool allowUnknownOptions = true;
      variables = new PO::variables_map;

      // Read options from command line:
      PO::store(PO::parse_command_line(argn,args,*descriptions),*variables);
      PO::notify(*variables);

      // Read options from environment variables:
      PO::store(PO::parse_environment(*descriptions,envPrefix),*variables);
      PO::notify(*variables);

      // Read options from runtime config file:
      if (fileNameRunConfig.size() > 0) {
	 fstream in(fileNameRunConfig.c_str(),fstream::in);
	 if (in.good() == true) {
	    PO::store(PO::parse_config_file(in,*descriptions,allowUnknownOptions),*variables);
	 } else {
	    cerr << "Could not open runtime config file '" << fileNameRunConfig << "'" << endl;
	    exit(1);
	 }
	 in.close();
      }

      // Read options from user-specific config file:
      if (fileNameUserConfig.size() > 0) {
	 fstream in(fileNameUserConfig.c_str(),fstream::in);
	 if (in.good() == true) {
	    PO::store(PO::parse_config_file(in,*descriptions,allowUnknownOptions),*variables);
	 } else {
	    cerr << "Could not open user config file '" << fileNameUserConfig << "'" << endl;
	    exit(1);
	 }
	 in.close();
      }

      // Read options from global config file:
      if (fileNameGlobalConfig.size() > 0) {
	 fstream in(fileNameGlobalConfig.c_str(),fstream::in);
	 if (in.good() == true) {
	    PO::store(PO::parse_config_file(in,*descriptions,allowUnknownOptions),*variables);
	 } else {
	    cerr << "Could not open global config file '" << fileNameGlobalConfig << "'" << endl;
	    exit(1);
	 }
	 in.close();
      }

      // Process values read from all config files:
      PO::notify(*variables);

      delete variables; variables = NULL;
   }

   // Master process counts the number of option-value pairs to broadcast:
   unsigned int N_broadcasts = 0;
   if (rank == masterRank) {
      for (map<string,bool>::const_iterator it=optionBroadcasted.begin(); it!=optionBroadcasted.end(); ++it) {
	 if (it->second == false) ++N_broadcasts;
      }
   }
   MPI_Bcast(&N_broadcasts,1,MPI_Type<unsigned int>(),masterRank,comm);

   // Master process broadcasts new option-value pairs:
   const unsigned int maxLength = 512;
   char optionName[maxLength];
   char optionValue[maxLength];
   if (rank == masterRank) {
      for (map<string,string>::const_iterator it=options.begin(); it!=options.end(); ++it) {
	 map<string,bool>::iterator jt=optionBroadcasted.find(it->first);
	 if (jt->second == true) continue;
	 strncpy(optionName,it->first.c_str(),maxLength-1);
	 strncpy(optionValue,it->second.c_str(),maxLength-1);
	 MPI_Bcast(optionName,maxLength,MPI_Type<char>(),masterRank,comm);
	 MPI_Bcast(optionValue,maxLength,MPI_Type<char>(),masterRank,comm);
	 jt->second = true;
      }
   } else {
      for (unsigned int n=0; n<N_broadcasts; ++n) {
	 MPI_Bcast(optionName,maxLength,MPI_Type<char>(),masterRank,comm);
	 MPI_Bcast(optionValue,maxLength,MPI_Type<char>(),masterRank,comm);
	 options.insert(make_pair(optionName,optionValue));
      }
   }

   // Master process counts number of vector options that have not been broadcasted:
   N_broadcasts = 0;
   if (rank == masterRank) {
      for (map<string,bool>::const_iterator it=vectorOptionBroadcasted.begin(); it!=vectorOptionBroadcasted.end(); ++it) {
	 if (it->second == true) continue;
	 ++N_broadcasts;
      }
   }
   MPI_Bcast(&N_broadcasts,1,MPI_Type<unsigned int>(),masterRank,comm);

   // Master process broadcasts new vector options:
   size_t vectorSize = 0;
   if (rank == masterRank) {
      for (map<string,vector<string> >::const_iterator it=vectorOptions.begin(); it!=vectorOptions.end(); ++it) {
	 map<string,bool>::iterator jt=vectorOptionBroadcasted.find(it->first);
	 if (jt->second == true) continue;

	 // Broadcast vector name:
	 strncpy(optionName,it->first.c_str(),maxLength-1);
	 MPI_Bcast(optionName,maxLength,MPI_Type<char>(),masterRank,comm);
	 // Broadcast vector size:
	 vectorSize = it->second.size();
	 MPI_Bcast(&vectorSize,1,MPI_Type<size_t>(),masterRank,comm);
	 // Broadcast vector contents:
	 for (size_t i=0; i<vectorSize; ++i) {
	    strncpy(optionValue,it->second[i].c_str(),maxLength-1);
	    MPI_Bcast(optionValue,maxLength,MPI_Type<char>(),masterRank,comm);
	 }
	 // Mark vector option as broadcasted:
	 jt->second = true;
      }
   } else {
      for (unsigned int i=0; i<N_broadcasts; ++i) {
	 // Receive option name:
	 MPI_Bcast(optionName,maxLength,MPI_Type<char>(),masterRank,comm);
	 // Receive vector size:
	 MPI_Bcast(&vectorSize,1,MPI_Type<size_t>(),masterRank,comm);
	 // Receive vector contents:
	 pair<map<string,vector<string> >::iterator,bool> result = 
	   vectorOptions.insert(make_pair(optionName,vector<string>()));
	 for (size_t i=0; i<vectorSize; ++i) {
	    MPI_Bcast(optionValue,maxLength,MPI_Type<char>(),masterRank,comm);
	    result.first->second.push_back(optionValue);
	 }
      }
   }
   return true;
}

bool ConfigReader::printHelpMessage() {
   if (initialized == false) return false;
   if (rank != masterRank) return true;
   if (descriptions != NULL) cout << *descriptions << endl;
   return true;
}
