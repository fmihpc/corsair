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
#include <map>
#include <vector>
#include <mpi.h>
#include <utility>
#include <iomanip>
#include <limits>
#include <cmath>

#include <mpiconversion.h>
#include "profiler.h"

using namespace std;

namespace profile {
   #ifdef PROFILE

   void sendValues();

   struct Leaf {
      int seqID;
      Leaf* parent;
      Leaf* address;
      map<int,Leaf> daughters;
      double t_start;
      double t_total;
      unsigned int calls;

      Leaf(): seqID(-1),parent(this),address(this),t_start(0.0),t_total(0.0),calls(0) { }

      Leaf(const Leaf& l) {
	 seqID = l.seqID;
	 parent = l.parent;
	 address = this;
	 daughters = l.daughters;
	 t_start = l.t_start;
	 t_total = l.t_total;
	 calls = l.calls;
      }
   };

   static MPI_Comm comm;
   static int N_processes;
   static int myrank;
   static int masterRank;

   static int nameWidth = 40;
   static int numberWidth = 11;
   static int rankWidth = 6;

   static int seqCounter = 0;
   static vector<string> names;
   static Leaf callTree;
   static Leaf* current = &callTree;
   static map<int,Leaf>::iterator leafIt;

   void print(Leaf& leaf,int nameID,int indent) {
      // Bcast number of characters in next item:
      string name = names[nameID];
      int N_characters = name.size();
      MPI_Bcast(&N_characters,1,MPI_Type<int>(),masterRank,comm);

      // Bcast name of next item in sequence:
      char* arr = new char[N_characters+1];
      arr[N_characters] = '\0';
      name.copy(arr,name.size(),0);
      MPI_Bcast(arr,N_characters,MPI_Type<char>(),masterRank,comm);
      name = arr;
      delete [] arr; arr = NULL;

      // Receive data from other processes:
      double* timeArray       = new double[N_processes];
      unsigned int* callArray = new unsigned int[N_processes];

      MPI_Gather(&leaf.t_total,1,MPI_Type<double>(),timeArray,1,MPI_Type<double>(),masterRank,comm);
      MPI_Gather(&leaf.calls,1,MPI_Type<unsigned int>(),callArray,1,MPI_Type<unsigned int>(),masterRank,comm);

      // Calculate statistics:
      double average = 0.0;
      double minValue = numeric_limits<double>::max();
      double maxValue = -numeric_limits<double>::max();
      double deviation = 0.0;

      double callAverage = 0.0;
      unsigned int callMin = numeric_limits<unsigned int>::max();
      unsigned int callMax = -numeric_limits<unsigned int>::max();

      int minProcRank = 0;
      int maxProcRank = 0;

      size_t callSum = 0;
      for (int i=0; i<N_processes; ++i) {
	 callSum += callArray[i];
	 if (callArray[i] < callMin) callMin = callArray[i];
	 if (callArray[i] > callMax) callMax = callArray[i];
      }
      callAverage = callSum / N_processes;

      for (int i=0; i<N_processes; ++i) {
	 average += timeArray[i];
	 if (timeArray[i] < minValue) {minValue = timeArray[i]; minProcRank = i;}
	 if (timeArray[i] > maxValue) {maxValue = timeArray[i]; maxProcRank = i;}
      }
      average /= N_processes;

      for (int i=0; i<N_processes; ++i) {
	 deviation += (timeArray[i]-average)*(timeArray[i]-average);
      }
      deviation = sqrt(deviation / N_processes);

      delete [] timeArray; timeArray = NULL;
      delete [] callArray; callArray = NULL;

      // Insert indentation and counter name
      string s;
      for (int i=0; i<indent; ++i) s += "  ";
      s += names[nameID];

      cout << setw(nameWidth) << s;

      // Write call stats
      cout << setw(numberWidth-1) << callAverage << ' ';
      //cout << setw(numberWidth-1) << callMin << ' ';
      //cout << setw(numberWidth-1) << callMax << ' ';
      //cout << setw(numberWidth-1) << callDeviation << ' ';

      // Write time stats
      cout << setw(numberWidth-1) << average << ' ';
      cout << setw(numberWidth-1) << minValue << ' ' << setw(rankWidth-1) << minProcRank << ' ';
      cout << setw(numberWidth-1) << maxValue << ' ' << setw(rankWidth-1) << maxProcRank << ' ';
      cout << setw(numberWidth-1) << deviation << ' ';
      cout << endl;

      for (map<int,Leaf>::iterator it=leaf.daughters.begin(); it!=leaf.daughters.end(); ++it) {
	 print(it->second,it->first,indent+1);
      }

      // Send "go up one level" command
      N_characters = 0;
      MPI_Bcast(&N_characters,1,MPI_Type<int>(),masterRank,comm);
   }

   void print(MPI_Comm comm,int masterRank) {
      MPI_Comm_size(comm,&N_processes);
      MPI_Comm_rank(comm,&myrank);
      profile::masterRank = masterRank;
      profile::comm = comm;

      if (myrank != masterRank) {
	 return sendValues();
      }

      cout << left;
      cout << resetiosflags(ios::floatfield);
      cout << setprecision(numberWidth-6);

      cout << setw(nameWidth) << "Name";
      cout << setw(numberWidth) << "Calls";
      //cout << setw(4*numberWidth) << "Total time (s)";
      cout << setw(6*numberWidth) << "Total time (s)";
      //cout << setw(4*numberWidth) << "Workunits / s";
      cout << endl;

      cout << setw(nameWidth) << " ";
      cout << setw(numberWidth) << "Total";
      for (int i=0; i<1; ++i) {
	 cout << setw(numberWidth) << "Average";
	 cout << setw(numberWidth) << "Min";
	 cout << setw(rankWidth) << "Rank";
	 cout << setw(numberWidth) << "Max";
	 cout << setw(rankWidth) << "Rank";
	 cout << setw(numberWidth) << "Std.dev.";
      }
      cout << endl;

      current = &callTree;      
      for (leafIt = current->daughters.begin(); leafIt != current->daughters.end(); ++leafIt) {
	 print(leafIt->second,leafIt->first,0);
      }

      int exitValue = -1;
      MPI_Bcast(&exitValue,1,MPI_Type<int>(),masterRank,comm);
   }

   void sendValues() {
      string name;
      int N_characters = 0;
      current = &callTree;      

      while (N_characters >= 0) {
	 // Receive number of characters from master, also used 
	 // to send exit / ascend commands:
	 MPI_Bcast(&N_characters,1,MPI_Type<int>(),masterRank,comm);
	 if (N_characters < 0) break;

	 // Ascend one level in sequence:
	 if (N_characters == 0) {
	    //cerr << "Received go up one level command" << endl;
	    current = current->parent;
	    continue;
	 }

	 // Master bcasts the name of the next item in sequence:
	 char* arr = new char[N_characters+1];
	 arr[N_characters] = '\0';
	 MPI_Bcast(arr,N_characters,MPI_Type<char>(),masterRank,comm);
	 name = arr;
	 delete [] arr; arr = NULL;

	 // Find all indices corresponding to the name of next item from local vector names:
	 vector<int> indices;
	 for (size_t i=0; i<names.size(); ++i) {
	    if (names[i] == name) indices.push_back(i);
	 }
	 // If name was not found, insert it:
	 if (indices.size() == 0) {
	    names.push_back(name);
	    indices.push_back(names.size()-1);
	 }

	 // Search nameID from daughters of current Leaf. Search will stop to first match:
	 int index = -1;
	 map<int,Leaf>::iterator it = current->daughters.end();
	 for (size_t i=0; i<indices.size(); ++i) {
	    it = current->daughters.find(indices[i]);
	    if (it != current->daughters.end()) {
	       index = indices[i];
	       break;
	    }
	 }

	 // Descend one level in sequence. If this process did not have 
	 // the next requested item, create it:
	 //map<int,Leaf>::iterator it = current->daughters.find(index);
	 if (it == current->daughters.end()) {
	    current->daughters[index].parent = current->address;
	    current->daughters[index].seqID = index;
	    it = current->daughters.find(index);
	 }
	 current = it->second.address;

	 // Gather data from current counter to master:
	 MPI_Gather(&(current->t_total),1,MPI_Type<double>(),NULL,0,MPI_Type<double>(),masterRank,comm);
	 MPI_Gather(&(current->calls),1,MPI_Type<unsigned int>(),NULL,0,MPI_Type<double>(),masterRank,comm);
      }
   }

   void start(const std::string& name,int& nameID) {
      #ifndef NDEBUG
         if (nameID >= (int)names.size()) {
	    cerr << "(PROFILER) ERROR: nameID has illegal value " << nameID << " should be < " << names.size() << endl;
	    exit(1);
	 }
      #endif

      // If name has not been labeled, assign a unique ID to it:
      if (nameID < 0) {
	 names.push_back(name);
	 nameID = names.size()-1;
      }

      // If this is the first time this sequence is called, create an entry for it:
      leafIt = current->daughters.find(nameID);
      if (leafIt == current->daughters.end()) {
	 current->daughters[nameID].parent = current->address;
	 current->daughters[nameID].seqID = seqCounter;
	 leafIt = current->daughters.find(nameID);
	 ++seqCounter;
      }

      current = leafIt->second.address;
      ++current->calls;
      current->t_start = MPI_Wtime();
   }

   void stop() {
      current->t_total += (MPI_Wtime() - current->t_start);
      current = current->parent;
   }
   #endif
}
