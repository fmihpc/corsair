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

#ifndef SIMULATION_H
#define SIMULATION_H

#include <mpi.h>
#include <definitions.h>
#include <profiler.h>
#include <pargrid_definitions.h>
#include <vlsv_common.h>

namespace block {
   const int WIDTH_X = 1;
   const int WIDTH_Y = 1;
   const int WIDTH_Z = 1;
   const int SIZE  = WIDTH_X*WIDTH_Y*WIDTH_Z;
   const int SIZE_PLUS_ONE_LAYER = (WIDTH_X+2)*(WIDTH_Y+2)*(WIDTH_Z+2);   
}

/** Container for generic simulation control variables. 
 * Only primitive datatypes are allowed here so that this 
 * header can be included by other header files.*/
struct Simulation {
   Simulation();
   ~Simulation();

   bool finalize();
   bool initialize(int argn,char** args);
   
   int argn;                               /**< Number of command line arguments.*/
   char** args;                            /**< Command line arguments.*/
   bool runTests;
   bool initializing;                      /**< If true, simulation is still initializing, i.e., the main 
					    * propagation loop hasn't started yet.*/
   
   // ******************************************** //
   // ***** GENERIC VARIABLES RELATED TO MPI ***** //
   // ******************************************** //

   MPI_Comm comm;                          /**< MPI communicator for the whole simulation, typically MPI_COMM_WORLD.*/
   const int MASTER_RANK;                  /**< MPI rank of master process in communicator comm.*/
   int mpiProcesses;                       /**< Number of MPI processes in communicator comm.*/
   int mpiRank;                            /**< MPI rank of the process in communicator comm.*/

   // ********************************************* //
   // ***** VARIABLES RELATED TO DATA WRITING ***** //
   // ********************************************* //

   bool atDataSaveStep;                    /**< If true, simulation data is saved on this time step.*/
   bool dataIntervalIsTime;
   unsigned int dataIntervalInteger;
   Real dataIntervalFloat;
   Real t_previousDataSave;
   std::string meshFileName;               /**< Name of output file that contains currently valid mesh.*/
   bool meshAlwaysWritten;                 /**< If true, mesh is written to every VLSV file. Default behaviour
					    * is to write the mesh only if it has changed, i.e. after repartitioning.*/

   // *************************************************** //
   // ***** VARIABLES RELATED TO SIMULATION RESTART ***** //
   // *************************************************** //
   
   bool restarted;                         /**< If true, simulation was restarted.*/
   unsigned int restartTimestep;           /**< Time step when simulation was restarted (if applicable).*/
   unsigned int restartWriteInterval;      /**< Interval, in time steps, between restart file writes.*/
   std::string restartFilenamePrefix;      /**< Restart filename prefix.*/
   std::size_t restartMajorInterval;       /**< Major interval, measured in timesteps, at which restart files are always kept.
					    * If this equals to zero value then all restart files are tested for deletion.*/
   std::size_t restartMinorFileAmount;     /**< Number of most recent restart files that are always kept.*/
   
   // *********************************************** //
   // ***** VARIABLES RELATED TO BALANCING LOAD ***** //
   // *********************************************** //
   
   bool countPropagTime;                   /**< If true, propagators should measure propagation time of each local cell.*/
   int repartitionCheckInterval;           /**< Interval, in number of time steps, when repartitioning is done.*/
   bool meshRepartitioned;                 /**< If true, mesh has been repartitioned since last time step 
					    * and variables depending on partitioning need to be recalculated.*/
   Real maximumLoadImbalance;
   unsigned int meshChangedStep;           /**< Time step when mesh changed last time (=repartitioned).*/
   unsigned int meshWrittenStep;           /**< Time step when mesh was written last time. 
					    * If meshChangedStep is larger than meshWrittenStep, the mesh 
					    * needs to be rewritten to output file.*/

   // ************************************************ //
   // ***** VARIABLES RELATED TO SIMULATION MESH ***** //
   // ************************************************ //
   
   vlsv::geometry::type meshGeometry;      /**< Used coordinate system. Set by GridBuilder.*/
   unsigned int x_blocks;                  /**< Number of mesh blocks in x-direction. Value is set by a GridBuilder class.*/
   unsigned int y_blocks;                  /**< Number of mesh blocks in y-direction. Value is set by a GridBuilder class.*/
   unsigned int z_blocks;                  /**< Number of mesh blocks in z-direction. Value is set by a GridBuilder class.*/
   
   bool dx_uniform;                        /**< If true, cell spacing in x-direction is constant.
					    Set to 'true' in constructor. Value may be re-set by a GridBuilder class.*/
   bool dy_uniform;                        /**< If true, cell spacing in y-direction is constant.
					    Set to 'true' in constructor. Value may be re-set by a GridBuilder class.*/
   bool dz_uniform;                        /**< If true, cell spacing in z-direction is constant.
					    Set to 'true' in constructor. Value may be re-set by a GridBuilder class.*/

   bool x_periodic;                        /**< If true, grid is periodic in x-direction.
					    Set to 'false' in constructor. Value may be re-set by a GridBuilder class.*/
   bool y_periodic;                        /**< If true, grid is periodic in y-direction.
					    Set to 'false' in constructor. Value may be re-set by a GridBuilder class.*/
   bool z_periodic;                        /**< If true, grid is periodic in z-direction..
					    Set to 'false' in constructor. Value may be re-set by a GridBuilder class.*/
   Real x_min,x_max,y_min,y_max,z_min,z_max; /** Simulation domain boundary coordinates
					    Set as zero in constructor. Value may be re-set by a GridBuilder class.*/

   Real* RESTRICT x_crds_node;             /**< Pointer to array containing node x-coordinates. Number of x-nodes
					    * is equal to number of x-cells plus one.
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   Real* RESTRICT y_crds_node;             /**< Pointer to array containing node y-coordinates. Number of y-nodes
					    * is equal to number of y-cells plus one.
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   Real* RESTRICT z_crds_node;             /**< Pointer to array containing node z-coordinates. Number of z-nodes
					    * is equal to number of z-cells plus one.
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   Real* RESTRICT dx_block;                /**< Pointer to array containing block sizes in x-direction.
					    * The size of array equals the number of mesh blocks in x-direction.
					    * Block size is equal to sum of its cell sizes, i.e.,
					    * dx_block = sum dx_cell(i), where 0 < i < block::WIDTH_X. 
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   Real* RESTRICT dy_block;                /**< Pointer to array containing block sizes in y-direction.
					    * The size of array equals the number of mesh blocks in y-direction.
					    * Block size is equal to sum of its cell sizes, i.e.,
					    * dy_block = sum dy_cell(i), where 0 < i < block::WIDTH_Y.
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   Real* RESTRICT dz_block;                /**< Pointer to array containing block sizes in z-direction.
					    * The size of array equals the number of mesh blocks in z-direction.
					    * Block size is equal to sum of its cell sizes, i.e.,
					    * dz_block = sum dz_cell(i), where 0 < i < block::WIDTH_Z.
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   Real* RESTRICT dx_cell;                 /**< Pointer to array containing cell sizes in x-direction.
					    * If dx_uniform equals true only dx_cell[0] is valid, otherwise the 
					    * size of array equals the number of cells (not blocks) in x-direction.
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   Real* RESTRICT dy_cell;                 /**< Pointer to array containing cell sizes in y-direction.
					    * If dy_uniform equals true only dy_cell[0] is valid, otherwise the 
					    * size of array equals the number of cells (not blocks) in y-direction.
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   Real* RESTRICT dz_cell;                 /**< Pointer to array containing cell sizes in z-direction.
					    * If dz_uniform equals true only dz_cell[0] is valid, otherwise the 
					    * size of array equals the number of cells (not blocks) in z-direction.
					    * This array is created by a GridBuilder class and deleted in destructor.*/
   
   Real dt;                                /**< Current global simulation time step.*/
   Real t;                                 /**< Current simulation time.*/
   uint32_t timestep;                      /**< Current simulation timestep.*/
   Real maximumTime;                       /**< Simulation exits if time reached this value. A negative value
					    * negative value indicates that maximumTime is not used as an end condition.
					    * maximumTime and maximumTimesteps cannot both have negative value(s).*/					   
   uint32_t maximumTimesteps;              /**< Simulation exits after this many timesteps have been taken. A negative value
					    * indicates that maximumTimesteps is not used as an end condition.
					    * maximumTime and maximumTimesteps cannot both have negative value(s).*/
   
   static pargrid::DataID crdsDataID;             /**< ParGrid data ID associated with cell coordinate array.*/
   static pargrid::StencilID crdsStencilID;       /**< ID of ParGrid Stencil that is used to sync cell coordinate array.*/

   pargrid::StencilID defaultStencilID;           /**< ID of ParGrid default Stencil, equals to pargrid::DEFAULT_STENCIL.*/
   pargrid::StencilID inverseStencilID;           /**< ID of inverse of ParGrid default Stencil. This Stencil is 
						   * allocated in Corsair main.cpp and is used by particle lists.*/
   
 private:
   Simulation(const Simulation& s);
};

#endif
