#ifndef MAIN_H
#define MAIN_H

#include <vector>

#include <configreader.h>
#include <dataoperatorcontainer.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <object_factories.h>
#include <particle_list_base.h>

namespace corsair {

   struct ObjectWrapper {
      ConfigReader configReader;
      DataOperatorContainer dataOperatorContainer;
      ObjectFactories objectFactories;
      std::vector<ParticleListBase*> particleLists;
      Simulation sim;
      SimulationClasses simClasses;
   };

   ObjectWrapper& getObjectWrapper();

   namespace compuload {
      enum Items {
	 SECONDS_PER_CELL,
	 N_PARTICLES,
	 SIZE
      };
   }

   struct ComputationalLoad {
      Real loadAvg[compuload::SIZE];
      Real loadMin[compuload::SIZE];
      Real loadMax[compuload::SIZE];
      Real ratios[compuload::SIZE];
   };

   void evaluateComputationalLoad(corsair::ObjectWrapper& owrapper,corsair::ComputationalLoad& load);

} // namespace kernel

#endif
