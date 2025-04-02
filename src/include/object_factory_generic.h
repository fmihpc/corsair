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

#ifndef OBJECT_FACTORY_GENERIC_H
#define OBJECT_FACTORY_GENERIC_H

#include <map>
#include <string>

/** Generic Object factory. Object factory is a class that, upon request, 
 * creates new objects. Object factory does not now what kind of objects it 
 * can create, user must provide necessary blueprints by registering maker functions
 * at runtime. Makers are actually wrappers that simply call class constructors.
 * In reality an object factory can make objects of single type, which is given 
 * as template parameter.
 * 
 * In Corsair object factories are used to create grid builders, particle accumulators, 
 * injectors, etc. Object factories are stored in object_factories.h file.
 * 
 * Note that one cannot use named variables in declarations, i.e., it is invalid 
 * to declare an std::map as "std::map<std::string name,bool value>". Correct way
 * to declare it is "std::map<std::string,bool>". In a similar manner a map 
 * containing a function pointer has to be declared as "std::map<std::string,PRODUCT* (*)()>".
 * This is synonymous to first defining an alias for function pointer "typedef PRODUCT* (*ProductMaker)()",
 * and then declaring "std::map<std::string,ProductMaker>". However, the typedef below 
 * depends on a template parameter and thus cannot be defined outside ObjectFactoryGeneric, 
 * so we stick to the long version of map declarations.
 */
template<class PRODUCT>
class ObjectFactoryGeneric {
 public:
   ObjectFactoryGeneric();
   ~ObjectFactoryGeneric();

   PRODUCT* create(const std::string& name) const;
   bool exists(const std::string& name) const;
   bool registerMaker(const std::string& name,PRODUCT* (*ProductMaker)());

 private:

   std::map<std::string,PRODUCT* (*)()> registeredMakers; /**< Container for all registered PRODUCTs. In reality function 
							   * pointers that return a new PRODUCT are stored here.*/							  

};

/** Default constructor (empty).*/
template<class PRODUCT>
ObjectFactoryGeneric<PRODUCT>::ObjectFactoryGeneric() { }

/** Default destructor (empty).*/
template<class PRODUCT>
ObjectFactoryGeneric<PRODUCT>::~ObjectFactoryGeneric() { }

/** Return a new PRODUCT of the requested type, if its maker function has been registered.
 * @param name Name of the requested object. Must correspond to one of the names 
 * that were given in ObjectFactoryGeneric::registerMaker function call.
 * @return If requested object could be manufactured, pointer to new object is returned.
 * Otherwise a NULL pointer is returned.*/
template<class PRODUCT>
PRODUCT* ObjectFactoryGeneric<PRODUCT>::create(const std::string& name) const {
   typename std::map<std::string,PRODUCT* (*)()>::const_iterator it = registeredMakers.find(name);
   if (it != registeredMakers.end()) return (*it->second)();
   return NULL;
}

/** Test if this factory can manufacture given object.
 * @param name Name of the requested object. Must correspond to one of the names 
 * that were given in ObjectFactoryGeneric::registerMaker function call.
 * @return If true, this factory can manufacture requested objects.*/
template<class PRODUCT>
bool ObjectFactoryGeneric<PRODUCT>::exists(const std::string& name) const {
   // Search for given PRODUCT Maker:
   typename std::map<std::string,PRODUCT* (*)()>::const_iterator it = registeredMakers.find(name);
   if (it != registeredMakers.end()) return true;
   return false;
}

/** Add new Maker to factory.
 * @name Unique name for the new object type.
 * @productMaker Function that creates a new object of class PRODUCT.
 * @return If true, given Maker function was successfully registered. If false, most probably 
 * a Maker function with the same name already exists in this factory.*/
template<class PRODUCT>
bool ObjectFactoryGeneric<PRODUCT>::registerMaker(const std::string& name,PRODUCT* (*productMaker)()) {
   // Check that a PRODUCT Maker with the same name doesn't already exist:
   typename std::map<std::string,PRODUCT* (*)()>::const_iterator it = registeredMakers.find(name);
   if (it != registeredMakers.end()) return false;

   // Add maker to list of known makers:
   registeredMakers[name] = productMaker;
   return true;
}

#endif
