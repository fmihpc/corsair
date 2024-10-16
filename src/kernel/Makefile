include ../../Makefile.${ARCH}

# Names of files that should be compiled and 
# inserted to libsimkernel.a:

OBJS = 	configreader.o constants.o datawriter.o main.o\
	mpifile.o mpilogger.o\
	profiler.o restart_writer.o simulation.o\
	randomnumber.o

OBJS2 = particle_list_base.o base_class_particle_accumulator.o\
	base_class_particle_boundary_condition.o\
	base_class_particle_injector.o base_class_particle_propagator.o\
	root_finding.o triangulated_plane.o triangulated_sphere.o\
	ucd_mesh.o

# Make targets

default:
	${MAKE} libsimkernel.a
	${MAKE} libplistbase.a

clean:
	rm -rf *.o *.a *~
	rm -f ../lib/libsimkernel.a

libplistbase.a: ${OBJS2}
	${AR} r libplistbase.a ${OBJS2}
	ln -f -s ${CURDIR}/libplistbase.a -t ../lib

libsimkernel.a: ${OBJS}
	${AR} r libsimkernel.a ${OBJS}
	ln -f -s ${CURDIR}/libsimkernel.a -t ../lib

# Adjust include paths

override INCS+=${INC} ${INC_BOOST} ${INC_PARGRID} ${INC_VLSV} ${INC_ZOLTAN}
override INCS+=-I../include

# Dependencies

DEPS =../include/simulation.h ../include/constants.h ../include/randomnumber.h ../include/mpilogger.h\
     ../include/definitions.h ../include/cell.h ../include/profiler.h ../include/particle_list_base.h\
	particle_list_base.cpp

DEPS_BASE_PACCUMULATOR = ${DEPS} ../include/base_class_particle_accumulator.h base_class_particle_accumulator.cpp
DEPS_BASE_PBCOND = ${DEPS} ../include/base_class_particle_boundary_condition.h base_class_particle_boundary_condition.cpp 
DEPS_BASE_PINJECTOR = ${DEPS} ../include/base_class_particle_injector.h base_class_particle_injector.cpp
DEPS_BASE_PPROPAGATOR = ${DEPS} ../include/base_class_particle_propagator.h base_class_particle_propagator.cpp
DEPS_CONFIGREADER = ${DEPS} ../include/configreader.h configreader.cpp
DEPS_CONSTANTS = ../include/constants.h
DEPS_DATAWRITER = ${DEPS} ../include/datawriter.h datawriter.cpp
DEPS_MAIN = ${DEPS} ../include/main.h ../include/object_factories.h main.cpp
DEPS_MPIFILE = ${DEPS} ../include/mpifile.h mpifile.cpp
DEPS_MPILOGGER = ${DEPS} ../include/mpifile.h ../include/mpilogger.h mpilogger.cpp
DEPS_PLIST_BASE=${DEPS} ../include/particle_list_base.h particle_list_base.cpp
DEPS_PROFILER = ../include/profiler.h profiler.cpp
DEPS_RANDOM = ../include/randomnumber.h randomnumber_default.cpp
DEPS_ROOT_FIND = ../include/root_finding.h root_finding.cpp
DEPS_RESTART =../include/simulation.h ../include/simulationclasses.h ../include/restart_writer.h restart_writer.cpp
DEPS_SIMULATION = ../include/simulation.h ../include/definitions.h ../include/profiler.h
DEPS_TRIANGULATED_PLANE =  ${DEPS_UCD_MESH} ../include/triangulated_plane.h triangulated_plane.cpp
DEPS_TRIANGULATED_SPHERE = ${DEPS_UCD_MESH} ../include/triangulated_sphere.h triangulated_sphere.cpp
DEPS_UCD_MESH = ../include/definitions.h ../include/ucd_mesh.h ucd_mesh.cpp

# Compilation rules

.cpp.o:
	${CMP} ${CXXFLAGS} ${FLAGS} -c $< ${INCS}

base_class_particle_accumulator.o: ${DEPS_BASE_PACCUMULATOR}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c base_class_particle_accumulator.cpp ${INC_BOOST} ${INCS}

base_class_particle_boundary_condition.o: ${DEPS_BASE_PBCOND}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c base_class_particle_boundary_condition.cpp ${INC_BOOST} ${INCS}

base_class_particle_injector.o: ${DEPS_BASE_PINJECTOR}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c base_class_particle_injector.cpp ${INC_BOOST} ${INCS}

base_class_particle_propagator.o: ${DEPS_BASE_PINJECTOR}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c base_class_particle_propagator.cpp ${INC_BOOST} ${INCS}

configreader.o: ${DEPS_CONFIGREADER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c configreader.cpp ${INC_BOOST} ${INCS}

constants.o: ${DEPS_CONSTANTS}
	${CMP} ${CXXFLAGS} ${FLAGS} -c constants.cpp -I../include

datawriter.o: ${DEPS_DATAWRITER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c datawriter.cpp ${INCS}

main.o: ${DEPS_MAIN}
	${CMP} ${CXXFLAGS} ${FLAGS} -c main.cpp ${INCS} -I../gridbuilder

mpifile.o: ${DEPS_MPIFILE}
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpifile.cpp ${INCS}

mpilogger.o: ${DEPS_MPILOGGER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpilogger.cpp ${INCS}

particle_list_base.o: ${DEPS_PLIST_BASE}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_list_base.cpp ${INCS}

profiler.o: ${DEPS_PROFILER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c profiler.cpp ${INC_PARGRID} -I../include ${INC}

randomnumber.o: ${DEPS_RANDOM}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -o randomnumber.o -c randomnumber_default.cpp ${INCS}

restart_writer.o: ${DEPS_RESTART}
	${CMP} ${CXXFLAGS} ${FLAGS} -c restart_writer.cpp ${INCS}

root_finding.o: ${DEPS_ROOT_FIND}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c root_finding.cpp ${INCS}

simulation.o: ${DEPS_SIMULATION}
	${CMP} ${CXXFLAGS} ${FLAGS} -c simulation.cpp ${INCS}

triangulated_plane.o: ${DEPS_TRIANGULATED_PLANE}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c triangulated_plane.cpp ${INCS}

triangulated_sphere.o: ${DEPS_TRIANGULATED_SPHERE}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c triangulated_sphere.cpp ${INCS}

ucd_mesh.o: ${DEPS_UCD_MESH}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c ucd_mesh.cpp ${INCS}
