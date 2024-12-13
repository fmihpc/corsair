####################################################
### SYSTEM-INDEPENDENT MAKEFILE, DO NOT EDIT !!! ###
####################################################

# Set default architecture, can be overwritten by command line argument
ARCH ?= arto
include Makefile.${ARCH}

# Set default archiver, can also be set in Makefile.arch
AR ?= ar

# Default project, can be set in Makefile.arch or overwritten
# by command line argument
SIM ?= example_particle

# Add src/user/SIM directory into include paths
INCS=-I${CURDIR}/src/user/${SIM}

# Distribution package name
DIR=corsair
DIST=corsair_v01_000.tar

# Build targets

default:
	${MAKE} "ARCH=${ARCH}" "SIM=${SIM}" "INCS=${INCS}" -C src

clean:
	rm -rf *~ corsair_* *.vlsv *.silo logfile.txt
	${MAKE} clean "ARCH=${ARCH}" "SIM=${SIM}" -C src

simclean:
	rm -rf *~ corsair_${SIM} *.vlsv *.silo logfile.txt
	${MAKE} simclean "ARCH=${ARCH}" "SIM=${SIM}" -C src

dist:
	ln -s ${CURDIR} ${DIR}
	tar -rf ${DIST} ${DIR}/Makefile ${DIR}/Makefile.arch
	tar -rf ${DIST} ${DIR}/COPYING ${DIR}/Doxyfile
	tar -rf ${DIST} ${DIR}/src/lib
	tar -rf ${DIST} ${DIR}/src/Makefile
	tar -rf ${DIST} ${DIR}/src/dataoperator/*
	tar -rf ${DIST} ${DIR}/src/gridbuilder/*
	tar -rf ${DIST} ${DIR}/src/include/*
	tar -rf ${DIST} ${DIR}/src/kernel/*
	tar -rf ${DIST} ${DIR}/src/particleinjector/*
	tar -rf ${DIST} ${DIR}/src/particlepropagator/*
	tar -rf ${DIST} ${DIR}/src/user/rectcuboid/*
	tar -rf ${DIST} ${DIR}/src/user/example_advection/*
	tar -rf ${DIST} ${DIR}/src/user/merka/*
	tar -rf ${DIST} ${DIR}/doc/*.pdf
	gzip -9 ${DIST}
	rm ${DIR}
