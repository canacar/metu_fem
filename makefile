# $Id: makefile,v 1.3 2008/07/09 05:50:18 canacar Exp $
# make file based on minimum given in the PetSc manual

ALL: forward
CFLAGS = -g -Wall -Wshadow 
FFLAGS =
CPPFLAGS =
FPPFLAGS =

FORWARD_SRC = forward.cc shape.cc shape8.cc shape20.cc dmatrix.cc \
              point3.cc femmesh3.cc hptimer.cc engine.cc meshutil.cc \
	      shape4.cc shape10.cc scache.cc sphere3.cc

FORWARD_OBJ = ${FORWARD_SRC:cc=o}

#include ${PETSC_DIR}/bmake/common/base
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

forward: ${FORWARD_OBJ} chkopts
	${CLINKER}  -g -L/usr/local/lib -o $@ ${FORWARD_OBJ} ${PETSC_LIB} -lz
