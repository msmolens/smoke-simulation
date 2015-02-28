######################################################################
# Particle System API
#
# Copyright 1998 by David K. McAllister.
#
######################################################################

C++ = g++

COPT = -O3

CFLAGS = $(COPT) $(COMPFLAGS) -I. -I..

POBJS =action_api.o actions.o opengl.o system.o

ALL = libparticle.a 

all: $(ALL)

.cpp.o:
	$(C++) $(CFLAGS) -c $<

libparticle.a: $(POBJS)
	rm -f $@
	ar clq $@ $(POBJS)

clean:
	rm -f libpar*
	rm -f *.o
	rm -rf *ii_files
