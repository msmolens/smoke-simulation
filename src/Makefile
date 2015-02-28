CC = g++

PLIB =../

#COPT =-O2 -g #-DDEPTH_TEST
COPT =-O3 #-DDEPTH_TEST
CFLAGS = $(COPT) $(COMPFLAGS) -I$(PLIB) -I../../glui/include -I$(HOME)/include

LIBDIR = -L$(PLIB)/particle -L../../glui/lib
LIBSG = $(LIBDIR) `Magick++-config --cxxflags --cppflags --ldflags --libs`  -lparticle -lglut -lGL -lGLU -lglui -lm
LDFLAGS =-IPA

ALL = smoke

all: $(ALL)

.cpp.o:
	$(CC) $(CFLAGS) -c $<

pspray: PSpray.o
	$(CC) $(LDFLAGS) -o $@ PSpray.o $(LIBSG)

smoke: smoke.o Noise.o turbulence.o ImageWriter.o
	$(CC) $(LDFLAGS) -o $@ *.o $(LIBSG)

Noise.o:
	$(CC) $(CFLAGS) -c Noise.cpp

turbulence.o: Noise.o
	$(CC) $(CFLAGS) -c turbulence.cpp

ImageWriter.o:
	$(CC) $(CFLAGS) -c ImageWriter.cpp

clean:
	rm -f $(ALL)
	rm -f *.o

movie:
	./convert
	rm -f *.jpg
	mencoder -mf on:type=sgi -ovc lavc -lavcopts vcodec=mpeg4 -o output_lavc.avi \*.sgi
	rm -f *.sgi

movieclean:
	rm -f *.jpg
	rm -f *.avi


