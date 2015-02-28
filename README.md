Smoke and Steam Simulation at Interactive Rates
===============================================

Max Smolens  
UNC-Chapel Hill  
COMP 259: Physically-Based Modeling, Simulation and Animation  
Final Project: Simulating Smoke and Steam at Interactive Rates  
Spring 2003  

![Smoke simulation screenshot](/doc/images/darksmoke.jpg?raw=true)

## Overview

This project implements much of the smoke simulation system described
in *Holtkamper, T. Real-time Gaseous Phenomena—A Phenomenological
Approach to Interactive Smoke and Stream, 2003*.  A concentration is
made on the visual effect of the smoke, rather than a
physically-correct simulation.

See [doc/index.html](doc/index.html) for details on the technique.

## Compiling

The program compiles on Linux using the provided Makefile.  I used
David McAllister's Particle Systems API, but had to modify the
API to allow finer-grained control over individual particles.  The
modified API code is provided in the lib/particle/ directory.  The program
also requires the ImageMagick and GLUI libraries.

## Controls

The simulation and rendering parameters can be changed via the GUI.
The W, A, S and D keys allow the view to be translated and rotated,
always facing the smoke.  The T and G keys translate in the
z-direction.

## Movies

The program can output screendumps as it runs.  The provided Makefile
can generate a movie from the screendumps, given that mencoder
(http://www.mplayerhq.hu) is installed by using 'make movie'.
