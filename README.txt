This is a repository created by Paolo Sylos Labini in order to test the Glenn Lawyer Expected Force algorithm on SNAP graphs in C++ .
exffunction.cpp copyright Glenn Lawyer, 2013.
--------------------------------------------------------------

USAGE via Makefile
make all compile and test

make compile will compile the code and create executable named ExfForce. 

make run_test run test


CONTENTS
------------------------------------------------------
exffunction.cpp is the Glenn Lawyer original function. Calculates the expected force of a node.

main.cpp loads a graph from a (number of) text file(s) and calculate the expected force of the nodes. Files must contain a full, sorted edgelists such as
0  2
1  2
2  0
2  1

stdafx.h is an header for various useful libraries and the exfccp function.

fb_full.txt is a test graph. 
