This is a repository created by Paolo Sylos Labini in order to test the Glenn Lawyer "Expected Force" algorithm on SNAP graphs in C++ .
"exffunction.cpp" copyright Glenn Lawyer, 2013.
--------------------------------------------------------------

USAGE (Visual Studio)
Download latest SNAP  distribution from https://snap.stanford.edu/snap/download.html

Configure Visual Studio for use with SNAP, which requires that the character set is set to Multi-byte and that the locations of SNAP include directories are specified.

To set the character set to Multi-byte 
go to Properties --> Configuration Properties --> General --> Projects Defaults --> Character Set --> Select "Use Multi-Byte Character Set".

To specify the location of SNAP include directories, 
go to Options --> Preferences --> C++ Directories --> Include Directories and add the paths on your system to SNAP folders glib-core, snap-core and snap-adv.

CONTENTS
------------------------------------------------------
"exffunction.cpp" is the Glenn Lawyer original function.

"main.cpp" loads a SNAP graph and calculate the expected force of its nodes using the function exfccp defined in "exffunction.cpp"

"stdafx.h" is an header for various useful libraries and the exfccp function.

"test.graph" is a 100-nodes SNAP graph saved in binary. 

Everything else correctly implements SNAP.

