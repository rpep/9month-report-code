# 9month-report-code

To build the code requires a compiler with support for the OpenMP library, an installation of Python, and Cython.
This code has only been tested on Linux and OS X. To build on these platforms, set the environment variable CC to a compiler satisfying the requirements, and then type

> make
 
from the build folder.

To adjust the number of particles, adjust the variable N in the file pyfmm.pyx and rerun the make command.

To run the code, from a Python prompt, run
> import pyfmm
