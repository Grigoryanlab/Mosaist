# Mosaist

Mosaist is a C++ library for working with and designing protein structures and sequences. It is publicly available under the terms of the
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).


# Installation

Running `make` (after appropriately modifying `makefile`, as necessary) will compile all library components, tests, and programs.
Running programs without command line arguments will display usage statements. Programs `search` and `design` may be of particular interest to users:
* `search` demonstrates the use of our efficient protein sub-structure search algorithm `FASST`. See `programs/search.cpp` for how to use `FASST`, with the relevant library code found in `srs/mstfasst.cpp` and `include/mstfasst.h`. Use program `fasstDB` to build a searchable database out of an arbitrary set of structures.
* `design` uses our newest efficient implementation of the dTERMen design method (see [Zhou et al.](http://dx.doi.org/10.1073/pnas.1908723117)) to enable a variety of structure-based design and optimization tasks.


# Python interface

A minimal python interface to Mosaist, which relies on the `Boost.Python` library, is implemented in `src/mstpython.cpp`. At this time, this is not
the primary means by which we (in the Grigoryan lab) use Mosaist ourselves, so the interface is provided mainly as a starting point for folks who might find
such an interface useful. Running `make libs python` (after installing `Boost.Python` for your version of python) should produce a python library
under `libs/mstpython.so`, which can then be accessed in python (e.g., `impost mstpython as mst`), provided the `.so` file is placed in the library path.
However, we have found that compiling `Boost.Python` interfaces tends to be very brittle and compilation details often need to be adjusted based on
machine and operating system. At the very least, you will likely need to modify the `pythonExec` variable within `makefile` to point to your python interpreter.
