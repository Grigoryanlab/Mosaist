# Mosaist

Mosaist is a C++ library for working with and designing protein structures and sequences. It is publicly available under the terms of the
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).

# Installation

Running `make` (after appropriately modifying `makefile`, as necessary), will compile all library components, tests, and programs.
Programs `search` and `design` may be of particular interest to users:
* `search` demonstrates the use of our efficient protein sub-structure search algorithm `FASST`. See `programs/search.cpp` for how to use `FASST`, with the relevant library code found in `srs/mstfasst.cpp` and `include/mstfasst.h`.
* `design` uses our newest efficient implementation of the dTERMen design method (see [Zhou et al.](http://dx.doi.org/10.1073/pnas.1908723117)) to enable a variety of design and optimization tasks.
