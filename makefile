CC = g++
CPPFLAGS=-O3 -std=c++11 -fPIC
ODIR = objs
SDIR = src
LDIR = lib
INCDIR = include

# python.boost stuff
uname := $(shell uname -s)
ifeq ($(uname),Linux)
    PYLIBPATH = $(shell python-config --exec-prefix)/lib64
else
    PYLIBPATH = $(shell python-config --exec-prefix)/lib
endif
PYLIB = -L$(PYLIBPATH) -L$(LDIR) $(shell python-config --libs) -lboost_python -Wl,-whole-archive -lmst -Wl,-no-whole-archive
PYOPTS = $(shell python-config --includes) -O2 -fPIC -std=c++11 -I$(INCDIR)

SOURCE  = mstfuser mstoptim mstlinalg mstoptions msttypes msttransforms mstrotlib mstmagic mstcondeg mstsequence mstsystem mstfasst
OBJECTS = $(patsubst %,$(ODIR)/%.o, $(SOURCE))

all: setup $(OBJECTS) libs tests
objs: $(OBJECTS)

setup:
	mkdir -p $(ODIR)
	mkdir -p bin
	mkdir -p lib

$(OBJECTS): $(ODIR)/%.o : $(SDIR)/%.cpp $(INCDIR)/%.h
	$(CC) $(CPPFLAGS) -I$(INCDIR) -c -o $@ $<

tests: $(OBJECTS)
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/mstoptions.o $(ODIR)/msttransforms.o $(ODIR)/mstfasst.o $(ODIR)/mstsequence.o tests/testFASST.cpp -o bin/testFASST
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o tests/test.cpp -o bin/test
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/msttransforms.o $(ODIR)/mstlinalg.o tests/testTransforms.cpp -o bin/testTransforms
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/msttransforms.o $(ODIR)/mstrotlib.o tests/testRotlib.cpp -o bin/testRotlib
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/mstmagic.o tests/testTERMUtils.cpp -o bin/testTERMUtils
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/mstrotlib.o $(ODIR)/msttransforms.o $(ODIR)/mstcondeg.o $(ODIR)/mstoptions.o $(ODIR)/mstsystem.o tests/testConFind.cpp -o bin/testConFind
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/mstrotlib.o $(ODIR)/msttransforms.o $(ODIR)/mstcondeg.o $(ODIR)/mstsystem.o tests/findBestFreedom.cpp -o bin/findBestFreedom
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/msttransforms.o $(ODIR)/mstoptim.o $(ODIR)/mstfuser.o $(ODIR)/mstlinalg.o tests/testFuser.cpp -o bin/testFuser
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/msttransforms.o $(ODIR)/mstoptim.o $(ODIR)/mstfuser.o $(ODIR)/mstlinalg.o tests/testAutofuser.cpp -o bin/testAutofuser
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/mstoptions.o tests/testClusterer.cpp -o bin/testClusterer
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/mstsystem.o programs/renumber.cpp -o bin/renumber
	$(CC) $(CPPFLAGS) -I./$(INCDIR) $(ODIR)/msttypes.o $(ODIR)/mstfasst.o $(ODIR)/mstoptions.o $(ODIR)/msttransforms.o $(ODIR)/mstsequence.o programs/findTERMs.cpp -o bin/findTERMs

libs: $(OBJECTS)
	ar rs lib/libmst.a $(ODIR)/mstoptions.o $(ODIR)/msttypes.o $(ODIR)/mstsequence.o $(ODIR)/mstsystem.o
	ar rs lib/libmsttrans.a $(ODIR)/msttransforms.o
	ar rs lib/libmstmagic.a $(ODIR)/msttypes.o $(ODIR)/mstmagic.o
	ar rs lib/libmstcondeg.a $(ODIR)/mstcondeg.o $(ODIR)/mstrotlib.o $(ODIR)/msttransforms.o
	ar rs lib/libmstlinalg.a $(ODIR)/mstlinalg.o
	ar rs lib/libmstoptim.a $(ODIR)/mstoptim.o $(ODIR)/mstlinalg.o
	ar rs lib/libmstfuser.a $(ODIR)/mstfuser.o $(ODIR)/msttypes.o $(ODIR)/mstoptim.o $(ODIR)/msttransforms.o $(ODIR)/mstlinalg.o
	ar rs lib/libmstfasst.a $(ODIR)/mstfasst.o $(ODIR)/msttypes.o $(ODIR)/mstsequence.o $(ODIR)/msttransforms.o


python: $(LDIR)/mstpython.so

$(LDIR)/mstpython.so: $(ODIR)/mstpython.o
	$(CC) $(PYLIB) -Wl,-rpath,$(PYLIBPATH) -shared $< -o $@

$(ODIR)/mstpython.o: $(SDIR)/mstpython.cpp makefile
	$(CC) $(PYOPTS) -c $< -o $@



.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so lib/*.a
