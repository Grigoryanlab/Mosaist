CC = g++
CPPFLAGS=-O3 -std=c++11
ODIR = objs
SDIR = src
INCDIR = include

SOURCE  = mstfuser mstoptim mstlinalg mstoptions msttypes msttransforms mstrotlib mstmagic mstcondeg mstsequence mstsystem
OBJECTS = $(patsubst %,$(ODIR)/%.o, $(SOURCE))
PROGRAMS = master createPDS parsePDS

all: setup $(OBJECTS) libs tests
objs: $(OBJECTS)

setup:
	mkdir -p $(ODIR)
	mkdir -p bin
	mkdir -p lib

$(OBJECTS): objs/%.o : src/%.cpp include/%.h
	$(CC) $(CPPFLAGS) -I$(INCDIR) -c -o $@ $<

tests: $(OBJECTS)
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o tests/test.cpp -o bin/test
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o objs/mstlinalg.o tests/testTransforms.cpp -o bin/testTransforms
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o objs/mstrotlib.o tests/testRotlib.cpp -o bin/testRotlib
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o objs/mstmagic.o tests/testTERMUtils.cpp -o bin/testTERMUtils
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o objs/mstrotlib.o objs/msttransforms.o objs/mstcondeg.o objs/mstoptions.o objs/mstsystem.o tests/testConFind.cpp -o bin/testConFind
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o objs/mstrotlib.o objs/msttransforms.o objs/mstcondeg.o objs/mstsystem.o tests/findBestFreedom.cpp -o bin/findBestFreedom
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o objs/mstoptim.o objs/mstfuser.o objs/mstlinalg.o tests/testFuser.cpp -o bin/testFuser
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o objs/mstoptions.o tests/testClusterer.cpp -o bin/testClusterer

libs: $(OBJECTS)
	ar rs lib/libmst.a objs/mstoptions.o objs/msttypes.o objs/mstsequence.o objs/mstsystem.o
	ar rs lib/libmsttrans.a objs/msttransforms.o
	ar rs lib/libmstmagic.a objs/msttypes.o objs/mstmagic.o
	ar rs lib/libmstcondeg.a objs/mstcondeg.o objs/mstrotlib.o objs/msttransforms.o
	ar rs lib/libmstlinalg.a objs/mstlinalg.o
	ar rs lib/libmstoptim.a objs/mstoptim.o objs/mstlinalg.o
	ar rs lib/libmstfuser.a objs/mstfuser.o objs/msttypes.o objs/mstoptim.o
	$(CC) $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o objs/mstoptim.o objs/mstfuser.o objs/mstlinalg.o tests/fuserTest4.cpp -o bin/fuserTest4


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
