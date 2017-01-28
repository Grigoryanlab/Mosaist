CPPFLAGS=-O3 -std=c++11

all:
	mkdir -p objs; cd objs; g++ $(CPPFLAGS) -I../include -c ../src/msttypes.cpp ../src/msttransforms.cpp ../src/mstrotlib.cpp ../src/mstmagic.cpp
	mkdir -p bin; g++ $(CPPFLAGS) -I./include objs/msttypes.o tests/test.cpp -o bin/test
	g++ $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o tests/testTransforms.cpp -o bin/testTransforms
	g++ $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o objs/mstrotlib.o tests/testRotlib.cpp -o bin/testRotlib
	g++ $(CPPFLAGS) -I./include objs/msttypes.o objs/mstmagic.o tests/testTERMUtils.cpp -o bin/testTERMUtils
	mkdir -p lib; ar rs lib/libmst.a objs/msttypes.o; ar rs lib/libmstmagic.a objs/msttypes.o objs/mstmagic.o
