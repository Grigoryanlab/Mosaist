CPPFLAGS=-g -std=c++11

all:
	mkdir -p objs; cd objs; g++ $(CPPFLAGS) -I../include -c ../src/msttypes.cpp ../src/msttransforms.cpp ../src/mstrotlib.cpp ../src/mstmagic.cpp ../src/mstcondeg.cpp
	mkdir -p bin; g++ $(CPPFLAGS) -I./include objs/msttypes.o tests/test.cpp -o bin/test
	g++ $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o tests/testTransforms.cpp -o bin/testTransforms
	g++ $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o objs/mstrotlib.o tests/testRotlib.cpp -o bin/testRotlib
	g++ $(CPPFLAGS) -I./include objs/msttypes.o objs/mstmagic.o tests/testTERMUtils.cpp -o bin/testTERMUtils
	g++ $(CPPFLAGS) -I./include objs/msttypes.o objs/mstrotlib.o objs/msttransforms.o objs/mstcondeg.o tests/testConFind.cpp -o bin/testConFind
	mkdir -p lib
	ar rs lib/libmst.a objs/msttypes.o
	ar rs lib/libmstmagic.a objs/msttypes.o objs/mstmagic.o
	ar rs lib/libmstcondeg.a objs/mstcondeg.o objs/mstrotlib.o objs/msttransforms.o
