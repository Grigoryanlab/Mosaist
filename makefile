CPPFLAGS=-g

all:
	mkdir -p objs; cd objs; g++ $(CPPFLAGS) -I../include -c ../src/msttypes.cpp ../src/msttransforms.cpp
	mkdir -p bin; g++ $(CPPFLAGS) -I./include objs/msttypes.o tests/test.cpp -o bin/test
	g++ $(CPPFLAGS) -I./include objs/msttypes.o objs/msttransforms.o tests/testTransforms.cpp -o bin/testTransforms
	mkdir -p lib; ar rs lib/libmst.a objs/msttypes.o
