all:
	mkdir -p objs; cd objs; g++ -O3 -I../include -c ../src/msttypes.cpp ../src/msttransforms.cpp
	mkdir -p bin; g++ -O3 -I./include objs/msttypes.o tests/test.cpp -o bin/test
	mkdir -p lib; ar rs lib/libmst.a objs/msttypes.o
