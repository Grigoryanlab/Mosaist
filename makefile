all:
	cd objs; g++ -O3 -I../include -c ../src/msttypes.cpp
	g++ -O3 -I./include objs/msttypes.o tests/test.cpp -o bin/test
	ar rs lib/libmst.a objs/msttypes.o
