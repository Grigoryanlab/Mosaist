all:
	g++ -g -c msttypes.cpp
	g++ -g msttypes.o test.cpp -o test
	ar rs libmst.a msttypes.o
