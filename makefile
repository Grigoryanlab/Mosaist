all:
	g++ -g -c mstlib.cpp
	g++ -g mstlib.o test.cpp -o test
