all: interpLib.o
	g++ main.cpp interpLib.o -std=c++11 -g -Wall -o int  -lboost_system -lboost_filesystem

interpLib.o: splineInterp/interpLib.cpp
	g++ -c -std=c++11 -g -Wall splineInterp/interpLib.cpp
	
clean:
	rm int *.o
