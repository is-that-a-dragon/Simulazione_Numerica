CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex10.exe : random.o ex10.o city.o
	mpicxx random.o ex10.o city.o -o ex10.exe

ex10.o : ex10.cpp city.h
	mpicxx -c ex10.cpp -o ex10.o $(CFLAGS)

TSP_int.exe : random.o TSP_int.o city.o
	$(CC) random.o TSP_int.o city.o -o TSP_int.exe

TSP_int.o : TSP_int.cpp city.h
	$(CC) -c TSP_int.cpp -o TSP_int.o $(CFLAGS)

city.o : city.cpp city.h
	mpicxx -c city.cpp -o city.o $(CFLAGS)

random.o : random.cpp random.h
	mpicxx -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o prova.dat average_fit.dat