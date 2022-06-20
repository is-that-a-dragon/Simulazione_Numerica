CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex9.exe : random.o ex9.o city.o
	$(CC) random.o ex9.o city.o -o ex9.exe

ex9.o : ex9.cpp
	$(CC) -c ex9.cpp -o ex9.o $(CFLAGS)

city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o prova.dat