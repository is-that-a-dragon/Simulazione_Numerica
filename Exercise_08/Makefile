CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex8.exe : random.o integral.o ex8.o
	$(CC) random.o integral.o ex8.o -o ex8.exe

ex8.o : ex8.cpp integral.h
	$(CC) -c ex8.cpp -o ex8.o $(CFLAGS)

integral.o : integral.cpp integral.h
	$(CC) -c integral.cpp -o integral.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o *.exe