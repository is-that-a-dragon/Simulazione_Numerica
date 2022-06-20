CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex1.exe : random.o ex1.o
	$(CC) random.o ex1.o -o ex1.exe
ex2.exe : random.o ex2.o
	$(CC) random.o ex2.o -o ex2.exe
ex3.exe : random.o ex3.o
	$(CC) random.o ex3.o -o ex3.exe

ex3.o : ex3.cpp libre.h
	$(CC) -c ex3.cpp -o ex3.o $(CFLAGS)
ex2.o : ex2.cpp libre.h
	$(CC) -c ex2.cpp -o ex2.o $(CFLAGS)
ex1.o : ex1.cpp libre.h
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe