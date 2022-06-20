CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex1.exe : random.o IntegralMC.o ex1.o
	$(CC) random.o IntegralMC.o ex1.o -o ex1.exe
ex2.exe : random.o ex2.o
	$(CC) random.o ex2.o -o ex2.exe

ex1.o : ex1.cpp libre.h funzioni.h IntegralMC.h
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)
ex2.o : ex2.cpp libre.h
	$(CC) -c ex2.cpp -o ex2.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
IntegralMC.o : IntegralMC.cpp IntegralMC.h
	$(CC) -c IntegralMC.cpp -o IntegralMC.o $(CFLAGS)
clean :
	rm *.o *.exe