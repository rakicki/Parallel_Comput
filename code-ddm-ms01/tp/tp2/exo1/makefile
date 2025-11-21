GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/ 

tp2-1: exo1.o
	$(GCC) $(FLAGS) exo1.o -o tp2-1

exo1.o: main.cpp
	$(GCC) $(FLAGS) $(INC) -c main.cpp -o exo1.o

clean:
	rm *.o
