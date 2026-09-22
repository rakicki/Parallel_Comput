GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/ 

tp2-2: exo2.o
	$(GCC) $(FLAGS) exo2.o -o tp2-2

exo2.o: tp2-exo2.cpp
	$(GCC) $(FLAGS) $(INC) -c tp2-exo2.cpp -o exo2.o

clean:
	rm *.o
