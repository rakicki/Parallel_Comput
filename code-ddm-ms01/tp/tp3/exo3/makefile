GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/ 

tp3-exo3: exo3.o
	$(GCC) $(FLAGS) exo3.o -o tp3-exo3

exo3.o: tp3-exo3.cpp
	$(GCC) $(FLAGS) $(INC) -c tp3-exo3.cpp -o exo3.o
clean:
	rm *.o
