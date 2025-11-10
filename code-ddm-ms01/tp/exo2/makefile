GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../eigen/ -I../../femtool/

exo2: exo2.o
	$(GCC) $(FLAGS) exo2.o -o exo2

exo2.o: main.cpp
	$(GCC) $(FLAGS) $(INC) -c main.cpp -o exo2.o

clean:
	rm *.o
