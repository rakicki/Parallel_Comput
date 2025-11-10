GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../eigen/ -I../../femtool/

tp1: tp1.o
	$(GCC) $(FLAGS) tp1.o -o tp1

tp1.o: main.cpp
	$(GCC) $(FLAGS) $(INC) -c main.cpp -o tp1.o

clean:
	rm *.o
