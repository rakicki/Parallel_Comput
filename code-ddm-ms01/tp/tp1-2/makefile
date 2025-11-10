GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../eigen/ -I../../femtool/

tp1-2: tp1-2.o
	$(GCC) $(FLAGS) tp1-2.o -o tp1-2

tp1-2.o: main.cpp
	$(GCC) $(FLAGS) $(INC) -c main.cpp -o tp1-2.o

clean:
	rm *.o
