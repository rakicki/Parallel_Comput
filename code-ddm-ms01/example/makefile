GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../eigen/ -I../femtool/

example: example.o
	$(GCC) $(FLAGS) example.o -o example

example.o: example.cpp
	$(GCC) $(FLAGS) $(INC) -c example.cpp -o example.o

clean:
	rm *.o
