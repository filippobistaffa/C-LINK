.PHONY:

all:
	g++ -Wall -Ofast -march=native -funroll-loops -funsafe-loop-optimizations -falign-functions=16 -falign-loops=16 *.c -o clink

run:
	./clink
