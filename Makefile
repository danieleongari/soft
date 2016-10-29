HEADERS = SOFT.h SOFT.h

default: SOFT.x

SOFT.o: SOFT.c $(HEADERS)
	gcc -c SOFT.c -o SOFT.o

SOFT.x: SOFT.o
	gcc SOFT.o -o SOFT.x -lm

clean:
	-rm -f SOFT.o
	-rm -f SOFT.x
