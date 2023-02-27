CC = gcc
CFLAGS = -Wall

hello : main.c
	${CC} ${CFLAGS} main.c -o main

clean : 
	rm -f main

all :
	make hello
	./main
	make clean