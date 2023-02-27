CC = gcc
CFLAGS = -Wall
TARGET = run

hello : main.c src/functions.o
	${CC} ${CFLAGS} -o ${TARGET} $^

clean : 
	rm -f ${TARGET}
	rm -f src/*.o

cleanFiles : 
	rm -f output/*.txt

all :
	make hello
	./${TARGET}
	make clean