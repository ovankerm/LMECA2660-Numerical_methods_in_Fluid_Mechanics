CC = gcc
CFLAGS = -Wall
TARGET = run
SOURCES = main.c src/functions.o

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

hello : ${SOURCES}
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