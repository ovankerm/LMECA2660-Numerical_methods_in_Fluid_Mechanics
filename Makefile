CC = clang
CFLAGS = -O3 -Wall #-g -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls
TARGET = run
SOURCES = main.c src/functions.o

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

exec : ${SOURCES}
	${CC} ${CFLAGS} -o ${TARGET} $^

clean : 
	rm -f ${TARGET}
	rm -f src/*.o

cleanFiles : 
	rm -f output/*.txt

all :
	make exec
	./${TARGET}
	make clean