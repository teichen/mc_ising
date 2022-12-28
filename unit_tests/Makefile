SHELL = /bin/sh

# make test
# make clean

OBJS_TEST = model_test.o ../SimSpace.o
CFLAGS =
CC = clang++
INCLUDES =
LIBS =

all:test

test:${OBJS_TEST}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_TEST} ${LIBS}

clean:
	-rm -f *.o core *.core

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<


