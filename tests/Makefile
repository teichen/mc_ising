SHELL = /bin/sh

# make all
# make ising_test
# make glauber_test
# make clean

OBJS_MODEL_TEST = ising_test.o ../SimSpace.o ../IsingModel.o
OBJS_GLAUBER_TEST = glauber_test.o ../GlauberIsing.o ../SimSpace.o ../MonteCarlo.o ../IsingModel.o
CFLAGS = -g -O0
CC = clang++
INCLUDES =
LIBS = -L/usr/local/lib -L/usr/local/Cellar/gperftools/2.10/lib -lprofiler

all:ising_test glauber_test

ising_test:${OBJS_MODEL_TEST}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_MODEL_TEST} ${LIBS}

glauber_test:${OBJS_GLAUBER_TEST}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_GLAUBER_TEST} ${LIBS}

clean:
	-rm -f *.o core *.core *.dat

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<


