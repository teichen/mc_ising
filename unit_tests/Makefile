SHELL = /bin/sh

# make all
# make model_test
# make ising_test
# make clean

OBJS_MODEL_TEST = model_test.o ../SimSpace.o ../LGModel.o
OBJS_ISING_TEST = ising_test.o ../ising.o ../SimSpace.o ../MonteCarlo.o ../LGModel.o
CFLAGS =
CC = clang++
INCLUDES =
LIBS =

all:model_test ising_test

model_test:${OBJS_MODEL_TEST}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_MODEL_TEST} ${LIBS}

ising_test:${OBJS_ISING_TEST}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_ISING_TEST} ${LIBS}

clean:
	-rm -f *.o core *.core *.dat

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<


