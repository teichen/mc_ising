SHELL = /bin/sh

# make all
# make install
# make clean
# make isingSim
# make snapshot

OBJS_MC = main.o ising.o SimSpace.o MonteCarlo.o LGModel.o
OBJS_SNAP = snapshot.o Visualizer.o SimSpace.o MonteCarlo.o
CFLAGS = -g -openmp
CC = icpc
INCLUDES = 
LIBS = 

all:isingSim snapshot

isingSim:${OBJS_MC}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_MC} ${LIBS}

snapshot:${OBJS_SNAP}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_SNAP} ${LIBS}

clean:
	-rm -f *.o core *.core

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<


