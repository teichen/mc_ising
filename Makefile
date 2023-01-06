SHELL = /bin/sh

# make all
# make clean
# make glauber_ising
# make snapshot

OBJS_MC = main.o GlauberIsing.o SimSpace.o MonteCarlo.o IsingModel.o
OBJS_SNAP = snapshot.o Visualizer.o SimSpace.o MonteCarlo.o
CFLAGS = -g -O0
CC = clang++
INCLUDES = 
LIBS = -L/usr/local/Cellar/gperftools/2.10/lib -lprofiler

all:glauber_ising snapshot

glauber_ising:${OBJS_MC}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_MC} ${LIBS}

snapshot:${OBJS_SNAP}
	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS_SNAP} ${LIBS}

clean:
	-rm -f *.o core *.core *.dat

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<


