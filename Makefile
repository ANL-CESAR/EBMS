CC=mpicc
CFLAGS=-g -O3 -std=gnu99

OTHER_OBJS=comm.o runtime_parameters.o matrix.o process_input.o

all : ebmc-rget ebmc-iallgather
ebmc-rget : ebmc-rget.o $(OTHER_OBJS)
ebmc-iallgather : ebmc-iallgather.o $(OTHER_OBJS)

clean:
	rm -f *.o *~ ebmc-rget ebmc-iallgather
