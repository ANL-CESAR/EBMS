EXE=ebmc
CC=mpicc -g -O3 -std=gnu99
#CC=tau_cc.sh -std=gnu99
${EXE} : ${EXE}.o runtime_parameters.o matrix.o process_input.o
	${CC} -o ${EXE} ${EXE}.o runtime_parameters.o matrix.o process_input.o
matrix.o : matrix.c matrix.h
	${CC} -c matrix.c -I.
runtime_parameters.o : runtime_parameters.c runtime_parameters.h
	${CC} -c runtime_parameters.c -I.
process_input.o : process_input.c
	${CC} -c process_input.c -I.
${EXE}.o : ${EXE}.c
	${CC} -c $^ -I.
clean:
	rm -f *.o *~ ${EXE} *.cobaltlog #*output *error
