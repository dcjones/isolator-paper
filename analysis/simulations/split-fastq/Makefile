
CFLAGS=-O3
OBJ=common.o parse.o split-fastq.o

split-fastq : $(OBJ)
	gcc ${CFLAGS} -o $@ $^

clean :
	rm -f *.o split-fastq

