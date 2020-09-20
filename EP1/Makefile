CC      = gcc
CFLAGS  = -Wall -std=c99 -g -O0 -pedantic -Wno-unused-result -Wno-unused-function
RM      = rm
OBJS = util.o ep1.o
#------------------------------------------------

# chama o linker
ep1: $(OBJS)
	$(CC)  $(OBJS) -o ep1 -lm && rm *.o

# gera o objeto do EP1
ep1.o: ep1.c ep1.h util.o
	$(CC) $(CFLAGS) -c ep1.c

# gera o objeto das funcoes de uso geral
util.o: util.c util.h
	$(CC) $(CFLAGS) -c util.c

clean:
	$(RM) *.o ep1

