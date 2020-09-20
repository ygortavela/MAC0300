CC      = gcc
CFLAGS  = -Wall -std=c99 -g -O0 -pedantic -Wno-unused-result -Wno-unused-function
RM      = rm
OBJS    = util.o ep1.o

ep1: $(OBJS)
	$(CC) $(OBJS) -o ep1 -lm && rm *.o

ep1.o: ep1.c ep1.h util.o
	$(CC) $(CFLAGS) -c ep1.c

util.o: util.c util.h
	$(CC) $(CFLAGS) -c util.c

clean:
	$(RM) *.o ep1

