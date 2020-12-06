CC      = gcc
CFLAGS  = -Wall -std=c99 -g -O0 -pedantic -Wno-unused-result -Wno-unused-function
RM      = rm
OBJS    = util.o

all: $(OBJS)
	$(CC) $(OBJS) -DTIME_BENCHMARK -o ep3 -lm ep3.c

debug: $(OBJS)
	$(CC) $(OBJS) -DDEBUG -o ep3 -lm ep3.c

util.o: util.c util.h
	$(CC) $(CFLAGS) -c util.c

clean:
	$(RM) *.o ep3

