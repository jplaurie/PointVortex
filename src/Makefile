OBJS = main.o read.o

#This is for mac osx
CC = g++-mp-14 -O3 -fopenmp
DEBUG = 
CFLAGS = -c $(DEBUG) -I/opt/local/include
LFLAGS = $(DEBUG) -L/opt/local/lib -lm -lpthread 


point: $(OBJS)
	$(CC) -o point $(OBJS) $(LFLAGS)

main.o: main.cpp read.h vortex.h
	$(CC) $(CFLAGS) main.cpp

read.o: read.cpp read.h
	$(CC) $(CFLAGS) read.cpp

clean:
	rm -rfv *.o ./output/ *~ point

