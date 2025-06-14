OBJS = main.o 

#This is for mac osx
CC = g++-mp-14 -O3 -fopenmp
DEBUG = 
CFLAGS = -c $(DEBUG) -I/opt/local/include
LFLAGS = $(DEBUG) -L/opt/local/lib -lm -lpthread 


point: $(OBJS)
	$(CC) -o point $(OBJS) $(LFLAGS)

main.o: main.cpp 
	$(CC) $(CFLAGS) main.cpp

#read.o: read.cpp 
#	$(CC) $(CFLAGS) read.cpp

clean:
	rm -rfv *.o ./output/ *~ point

