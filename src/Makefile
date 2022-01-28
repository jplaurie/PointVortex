OBJS = main.o step.o initialize.o output.o hamiltonian.o remove.o

# to use openmp include -fopenmp
#This is for linux machines
#CC = g++ -O3 -std=c++14 -fopenmp
#DEBUG = 
#CFLAGS = -c $(DEBUG) -I/usr/local/include
#LFLAGS = $(DEBUG) -L/usr/local/lib -lm -lpthread -lopenblas -llapack -larmadillo

#CC = g++ -O3 -std=c++14 -fopenmp
#DEBUG =
#CFLAGS = -c $(DEBUG) -I/apps/lib64/openblas/include -I/home/lauriej/libs/include -I/home/lauriej/libs/usr/local/include
#LFLAGS = $(DEBUG) -L/apps/lib64/openblas/lib -L/home/lauriej/libs/lib -L/home/lauriej/libs/usr/local/lib64 -lm -lpthread -larmadillo -lopenblas -llapack



#This is for mac osx
CC = g++-mp-10 -O2 -std=c++14 -fopenmp
DEBUG = 
CFLAGS = -c $(DEBUG) -I/opt/local/include
LFLAGS = $(DEBUG) -L/opt/local/lib -lm -lpthread -larmadillo


point: $(OBJS)
	$(CC) -o point $(OBJS) $(LFLAGS)

main.o: main.cpp Const.h
	$(CC) $(CFLAGS) main.cpp

step.o: step.cpp Const.h
	$(CC) $(CFLAGS) step.cpp

initialize.o: initialize.cpp Const.h
	$(CC) $(CFLAGS) initialize.cpp

output.o: output.cpp Const.h
	$(CC) $(CFLAGS) output.cpp

hamiltonian.o: hamiltonian.cpp Const.h
	$(CC) $(CFLAGS) hamiltonian.cpp 

remove.o: remove.cpp Const.h
	$(CC) $(CFLAGS) remove.cpp

clean:
	rm -rfv *.o ./Output/ *~ point

