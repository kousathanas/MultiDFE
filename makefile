CC = gcc
CFLAGS = -g -O3 -lm -lgsl -lgslcblas -w

MultiDFE: 
	$(CC) -o MultiDFE *.c $(CFLAGS)

clean: 
		rm MultiDFE *.o
