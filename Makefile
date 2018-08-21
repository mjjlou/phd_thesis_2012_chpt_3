CC = gcc
#CFLAGS =  -O3
CFLAGS =  -g #debug  

# Libraries for linking and flags for the C preprocessor. 
LIBS = -lm 

default: gps_random

###################################################################
# Various commands
RM = rm -f

SRCS =	gasdev.c  expdev.c gps.c  nrutil.c  ran1.c gfsr8.c poidev.c
OBJS =	gasdev.o  expdev.o gps.o  nrutil.o  ran1.o gfsr8.o poidev.o
DEFS =  nrutil.h

gps_random:  $(OBJS) $(DEFS)
	$(CC) -o gps_random $(OBJS) $(LIBS)

clean: 
	-rm *.o 

.c.o: $(DEFS)
	$(CC) -c $(CFLAGS)  $< 

