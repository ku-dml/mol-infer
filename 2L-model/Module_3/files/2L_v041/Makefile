# Makefile
CC = g++
CFLAGS = -O2 -Wall -std=c++11
RM = rm -f

# cleaning rules.
clean:
	$(RM) *.o *~ \#*

### v041 ###
OBJ = fv_2L.o fv_common.o compute_fc.o
HEADER = fv_def.h

FV_2L_V041: $(OBJ) $(HEADER)
	$(CC) $(CFLAGS) -o FV_2L_V041 $(OBJ)





.cpp.o:
	$(CC) $(CFLAGS) -c $<
