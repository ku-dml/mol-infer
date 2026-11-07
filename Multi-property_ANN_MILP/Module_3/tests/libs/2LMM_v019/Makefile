# Makefile
CC = g++
CFLAGS = -O2 -Wall -std=c++20
RM = rm -f

# cleaning rules.
clean:
	$(RM) *.o *~ \#*

FV_2LMM_V019: fv_2LMM.cpp fv_common.hpp compute_fc.hpp fv_def.h
	$(CC) $(CFLAGS) -o FV_2LMM_V019 fv_2LMM.cpp

.cpp.o:
	$(CC) $(CFLAGS) -c $<
