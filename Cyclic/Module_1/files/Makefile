# Makefile

CC = g++
CFLAGS = -O2 -Wall -std=c++11
RM = rm -f

clean:
	$(RM) *.o *~ \#*

distclean:
	$(RM) *.o *~ \#* FV_ec FV_proj CHECKER

FV_ec: fv_ec.cpp
	$(CC) $(CFLAGS) -o FV_ec fv_ec.cpp

FV_proj: fv_proj.cpp
	$(CC) $(CFLAGS) -o FV_proj fv_proj.cpp

CHECKER: cycle_checker.cpp
	$(CC) $(CFLAGS) -o CHECKER cycle_checker.cpp
