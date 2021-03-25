# Makefile

CC = g++
CFLAGS = -O3 -std=c++11
RM = rm -f

clean:
	$(RM) *.o *~ \#*

distclean:
	$(RM) *.o *~ \#* generate_isomers generate_partition

generate_isomers: ./main/generate_isomers.cpp
	$(CC) $(CFLAGS) -o generate_isomers ./main/generate_isomers.cpp

generate_partition: ./main/generate_partition.cpp
	$(CC) $(CFLAGS) -o generate_partition ./main/generate_partition.cpp

