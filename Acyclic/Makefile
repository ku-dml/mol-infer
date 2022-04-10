CC=g++
CFLAGS=-O3 -std=c++11

all: bin/fv4_in_ex bin/2-branches bin/3-branches

bin/fv4_in_ex: ./Module_1/fv4_in_ex.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/2-branches: ./Module_4/2-branches/2-branches.cpp \
                ./Module_4/2-branches/include/cross_timer.h \
                ./Module_4/2-branches/include/data_structures.hpp \
                ./Module_4/2-branches/include/debug.h \
                ./Module_4/2-branches/include/fringe_tree.hpp \
                ./Module_4/2-branches/include/tools.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/3-branches: ./Module_4/3-branches/3-branches.cpp \
                ./Module_4/3-branches/include/cross_timer.h \
                ./Module_4/3-branches/include/data_structures.hpp \
                ./Module_4/3-branches/include/debug.h \
                ./Module_4/3-branches/include/fringe_tree.hpp \
                ./Module_4/3-branches/include/tools.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -rf bin