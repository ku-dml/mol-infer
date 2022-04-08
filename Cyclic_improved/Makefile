CC=g++
CFLAGS=-O3 -std=c++11

all: bin/CHECKER bin/FV_ec bin/FV_proj bin/fv4_cyclic_stdout \
     bin/generate_isomers bin/generate_partition

bin/CHECKER: ./Module_1/files/cycle_checker.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/FV_ec: ./Module_1/files/fv_ec.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/FV_proj: ./Module_1/files/fv_proj.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/fv4_cyclic_stdout: ./Module_3/files/fv4_cyclic_stdout.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/generate_isomers: ./Module_4/files/main/generate_isomers.cpp \
                      ./Module_4/files/include/chemical_graph.hpp \
                      ./Module_4/files/include/cross_timer.h \
                      ./Module_4/files/include/data_structures.hpp \
                      ./Module_4/files/include/fringe_tree.hpp \
                      ./Module_4/files/include/tools.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/generate_partition: ./Module_4/files/main/gen_partition.cpp \
                        ./Module_4/files/include/chemical_graph.hpp \
                        ./Module_4/files/include/cross_timer.h \
                        ./Module_4/files/include/data_structures.hpp \
                        ./Module_4/files/include/fringe_tree.hpp \
                        ./Module_4/files/include/tools.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -rf bin