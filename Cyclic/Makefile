CC=g++
CFLAGS=-O3 -std=c++11

all: bin/CHECKER bin/FV_ec bin/FV_proj \
     bin/generate_isomers bin/generate_partition

bin/CHECKER: ./src/Module_1/cycle_checker.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/FV_ec: ./src/Module_1/fv_ec.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/FV_proj: ./src/Module_1/fv_proj.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/generate_isomers: ./src/Module_4/main/generate_isomers.cpp \
                      ./src/Module_4/include/chemical_graph.hpp \
                      ./src/Module_4/include/cross_timer.h \
                      ./src/Module_4/include/data_structures.hpp \
                      ./src/Module_4/include/fringe_tree.hpp \
                      ./src/Module_4/include/tools.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/generate_partition: ./src/Module_4/main/generate_partition.cpp \
                        ./src/Module_4/include/chemical_graph.hpp \
                        ./src/Module_4/include/cross_timer.h \
                        ./src/Module_4/include/data_structures.hpp \
                        ./src/Module_4/include/fringe_tree.hpp \
                        ./src/Module_4/include/tools.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -rf bin