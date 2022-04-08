CC=g++
CFLAGS=-O3 -std=c++11

all: bin/2L_FV bin/generate_isomers

bin/2L_FV: ./Module_1/files/2L_fv_main.cpp \
           ./Module_1/files/fringe/ChemicalGraph.hpp \
		   ./Module_1/files/fringe/commonused.hpp \
		   ./Module_1/files/fringe/cross_timer.h \
		   ./Module_1/files/fringe/RootedGraph.hpp \
		   ./Module_1/files/fringe/TopologyGraph.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/generate_isomers: ./Module_4/main/generate_isomers.cpp \
                      ./Module_4/include/chemical_graph.hpp \
					  ./Module_4/include/cross_timer.h \
					  ./Module_4/include/data_structures.hpp \
					  ./Module_4/include/fringe_tree.hpp \
					  ./Module_4/include/tools.hpp \
					  ./Module_4/include/tree_signature.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -rf bin