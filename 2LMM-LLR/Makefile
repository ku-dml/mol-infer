CC=g++
CFLAGS=-O3 -std=c++11

all: bin/FV_2LMM_V018 bin/generate_isomers

bin/FV_2LMM_V018: ./src/Module_1/fv_2LMM.cpp \
                  ./src/Module_1/compute_fc.hpp \
				  ./src/Module_1/fv_common.h \
				  ./src/Module_1/fv_common.hpp \
				  ./src/Module_1/fv_def.h \
				  ./src/Module_1/fringe/ChemicalGraph.hpp \
				  ./src/Module_1/fringe/commonused.hpp \
				  ./src/Module_1/fringe/cross_timer.h \
				  ./src/Module_1/fringe/RootedGraph.hpp \
				  ./src/Module_1/fringe/TopologyGraph.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

bin/generate_isomers: ./src/Module_4/main/generate_isomers.cpp \
                      ./src/Module_4/include/chemical_graph.hpp \
					  ./src/Module_4/include/cross_timer.h \
					  ./src/Module_4/include/data_structures.hpp \
					  ./src/Module_4/include/fringe_tree.hpp \
					  ./src/Module_4/include/tools.hpp \
					  ./src/Module_4/include/tree_signature.hpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -rf bin