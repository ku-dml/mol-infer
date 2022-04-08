all:
	$(MAKE) -C Acyclic
	$(MAKE) -C Cyclic
	$(MAKE) -C Cyclic_improved
	$(MAKE) -C 2L-model
	$(MAKE) -C 2LMM-LLR

.PHONY: clean
clean:
	$(MAKE) -C Acyclic clean
	$(MAKE) -C Cyclic clean
	$(MAKE) -C Cyclic_improved clean
	$(MAKE) -C 2L-model clean
	$(MAKE) -C 2LMM-LLR clean