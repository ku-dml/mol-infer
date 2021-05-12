# Makefile

CC = g++
CFLAGS = -O3 -std=c++11
RM = rm -f

# Determine the OS
# Starting with default Unix
OS = Unix
ifneq ($(shell uname -a | grep -i Darwin),)
	OS = osx
endif
ifneq ($(shell uname -a | grep -i Windows),)
	OS = windows
endif
ifneq ($(shell uname -a | grep -i Linux),)
	OS = linux
endif
ifneq ($(shell uname -a | grep -i Cygwin),)
	OS = cygwin
endif

# Executable files prefix
BIN = executable_programs

# Compiled Program folder
SUBF = $(BIN)/$(OS)


# Compile everything with a single make
all: directories generate_isomers generate_partition FV_ec FV_proj CHECKER clean
	
clean:
	$(RM) *.o *~ \#*

generate_isomers: ./source_code/Module_4/main/generate_isomers.cpp
	$(CC) $(CFLAGS) -o ./$(SUBF)/generate_isomers ./source_code/Module_4/main/generate_isomers.cpp

generate_partition: ./source_code/Module_4/main/generate_partition.cpp
	$(CC) $(CFLAGS) -o ./$(SUBF)/generate_partition ./source_code/Module_4/main/generate_partition.cpp

FV_ec: ./source_code/Module_1/fv_ec.cpp
	$(CC) $(CFLAGS) -o ./$(SUBF)/FV_ec ./source_code/Module_1/fv_ec.cpp

FV_proj: ./source_code/Module_1/fv_proj.cpp
	$(CC) $(CFLAGS) -o ./$(SUBF)/FV_proj ./source_code/Module_1/fv_proj.cpp

CHECKER: ./source_code/Module_1/cycle_checker.cpp
	$(CC) $(CFLAGS) -o ./$(SUBF)/CHECKER ./source_code/Module_1/cycle_checker.cpp

directories: $(SUBF)

$(SUBF):
	mkdir -p $(SUBF)
