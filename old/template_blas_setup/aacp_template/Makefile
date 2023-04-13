# Define compiler and default c++-flags
CC=g++-8
CPPFLAGS=-Wall -g1 -m64 -std=c++17
CPPFLAGS_EXTRA=
LDFLAGS=-lpthread
LDFLAGS_EXTRA=

# Define additional search path for include files
CPPFLAGS_EXTRA+=-I$(MKLROOT)/include/

# Define library search path
LDFLAGS_EXTRA+=-L$(MKLROOT)/lib/ -L$(MKLROOT)/../../compiler/2021.1.1/mac/compiler/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_avx2 -liomp5 -lm -ldl

# Define explicit library identifiers
LDFLAGS_EXTRA+=

# get directory Makefile is located in, we use it as default search path
CPPFLAGS_EXTRA+=-I$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

main : main.o
	$(CC) $(CPPFLAGS) $(CPPFLAGS_EXTRA) $@.o -o $@ $(LDFLAGS) $(LDFLAGS_EXTRA)

main.o : main.cpp la_wrapper.h la_base_obj.h la_operations.h
	$(CC) -c $(CPPFLAGS) $(CPPFLAGS_EXTRA) $(@:.o=.cpp) -o $@