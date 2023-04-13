#Compiler
CC=g++-8

#MKL
CPPFLAGS_EXTRA+=-I${MKLROOT}/include
LDFLAGS_EXTRA+=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 

#Debugging flags
CPPFLAGS_EXTRA+= -Wno-unused-parameter -Wno-useless-cast -Wno-unused-variable -Wno-pragmas -Wno-shadow  -Wno-duplicated-branches -Wno-old-style-cast

