# ========== SETUP SOURCE FILES ==========

EXCS = main
SRCS = $(wildcard *.cpp)

# ========== SETUP COMPILER FLAGS ==========

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

# ==========  SETUP DEPENDENCY FILE GENERATION ==========

# define dependency directory
DEP_DIR := .deps

# define flags for dependency file generation
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.d

# define target rule for object file that also generates dependency files
%.o : %.cpp
%.o : %.cpp $(DEP_DIR)/%.d | $(DEP_DIR)
	@echo building $@
	$(CC) $(DEPFLAGS) $(CPPFLAGS) $(CPPFLAGS_EXTRA) -c $(LDFLAGS) $(LDFLAGS_EXTRA) $< 

# define target rule for dependency directory
$(DEP_DIR):
	@mkdir -p $@

# and of course we also need a list of dependency files
DEPFILES := $(SRCS:%.cpp=$(DEP_DIR)/%.d)

# make every dependeny file a target
$(DEPFILES):

# and include the dependency files
include $(wildcard $(DEPFILES))

# ==========  BUILD TARGETS ==========

# target for actually building executables
$(EXCS): $(EXCS:%=%.o)
	@echo building $@ from $<
	$(CC) $(CPPFLAGS) $(CPPFLAGS_EXTRA) $@.o -o $@ $(LDFLAGS) $(LDFLAGS_EXTRA)

clean-excs:
	@echo removing
	@echo $(EXCS) | xargs rm -v

clean-deps:
	@echo removing
	@find $(DEP_DIR) -iname "*.d" -delete | xargs rm -v
	@rm -rv $(DEP_DIR)

clean-objects:
	@echo removing
	@find . -iname "*.o" -type f | xargs rm -v
	@find . -iname "*.dSYM" -type d | xargs rm -rv

.PHONY: build clean clean-excs clean-deps clean-objects

build : $(EXCS)
clean : clean-deps clean-objects clean-excs
