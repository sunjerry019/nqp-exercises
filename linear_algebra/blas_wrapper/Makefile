# ========== SETUP SOURCE FILES ==========

EXCS = main
SRCS = $(wildcard *.cpp)

# ========== SETUP COMPILER FLAGS ==========

# Define compiler and default c++-flags
CC=g++-8
CPPFLAGS=-Wall -g1 -m64 -std=c++17
CPPFLAGS_EXTRA=
LDFLAGS=
LDFLAGS_EXTRA=

make_inc=../make.inc

# include optional make.inc file
include $(make_inc)

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
	$(CC) $(DEPFLAGS) $(CPPFLAGS) $(CPPFLAGS_EXTRA) -c $< 

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
$(EXCS): $(SRCS:%.cpp=%.o)
	@echo building $@ from $<
	$(CC) $(CPPFLAGS) $(CPPFLAGS_EXTRA) $@.o -o $@ $(LDFLAGS) $(LDFLAGS_EXTRA)

clean-excs:
	@echo $(EXCS) | xargs rm -v || true

clean-deps:
	@rm -rv $(DEP_DIR) || true

clean-objects:
	@find . -iname "*.o" -type f | xargs rm -v || true

.PHONY: build clean clean-excs clean-deps clean-objects

build : $(EXCS)
clean : clean-deps clean-objects clean-excs
