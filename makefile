# this makefile uses snippets / patterns from:
#	https://stackoverflow.com/questions/2908057/can-i-compile-all-cpp-files-in-src-to-os-in-obj-then-link-to-binary-in
#	https://stackoverflow.com/questions/3261737/makefile-change-variable-value-depending-on-a-target
#	https://stackoverflow.com/questions/2205603/conditional-dependency-with-make-gmake
#	https://www.cmcrossroads.com/article/trouble-wildcard

# how to use this makefile:
#	1) run `make` or `make all` to compile and create everything
#	2) run `make clean` to delete all output files (objects, libraries, binaries), leaving just the output directories (objs, lib, bin)
#	3) run `make setup` to create the output directories (objs, lib, bin) and nothing else
#	4) run `make libs` to create all MST libraries
#	5) run `make python` to create the boost.python shared object for using MST in python
#	6) if $target is a target name (i.e. a recognized binary), run `make $target` to compile the target
#	7) if $library is a library name (i.e. a recognized MST library), run `make $library` to create the library
#	8) targets, libraries, and helper source files can be compiled directly via their pathnames if needed, e.g. `make bin/$target`, `make objs/$helper.o`, etc.

# how to maintain this makefile:
#	to add a target binary named $target
#		1) add its name to TESTS, PROGRAMS, or your own variable that gets added to TARGETS
#		2) specify its dependencies in the variable $target_DEPS (e.g. see test_DEPS below)
#	to add a helper source, i.e. a .cpp that does not have its own `main`, add its name to HELPERS
#	to add a library named $lib,
#		1) add its name to LIBRARIES
#		2) specify its dependencies in the variable $lib_DEPS (e.g. see libmst_DEPS below)
#	if any external headers or libraries become needed, add their directories to INC_DIRS and LIB_DIRS, respectively
#	note that it's assumed that target and library names are unique - if two targets, two libraries, or a target and a library have the same name, compilation will not work as expected
#	also note that it's assumed that no file in this directory has the same name as a phony target (all, clean, libs, python, setup)
#		if a file is created with one of these names, compilation will not work as expected
#	on a pedantic note, I have alphabetized the lists of targets, libraries, and dependencies to make it easier to find things
#		consider maintaining this so it's easier to determine whether a target, library, or dependency already exists!

# customizations
# define environmental variable INCLUDE_ARMA if you want to compile with Armadillo C++ linear algebra library (needed for some more complex things in mstlinalg)

# stuff meant to be regularly updated:

# flags
CC := g++
CPP_FLAGS := -std=c++11 -fPIC
DEBUG_FLAGS := -g3 #-g3 -rdynamic -gdwarf-3

# essential directories
INCD := include
OBJD := objs
LIBD := lib
BIND := bin

# source directories
SRCD		:= src
TESTD		:= tests
PROGRAMD	:= programs

# header and external library directories
INC_DIRS := $(INCD)
LIB_DIRS := 

# armadillo-dependent stuff
ifdef INCLUDE_ARMA
  CPP_FLAGS := $(CPP_FLAGS) -DARMA
  uname := $(shell uname -s)
  ifeq ($(uname),Linux)
    LDLIBS := -larmadillo -llapack -lcblas
  else
    LDLIBS := -larmadillo
  endif
  ARMA_PROGRAMS			:= chainGrow
  chainGrow_DEPS		:= msttypes mstfasst mstcondeg mstfuser mstrotlib msttransforms mstsequence mstoptim mstlinalg mstoptions mstmagic
endif

# targets and MST libraries
TESTS		:= findBestFreedom test testAutofuser testConFind testClusterer testSequence testStride testFASST testFuser testGrads testParsing testRestrictSiteAlphabet testRotlib testTERMUtils testTransforms testdTERMen testTermanal
PROGRAMS	:= findTERMs renumber TERMify subMatrix fasstDB bind analyzeLandscape extractSegments design pairEnergies search scoreStructure clusterStructs connect $(ARMA_PROGRAMS)
TARGETS		:= $(TESTS) $(PROGRAMS)
HELPERS		:= mstcondeg mstexternal mstfasst mstfuser mstlinalg mstmagic mstoptim mstoptions mstrotlib mstsequence mstsystem msttransforms msttypes msttermanal
LIBRARIES	:= libmst libmstcondeg libmstfasst libmstfasstcache libmstfuser libmstlinalg libmstmagic libmstoptim libmsttrans libdtermen

# target dependencies
findBestFreedom_DEPS	:= mstcondeg mstrotlib mstsystem msttransforms msttypes
test_DEPS			:= msttypes mstsystem
testAutofuser_DEPS		:= mstfuser mstlinalg mstoptim msttransforms msttypes
testConFind_DEPS		:= mstcondeg mstoptions mstrotlib mstsystem msttransforms msttypes
testClusterer_DEPS		:= mstoptions msttypes mstfasst msttransforms mstsequence
testSequence_DEPS		:= mstoptions msttypes mstsequence
testFASST_DEPS			:= mstfasst mstoptions mstsequence msttransforms msttypes mstsystem
testFuser_DEPS			:= mstfuser mstlinalg mstoptim msttransforms msttypes
testGrads_DEPS			:= msttypes
testParsing_DEPS		:= msttypes
testRestrictSiteAlphabet_DEPS   := msttypes mstfasst dtermen msttransforms mstsequence mstrotlib mstcondeg mstoptions mstmagic mstsystem
testRotlib_DEPS			:= mstrotlib msttransforms msttypes
testStride_DEPS			:= msttypes mstexternal mstsystem
testTERMUtils_DEPS		:= mstmagic msttypes mstcondeg mstrotlib msttransforms
testTransforms_DEPS		:= mstlinalg msttransforms msttypes
testTermanal_DEPS		:= msttermanal msttypes mstrotlib mstcondeg mstfasst mstoptions mstsequence msttransforms mstmagic
findTERMs_DEPS			:= mstfasst mstoptions mstsequence msttransforms msttypes
renumber_DEPS			:= mstsystem msttypes mstoptions
extractSegments_DEPS		:= msttypes msttransforms mstsequence mstoptions mstfasst dtermen mstcondeg mstrotlib mstmagic mstlinalg
TERMify_DEPS			:= msttypes mstfasst mstcondeg mstfuser mstrotlib msttransforms mstsequence mstoptim mstlinalg mstoptions mstmagic mstfasstcache mstsystem
bind_DEPS			:= msttypes mstfasst mstcondeg mstrotlib msttransforms mstsequence mstoptions mstmagic
connect_DEPS			:= msttypes mstfasst mstcondeg mstrotlib msttransforms mstsequence mstoptions
subMatrix_DEPS			:= msttypes mstfasst mstcondeg mstrotlib msttransforms mstsequence mstoptions
fasstDB_DEPS			:= msttypes mstfasst mstrotlib mstoptions msttransforms mstsequence mstsystem mstcondeg mstexternal
testdTERMen_DEPS		:= msttypes mstfasst dtermen msttransforms mstsequence mstrotlib mstcondeg mstoptions mstmagic mstsystem
design_DEPS			:= msttypes mstfasst dtermen msttransforms mstsequence mstrotlib mstcondeg mstoptions mstmagic mstsystem
pairEnergies_DEPS		:= msttypes mstfasst dtermen msttransforms mstsequence mstrotlib mstcondeg mstoptions mstmagic mstsystem
analyzeLandscape_DEPS		:= msttypes msttransforms mstsequence mstoptions mstfasst dtermen mstcondeg mstrotlib mstmagic
search_DEPS			:= mstfasst mstoptions mstsequence msttransforms msttypes mstsystem
scoreStructure_DEPS		:= msttermanal msttypes mstrotlib mstcondeg mstfasst mstoptions mstsequence msttransforms mstmagic
clusterStructs_DEPS		:= mstcondeg mstoptions mstrotlib mstsystem msttransforms msttypes mstrotlib mstsequence

# MST library dependencies
libmst_DEPS			:= mstoptions mstsequence mstsystem msttypes
libmstcondeg_DEPS		:= mstcondeg mstrotlib msttransforms
libmstfasst_DEPS		:= mstfasst mstsequence msttransforms msttypes
libmstfasstcache_DEPS	:= mstfasst mstfasstcache mstsequence msttransforms msttypes
libmstfuser_DEPS		:= mstfuser mstlinalg mstoptim msttransforms msttypes
libmstlinalg_DEPS		:= mstlinalg
libmstmagic_DEPS		:= msttypes mstmagic mstcondeg
libmstoptim_DEPS		:= mstoptim mstlinalg
libmsttrans_DEPS		:= msttransforms
libdtermen_DEPS			:= mstsequence msttypes dtermen

# stuff not meant to be regularly updated (that is, makefile magic that shouldn't have to be updated until a major change to the project structure is made):

# collect targets, helpers, dependencies, and associated files
TARGET_FILES := $(patsubst %, $(TARGETD)/%.cpp, $(TARGETS))
TARGET_OBJ_FILES := $(patsubst $(TARGETD)/%.cpp, $(OBJD)/%.o, $(TARGET_FILES))
HELPER_OBJ_FILES := $(patsubst %, $(OBJD)/%.o, $(HELPERS))

OBJ_FILES := $(TARGET_OBJ_FILES) $(HELPER_OBJ_FILES)
LIB_FILES := $(patsubst %, $(LIBD)/%.a, $(LIBRARIES))
BIN_FILES := $(patsubst %, $(BIND)/%, $(TARGETS))

TARGET_DEPS := $(foreach target, $(TARGETS), $($(target)_DEPS))
LIB_DEPS := $(foreach lib, $(LIBRARIES), $($(lib)_DEPS))
DEPS := $(TARGET_DEPS) $(LIB_DEPS)

# construct 'include' and 'library' flags from already-specified directories
INC := $(patsubst %, -I%, $(INC_DIRS))
LIB := $(patsubst %, -L%, $(LIB_DIRS))

# construct the flags that will be included
FLAGS := $(CPP_FLAGS)
ifdef DEBUG_FLAGS
	FLAGS += $(DEBUG_FLAGS)
endif

# variables to compile the boost.python shared object
uname := $(shell uname -s)
ifeq ($(uname),Linux)
	pythonExec := python
	PYLIB_PATH = $(shell $(pythonExec)-config --exec-prefix)/lib64
else
	pythonExec := python3.8
	PYTHON_SUFFIX = $(shell $(pythonExec) -c "import sys; print(''.join(map(str,sys.version_info[0:2])));")
	PYLIB_PATH = $(shell $(pythonExec) -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'));")
	PYLIB = -L$(PYLIB_PATH) -L$(LIBD) -ldl -framework CoreFoundation -undefined dynamic_lookup -lboost_python$(PYTHON_SUFFIX) -lboost_numpy$(PYTHON_SUFFIX) -ldtermen -lmst -lmstfasst -lmstcondeg -lmstoptim -lmstmagic -lmstfuser $(LDLIBS)
endif
PY_INCLUDES = $(shell $(pythonExec)-config --includes)
PY_SITE_INCLUDE_PARENT = $(shell $(pythonExec)-config --exec-prefix)
PYFLAGS = $(PY_INCLUDES) -I$(PY_SITE_INCLUDE_PARENT)/include -O3 -fPIC -std=c++11 $(INC) $(LIB) $(CONDA_INC)

# phony targets (targets that aren't files should be specified as phony so that they aren't remade each time `make` is run)
.PHONY: all clean libs python setup

# make everything that can be made
# note that the 'all' target should go first so that `make` is equivalent to `make all`
all: $(TARGETS) $(LIBRARIES)

# delete every output file
clean:
	rm -f $(OBJD)/* $(LIBD)/*.a $(BIND)/*

# make every library that can be made
libs: $(patsubst %, $(LIBD)/%.a, $(LIBRARIES))

# make the boost.python shared object
python: $(LIBD)/mstpython.so

# create the output directories that will be needed to make anything
setup: $(OBJD) $(LIBD) $(BIND)

# recipes to ensure each output directory exists when making anything that needs it even if `make setup` hasn't been run
MAKE_DIR = mkdir -p $@

$(OBJD):
	$(MAKE_DIR)
$(LIBD):
	$(MAKE_DIR)
$(BIND):
	$(MAKE_DIR)

$(OBJ_FILES): | $(OBJD)
$(LIB_FILES): | $(LIBD)
$(BIN_FILES): | $(BIND)

# recipes to compile object files (targets, helpers, and boost.python)
COMPILE_OBJ = $(CC) $(FLAGS) -c -o $@ $(INC) $(LIB) $<

$(OBJD)/%.o: $(SRCD)/%.cpp $(INCD)/%.h | $(OBJD)
	$(COMPILE_OBJ)
$(OBJD)/mstpython.o: $(SRCD)/mstpython.cpp | $(OBJD)
	$(CC) $(PYFLAGS) -c -o $@ $<

# secondary expansion is used to specify each target/library's dependencies
#	specifically, "$(LIBD)/%.a: $$($$*_DEPS)" specifies that the library named $lib depends on what's stored in a variable named $lib_DEPS
#	each $lib_DEPS contains the names of the dependencies of $lib
#	to get the files, the functions DEP_HEADER_FILE_MAP and DEP_OBJ_FILE_MAP (defined below) map each name to its corresponding header and object file, respectively
#	this kind of variable substitution requires two passes (since $* is not defined yet on the first pass through)
.SECONDEXPANSION:

# maps from each dependency name to its respective header or object file
DEP_HEADER_FILE_MAP = $(patsubst %, $(INCD)/%.h, $1)
DEP_OBJ_FILE_MAP = $(patsubst %, $(OBJD)/%.o, $1)

# recipes to compile objects
$(OBJD)/%.o: $(TESTD)/%.cpp $$(foreach dep, $$($$*_DEPS), $$(call DEP_HEADER_FILE_MAP, $$(dep))) | $(OBJD)
	$(COMPILE_OBJ)
$(OBJD)/%.o: $(PROGRAMD)/%.cpp $$(foreach dep, $$($$*_DEPS), $$(call DEP_HEADER_FILE_MAP, $$(dep))) | $(OBJD)
	$(COMPILE_OBJ)

# recipes to compile libraries
COMPILE_LIB = ar rs $(LIBD)/$*.a $^

$(LIBD)/%.a: $$(foreach dep, $$($$*_DEPS), $$(call DEP_OBJ_FILE_MAP, $$(dep))) | $(LIBD)
	$(COMPILE_LIB)
$(LIBD)/mstpython.so: $(OBJD)/mstpython.o
	$(CC) $(PYLIB) -Wl,-rpath,$(PYLIB_PATH) -shared -o $@ $<

# recipe to compile targets
COMPILE_BIN = $(CC) $(FLAGS) $(LDLIBS) -o $@ $(INC) $(LIB) $^

$(BIND)/%: $(OBJD)/%.o $$(foreach dep, $$($$*_DEPS), $$(call DEP_OBJ_FILE_MAP, $$(dep))) | $(BIND)
	$(COMPILE_BIN)

# recipes that map names to literal files
$(DEPS): $(OBJD)/$$@.o
$(LIBRARIES): $(LIBD)/$$@.a
$(TARGETS): $(BIND)/$$@
