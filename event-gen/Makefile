# Makefile for Jet Images
# Author: Luke de Oliveira (luke.deoliveira@yale.edu)

# --- set directories
BIN          := bin
SRC          := src
INC          := include
LIB          := $(CURDIR)/lib

#  set search path
vpath %.o    $(BIN)
vpath %.cc  $(SRC) 
vpath %.h   $(INC) 

# --- set compiler and flags

CXX          ?= g++
CXXFLAGS     := -Wall -fPIC -I$(INC) -g -std=c++11 
# LIBRARY      := $(LIB)/libdataframe.a

ifeq ($(CXX),clang++)
CXXFLAGS += -stdlib=libc++
endif

# --- HEP flags

FASTJETFLAGS := $(shell fastjet-config --cxxflags --plugins)
FASTJETLDFLAGS =
FASTJETLIBS = $(shell fastjet-config --libs --plugins)

ROOTFLAGS    := $(shell root-config --cflags)
ROOTLDFLAGS   := $(shell root-config --ldflags)
ROOTLIBS     := $(shell root-config --glibs)

PYTHIAFLAGS  := $(shell pythia8-config --cxxflags)
PYTHIALDFLAGS:= $(shell pythia8-config --ldflags)
PYTHIALIBS = $(shell pythia8-config --libs)

#HEPLIBS := -lSubjetJVF -lVertexJets 
HEPLIBS += $(FASTJETLIBS)
HEPLIBS += $(ROOTLIBS)
HEPLIBS += $(PYTHIALIBS)


# --- Add HEP flags to common stuff

CXXFLAGS += $(ROOTFLAGS)
CXXFLAGS += $(FASTJETFLAGS)
CXXFLAGS += $(PYTHIAFLAGS)

LIBS     += $(HEPLIBS)

LDFLAGS  = $(ROOTLDFLAGS) $(PYTHIALDFLAGS) $(FASTJETLDFLAGS)


# -- Boost Flags

BOOSTINCDIR = /usr/local/include
BOOSTLIBLOCATION = /usr/local/lib

CXXFLAGS += -I$(BOOSTINCDIR)
LIBS += -L$(BOOSTLIBLOCATION) -lboost_program_options 

# -- building excecutable

OBJ := MI.o MIAnalysis.o MITools.o

EXECUTABLE   := MI


all: $(EXECUTABLE)
	@echo "jet-images build sucessful."

$(EXECUTABLE): $(OBJ:%=$(BIN)/%)
	@echo "linking $^ --> $@"
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS) 

# $(LIBRARY): $(LIB_OBJ:%=$(BIN)/%)
# 	@mkdir -p $(LIB)
# 	@echo "linking objects to --> $@"
# 	@ar rc $@ $^ && ranlib $@

# --------------------------------------------------

# compile rule
$(BIN)/%.o: %.cc
	@echo compiling $<
	@mkdir -p $(BIN)
	$(CXX) -c $(CXXFLAGS) $< -o $@

# use auto dependency generation
DEP = $(BIN)

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),rmdep)
ifneq ($(MAKECMDGOALS),purge)
include  $(ALLOBJ:%.o=$(DEP)/%.d)
endif
endif
endif

DEPTARGSTR = -MT $(BIN)/$*.o -MT $(DEP)/$*.d
$(DEP)/%.d: %.cxx
	@echo making dependencies for $<
	@mkdir -p $(DEP)
	@$(CXX) -MM -MP $(DEPTARGSTR) $(CXXFLAGS) $< -o $@ 

# clean
.PHONY : clean rmdep
CLEANLIST     = *~ *.o *.o~ *.d core 

clean:
	rm -fr $(CLEANLIST) $(CLEANLIST:%=$(BIN)/%) $(CLEANLIST:%=$(DEP)/%)
	rm -fr $(BIN) 

purge:
	rm -fr $(CLEANLIST) $(CLEANLIST:%=$(BIN)/%) $(CLEANLIST:%=$(DEP)/%)
	rm -fr $(BIN) 
	rm -fr $(EXECUTABLE) $(LIB)

rmdep: 
	rm -f $(DEP)/*.d