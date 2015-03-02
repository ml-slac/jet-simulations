# --------------------------------------------- #
# Makefile for Recluster code                        #
# Pascal Nef, March 6th 2014                    #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS =   -O2 -Wall 

.PHONY: clean debug all

all: Recluster

Recluster:  Recluster.so ReclusterTools.so RC.so ReclusterAnalysis.so
	$(CXX) Recluster.so ReclusterTools.so ReclusterAnalysis.so RC.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs` -lEG -lEGPythia8 \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lSubjetJVF  -lVertexJets \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

Recluster.so: Recluster.C    ReclusterTools.so ReclusterAnalysis.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

RC.so: RC.cc RC.h RCConfiguration.h 
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

ReclusterTools.so : ReclusterTools.cc ReclusterTools.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags`

ReclusterAnalysis.so : ReclusterAnalysis.cc ReclusterAnalysis.h 
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs` 

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~
