# --------------------------------------------- #
# Makefile for Mutual Information code          #
# Benjamin Nachman, Aug 21, 2014                # 
# Originally by Pascal Nef, March 6th 2014      #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXX           ?= g++
CXXFLAGS =   -O2 -Wall 

.PHONY: clean debug all

all: MI

MI:  MI.so MITools.so MIAnalysis.so
	$(CXX) MI.so MITools.so MIAnalysis.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs` -lEG -lEGPythia8 \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lSubjetJVF  -lVertexJets \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

MI.so: MI.C    MITools.so MIAnalysis.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

MITools.so : MITools.cc MITools.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags`

MIAnalysis.so : MIAnalysis.cc MIAnalysis.h 
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
