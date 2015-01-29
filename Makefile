CXX = g++

ROOT_CFLAGS = `$(ROOTSYS)/bin/root-config --cflags --glibs`
ROOT_LIBS = `$(ROOTSYS)/bin/root-config --libs --glibs` -lTree

CURRENT_DIR=$(shell pwd)

SOURCES=$(CURRENT_DIR)/TData.cpp $(CURRENT_DIR)/DataManager.cpp  $(CURRENT_DIR)/Globals.cpp $(CURRENT_DIR)/Analyzer.cpp
OBJECTS=$(SOURCES:.cpp=.o)
HEADERS=$(SOURCES:.cpp=.h)



EXECUTABLE=Minimize

CXXFLAGS= -g -O2 -Wall $(ROOT_CFLAGS) -fPIC -I$(ROOTSYS)/include -I$(CURRENT_DIR)
LDFLAGS=$(ROOT_LIBS)


all: $(EXECUTABLE) libDataClasses.so 

$(EXECUTABLE): main.cpp $(OBJECTS) $(CURRENT_DIR)/libDataClasses.so
	$(CXX)  $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

#TInput.o: TInput.cpp $(HEADERS)
#	$(CXX) -c -o $@ $< $(CXXFLAGS) $(ROOT_FLAGS)

DataClassesDict.cpp: $(HEADERS) LinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c $^

$(CURRENT_DIR)/libDataClasses.so: $(OBJECTS) DataClassesDict.o 
	$(CXX) -Wl,--no-as-needed $(ROOT_LIBS) -shared -fPIC -o $@ $(shell root-config --ldflags) -I$(ROOTSYS)/include $^

libDataClasses.rootmap: libDataClasses.so LinkDef.h
	rlibmap -f -o $@ -l libDataClasses.so -c LinkDef.h



print:
	@echo "Objects $(OBJECTS)"
	@echo "Headers $(HEADERS)"

.PHONY: clean

clean:
	rm *.o 
	rm *~
