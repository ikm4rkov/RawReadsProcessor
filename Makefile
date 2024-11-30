CXX = g++
CXXFLAGS = -std=c++11

all: BridgeSplitter EndsProcessor

BridgeSplitter: V1.3.1.cpp
  $(CXX) $(CXXFLAGS) -o BridgeSplitter V1.3.1.cpp $(LDFLAGS)

EndsProcessor: EndsProcessor.cpp
  $(CXX) $(CXXFLAGS) -o EndsProcessor EndsProcessor.cpp $(LDFLAGS)

clean:
  rm -f BridgeSplitter EndsProcessor

.PHONY: all clean
