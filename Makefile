CXX = g++
CXXFLAGS = -std=c++17 -I inc -O3
LIBS = -lgsl -lgslcblas -lm -fopenmp
SRC = src/base.cpp

spectrum: src/base.cpp src/spectrum.cpp
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o bin/$@

tabulate: src/base.cpp src/tabulate.cpp
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o bin/$@

distribution: src/base.cpp src/distribution.cpp
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o bin/$@

clean:
	rm -f bin/*


