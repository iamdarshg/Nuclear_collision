CXX = g++ src/main.cpp
CXXFLAGS = -O3 -std=c++17 -I./apfel/include -I/usr/include/gsl
LDFLAGS = -L/usr/local/lib  -lgsl -lgslcblas -lm


SOURCES = src/main.cpp src/matrix_elements.hpp src/pdf_evolution.hpp src/parton_distribution.hpp src/parton_distribution_set.hpp src/pdf_evolution.cpp
OBJECTS = $(.cpp=.o)
EXECUTABLE = hg199_sim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f *.o $(EXECUTABLE)
