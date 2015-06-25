
CXX=g++
CFLAGS=-c -O2 -Wall
LDFLAGS=
SOURCES = ordering-annealing.cpp ordering-io.cpp ordering-tools.cpp
OBJECTS = $(SOURCES:.cpp=.o)

EXECUTABLES = ordering-kernel ordering-exhaustive

all: $(EXECUTABLES)

$(EXECUTABLES): $(OBJECTS)

.cpp.o:
	$(CXX) $(CFLAGS) $< -o $@

.cpp:
	$(CXX) $(LDFLAGS) $(OBJECTS) $< -o $@ #$(addprefix ../../bin/, $@)

clean:
	rm -rf $(OBJECTS) $(EXECUTABLES)
