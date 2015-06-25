
CXX=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES = greedy.cpp
OBJECTS = $(SOURCES:.cpp=.o)

EXECUTABLES = main_boxes_greedy_nbox_ic_wdiag_kernel

all: $(EXECUTABLES)

$(EXECUTABLES): $(OBJECTS)

.cpp.o:
	$(CXX) $(CFLAGS) $< -o $@

.cpp:
	$(CXX) $(LDFLAGS) $(OBJECTS) $< -o $@

clean:
	rm -rf $(OBJECTS) $(EXECUTABLES)