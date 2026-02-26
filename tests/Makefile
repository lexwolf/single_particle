CC = g++
CXXFLAGS += -Iinclude
CFLAGS = -Wall -I../src/headers -L../src/lib
LIBS = -lgsl -lgslcblas -lm -larmadillo

# List of test programs
TARGETS = time_behavior threshold steady_state

# Define how to build each target
all: $(TARGETS)

time_behavior: time_behavior.cxx
	$(CC) $(CFLAGS) $(CXXFLAGS) time_behavior.cxx -o time_behavior $(LIBS)

threshold: threshold.cxx
	$(CC) $(CFLAGS) $(CXXFLAGS) threshold.cxx -o threshold $(LIBS)

steady_state: steady_state.cxx
	$(CC) $(CFLAGS) $(CXXFLAGS) steady_state.cxx -o steady_state $(LIBS)

clean:
	rm -f $(TARGETS)
