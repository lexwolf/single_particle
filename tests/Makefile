CXX = g++
CPPFLAGS = -I../include -I/usr/include/ -I/usr/include/eigen3 -DCUP_BACKEND_QUASI_STATIC
CXXFLAGS = -Wall
LDFLAGS = -L/usr/local/lib
LIBS = -lgsl -lgslcblas -lm -larmadillo

# List of test programs
TARGETS = time_behavior threshold steady_state time_behavior_Q

# Define how to build each target
all: $(TARGETS)

time_behavior: time_behavior.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) time_behavior.cxx -o time_behavior $(LIBS)

threshold: threshold.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) threshold.cxx -o threshold $(LIBS)

steady_state: steady_state.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) steady_state.cxx -o steady_state $(LIBS)

time_behavior_Q: time_behavior_Q.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) time_behavior_Q.cxx -o time_behavior_Q $(LIBS)

clean:
	rm -f $(TARGETS)
