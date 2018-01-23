CXX = g++
CFLAGS = -Wall
EXEC = temp-plotting.o
EXEC_SRCS = temp-plotting.cc
SRCS := $(wildcard *.cc)
SRCS := $(filter-out $(EXEC_SRCS), $(SRCS))
DEP = $(SRCS:.cc=.o)
GNU = -lboost_iostreams -lboost_system -lboost_filesystem
all: $(DEP) $(EXEC)
clean:
	rm -rf *.o
%.o: %.cc
	$(CXX) $(CFLAGS) $(GNU) -c $< -o $@
$(EXEC): $(DEP) $(EXEC_SRCS)
	$(CXX) $(CFLAGS) $(DEP) -o $@ $(EXEC_SRCS)
