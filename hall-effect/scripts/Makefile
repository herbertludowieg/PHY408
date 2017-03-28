CXX = g++
CFLAGS = -Wall
EXEC_SRC = main.cc
EXEC = main.o
SRCS := $(wildcard *.cc)
SRCS := $(filter-out $(EXEC_SRC), $(SRCS))
DEP = $(SRCS:.cc=.o)
all: $(DEP) $(EXEC)
clean:
	rm -rf *.o
%.o: %.cc
	$(CXX) $(CFLAGS) -c $< -o $@
$(EXEC) : $(DEP) $(SRC)
	$(CXX) $(CFLAGS) $(DEP) -o $@ $(EXEC_SRC)
