BINDIR  = bin
HEADERS = Permutation.h State.h Trail.h TrailExtension.h Units.h MyTree.h Exception.h Tree.h
SOURCES = Permutation.cpp State.cpp Trail.cpp TrailExtension.cpp Units.cpp MyTree.cpp main.cpp
OBJECTS = $(addprefix $(BINDIR)/, $(notdir $(SOURCES:.cpp=.o)))
EXEC    = $(BINDIR)/cpmtrail
STANDALONE = $(BINDIR)/distribution $(BINDIR)/bruteforce

CXXFLAGS += -std=c++14 -O3 -Wall -Wextra -flto -march=native -funroll-loops
LDFLAGS  += -lpthread

.PHONY: all clean doc
all: $(EXEC) $(STANDALONE)

$(BINDIR):
	mkdir -p $(BINDIR)

$(STANDALONE): %: %.o $(BINDIR) 
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

$(EXEC): $(OBJECTS) $(BINDIR)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -o $@

$(BINDIR)/%.o: %.cpp $(BINDIR)
	$(CXX) -c $(CXXFLAGS) $< -o $(BINDIR)/$*.o

clean:
	rm -rf $(BINDIR)

doc:
	doxygen Doxyfile
