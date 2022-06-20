CXX = g++
EXEC = ./bin/FastRemap 
SRCDIR = src
SRCEXT = cpp 
OBJEXT = o 
BUILDDIR = build

INC = ./seqan2/include
INC_PARAMS=$(foreach d, $(INC), -I$d)

LIB = ./zlib
LIB_PARAMS=$(foreach d, $(LIB), -L$d)
LDFLAGS = $(LIB_PARAMS) -lz -lpthread

CFLAGS = $(INC_PARAMS) -c -O3 -std=c++2a -pthread -DNDEBUG -DSEQAN_HAS_ZLIB=1 -DSEQAN_DISABLE_VERSION_CHECK=YES -W -Wall -pedantic
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))

OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
OBJECTS = $(SOURCES:.cpp=.o) 

$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS) 

.cpp.o:
	$(CXX) $(CFLAGS) $< -o $@

clean:
	rm -f $(BUILDDIR)/*.o $(SRCDIR)/*.o $(EXEC) 

