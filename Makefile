PREFIX = .
CC = g++
INCLUDES= -I./seqan2/include
EXECUTABLE = ./bin/FastRemap 
SRCDIR = src
SRCEXT = cpp 
OBJEXT = o 
BUILDDIR = build
LDFLAGS = -lz -lpthread
CFLAGS = -c -O3 -std=c++2a -pthread -O3 -DNDEBUG -DSEQAN_HAS_ZLIB=1 -DSEQAN_DISABLE_VERSION_CHECK=YES -W -Wall -pedantic -L./zlib
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))

OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
OBJECTS = $(SOURCES:.cpp=.o) 

$(EXECUTABLE): $(OBJECTS)
	$(CC) -std=c++11 $(OBJECTS) -o $@ $(LDFLAGS) 

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -f $(BUILDDIR)/*.o $(SRCDIR)/*.o bin/FastRemap FastRemap 

