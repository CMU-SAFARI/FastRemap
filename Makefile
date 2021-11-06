CC = g++
SRCDIR = src
SRCEXT = cpp 
OBJEXT = o 
BUILDDIR = build
LDFLAGS = -L./zlib -lz -lpthread
CFLAGS = -c -O3 -Wall -std=c++2a -pthread -I ./seqan2/include -I./seqan2/include -I./zlib -DSEQAN_HAS_ZLIB=1 
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))


OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
OBJECTS = $(SOURCES:.cpp=.o) 

EXECUTABLE = bin/FastRemap 

$(EXECUTABLE): $(OBJECTS)
	$(CC) -std=c++11 $(OBJECTS) -o $@ $(LDFLAGS) 
	cp bin/FastRemap ./

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(BUILDDIR)/*.o $(SRCDIR)/*.o bin/FastRemap 

