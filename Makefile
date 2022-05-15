CC = g++
SRCDIR = src
SRCEXT = cpp 
OBJEXT = o 
BUILDDIR = build
LDFLAGS = -L./zlib -lz -lpthread
CFLAGS = -c -O3 -std=c++2a -pthread -I ./seqan2/include -I./seqan2/include -I./zlib -O3 -DNDEBUG -DSEQAN_HAS_ZLIB=1 -DSEQAN_DISABLE_VERSION_CHECK=YES -W -Wall -pedantic
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
	rm -f $(BUILDDIR)/*.o $(SRCDIR)/*.o bin/FastRemap FastRemap 

