CXX      ?= g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -g -O2
LIBS      = -lm -lz -lpthread

PREFIX    = $(DESTDIR)/usr/local
BINDIR    = $(PREFIX)/bin

ofiles    = src/main.cpp.o
hfiles    = $(wildcard src/*.h)

.PHONY: all clean install

all: disty

install: disty
	install disty $(BINDIR)/disty

disty: $(ofiles)
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(ofiles) -o $@ -L. $(LIBS)

src/%.cpp.o: src/%.cpp $(hfiles)
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $< -o $@

clean:
	rm -f src/*.o
	rm -f disty

