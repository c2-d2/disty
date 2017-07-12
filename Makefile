CXX      ?= g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -g -O2
LIBS      = -lm -lz -lpthread

PREFIX    = $(DESTDIR)/usr/local
BINDIR    = $(PREFIX)/bin

ofiles    = distmat/main.cpp.o
hfiles    = $(wildcard distmat/*.h)

.PHONY: all clean install

all: distmat/distmat

install: distmat/distmat
	install distmap/distmat $(BINDIR)/distmat

distmat/distmat: $(ofiles)
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(ofiles) -o $@ -L. $(LIBS)

distmat/%.cpp.o: distmat/%.cpp $(hfiles)
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $< -o $@

clean:
	rm -f distmat/*.o
	rm -f distmat/distmat

