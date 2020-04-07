DIR=$(shell pwd)
HTSLIB=$(DIR)/htslib
SEQAN=$(DIR)/seqan

CXXFLAGS+=-I.
CXXFLAGS+=-isystem $(SEQAN)/include
CXXFLAGS+=-I$(HTSLIB)
CXXFLAGS+=-pthread 
CXXFLAGS+=-Wfatal-errors

LDFLAGS=-g htslib/libhts.a -lz -lbz2 -llzma -lboost_iostreams

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1

#Debug build
#CXXFLAGS+=-O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_HAS_ZLIB=1

# set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x

all: read_haps

read_haps: read_haps.cc htslib/libhts.a seqan/include 
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

seqan/include:
	tar xf seqan.tgz

htslib/libhts.a:
	tar xf htslib.tgz
	$(MAKE) -C htslib libhts.a
