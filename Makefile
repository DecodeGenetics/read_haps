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
CXXFLAGS+=-std=c++0x -Wall

all: read_haps

read_haps: read_haps.cc htslib/libhts.a seqan/.done
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

seqan/.done:
	tar xf seqan.tgz
	touch $@

htslib/libhts.a:
	tar xf htslib.tgz
	$(MAKE) -C htslib libhts.a

test: read_haps test/small_expected
	./read_haps -fa test/small_genome.fa \
	  test/HG002_chr1_10_20MB.bam test/small_hq_markers \
    test/HG002_chr1_10_20MB.vcf.gz > test/small_test
	cmp --silent test/small_test test/small_expected && echo "Test was succesful!"\
	  || echo "Test failed!" && diff test/small_test test/small_expected

.PHONY: all test
