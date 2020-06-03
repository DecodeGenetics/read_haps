# read_haps
Program to estimate human contamination in Illumina WGS files

## Installation
### Get binary
Fetch a static binary (Linux 64bit) using

```sh
wget https://github.com/DecodeGenetics/read_haps/releases/download/v0.1.1/read_haps
chmod a+x read_haps
```

### Compilation
You can also compile from source if you have the following:
 * C++ compiler with C++11 support
 * Boost >=1.52
 * zlib
 * bz2

Download or clone this repository and type `make` to build.

### Testing
You can test your `read_haps` compilation using

```sh
make test
```

## Usage
 read_haps [OPTIONS] "BAMFILE" "RELIABLE_SNP_FILE" "VCF_FILE"

SYNOPSIS
    read_haps [OPTIONS] "BAMFILE" "RELIABLE_SNP_FILE" "VCF_FILE"

DESCRIPTION
    Determines evidence for DNA contamination in a bam file
    in a diploid individual (human) using evidence of three haplotypes
    between sets read-pair adjacent SNPs using a set of SNPs that
    have been shown to give reliable genotypes on a population scale and a VCF file

    BAMFILE - The (sorted) bam file being consider, along with bai
    RELIABLE_SNP_FILE - List of chromosome and position with positions ordered (Provided file high_quality_markers_deCODE_2015.txt), only biallelic SNPs are used.
    VCF_FILE - A file containing genotype calls

    -h, --help
          Display the help message.
    -q, --qual INT
          Minimum bp quality Default: 30.
    -pl, --phred INT
          Minimum phred likelihood Default: 40.
    -c, --cigarlen INT
          Max cigar len Default: 1.
    -w, --window INT
          Read pair window Default: 1000.
    -mq, --mapq INT
          Min mapping quality Default: 30.
    -np, --npairs INT
          Min number of SNP pairs Default: 10000.
    -e, --error DOUBLE
          Max double error pair fraction Default: 0.002.
    -fa, --fa FA
          Fasta file Default: genome.fa.
    -v, --verbose
          verbose
    -p, --pairs
          Output SNP pairs
    -i, --indels
          Use indels

DETAILS
  read_haps takes as input a bam file, a set of reliable markers and a VCF file. It starts by finding all heterozygous reliable markers in the VCF file.  Only biallelic
  markers are considered and  marker is assumed to have alleles 0 and 1. Two heterozygous markers can have phase 00-11 (parity) or 01-10 (non-parity), i.e. for parity 
  phasing the individual has two haplotypes one with allele 0 on both markers and the other with allele 1 on both markers.  Observing both parity and non-parity 
  suggests the presence of three haplotypes between the two markers, but can also be explained by genotyping error and sequencing error.  To limit genotyping error
  the marker set is limited to markers that have been extensively verified to not have genotyping errors at deCODE genetics and markers that the genotyping software 
  has made high confidence calls.  To limit sequencing errors only reads with high mapQ and bases with high bp quality.  Three haplotypes should generally not be 
  observed for a sample from a single diploid individual (human), but may occur in regions of structural variation.  Three haplotypes will be commonly observed in 
  a sample that contains the mixture of the DNA of two individuals. Pairs of heterozygous markers overlapped by two read pairs in parity and two reads pairs in 
  non-parity are considered evidence for three haplotypes at the marker pair.  Samples with multiple such pairs are considered contaminated.


# Example output

SNP_PAIRS ERROR_PAIRS DOUBLE_ERROR_PAIR_COUNT DOUBLE_ERROR_FRACTION REL_ERROR_FRACTION NONSENSE_FRACTION PASS_FAIL REASON

905879 955 18 1.98702e-05 0.000490452 0.00196273 PASS -

SNP_PAIRS - Number of marker pairs that are connnected by at least one haplotype

ERROR_PAIRS - Number of markers pairs where there is evidence of both parity and non-parity haplotypes

DOUBLE_ERROR_PAIR_COUNT - Number of marker pairs where there are two read pairs evidencing both parity and non-parity

DOUBLE_ERROR_FRACTION - Fraction of marker pairs with two read pairs evidencing parity and non-parity

REL_ERROR_FRACTION - Average fraction of reads from the smaller group (parity, non-parity)

NONSENSE_FRACTION - Fraction of read pairs reporting neither of the two alleles at the markers considered.

PASS_FAIL - PASS or FAIL

REASON - "-" for PASS, FAIL can be READ_PAIR_COUNT, which will occur if a low coverage sample is used or CONTAMINATION

to further investigate the contamination the -pu option can be used, reporting the location of all read pairs showing evidence of three haplotypes.
