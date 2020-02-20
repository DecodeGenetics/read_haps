#define SEQAN_HAS_ZLIB 1

#include <boost/lexical_cast.hpp>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdlib.h>
#include <set>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/store.h>
#include <seqan/vcf_io.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

using namespace seqan;
using namespace std;

#define MAX_READS_STORED 10000

struct readhaps_options{
public:
  CharString bam;
  CharString variants;
  CharString faFile;
  int maxCigarLen;
  bool verbose;
  bool outputSNPpairs;
  int qual;
  int mq;
  int readPairWindow;
  String<Dna5> chrSeq;
  CharString positionsUsedFile;
  float maxDoubleError;
  int nPairsMin;
  int minPL;
  bool seqanAlign;
  int maxSNPsInReadPair; //For vector allocation, the number of SNPs in the interval from begin first to end second, might be increased for mate pair/long read data
  readhaps_options(): faFile("genome.fa"), maxCigarLen(1), verbose(false), outputSNPpairs(false), qual( 30), mq(30), readPairWindow( 1000 ),  maxDoubleError(0.002), nPairsMin(10000), minPL(40),maxSNPsInReadPair(500), seqanAlign(false){};
};

readhaps_options O;


class markerPair{
public:
    int parity;
    int nonparity;
    int errorpair;
};

class markerInfo{
public:
  int pos;
  string a1,a2;
  CharString name;
  markerInfo(){
    pos = 1000000000;
  }
};

class haploSummary{
public:
  int nPairs;
  int nSingleErrorPairs;
  int nDoubleErrorPairs;
  int nBSPairs;
  float errorFractionCount;
  void update( map< int, markerPair> sps ){
    for( auto it = sps.begin(); it != sps.end(); it++ ){
      nPairs++;
      if( (*it).second.errorpair > 0 ) nBSPairs++;
      if( (*it).second.parity > 0 and (*it).second.nonparity > 0 ) nSingleErrorPairs++;
      if( (*it).second.parity > 1 and (*it).second.nonparity > 1 ) nDoubleErrorPairs++;
      errorFractionCount += (1.0*(*it).second.errorpair)/(1.0*((*it).second.parity+(*it).second.nonparity+(*it).second.errorpair));
    }
  }
  void report(){
    cout << "SNP_PAIRS ERROR_PAIRS DOUBLE_ERROR_PAIR_COUNT DOUBLE_ERROR_FRACTION REL_ERROR_FRACTION NONSENSE_FRACTION PASS_FAIL REASON" << endl;
    if( nPairs == 0 ){
      cout << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " FAIL READ_PAIR_COUNT" << endl;
    }else{
      float doubleErrorFrac = (1.0*nDoubleErrorPairs)/(1.0*nPairs);
      float errorFraction = errorFractionCount/(1.0*nPairs);
      float nonsenseFraction = (1.0*nBSPairs)/(1.0*nPairs);
      cout << nPairs << " " << nSingleErrorPairs << " " << nDoubleErrorPairs << " " << doubleErrorFrac << " " << errorFraction << " " << nonsenseFraction << " ";
      if( doubleErrorFrac < O.maxDoubleError ){
	if( nPairs > O.nPairsMin ){
	  cout << "PASS -" << endl;
	}else{
	  cout << "FAIL READ_PAIR_COUNT" << endl;
	}
      }else{
	cout << "FAIL CONTAMINATION" << endl;
      }
    }
  }
  haploSummary(): nPairs(0), nSingleErrorPairs(0), nDoubleErrorPairs(0), nBSPairs(0), errorFractionCount(0.0){};
};




void outputSNPpairs( vector< markerInfo>& mI, int orig, map< int, markerPair> sps ){
    if( O.verbose ) cout << "outputSNPpairs " << orig << " " << mI[orig].name << endl;
    for( auto it = sps.begin(); it != sps.end(); it++ ){
        cout << mI[orig].name << " " << mI[(*it).first].name  << " " << (*it).second.parity << " " << (*it).second.nonparity << " " << (*it).second.errorpair << endl;
    }
}

int addReadToVectors( vector< markerInfo>& mI, map< int, int>& posHash, BamAlignmentRecord& bar, int& vectorI, vector< int>& snpIDs, vector< bool>& refCarrier, vector< bool>& altCarrier ){
  int snpBeg;
  if( posHash.count( bar.beginPos/100 ) != 0 ){
    snpBeg = posHash[bar.beginPos/100];
  }else if( posHash.count( bar.beginPos/100 +1 ) != 0 ){
    snpBeg = posHash[bar.beginPos/100+1];
  }else if( posHash.count( bar.beginPos/100 +2 ) != 0 ){
    snpBeg = posHash[bar.beginPos/100+2];
  }else{  
    return 1;
  }
  while( snpBeg < (int) mI.size() and mI[snpBeg].pos < bar.beginPos )
    snpBeg++;
  if( snpBeg == (int) mI.size() )
    return 1;
  int alnLen = (int) length( bar.seq );
  if( O.seqanAlign ){
    alnLen = (int) getAlignmentLengthInRef( bar );
  }
  if( mI[snpBeg].pos > bar.beginPos + alnLen )
    return 1;    
  
  Align<Dna5String> align;
  if( O.seqanAlign )
    bamRecordToAlignment(align, O.chrSeq, bar);
  while( snpBeg < (int) mI.size() and mI[snpBeg].pos < bar.beginPos + alnLen ){
    int barPos = 0;
    if( O.seqanAlign ){
      toSourcePosition( row( align, 1 ), toViewPosition( row( align, 0), mI[snpBeg].pos)  );
    }else{
      barPos = mI[snpBeg].pos - bar.beginPos;
    }
    if( O.verbose ) cout << "barPos " << barPos << endl;
    if( bar.qual[barPos] - 33 >= O.qual ){
      if( vectorI < (int) snpIDs.size() ){
	snpIDs[vectorI] = snpBeg;
	refCarrier[vectorI] = mI[snpBeg].a1[0] == bar.seq[barPos];
	altCarrier[vectorI] = mI[snpBeg].a2[0] == bar.seq[barPos];
	vectorI++;
      }
    }
    snpBeg++;
  }
  return 0;
}


void addReadPairToSNPpairs( vector< markerInfo>& mI, map< int, int>& posHash, map< int, map<int, markerPair> >& snpPairs, std::list<int>& snpList,  BamAlignmentRecord& bar1, BamAlignmentRecord& bar2 ){
  vector< int> snpIDs;
  snpIDs.resize( O.maxSNPsInReadPair );
  vector< bool> refCarrier;
  refCarrier.resize( O.maxSNPsInReadPair );
  vector< bool> altCarrier;
  altCarrier.resize( O.maxSNPsInReadPair );
  int nSNPsInReadPair = 0;
  if( (int) length( bar1.cigar ) <= O.maxCigarLen and (int) length( bar2.cigar ) <= O.maxCigarLen and bar1.mapQ >= O.mq and bar2.mapQ >= O.mq and (not hasFlagUnmapped( bar1 )) and 
      (not hasFlagUnmapped( bar2 )) and (not hasFlagQCNoPass( bar1 )) and (not hasFlagQCNoPass( bar2 )) and (not hasFlagDuplicate( bar1 )) and (not hasFlagDuplicate( bar2 )) and  
      hasFlagAllProper(bar1) and  hasFlagAllProper(bar2) and (not hasFlagSecondary(bar1)) and (not hasFlagSecondary(bar2)) ){
    addReadToVectors( mI, posHash, bar1, nSNPsInReadPair, snpIDs, refCarrier, altCarrier );
    addReadToVectors( mI, posHash, bar2, nSNPsInReadPair, snpIDs, refCarrier, altCarrier );
  }
  if( O.verbose ) cout << "Adding read " << nSNPsInReadPair << endl;
  for( int i = 0; (i < nSNPsInReadPair) and (i < O.maxSNPsInReadPair); i++ ){
    for( int j = i+1; (j < nSNPsInReadPair) and (j < O.maxSNPsInReadPair); j++ ){
      int first = snpIDs[i]; int second = snpIDs[j];
      if( snpIDs[i] >= snpIDs[j] ){  // Markers on the second read could also be on the first read, we don't want to double count those
	first = snpIDs[j]; second = snpIDs[i];
      }else{
	if( snpPairs.count( first ) == 0 ){
	  snpPairs[first] ={};
	  snpList.emplace_front( first );
	}
	if( snpPairs[first].count( second ) == 0 ){
	  snpPairs[first][second].parity = 0;
	  snpPairs[first][second].nonparity = 0;
	  snpPairs[first][second].errorpair = 0;
	}
        
	if( (refCarrier[i] and refCarrier[j]) || (altCarrier[i] and altCarrier[j]) )
	  snpPairs[first][second].parity++;
	else if( (refCarrier[i] and altCarrier[j]) || (altCarrier[i] and refCarrier[j]))
	  snpPairs[first][second].nonparity++;
	else
	  snpPairs[first][second].errorpair++;
      }
    }
  }
}

int processBamRegion( vector< markerInfo>& mI,  map< int, int>& posHash, CharString& cChrom, HtsFile& bamS, haploSummary& hs ){
    
    
    //    BamHeader header;
    //readHeader(header, bamS );

  if( O.verbose ) cout << "reading Bam region " << cChrom << " " << 1 << " " << length(O.chrSeq ) << endl;
    
    /*   unsigned int rID = 0;
    if (!getIdByName( rID, contigNamesCache( context( bamS)   ), cChrom )){
      cerr << "ERROR: Reference sequence named " << cChrom << " not known.\n";
      return 1;
    }
    
    if( O.verbose ) cerr << "GotID " << cChrom << " " << rID <<  " " << length(O.chrSeq ) << endl; */
    // Jump the BGZF stream to this position.
    if (!setRegion(bamS, toCString(cChrom), 1, length(O.chrSeq)))
    {
        cerr << "ERROR: Could not jump to " << cChrom << "\n";
        return 1;
    }
    
    map< CharString, BamAlignmentRecord > bars;
    std::list< CharString> barList;
    map< int, map< int, markerPair> > snpPairs;
    std::list< int> snpList;

    BamAlignmentRecord record;

    unsigned lrID = 0;
    bool rIDunknown = true;
    bool breakWhile = false;

    while (readRegion(record, bamS) and (not breakWhile))
    {
      if( rIDunknown ){  // Silly hack 
	rIDunknown = false;
	lrID = record.rID;
      }
      if( bars.count( record.qName ) == 0 ){
	if( not hasFlagSecondary( record ) ){
	  bars[record.qName] = record;
	  barList.emplace_front( record.qName );
	  if( O.verbose ) cout << "added " << record.qName << endl;
	}
      }else{
	if( record.beginPos < bars[record.qName].beginPos )
	  cout << "WTF " << record.qName << " " << record.beginPos << " " << bars[record.qName].beginPos << endl;
	if( (record.beginPos <= (bars[record.qName].beginPos + O.readPairWindow)) and 
	    (hasFlagRC(record) xor hasFlagRC(bars[record.qName])) )
	  addReadPairToSNPpairs( mI, posHash, snpPairs, snpList, bars[record.qName], record );  // Changed order here, old order does not make sense
	barList.remove(record.qName);
	bars.erase( record.qName);
	if( O.verbose ) cout << "removed " << record.qName << endl;
      }
      while( (not barList.empty()) and ((bars[barList.back()].beginPos < record.beginPos - O.readPairWindow) or (bars.size() > MAX_READS_STORED)) ){
	if( O.verbose ) cout << barList.empty() << " barList empty" << endl;
	if( O.verbose ) cout << barList.back() << endl;
	bars.erase(barList.back());
	barList.pop_back();
      }
      
      while( (not snpList.empty()) and mI[snpList.back()].pos < record.beginPos - O.readPairWindow ){
	if( O.outputSNPpairs) outputSNPpairs( mI, snpList.back(), snpPairs[snpList.back()] );
	hs.update( snpPairs[snpList.back()] );
	snpPairs.erase( snpList.back() );
	snpList.pop_back();
      }
    }
    while( not snpList.empty() ){
      if( O.outputSNPpairs ) outputSNPpairs( mI, snpList.back(), snpPairs[snpList.back()] );
      snpPairs.erase( snpList.back() );
      snpList.pop_back();
    }
    return 0;
}


int parseReadHapArguments( int argc, char const ** argv ){
    
  ArgumentParser parser("read_haps");
  addUsageLine(parser, "[\\fIOPTIONS\\fP]  \"\\fIBAMFILE\\fP\" \"\\fIRELIABLE_SNP_FILE\\fP\"  \"\\fIVCF_FILE\\fP\" ");
  addDescription(parser, "Determines evidence for DNA contamination in a bam file ");
  addDescription(parser, "in a diploid individual (human) using evidence of three haplotypes");
  addDescription(parser, "between sets read-pair adjacent SNPs using a set of SNPs that ");
  addDescription(parser, "have been shown to give reliable genotypes on a population scale and a VCF file");
  addDescription(parser, "BAMFILE - The bam file being consider, along with bai");
  addDescription(parser, "RELIABLE_SNP_FILE - List of chromosome and position with positions ordered");
  addDescription(parser, "VCF_FILE - A file containing genotype calls");
  
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "bam"));
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "pos_used"));
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "variants"));
  addOption(parser, ArgParseOption("q", "qual", "Minimum bp quality ", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("pl", "phred", "Minimum phred likelihood ", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("c", "cigarlen", "Max cigar len ", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("w", "window", "Read pair window ", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("mq", "mapq", "Min mapping quality ", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("np", "npairs", "Min number of SNP pairs ", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("e", "error", "Max double error pair fraction ", ArgParseArgument::DOUBLE, "DOUBLE"));
  addOption(parser, ArgParseOption("fa", "fa", "Fasta file ", ArgParseArgument::STRING, "FA"));
  addOption(parser, ArgParseOption("v", "verbose", "verbose"));
  addOption(parser, ArgParseOption("p", "pairs", "Output SNP pairs"));
  setDefaultValue( parser, "qual", O.qual );
  setDefaultValue( parser, "phred", O.minPL );
  setDefaultValue( parser, "mq", O.mq );
  setDefaultValue( parser, "np", O.nPairsMin );
  setDefaultValue( parser, "e", O.maxDoubleError );
  setDefaultValue( parser, "cigarlen", O.maxCigarLen );
  setDefaultValue( parser, "window", O.readPairWindow );
  setDefaultValue( parser, "fa", O.faFile );
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // Only extract  options if the program will continue after parseReadHapArguments()
  if (res != ArgumentParser::PARSE_OK){
    return res;
  }
    
  //  cout << "Get arg" << endl;
  getArgumentValue( O.bam, parser, 0);
  getArgumentValue( O.positionsUsedFile, parser, 1);
  getArgumentValue( O.variants, parser, 2);

  //cout << "Got arg" << endl;
  O.verbose = isSet(parser, "verbose");
  O.outputSNPpairs = isSet( parser, "pairs" );
  if( isSet( parser, "qual" ) )
    getOptionValue( O.qual, parser, "qual");
  if( isSet( parser, "phred" ) )
    getOptionValue( O.minPL, parser, "phred");
  if( isSet( parser, "mq" ) )
    getOptionValue( O.mq, parser, "mq");
  if( isSet( parser, "np" ) )
    getOptionValue( O.nPairsMin, parser, "np");
  if( isSet( parser, "e" ) )
    getOptionValue( O.maxDoubleError, parser, "e");
  if( isSet( parser, "fa" ) )
    getOptionValue( O.faFile, parser, "fa");
  if( isSet( parser, "cigarlen" ) ){
    getOptionValue( O.maxCigarLen, parser, "cigarlen");
    if( O.maxCigarLen != 1 ) O.seqanAlign = true;
  }
  if( isSet( parser, "window" ) )
    getOptionValue( O.readPairWindow, parser, "window");
  
  return 0;
}


int main( int argc, char const ** argv )
{
  
  if( parseReadHapArguments( argc, argv ) != 0 )
    return 1;

  if( O.verbose ) cerr << "parsed args " << endl;
  vector< int> posPointers;
  map< CharString, BamAlignmentRecord> reads;
  list< CharString> readList;

  //  char buffer[1024];

  //  int nMarkers = 0;
  //string name;
  map<CharString, vector< markerInfo>> mI;
  map<CharString, map< int, int>> posHash;
  
  map<CharString, map< int, int>> posUsed;
  ifstream f;
  int markerP;
  string sChrom;
  f.open( toCString( O.positionsUsedFile ) );  // Need to catch error
  f >> sChrom;
  f >> markerP;
  
  while( f ){
    CharString chrom(sChrom.c_str());
    posUsed[chrom][markerP] = 1;
    f >> sChrom;
    f >> markerP;
  }
 
  if( O.verbose ) cerr << "read pos used file " << endl;
  // Open the input VCF file and prepare output vcf stream.
  //  VcfFileIn vcfIn(toCString(O.variants));
  //VcfHeader header;
  //readHeader(header, vcfIn);

  FaiIndex faiIndex;
  if( ! open( faiIndex,  toCString( O.faFile ) )  )
    std::cerr << "ERROR: Could not load FAI index path/to/index.fai\n";
  unsigned idx = 0;

  seqan::Tabix vcfIndex;
  seqan::open(vcfIndex, toCString(O.variants));

  //  seqan::setRegion(index, "chrX:A-B"); // "chrX" and "chrX:A" also supported
  seqan::VcfRecord record;

  
  if( O.verbose ) cerr << "read VCF header " << endl;
  map< CharString, int> mUsed;
  for( auto i = posUsed.begin(); i != posUsed.end(); i++ ){
    CharString cChrom = (*i).first;
    mI[cChrom].resize( 100 );
    mUsed[cChrom] = 0;
    unsigned lrID = 0;
    bool rIDunknown = true;
    bool breakWhile = false;

    setRegion( vcfIndex, toCString(cChrom));
    while ( readRecord(record, vcfIndex ) and (not breakWhile) ){
      if( rIDunknown ){  // Silly hack 
	rIDunknown = false;
	lrID = record.rID;
      }
      if( record.rID != lrID )
	breakWhile = true;
      StringSet< CharString> altSet;
      strSplit(altSet,record.alt,EqualsChar<','>());
      if( length( record.ref ) == 1 and length( altSet[0] ) == 1 ){
	bool heterozygote = false;
	// This is not robust, the GT field needs to be first, should ask seqan for the gt field
	unsigned int gtID = 0, plID;
	StringSet< CharString> formatIs;
	strSplit(formatIs,record.format,EqualsChar<':'>());
	CharString gt2 = concat( record.genotypeInfos );
	StringSet< CharString> gtIs;
	strSplit(gtIs,gt2,EqualsChar<':'>());
	getIdByName( gtID, formatIs, "GT" );
	getIdByName( plID, formatIs, "PL" );
	CharString cGT = getValueById( record.genotypeInfos, gtID );
	if( gtIs[gtID][0] == '0' and gtIs[gtID][2] == '1' )
	  heterozygote = true;
	StringSet< CharString> PLset;
	strSplit(PLset,gtIs[plID],EqualsChar<','>());
	bool highPL = false;
	try{  // The PL might be set as "-"
	  if( atoi( toCString( PLset[0] ) ) >= O.minPL and atoi( toCString( PLset[2] ) ) >= O.minPL ){
	    highPL = true;
	  }
	}catch(std::exception const & ex){}

	if( ( highPL and heterozygote and posUsed[cChrom].count( record.beginPos + 1 ) != 0) ){
	  mI[cChrom][mUsed[cChrom]].pos = record.beginPos;
	  mI[cChrom][mUsed[cChrom]].name = record.id;
	  if( mI[cChrom][mUsed[cChrom]].name  == "." )  // Use marker name if given, else use chr:pos
	    mI[cChrom][mUsed[cChrom]].name  = string( toCString( cChrom ))+string(":")+to_string( (long long unsigned int) (record.beginPos + 1) );
	  mI[cChrom][mUsed[cChrom]].a1 = string( toCString( record.ref ));
	  mI[cChrom][mUsed[cChrom]].a2 = string( toCString( altSet[0] ));
	  if( posHash[cChrom].count( mI[cChrom][mUsed[cChrom]].pos/100 ) == 0 ) posHash[cChrom][mI[cChrom][mUsed[cChrom]].pos/100] = mUsed[cChrom];
	  if( O.verbose ) cout << mI[cChrom][mUsed[cChrom]].pos << " " << mI[cChrom][mUsed[cChrom]].name << " " << mI[cChrom][mUsed[cChrom]].a1 << " " << mI[cChrom][mUsed[cChrom]].a2 << endl;
	  mUsed[cChrom]++;
	  if( mUsed[cChrom] >= (int) mI[cChrom].size() ) mI[cChrom].resize( 2*mI[cChrom].size() );
	}
      }
    }
    /*    for( auto h = posHash[cChrom].begin(); h != posHash[cChrom].end(); h++ )
      cout << "PH " << cChrom << " " << (*h).first << " " << (*h).second << endl;
      return 1;*/
  }
  for( auto i = mI.begin(); i != mI.end(); i++ ){
    mI[(*i).first].resize( mUsed[(*i).first] );
    if( O.verbose ) cout << (*i).first << " " << mI[(*i).first].size() << endl;
  }

  //  BamIndex<Bai> baiI;
  HtsFile bamS( toCString( O.bam), "r");
  if (!loadIndex(bamS)){
      // Build it if we cannot find it
      buildIndex(bamS);
      loadIndex(bamS);
  }
  //  if( O.verbose ) cerr << O.bam << " bam file " << endl;
  //if (initializeBam(toCString(O.bam), baiI, bamS) != 0){
  // cerr << "ERROR: Could not open " << O.bam << " for reading.\n";
  // return 1;
  // }

  haploSummary hs;
  for( auto i = mI.begin(); i != mI.end(); i++ ){
    CharString tChrom = (*i).first;
    if (!getIdByName(idx, faiIndex, tChrom ))
      std::cerr << "ERROR: FAI index has no entry for " << (*i).first << std::endl;
    readSequence(O.chrSeq, faiIndex, idx);
    processBamRegion( mI[(*i).first], posHash[(*i).first], tChrom, bamS, hs );
  }
  hs.report();
  return 0;
}
