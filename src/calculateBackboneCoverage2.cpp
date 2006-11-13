/*******************************************************************************
 * $Id: calculateBackboneCoverage.cpp,v 1.5 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/IntervalList.h"
#include "libMems/Islands.h"

using namespace std;
using namespace genome;
using namespace mems;

void print_usage( const char* pname ){
	cerr << "Usage: " << pname << " <XMFA alignment> <min bb sequence length> <max bb gap size> \n";
}


int main( int argc, const char* argv[] ){

try{
	if( argc <= 0 ){
		print_usage( "extractBackbone" );
		return -1;
	}
	if( argc < 4 ){
		print_usage( argv[0] );
		return -1;
	}
	
	string alignment_fname = argv[1];
	int64 min_bb_length = atol( argv[2] );
	int64 max_gap_length = atol( argv[3] );
	vector< string > sequence_fname;
	vector< gnSequence* > source_seqs;
	
	ifstream alignment_in;
	alignment_in.open( alignment_fname.c_str() );
	if( !alignment_in.is_open() ){
		cerr << "Error opening " << alignment_fname << endl;
		return -1;
	}
	
	cout << "Loading alignment...\n";
	IntervalList aligned_ivs;
	aligned_ivs.ReadStandardAlignment( alignment_in );
	
	MatchList ml;
	ml.seq_filename = aligned_ivs.seq_filename;
	ml.LoadSequences(&cout);
	aligned_ivs.seq_table = ml.seq_table;
	source_seqs = ml.seq_table;

	// add the sequence data to the interval list
	uint seq_count = source_seqs.size();
	cout << "Extracting backbone..." << endl;
	vector< GappedAlignment > backbone_data;
	simpleFindBackbone( aligned_ivs, min_bb_length, max_gap_length, backbone_data );

	IntervalList backbone_ivs;
	backbone_ivs.seq_table = aligned_ivs.seq_table;
	
	cout << "There are " << backbone_data.size() << " backbone segments\n";

	// count up the total length of backbone in each genome
	cout << "Averaging backbone lengths..." << endl;
	vector< gnSeqI > total_bb( seq_count, 0 );
	NumericMatrix< double > overall_identity;
	for( uint bbI = 0; bbI < backbone_data.size(); bbI++ ){
		vector<AbstractMatch*> tmp(1, &backbone_data[ bbI ]);
		backbone_ivs.push_back( Interval(tmp.begin(), tmp.end()) );
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			total_bb[ seqI ] += backbone_data[ bbI ].Length( seqI );
		}
	}
	vector< AbstractMatch* > bbivs;
	for( uint bbI = 0; bbI < backbone_ivs.size(); bbI++ )
		bbivs.push_back( &backbone_ivs[bbI] );
	BackboneIdentityMatrix( bbivs, aligned_ivs.seq_table, overall_identity );
		
	gnSeqI avg_bb = 0;
	for( uint seqI = 0; seqI < aligned_ivs.seq_table.size(); seqI++ ){
		cout << "seq " << seqI << " backbone: " << total_bb[ seqI ] << endl;
		avg_bb += total_bb[ seqI ];
	}
	avg_bb /= aligned_ivs.seq_table.size();
	cout << "Average: " << avg_bb << endl;
	
	// output the identity matrix
	cout << "Identity matrix: " << endl;
	overall_identity.print( cout );
	cout << endl;
	
}catch( gnException& gne ){
	cerr << gne << endl;
}catch( exception& e ){
	cerr << e.what() << endl;
}catch(...){

}
	return 0;
}
