#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "libGenome/gnFilter.h"
#include "libMems/IntervalList.h"
#include "libMems/MemSubsets.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/Matrix.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/Aligner.h"
#include "libGenome/gnFASSource.h"

using namespace std;
using namespace genome;
using namespace mems;



//
// baaad:  copied functions from ProgressiveAligner.h
//

template< typename PairType >
class LabelSort 
{
public:
	LabelSort( uint seqI ) : ssc( seqI ) {};
	bool operator()( const PairType& pt1, const PairType& pt2 )
	{
		return ssc( pt1.first, pt2.first );
	}
private:
	LabelSort();
	mems::SSC<mems::AbstractMatch> ssc;
};

template<class MatchVector>
void IdentifyBreakpoints( MatchVector& mlist, std::vector<gnSeqI>& breakpoints )
{
	if( mlist.size() == 0 )
		return;
	breakpoints = std::vector<gnSeqI>(1, mlist.size()-1);

	mems::SSC<mems::AbstractMatch> ssc(0);
	std::sort( mlist.begin(), mlist.end(), ssc );
	typedef typename MatchVector::value_type value_type;
	typedef std::pair< value_type, size_t > LabelPairType;
	std::vector< LabelPairType > label_list;
	typename MatchVector::iterator cur = mlist.begin();
	typename MatchVector::iterator end = mlist.end();
	size_t i = 0;
	for( ;cur != end; ++cur )
	{
		label_list.push_back( std::make_pair( *cur, i ) );
		++i;
	}

	uint seq_count = mlist[0]->SeqCount();
	// check for breakpoints in each sequence
	for( uint seqI = 1; seqI < seq_count; seqI++ )
	{
		LabelSort< LabelPairType > ls(seqI); 
		std::sort( label_list.begin(), label_list.end(), ls );

		typename std::vector< LabelPairType >::const_iterator prev = label_list.begin();
		typename std::vector< std::pair< typename MatchVector::value_type, size_t > >::const_iterator iter = label_list.begin();
		typename std::vector< std::pair< typename MatchVector::value_type, size_t > >::const_iterator lab_end = label_list.end();

		bool prev_orient = (*prev).first->Orientation(seqI) == (*prev).first->Orientation(0);
		if( !prev_orient )	// if we start in a different orientation than the ref seq there's a bp here
			breakpoints.push_back(prev->second);

		for( ++iter; iter != lab_end; ++iter )
		{
			bool cur_orient = (*iter).first->Orientation(seqI) == (*iter).first->Orientation(0);
			if( prev_orient == cur_orient &&
				( ( prev_orient && (*prev).second + 1 == (*iter).second) ||
				  ( !prev_orient && (*prev).second - 1 == (*iter).second) 
				)
			  )
			{
				prev_orient = cur_orient;
				++prev;
				continue;	// no breakpoint here
			}

			// always add the last match in a new block (scanning from left to right in seq 0)
			if( prev_orient )
				breakpoints.push_back( prev->second );
			if( !cur_orient )
				breakpoints.push_back( iter->second );

			prev_orient = cur_orient;
			++prev;
		}
		if( prev_orient )
			breakpoints.push_back( prev->second );
	}
	std::sort( breakpoints.begin(), breakpoints.end() );
	std::vector<gnSeqI>::iterator uni = std::unique( breakpoints.begin(), breakpoints.end() );
	breakpoints.erase( uni, breakpoints.end() );
}


template< class MatchVector >
void ComputeLCBs_v2( const MatchVector& meml, const std::vector<gnSeqI>& breakpoints, std::vector< MatchVector >& lcb_list )
{
	// there must be at least one end of a block defined
	if( breakpoints.size() < 1 )
		return;
		
	lcb_list.clear();
	
	// organize the LCBs into different MatchVector instances
	std::vector<gnSeqI>::const_iterator break_iter = breakpoints.begin();
	uint prev_break = 0;	// prev_break is the first match in the current block
	MatchVector lcb;
	for( ; break_iter != breakpoints.end(); ++break_iter ){
		// add the new MatchList to the set if it made the cut
		lcb_list.push_back( lcb );
		lcb_list.back().insert( lcb_list.back().end(), meml.begin() + prev_break, meml.begin() + *break_iter + 1 );
		prev_break = *break_iter + 1;
	}
}


template <class MatchVector>
void computeLCBAdjacencies_v3( const std::vector< MatchVector >& lcb_list, std::vector< double >& weights, std::vector< mems::LCB >& adjacencies )
{
	adjacencies.clear(); // start with no LCB adjacencies
	if( lcb_list.size() == 0 )
		return;	// there aren't any LCBs so there aren't any adjacencies!

	uint seq_count = lcb_list.front().front()->SeqCount();
	uint seqI;
	uint lcbI;
	for( lcbI = 0; lcbI < lcb_list.size(); ++lcbI ){
		mems::LCB lcb;
		std::vector<gnSeqI> left_end;
		std::vector<gnSeqI> length;
		std::vector<bool> orientation;
		FindBoundaries( lcb_list[lcbI], left_end, length, orientation );

		lcb.left_adjacency = std::vector<uint>( left_end.size(), -1 );
		lcb.right_adjacency = std::vector<uint>( left_end.size(), -1 );
		lcb.left_end = std::vector<int64>( left_end.size(), 0 );
		lcb.right_end = std::vector<int64>( left_end.size(), 0 );

		for( seqI = 0; seqI < seq_count; seqI++ ){
			// support "ragged edges" on the ends of LCBs
			if( left_end[seqI] == mems::NO_MATCH )
				continue;
			lcb.left_end[seqI] = left_end[seqI];
			lcb.right_end[seqI] = left_end[seqI] + length[seqI];
			if( !orientation[seqI] )
			{
				lcb.left_end[seqI] = -lcb.left_end[seqI];
				lcb.right_end[seqI] = -lcb.right_end[seqI];
			}
		}
		lcb.lcb_id = adjacencies.size();
		lcb.weight = weights[ lcbI ];
		lcb.to_be_deleted = false;
		adjacencies.push_back( lcb );
	}

	for( seqI = 0; seqI < seq_count; seqI++ ){
		mems::LCBLeftComparator llc( seqI );
		std::sort( adjacencies.begin(), adjacencies.end(), llc );
		for( lcbI = 1; lcbI + 1 < lcb_list.size(); lcbI++ ){
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
			adjacencies[ lcbI ].right_adjacency[ seqI ] = adjacencies[ lcbI + 1 ].lcb_id;
		}
		if( lcbI == lcb_list.size() )
			lcbI--;	// need to decrement when there is only a single LCB

		// set first and last lcb adjacencies to -1
		adjacencies[ 0 ].left_adjacency[ seqI ] = (uint)-1;
		adjacencies[ lcbI ].right_adjacency[ seqI ] = (uint)-1;
		if( lcbI > 0 ){
			adjacencies[ 0 ].right_adjacency[ seqI ] = adjacencies[ 1 ].lcb_id;
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
		}
	}
	mems::LCBIDComparator lic;
	std::sort( adjacencies.begin(), adjacencies.end(), lic );

}


void projectIntervalList( IntervalList& iv_list, vector< uint >& projection, vector< vector< MatchProjectionAdapter* > >& LCB_list, vector< LCB >& projected_adjs )
{
	vector< MatchProjectionAdapter* > mpa_list;
	// construct pairwise Interval projections
	for( size_t corI = 0; corI < iv_list.size(); corI++ )
	{
		size_t projI = 0;
		for( ; projI < projection.size(); ++projI )
			if( iv_list[corI].LeftEnd(projection[projI]) == NO_MATCH )
				break;
		if( projI != projection.size() )
			continue;
		MatchProjectionAdapter mpa_tmp( &iv_list[corI], projection );
		mpa_list.push_back( mpa_tmp.Copy() );
		if( mpa_list.back()->Orientation(0) == AbstractMatch::reverse )
			mpa_list.back()->Invert();
	}
	vector< gnSeqI > breakpoints;
	IdentifyBreakpoints( mpa_list, breakpoints );
	ComputeLCBs_v2( mpa_list, breakpoints, LCB_list );
	vector< double > lcb_scores( LCB_list.size(), 0 );
	computeLCBAdjacencies_v3( LCB_list, lcb_scores, projected_adjs );
}

//
// end baaad
//



int main( int argc, char* argv[] )
{
	if( argc < 6 )
	{
		cerr << "Usage: alignmentProjector <input xmfa> <output xmfa> <mfa seq input> <mfa seq output> <list of seqs to include, starting at 0>\n";
		return -1;
	}
	ifstream aln_in;
	aln_in.open( argv[1] );
	if( !aln_in.is_open() ){
		cerr << "Error opening " << argv[1] << endl;
		return -1;
	}
	ofstream aln_out;
	aln_out.open( argv[2] );
	if( !aln_out.is_open() ){
		cerr << "Error writing to " << argv[2] << endl;
		return -1;
	}
	string mfa_seqs = argv[3];
	string mfa_output = argv[4];
	
	try{
		IntervalList input_ivs;
		input_ivs.ReadStandardAlignment( aln_in );
		aln_in.close();

		MatchList ml;
		ml.seq_filename = input_ivs.seq_filename;
		ml.LoadMFASequences( mfa_seqs, 7, NULL, false );
		input_ivs.seq_table = ml.seq_table;

		// create a projection list
		vector< uint > projection;
		IntervalList proj_ivs;
		for( int i = 5; i < argc; ++i )
		{
			projection.push_back( atoi( argv[i] ) );
			proj_ivs.seq_filename.push_back( mfa_seqs );
			proj_ivs.seq_table.push_back( input_ivs.seq_table[projection.back()] );
		}

		vector< vector< MatchProjectionAdapter* > > LCB_list;
		vector< LCB > projected_adjs;
		projectIntervalList( input_ivs, projection, LCB_list, projected_adjs );

		cout << "projection has " << LCB_list.size() << " LCBs\n";
		proj_ivs.resize( LCB_list.size() );
		for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
			proj_ivs[lcbI].SetMatches( LCB_list[lcbI] );

		proj_ivs.WriteStandardAlignment( aln_out );

		gnSequence seq;
		seq.LoadSource( mfa_seqs );
		ofstream seq_out( mfa_output.c_str() );
		gnSequence proj_seq;
		for( size_t projI = 0; projI < projection.size(); ++projI )
			proj_seq += seq.contig(projection[projI]);
		gnFASSource::Write(proj_seq,seq_out,false,false);

	}catch( gnException& gne ){
		cerr << gne << endl;
		return -1;
	}catch( exception& e ){
		cerr << e.what() << endl;
		return -2;
	}catch( char const* c ){
		cerr << c << endl;
		return -3;
	}catch(...){
		cerr << "Unhandled exception" << endl;
		return -4;
	}
}

