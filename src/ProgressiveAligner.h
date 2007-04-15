/*******************************************************************************
 * $Id: ProgressiveAligner.h,v 1.23 2004/04/19 23:10:13 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifndef _ProgressiveAligner_h_
#define _ProgressiveAligner_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/Aligner.h"
#include "libMems/PhyloTree.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/Islands.h"
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/multi_array.hpp>
#include "SeedOccurrenceList.h"
#include "libMems/SubstitutionMatrix.h"
#include "libMems/MatchProjectionAdapter.h"

// this is a 99.9% score threshold derived from the EVD of
// simulations of homolgous sequence diverged to .75 substitutions per site and .05 indels per site
const mems::score_t DEFAULT_ISLAND_SCORE_THRESHOLD = 2727;

class SuperInterval
{
public:

	SuperInterval();
	/**
	 * Creates a new SuperInterval.
	 */
	SuperInterval( const mems::Interval& reference_iv );
	SuperInterval(const SuperInterval& siv);
	SuperInterval& operator=(const SuperInterval& siv);
	~SuperInterval(){}
	
	/** Returns the length */
	virtual gnSeqI Length() const { return length; }

	/** Sets the length to @param len */
	virtual void SetLength( gnSeqI len );

	virtual int64 LeftEnd() const { return left_end; }

	virtual void SetLeftEnd( const int64& left_end ) { this->left_end = left_end; }

	mems::Interval reference_iv;

	/** the index of the SuperInterval this is aligned to in c1 */
	size_t c1_siv;
	/** the index of the SuperInterval this is aligned to in c2 */
	size_t c2_siv;
	/** the index of the SuperInterval this is aligned to in the parent */
	size_t parent_siv;

	void CropLeft( gnSeqI amount );
	void CropRight( gnSeqI amount );

	bool operator<( const SuperInterval& si ) const{ return left_end < si.left_end; }

	void ValidateSelf() const;
protected:
	int64 left_end;
	int64 length;
};


class AlignmentTreeNode : public TreeNode
{
public:
	std::vector< SuperInterval > ordering;
	std::vector< boolean > parents_aligned;
	std::vector< boolean > children_aligned;
	genome::gnSequence* sequence;	/**< The sequence associated with this node, NULL for ancestral nodes */
	bool refined;	/**< true if iterative refinement has been applied to the alignment at this node */
};

template <class MatchType>
class LcbTrackingMatch
{ 
public:
	MatchType original_match;
	MatchType node_match;
	size_t match_id;	// used to index into global arrays of lcb_id and score
//	boost::multi_array< size_t, 2 > lcb_id;
//	boost::multi_array< double, 2 > score;
};
typedef LcbTrackingMatch< mems::AbstractMatch* > TrackingMatch;

/** 
 * This class is used to track relationships between LCBs during the LCB determination process.
 */
template <class MatchType>
class TrackingLCB
{
public:
	TrackingLCB(){}
	TrackingLCB( const TrackingLCB& l ){ *this = l; }
	/** Constructs a TrackingLCB from a pairwise LCB */
	TrackingLCB( const mems::LCB& l ){ *this = l; }
	TrackingLCB& operator=( const mems::LCB& l )
	{
		left_end[0] = l.left_end[0];
		left_end[1] = l.left_end[1];
		right_end[0] = l.right_end[0];
		right_end[1] = l.right_end[1];
		left_adjacency[0] = l.left_adjacency[0];
		left_adjacency[1] = l.left_adjacency[1];
		right_adjacency[0] = l.right_adjacency[0];
		right_adjacency[1] = l.right_adjacency[1];
		lcb_id = l.lcb_id;
		weight = l.weight;
		to_be_deleted = false;
		return *this;
	}
	int64 left_end[2];	/**< The left end position of the LCB in each sequence */
	int64 right_end[2];  /**< The right end position of the LCB in each sequence */
	uint left_adjacency[2];	/**< 'Pointers' (actually IDs) to the LCBs on the left in each sequence */
	uint right_adjacency[2];	/**< 'Pointers' (actually IDs) to the LCBs on the right in each sequence */
	double weight;		/**< The weight (or coverage) of this LCB */
	std::vector< MatchType > matches;
	int lcb_id;			/**< A numerical ID that can be assigned to this LCB */
	bool to_be_deleted;
};


struct AlnProgressTracker
{
	gnSeqI total_len;
	gnSeqI cur_leftend;
	double prev_progress;
};


double getDefaultBreakpointPenalty( std::vector< genome::gnSequence* >& sequences );

/**
 * Extends the Aligner class to compute alignments in parallel using MPI
 */
class ProgressiveAligner : public mems::Aligner
{
public:
	/** 
	 * Constructs an aligner for the specified number of sequences.
	 * @param seq_count 	The number of sequences that will be aligned with this Aligner
	 */
	ProgressiveAligner( uint seq_count );
	ProgressiveAligner( const ProgressiveAligner& al );
	ProgressiveAligner& operator=( const ProgressiveAligner& al );
	~ProgressiveAligner();

	/** sets the breakpoint penalty */
	void setBreakpointPenalty( double bp_penalty ){ breakpoint_penalty = bp_penalty; }
	/** assume all genomes are collinear when set to true */
	void setCollinear( boolean collinear ){ this->collinear_genomes = collinear; }
	/** use a list of precomputed matches instead of computing them */
	void setPairwiseMatches( mems::MatchList& pair_ml );
	/** use a precomputed guide tree stored in the given file */
	void setInputGuideTreeFileName( std::string& fname ){ this->input_guide_tree_fname = fname; }
	/** write the guide tree stored to the given file */
	void setOutputGuideTreeFileName( std::string& fname ){ this->output_guide_tree_fname = fname; }
	/** set the max length (in columns) of alignments passed to MUSCLE */
	void SetMaxGappedAlignmentLength( size_t len );

	/** Set whether iterative refinement using MUSCLE should be performed (true/false) */
	void setRefinement( bool refine ){ this->refine = refine; }
	/** Set whether iterative refinement using MUSCLE should be performed (true/false) */
	void setGappedAlignment( bool do_gapped_alignment ){ this->gapped_alignment = do_gapped_alignment; }

	void setPairwiseScoringScheme( const mems::PairwiseScoringScheme& pss ){ this->subst_scoring = pss; }

	enum LcbScoringScheme
	{
		AncestralScoring,
		AncestralSumOfPairsScoring,
		ExtantSumOfPairsScoring
	};

	/** set LCB the scoring scheme */
	void setLcbScoringScheme( LcbScoringScheme scheme ){ scoring_scheme = scheme; }
	LcbScoringScheme getLcbScoringScheme(void){ return scoring_scheme; }

	void setUseLcbWeightScaling( bool use_weight_scaling ){ this->use_weight_scaling = use_weight_scaling; }
	bool getUseLcbWeightScaling(void){ return this->use_weight_scaling; }

	void setBreakpointDistanceScale( double bp_dist_scale ){ this->bp_dist_scale = bp_dist_scale; }
	double getBreakpointDistanceScale(void){ return this->bp_dist_scale; }

	void setConservationDistanceScale( double conservation_dist_scale ){ this->conservation_dist_scale = conservation_dist_scale; }
	double getConservationDistanceScale(void){ return this->conservation_dist_scale; }

	void setBpDistEstimateMinScore( double min_score ){ this->bp_dist_estimate_score = min_score; }
	double getBpDistEstimateMinScore(void){ return this->bp_dist_estimate_score; }

	/** determine which extant sequences have been aligned at a given node */
	void getAlignedChildren( node_id_t node, std::vector< node_id_t >& descendants );

	/** chooses an ordering for aligned intervals at an ancestor node */
	void createAncestralOrdering( std::vector< mems::Interval* >& interval_list, std::vector< SuperInterval >& ancestral_sequence );

	/** constructs an alignment of node1 and node2 at their ancestor */
	void alignProfileToProfile( node_id_t node1, node_id_t node2, node_id_t ancestor );

	/** align the sequences at the designated pair of alignment tree nodes */
	void alignNodes( node_id_t node1, node_id_t node2, node_id_t ancestor );


	/** Given a set of sequences, construct and output an alignment as an IntervalList */
	void align( std::vector< genome::gnSequence* >& seq_table, mems::IntervalList& interval_list );

	void getPath( node_id_t first_n, node_id_t last_n, std::vector< node_id_t >& path );
	template<class MatchType>
	void propagateDescendantBreakpoints( node_id_t node1, uint seqI, std::vector< MatchType* >& iv_list );
	void linkSuperIntervals( node_id_t node1, uint seqI, node_id_t ancestor );
	void recursiveApplyAncestralBreakpoints( node_id_t ancestor );
	void extractAlignment( node_id_t ancestor, size_t super_iv, mems::GappedAlignment& gal );
	void extractAlignment( node_id_t ancestor, size_t super_iv, mems::CompactGappedAlignment<>& cga );

	void getPairwiseMatches( const std::vector< node_id_t >& node1_seqs, const std::vector< node_id_t >& node2_seqs, Matrix<mems::MatchList>& pairwise_matches );
	void getAncestralMatches( const std::vector< node_id_t > node1_seqs, const std::vector< node_id_t > node2_seqs, node_id_t node1, node_id_t node2, node_id_t ancestor, std::vector< mems::AbstractMatch* >& ancestral_matches );
	void getRepresentativeAncestralMatches( const std::vector< node_id_t > node1_seqs, const std::vector< node_id_t > node2_seqs, node_id_t node1, node_id_t node2, node_id_t ancestor, std::vector< mems::AbstractMatch* >& ancestral_matches );
	
	// functions for recursive anchor search
	
	template<class GappedAlignmentType>
	void recurseOnPairs( const std::vector<node_id_t>& node1_seqs, 
		const std::vector<node_id_t>& node2_seqs, const GappedAlignmentType& iv, 
		Matrix<mems::MatchList>& matches, Matrix< std::vector< mems::search_cache_t > >& search_cache_db, 
		Matrix< std::vector< mems::search_cache_t > >& new_cache_db );
	void pairwiseAnchorSearch( mems::MatchList& r_list, mems::Match* r_begin, mems::Match* r_end );

	void translateGappedCoordinates( std::vector<mems::AbstractMatch*>& ml, uint seqI, node_id_t extant, node_id_t ancestor );

	void doGappedAlignment( node_id_t ancestor, bool profile_aln );
	void refineAlignment( mems::GappedAlignment& gal, node_id_t ancestor, bool profile_aln, AlnProgressTracker& apt );
	void FixLeftEnds( node_id_t ancestor );
	void ConstructSuperIntervalFromMSA( node_id_t ancestor, size_t ans_siv, mems::GappedAlignment& gal );

	// determines LCBs among each pair of genomes using a somewhat stringent homology 
	// criteria.  fills the distance matrix with the number of breakpoints between each pair
	void CreatePairwiseBPDistance( boost::multi_array<double, 2>& bp_distmat );

	void constructLcbTrackingMatches( node_id_t ancestral_node, std::vector< mems::AbstractMatch* >& ancestral_matches, std::vector< LcbTrackingMatch< mems::AbstractMatch* > >& tracking_matches );

	void pairwiseScoreTrackingMatches( 
						std::vector< TrackingMatch >& tracking_matches, 
						std::vector<node_id_t>& node1_descendants, 
						std::vector<node_id_t>& node2_descendants,
						boost::multi_array< double, 3 >& tm_score_array
						);

	void computeAvgAncestralMatchScores( 
						std::vector< TrackingMatch >& tracking_matches, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants,
						boost::multi_array< double, 3 >& tm_score_array
						);

	void computeInternalNodeDistances( 
						boost::multi_array<double, 2>& bp_dist_mat, 
						boost::multi_array<double, 2>& cons_dist_mat, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants);

	bool validateSuperIntervals(node_id_t node1, node_id_t node2, node_id_t ancestor);
	bool validatePairwiseIntervals(node_id_t node1, node_id_t node2, std::vector<mems::Interval*>& pair_iv);


	void alignPP(mems::IntervalList& prof1, mems::IntervalList& prof2, mems::IntervalList& interval_list );

protected:
	void getAlignment( mems::IntervalList& interval_list );

	mems::MatchList original_ml;	/**< The list of matches calculated among all sequences.  Also contains the full sequences and sorted mer lists */
	PhyloTree< AlignmentTreeNode > alignment_tree;
	std::vector< uint > node_sequence_map;
	double breakpoint_penalty;
	std::string input_guide_tree_fname;
	std::string output_guide_tree_fname;
	boolean debug;
	boolean refine;

	std::vector< SeedOccurrenceList > sol_list;
	boost::multi_array<double, 2> bp_distance;	/**< pairwise breakpoint distances.  dims will be [seq_count][seq_count] */
	boost::multi_array<double, 2> conservation_distance;	/**< pairwise genome conservation distances.  dims will be [seq_count][seq_count] */

	LcbScoringScheme scoring_scheme;
	bool use_weight_scaling;

	double bp_dist_scale;
	double conservation_dist_scale;

	double bp_dist_estimate_score;	/**< the minimum LCB score to use when estimating BP distance.  should be conservative (high) */

	size_t max_gapped_alignment_length;

	mems::PairwiseScoringScheme subst_scoring;
};

extern bool debug_aligner;

typedef mems::UngappedLocalAlignment< mems::HybridAbstractMatch<> > ULA;

typedef std::vector< std::vector< ULA* > > backbone_list_t;
void applyIslands( mems::IntervalList& iv_list, backbone_list_t& bb_list, const mems::PairwiseScoringScheme& subst_scoring, mems::score_t score_threshold = DEFAULT_ISLAND_SCORE_THRESHOLD );
void writeBackboneSeqCoordinates( backbone_list_t& bb_list, mems::IntervalList& iv_list, std::ostream& bb_out );
void writeBackboneColumns( std::ostream& bb_out, backbone_list_t& bb_list );

	/** Select the next pair of nodes to align
	 *  The chosen pair will either be unaligned extant sequences or unaligned
	 *  ancestral sequences whose descendants have all been aligned.  The chosen pair has
	 *  the shortest path on the tree
	 *  When no sequences remain to be aligned, returns node1 == node2
	 */
void chooseNextAlignmentPair( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t& node1, node_id_t& node2, node_id_t& ancestor );

void markAligned( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t subject_node, node_id_t neighbor );

node_id_t createAlignmentTreeRoot( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t node1, node_id_t node2 );

// homogenizes an alignment tree and ordering to prepare for alignment
void prepareAlignmentTree( PhyloTree< AlignmentTreeNode >& alignment_tree );

inline
ProgressiveAligner::~ProgressiveAligner()
{
	for( size_t mI = 0; mI < original_ml.size(); mI++ )
		original_ml[mI]->Free();
}

template<class T>
class AbsolutComparator
{
public:
	boolean operator()(const T& a, const T& b) const
	{
		return (genome::absolut(a) < genome::absolut(b));
	}
};



template <class MatchVector>
void processNewMatch( uint seqI, MatchVector& new_matches, typename MatchVector::value_type& new_match )
{
	new_match->SetStart( seqI, 0 );
	if( new_match->Multiplicity() > 1 && new_match->Length(seqI) > 0 )
		new_matches.push_back( new_match );
	else
	{
		new_match->Free();
		new_match = NULL;
	}
}

/**
 * Delete overlapping regions in favor of the larger match.
 * This code isn't perfect, it can delete too many base pairs in some cases
 * @param	ml	The vector of matches
 * @param	seq_ids	The indexes of sequences in which overlaps should be eliminated
 * @param	eliminate_both	Delete both of the overlapping matches, instead of leaving one remaining
 */
template <class MatchVector>
void EliminateOverlaps_v2( MatchVector& ml, const std::vector< uint >& seq_ids, bool eliminate_both = false ){
	if( ml.size() < 2 )
		return;
	uint seq_count = ml[0]->SeqCount();
	for( uint sidI = 0; sidI < seq_ids.size(); sidI++ ){
		uint seqI = seq_ids[ sidI ];
		mems::SingleStartComparator<mems::AbstractMatch> msc( seqI );
		std::sort( ml.begin(), ml.end(), msc );
		int64 matchI = 0;
		int64 nextI = 0;
		int64 deleted_count = 0;
		MatchVector new_matches;

		// scan forward to first defined match
		for(; matchI != ml.size(); matchI++ )
			if( ml[ matchI ]->Start( seqI ) != mems::NO_MATCH )
				break;

		for(; matchI < ml.size(); matchI++ ){
			if( ml[ matchI ] == NULL )
				continue;
			
			for( nextI = matchI + 1; nextI < ml.size(); nextI++ ){
				if( ml[ nextI ] == NULL )
					continue;

				boolean deleted_matchI = false;
				// check for overlaps
				int64 startI = ml[ matchI ]->Start( seqI );
				int64 lenI = ml[ matchI ]->Length( seqI );
				int64 startJ = ml[ nextI ]->Start( seqI );
				int64 diff =  genome::absolut( startJ ) - genome::absolut( startI ) - lenI;

				if( diff >= 0 )
					break;	// there are no more overlaps

				diff = -diff;
				typename MatchVector::value_type new_match;
				bool mem_iter_smaller = ( ml[ nextI ]->Multiplicity() > ml[ matchI ]->Multiplicity() ) ||
					( ml[ nextI ]->Multiplicity() == ml[ matchI ]->Multiplicity() && ml[ nextI ]->Length(seqI) > ml[ matchI ]->Length(seqI) );

				// delete bases from the smaller match
				if( eliminate_both || mem_iter_smaller )
				{
					// mem_iter is smaller
					new_match = ml[matchI]->Copy();
					// erase base pairs from new_match
					if( diff >= lenI ){
//							cerr << "Deleting " << **mem_iter << " at the hands of\n" << **next_iter << endl;
						ml[ matchI ]->Free();
						ml[ matchI ] = NULL;
						matchI--;
						deleted_matchI = true;
						deleted_count++;
					}else{
						ml[ matchI ]->CropRight( diff, seqI );
						new_match->CropLeft( new_match->Length(seqI) - diff, seqI );
					}
					processNewMatch( seqI, new_matches, new_match );
				}
				if( eliminate_both || !mem_iter_smaller )
				{
					// match_iter is smaller
					new_match = ml[nextI]->Copy();
					// erase base pairs from new_match
					if( diff >= ml[ nextI ]->Length(seqI) ){
//							cerr << "Deleting " << **next_iter << " at the hands of\n" << **mem_iter << endl;
						ml[ nextI ]->Free();
						ml[ nextI ] = NULL;
						deleted_count++;
					}else{
						ml[ nextI ]->CropLeft( diff, seqI );
						new_match->CropRight( new_match->Length(seqI) - diff, seqI );
					}
					processNewMatch( seqI, new_matches, new_match );
				}
				if( deleted_matchI )
					break;
			}
		}

		if( deleted_count > 0 ){
			size_t cur = 0;
			for( size_t mI = 0; mI < ml.size(); ++mI )
				if( ml[mI] != NULL )
					ml[cur++] = ml[mI];
			ml.erase( ml.begin() + cur, ml.end() );
		}
		ml.insert( ml.end(), new_matches.begin(), new_matches.end() );
		new_matches.clear();
	}
}

template <class MatchVector>
void EliminateOverlaps_v2( MatchVector& ml )
{
	if( ml.size() < 2 )
		return;	// can't eliminate overlaps between fewer than 2 matches
	uint seq_count = ml[0]->SeqCount();
	std::vector< uint > seq_ids( seq_count );
	for( uint i = 0; i < seq_count; ++i )
		seq_ids[i] = i;
	EliminateOverlaps_v2( ml, seq_ids );
};

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

/**
 *  Redesign to be more intuitive.  left_adjacency is always left, regardless of LCB orientation
 */
inline
void computeLCBAdjacencies_v3( mems::IntervalList& iv_list, std::vector< double >& weights, std::vector< mems::LCB >& adjacencies ){
	std::vector< std::vector< mems::Interval* > > nivs;
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
		nivs.push_back( std::vector< mems::Interval* >( 1, &iv_list[ivI] ) );
	computeLCBAdjacencies_v3( nivs, weights, adjacencies );
}

/**
 * Takes a set of filtered LCB adjacencies and an unfiltered set of matches as input
 * returns a filtered set of matches that reflects the LCBs found
 */
template< class MatchVector >
void filterMatches_v2( std::vector< mems::LCB >& adjacencies, std::vector< MatchVector >& lcb_list, std::vector< double >& weights, MatchVector& deleted_matches ){
	if( lcb_list.size() < 1 )
		return;
	MatchVector lcb_tmp = lcb_list[ 0 ];
	lcb_tmp.clear();
	std::vector< MatchVector > filtered_lcbs( lcb_list.size(), lcb_tmp );
	uint lcbI;
	for( lcbI = 0; lcbI < adjacencies.size(); lcbI++ ){
		if( adjacencies[ lcbI ].lcb_id == lcbI ){
			filtered_lcbs[ lcbI ].insert( filtered_lcbs[ lcbI ].end(), lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end() );
			continue;
		}
		if( adjacencies[ lcbI ].lcb_id == -1 ){
			std::cerr << "weird";
			continue; 	// this one was removed
		}
		if( adjacencies[ lcbI ].lcb_id == -2 )
		{
			deleted_matches.insert( deleted_matches.end(), lcb_list[lcbI].begin(), lcb_list[lcbI].end() );
			continue; 	// this one was removed
		}

		// this one points elsewhere
		// search and update the union/find structure for the target
		std::stack< uint > visited_lcbs;
		visited_lcbs.push( lcbI );
		uint cur_lcb = adjacencies[ lcbI ].lcb_id;
		while( adjacencies[ cur_lcb ].lcb_id != cur_lcb ){
			visited_lcbs.push( cur_lcb );
			cur_lcb = adjacencies[ cur_lcb ].lcb_id;
			if( cur_lcb == -1 || cur_lcb == -2 ){
//				std::cerr << "improper hoodidge\n";
				break;	// this one points to an LCB that got deleted
			}
		}
		while( visited_lcbs.size() > 0 ){
			adjacencies[ visited_lcbs.top() ].lcb_id = cur_lcb;
			visited_lcbs.pop();
		}
		// add this LCB's matches to the target LCB.
		if( cur_lcb != -1 && cur_lcb != -2 )
			filtered_lcbs[ cur_lcb ].insert( filtered_lcbs[ cur_lcb ].end(), lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end() );
		else
			deleted_matches.insert( deleted_matches.end(), lcb_list[lcbI].begin(), lcb_list[lcbI].end() );
	}


	lcb_list.clear();
	std::vector< double > new_weights;
	for( lcbI = 0; lcbI < filtered_lcbs.size(); lcbI++ ){
		if( filtered_lcbs[ lcbI ].size() > 0 ){
			lcb_list.push_back( filtered_lcbs[ lcbI ] );
			new_weights.push_back( weights[lcbI] );
		}
	}

	// sort the matches inside consolidated LCBs
	mems::MatchStartComparator<mems::AbstractMatch> msc( 0 );
	for( lcbI = 0; lcbI < lcb_list.size(); lcbI++ ){
		std::sort( lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end(), msc );
	}

	// calculate the LCB adjacencies
	weights = new_weights;
	computeLCBAdjacencies_v3( lcb_list, weights, adjacencies );

}

template< class MatchVector >
uint64 SimpleGetLCBCoverage( MatchVector& lcb ){
	typename MatchVector::iterator match_iter = lcb.begin();
	uint64 coverage = 0;
	bool debug = true;
	for( ; match_iter != lcb.end(); ++match_iter ){
		double maxlen = 0;
		double minlen = 0;
		for( uint seqI = 0; seqI < (*match_iter)->SeqCount(); seqI++ )
		{
			if( (*match_iter)->LeftEnd(seqI) != mems::NO_MATCH )
			{
				maxlen += (*match_iter)->Length(seqI);
				if( (*match_iter)->Length(seqI) > minlen )
					minlen = (*match_iter)->Length(seqI);
			}
		}
		double score = exp( ((*match_iter)->AlignmentLength() - minlen) / (maxlen - minlen) );
		score *= maxlen;
		coverage += score;
	}
	return coverage;
}

template< class MatchVectorType >
void addUnalignedIntervals_v2( MatchVectorType& iv_list, std::set< uint > seq_set, std::vector<gnSeqI> seq_lengths )
{
	std::vector< mems::LCB > adjacencies;
	uint lcbI;
	uint seqI;
	uint seq_count = seq_lengths.size();


	if( seq_set.size() == 0 )
	{
		// if an empty seq set was passed then assume all seqs
		// should be processed
		for( seqI = 0; seqI < seq_count; seqI++ )
			seq_set.insert( seqI );
	}
	std::vector< std::vector< typename MatchVectorType::value_type > > ymmv;
	for( size_t ivI = 0; ivI < iv_list.size(); ++ivI )
		ymmv.push_back( std::vector< typename MatchVectorType::value_type >( 1, iv_list[ivI] ) );

	std::vector< double > scores( iv_list.size(), 0 );
	computeLCBAdjacencies_v3( ymmv, scores, adjacencies );

	std::vector< int > rightmost;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		rightmost.push_back( -1 );
	}

	for( lcbI = 0; lcbI <= adjacencies.size(); lcbI++ ){
		std::set< uint >::iterator seq_set_iterator = seq_set.begin();
		for( ; seq_set_iterator != seq_set.end(); seq_set_iterator++ ){
			seqI = *seq_set_iterator;
			// scan left
			int leftI;
			if( lcbI < adjacencies.size() ){
// left is always to the left!!
				leftI = adjacencies[ lcbI ].left_adjacency[ seqI ];
			}else
				leftI = rightmost[ seqI ];

			int rightI = lcbI < adjacencies.size() ? lcbI : -1;
// right is always to the right!!
			if( lcbI < adjacencies.size() )
				if( adjacencies[ lcbI ].right_adjacency[ seqI ] == -1 )
					rightmost[ seqI ] = lcbI;
			
			int64 left_start, right_start;
			mems::getGapBounds( seq_lengths, adjacencies, seqI, leftI, rightI, left_start, right_start );
			int64 gap_len =  genome::absolut( right_start ) - genome::absolut( left_start );
			if( gap_len > 0 ){
				mems::Match mm( seq_count );
				mems::Match* m = mm.Copy();
				for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
					m->SetStart( seqJ, 0 );
				}
				m->SetStart( seqI, left_start );
				m->SetLength( gap_len );
				mems::Interval iv;
				std::vector< mems::AbstractMatch* > tmpvec(1, m);
				iv.SetMatches( tmpvec );
				iv_list.push_back( iv.Copy() );
			}
		}
	}
}

inline
void projectIntervalList( mems::IntervalList& iv_list, std::vector< uint >& projection, std::vector< std::vector< mems::MatchProjectionAdapter* > >& LCB_list, std::vector< mems::LCB >& projected_adjs )
{
	std::vector< size_t > proj(projection.size());
	for( size_t i = 0; i < projection.size(); ++i )
		proj[i] = projection[i];
	std::vector< mems::MatchProjectionAdapter* > mpa_list;
	// construct pairwise Interval projections
	for( size_t corI = 0; corI < iv_list.size(); corI++ )
	{
		size_t projI = 0;
		for( ; projI < projection.size(); ++projI )
			if( iv_list[corI].LeftEnd(projection[projI]) == mems::NO_MATCH )
				break;
		if( projI != projection.size() )
			continue;
		mems::MatchProjectionAdapter mpa_tmp( &iv_list[corI], proj );
		mpa_list.push_back( mpa_tmp.Copy() );
		if( mpa_list.back()->Orientation(0) == mems::AbstractMatch::reverse )
			mpa_list.back()->Invert();
	}
	std::vector< gnSeqI > breakpoints;
	IdentifyBreakpoints( mpa_list, breakpoints );
	ComputeLCBs_v2( mpa_list, breakpoints, LCB_list );
	std::vector< double > lcb_scores( LCB_list.size(), 0 );
	computeLCBAdjacencies_v3( LCB_list, lcb_scores, projected_adjs );
}



#endif // _ProgressiveAligner_h_
