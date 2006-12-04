/*******************************************************************************
 * $Id: progressiveAligner.cpp,v 1.47 2004/04/19 23:10:30 darling Exp $
 * BEWARE!!
 * This code was created in the likeness of the flying spaghetti monster
 *
 * dedicated to Loren...
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "ProgressiveAligner.h"
#include "libMems/Aligner.h"
#include "libMems/MemSubsets.h"
#include "libMems/Islands.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MuscleInterface.h"	// it's the default gapped aligner
#include "libGenome/gnRAWSource.h"
#include "libMems/gnAlignedSequences.h"
#include "libMems/ClustalInterface.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/MatchProjectionAdapter.h"
#include "PairwiseMatchFinder.h"
#include "TreeUtilities.h"
#include "PairwiseMatchAdapter.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>

#include <map>
#include <fstream>	// for debugging
#include <sstream>
#include <stack>
#include <algorithm>
#include <limits>
#include <iomanip>

using namespace std;
using namespace genome;
using namespace mems;


bool debug_aligner = false;
bool progress_msgs = false;

bool debug_me = false;
static int dbg_count = 0; 	 

const uint LCB_UNASSIGNED = (std::numeric_limits<uint>::max)();

double min_window_size = 200;
double max_window_size = 20000;  // don't feed MUSCLE anything bigger than this
double min_density = .5;
double max_density = .9;
size_t max_gap_length = 3000;
size_t lcb_hangover = 300;


void mergeUnalignedIntervals( uint seqI, vector< Interval* >& iv_list, vector< Interval* >& new_list );

/**
 * Test code to ensure that an individual LCB is truly collinear
 * @return	true if the LCB is good
 */
boolean my_validateLCB( MatchList& lcb ){
	vector< Match* >::iterator lcb_iter = lcb.begin();
	if( lcb.size() == 0 )
		return true;
	uint seq_count = (*lcb_iter)->SeqCount();
	uint seqI = 0;
	boolean complain = false;
	for(; seqI < seq_count; seqI++ ){
		lcb_iter = lcb.begin();
		int64 prev_coord = 0;
		for(; lcb_iter != lcb.end(); ++lcb_iter ){
			if( (*lcb_iter)->Start( seqI ) == NO_MATCH )
				continue;
			else if( prev_coord != 0 && (*lcb_iter)->Start( seqI ) < prev_coord ){
				complain = true;
			}
			prev_coord = (*lcb_iter)->Start( seqI );
		}
	}
	return !complain;
}

template< class BoostMatType >
void print2d_matrix( BoostMatType& mat, std::ostream& os )
{
	for( size_t i = 0; i < mat.shape()[0]; ++i )
	{
		for( size_t j = 0; j < mat.shape()[1]; ++j )
		{
			if( j > 0 )
				os << "\t";
			os << mat[i][j];
		}
		os << endl;
	}
}

void printProgress( uint prev_prog, uint cur_prog, ostream& os )
{
	if( prev_prog != cur_prog )
	{
		if( cur_prog / 10 != prev_prog / 10 )
			os << endl;
		os << cur_prog << "%..";
		os.flush();
	}
}

double getDefaultBreakpointPenalty( std::vector< gnSequence* >& sequences )
{
	uint default_mer_size = MatchList::GetDefaultMerSize( sequences );
	double avg_seq_len = 0;
	for( size_t seqI = 0; seqI < sequences.size(); ++seqI )
		avg_seq_len += (double)sequences[seqI]->length();
	avg_seq_len /= (double)sequences.size();
	avg_seq_len = log( avg_seq_len ) / log( 2.0 );
	return avg_seq_len * 1500;	  // seems to work reasonably well
}


double getDefaultBpDistEstimateMinScore( std::vector< gnSequence* >& sequences )
{
	return 3.0 * getDefaultBreakpointPenalty(sequences);
}

template< class MatchVectorType >
double GetLCBUniquenessScore( MatchVectorType& mlist, vector< SeedOccurrenceList* >& sol_list, uint seed_size )
{
	double uni_score = 0;
	for( size_t mI = 0; mI < mlist.size(); mI++ )
	{
		for( size_t seqI = 0; seqI < sol_list.size(); seqI++ )
		{
			for( size_t merI = 0; merI < mlist[mI]->Length(seqI); merI++ )
			{
				double uni = (double)sol_list[seqI]->getFrequency(mlist[mI]->LeftEnd(seqI) + merI - 1);
				uni_score += 1 / uni;
				// FIXME:  is this a good idea or not?
//				if( merI + seed_size >= mlist[mI]->Length(seqI) )
//					break;
			}
		}
	}
	return uni_score;
}

typedef boost::multi_array< std::vector< TrackingLCB< TrackingMatch* > >, 2 > PairwiseLCBMatrix;

void getPairwiseLCBs( 
	uint nI, 
	uint nJ, 
	uint dI, 
	uint dJ, 
	vector< TrackingMatch* >& tracking_matches, 
	vector< TrackingLCB<TrackingMatch*> >& t_lcbs,
	boost::multi_array< double, 3 >& tm_score_array,
	boost::multi_array< size_t, 3 >& tm_lcb_id_array )
{
	// make a set of projection matches
	vector< AbstractMatch* > pair_matches;
	for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
	{
		if( tracking_matches[mI]->node_match->LeftEnd(nI) == NO_MATCH ||
			tracking_matches[mI]->node_match->LeftEnd(nJ) == NO_MATCH )
			continue;
		PairwiseMatchAdapter pma(tracking_matches[mI]->node_match, nI, nJ );
		pma.tm = tracking_matches[mI];
		if( pma.Orientation(0) == AbstractMatch::reverse )
			pma.Invert();
		pair_matches.push_back(pma.Copy());
	}
	// find LCBs...
	vector< gnSeqI > breakpoints;
	IdentifyBreakpoints( pair_matches, breakpoints );

	vector< vector< AbstractMatch* > > LCB_list;
	ComputeLCBs_v2( pair_matches, breakpoints, LCB_list );

	//
	// now compute scores on them
	//
	vector< double > lcb_scores(LCB_list.size());
	for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
	{
		double lcb_score = 0;
		for( size_t mI = 0; mI < LCB_list[lcbI].size(); ++mI )
		{
			PairwiseMatchAdapter* pma = (PairwiseMatchAdapter*)LCB_list[lcbI][mI];
			lcb_score += tm_score_array[pma->tm->match_id][dI][dJ];
		}
		lcb_scores[lcbI] = lcb_score;
	}

	// and build the pairwise adjacency list
	vector< LCB > adjacencies;
	computeLCBAdjacencies_v3( LCB_list, lcb_scores, adjacencies );

	t_lcbs.resize(adjacencies.size());
	for( size_t lcbI = 0; lcbI < adjacencies.size(); ++lcbI )
	{
		t_lcbs[lcbI] = adjacencies[lcbI];
		for( size_t mI = 0; mI < LCB_list[lcbI].size(); ++mI )
			t_lcbs[lcbI].matches.push_back( ((PairwiseMatchAdapter*)LCB_list[lcbI][mI])->tm );
		// sort them by ptr
		sort( t_lcbs[lcbI].matches.begin(), t_lcbs[lcbI].matches.end() );

		// set the match LCB ids appropriately
		for( size_t mI = 0; mI < t_lcbs[lcbI].matches.size(); ++mI )
			tm_lcb_id_array[t_lcbs[lcbI].matches[mI]->match_id][dI][dJ] = lcbI;
	}

	// free the memory used by pairwise matches
	for( size_t mI = 0; mI < pair_matches.size(); ++mI )
		pair_matches[mI]->Free();
}


/** removes an LCB from an LCB list and coalesces surrounding LCBs.  Returns the number of LCBs removed 
 *  After LCBs are removed, the adjacency list should be processed with filterLCBs()
 *  @param	id_remaps	This is populated with a list of LCB ids that were deleted or coalesced and now have a new LCB id
 *                      for each coalesced LCB, an entry of the form <old id, new id> is added, deleted LCBs have
 *						entries of the form <deleted, -1>.  Entries appear in the order operations were performed
 *						and the function undoLcbRemoval() can undo these operations in reverse order
 */
template< class LcbVector >
uint RemoveLCBandCoalesce( size_t lcbI, LcbVector& adjacencies, std::vector< double >& scores, std::vector< std::pair< uint, uint > >& id_remaps, std::vector< uint >& impact_list )
{
	uint seq_count = adjacencies[0].left_end.size();
	uint removed_count = 0;
	vector< uint > imp_tmp(seq_count * 10, LCB_UNASSIGNED);
	swap(impact_list, imp_tmp);
	size_t impactI = 0;
	id_remaps.clear();

	adjacencies[ lcbI ].lcb_id = -2;
	
	// update adjacencies
	uint seqI;
	uint left_adj;
	uint right_adj;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		left_adj = adjacencies[ lcbI ].left_adjacency[ seqI ];
		right_adj = adjacencies[ lcbI ].right_adjacency[ seqI ];
		if( left_adj != -1 )
			adjacencies[ left_adj ].right_adjacency[ seqI ] = right_adj;
		if( right_adj != -1 && right_adj != adjacencies.size() )
			adjacencies[ right_adj ].left_adjacency[ seqI ] = left_adj;
	}

	// populate the impact list -- LCBs whose removal scores may change due to this one's removal
	for( seqI = 0; seqI < seq_count; seqI++ ){
		left_adj = adjacencies[ lcbI ].left_adjacency[ seqI ];
		right_adj = adjacencies[ lcbI ].right_adjacency[ seqI ];
		impact_list[impactI++] = left_adj;
		impact_list[impactI++] = right_adj;
		for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
			if( left_adj != -1 )
			{
				impact_list[impactI++] = adjacencies[ left_adj ].left_adjacency[ seqJ ];
				impact_list[impactI++] = adjacencies[ left_adj ].right_adjacency[ seqJ ];
			}
			if( right_adj != -1 )
			{
				impact_list[impactI++] = adjacencies[ right_adj ].left_adjacency[ seqJ ];
				impact_list[impactI++] = adjacencies[ right_adj ].right_adjacency[ seqJ ];
			}
		}
	}

	// just deleted an lcb...
	id_remaps.push_back( make_pair( lcbI, -1 ) );
	removed_count++;

	// check for collapse
	for( seqI = 0; seqI < seq_count; seqI++ ){
		left_adj = adjacencies[ lcbI ].left_adjacency[ seqI ];
		right_adj = adjacencies[ lcbI ].right_adjacency[ seqI ];
		// find the real slim shady
		while( left_adj != -1 && adjacencies[ left_adj ].lcb_id != left_adj )
			left_adj = adjacencies[ left_adj ].left_adjacency[ seqI ];
		while( right_adj != -1 && adjacencies[ right_adj ].lcb_id != right_adj )
			right_adj = adjacencies[ right_adj ].right_adjacency[ seqI ];
		if( left_adj == -1 || right_adj == -1 )
			continue;	// can't collapse with a non-existant LCB!
		if( adjacencies[ left_adj ].lcb_id != left_adj ||
			adjacencies[ right_adj ].lcb_id != right_adj )
			if( seqI > 0 )
				continue;	// already coalesced
			else
				cerr << "trouble on down street\n";

		// check whether the two LCBs are adjacent in each sequence
		boolean orientation = adjacencies[ left_adj ].left_end[ seqI ] > 0 ? true : false;
		uint seqJ;
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			boolean j_orientation = adjacencies[ left_adj ].left_end[ seqJ ] > 0;
			if( j_orientation == orientation &&
				adjacencies[ left_adj ].right_adjacency[ seqJ ] != right_adj )
				break;
			if( j_orientation != orientation &&
				adjacencies[ left_adj ].left_adjacency[ seqJ ] != right_adj )
				break;
			// check that they are both in the same orientation
			if( adjacencies[ right_adj ].left_end[ seqJ ] > 0 != j_orientation )
				break;
		}

		if( seqJ != seq_count ||
			adjacencies[ left_adj ].to_be_deleted ||
			adjacencies[ right_adj ].to_be_deleted )
			continue;	// if these two aren't collinear, or one or both will get deleted, then don't coalesce
		

		// these two can be coalesced
		// do it.  do it now.
		id_remaps.push_back( make_pair( adjacencies[ right_adj ].lcb_id, left_adj ) );
		adjacencies[ right_adj ].lcb_id = left_adj;
		scores[ left_adj ] += scores[ right_adj ];
		adjacencies[ left_adj ].weight += adjacencies[ right_adj ].weight;

		// unlink right_adj from the adjacency list and
		// update left and right ends of left_adj
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			boolean j_orientation = adjacencies[ left_adj ].left_end[ seqJ ] > 0;
			uint rr_adj = adjacencies[ right_adj ].right_adjacency[ seqJ ];
			uint rl_adj = adjacencies[ right_adj ].left_adjacency[ seqJ ];
			if( j_orientation == orientation ){
				adjacencies[ left_adj ].right_end[ seqJ ] = adjacencies[ right_adj ].right_end[ seqJ ];
				adjacencies[ left_adj ].right_adjacency[ seqJ ] = rr_adj;
				if( rr_adj != -1 )
					adjacencies[ rr_adj ].left_adjacency[ seqJ ] = left_adj;
			}else{
				adjacencies[ left_adj ].left_end[ seqJ ] = adjacencies[ right_adj ].left_end[ seqJ ];
				adjacencies[ left_adj ].left_adjacency[ seqJ ] = rl_adj;
				if( rl_adj != -1 )
					adjacencies[ rl_adj ].right_adjacency[ seqJ ] = left_adj;
			}
		}
		// just coalesced two LCBs...
		removed_count++;
	}
	// uniquify the impact list and get rid of empty entries
	std::sort( impact_list.begin(), impact_list.end() );
	vector< uint >::iterator imp_end = std::unique( impact_list.begin(), impact_list.end() );
	vector< uint >::iterator imp_preend = std::lower_bound( impact_list.begin(), imp_end, LCB_UNASSIGNED );
	impact_list.erase( imp_preend, impact_list.end() );

	return removed_count;
}


template< class LcbVector >
void undoLcbRemoval( LcbVector& adjs, std::vector< std::pair< uint, uint > >& id_remaps )
{
	for( size_t rI = id_remaps.size(); rI > 0; --rI )
	{
		if( id_remaps[rI-1].second == -1 )
		{
			// this one was deleted
			// revert adjacencies
			uint lcbI = id_remaps[rI-1].first;
			uint seq_count = adjs[ lcbI ].left_adjacency.size();
			for( uint seqI = 0; seqI < seq_count; seqI++ )
			{
				uint left_adj = adjs[ lcbI ].left_adjacency[ seqI ];
				uint right_adj = adjs[ lcbI ].right_adjacency[ seqI ];
				if( left_adj != -1 )
					adjs[ left_adj ].right_adjacency[ seqI ] = lcbI;
				if( right_adj != -1 && right_adj != adjs.size() )
					adjs[ right_adj ].left_adjacency[ seqI ] = lcbI;
			}
			adjs[lcbI].lcb_id = lcbI;	// reset the lcb id
			adjs[lcbI].to_be_deleted = false;	// no longer TBD
		}else{
			// this one was coalesced
			// uncoalesce it
			uint lcbI = id_remaps[rI-1].first;
			uint lcbJ = id_remaps[rI-1].second;
			adjs[lcbI].lcb_id = lcbI;
			adjs[lcbJ].weight -= adjs[lcbI].weight;
			// link lcbI back in
			// TODO: fix right end and left end coordinates
			uint seq_count = adjs[ lcbI ].left_adjacency.size();
			for( uint seqI = 0; seqI < seq_count; ++seqI )
			{
				uint ladj = adjs[lcbI].left_adjacency[seqI];
				uint radj = adjs[lcbI].right_adjacency[seqI];
				if(  ladj == lcbJ )
				{
					adjs[lcbJ].right_adjacency[seqI] = lcbI;
					if( radj != -1 && radj != adjs.size())
						adjs[radj].left_adjacency[seqI] = lcbI;
				}else
				if(  radj == lcbJ )
				{
					adjs[lcbJ].left_adjacency[seqI] = lcbI;
					if( ladj != -1 && ladj != adjs.size())
						adjs[ladj].right_adjacency[seqI] = lcbI;
				}
			}
		}
	}
}


void printMatch( AbstractMatch* m, ostream& os )
{
	for( size_t ii = 0; ii < m->SeqCount(); ++ii )
	{
		if( ii > 0 )
			os << '\t';
		os << "(" << m->Start(ii) << "," << m->RightEnd(ii) << ")";
	}
}


/**
 * TODO: refactor this so the function calls and names make sense
 */
class EvenFasterSumOfPairsBreakpointScorer
{
public:
	EvenFasterSumOfPairsBreakpointScorer( 
		double breakpoint_penalty,
		boost::multi_array<double,2> bp_weight_matrix, 
		boost::multi_array<double,2> conservation_weight_matrix,
		vector< TrackingMatch* > tracking_match,
		PairwiseLCBMatrix& pairwise_adjacency_matrix,
		vector<node_id_t>& n1_descendants,
		vector<node_id_t>& n2_descendants,
		boost::multi_array< double, 3 >& tm_score_array,
		boost::multi_array< size_t, 3 >& tm_lcb_id_array,
		size_t seqI_begin,
		size_t seqI_end,
		size_t seqJ_begin,
		size_t seqJ_end
		) :
	  bp_penalty( breakpoint_penalty ),
	  bp_weights( bp_weight_matrix ), 
	  conservation_weights( conservation_weight_matrix ),
	  tracking_matches( tracking_match ),
	  pairwise_adjacencies( pairwise_adjacency_matrix ),
	  n1_des(n1_descendants),
	  n2_des(n2_descendants),
      tm_score_array(tm_score_array),
      tm_lcb_id_array(tm_lcb_id_array),
	  seqI_count(pairwise_adjacencies.shape()[0]),
	  seqJ_count(pairwise_adjacencies.shape()[1]),
	  seqI_first(seqI_begin),
	  seqI_last(seqI_end),
	  seqJ_first(seqJ_begin),
	  seqJ_last(seqJ_end)
	  {
		  std::sort(tracking_matches.begin(), tracking_matches.end());
		  pairwise_lcb_count.resize( boost::extents[pairwise_adjacencies.shape()[0]][pairwise_adjacencies.shape()[1]] );
		pairwise_lcb_score.resize( boost::extents[pairwise_adjacencies.shape()[0]][pairwise_adjacencies.shape()[1]] );;
		all_id_remaps.resize( boost::extents[pairwise_lcb_count.shape()[0]][pairwise_lcb_count.shape()[1]] );
		full_impact_list.resize( boost::extents[pairwise_lcb_count.shape()[0]][pairwise_lcb_count.shape()[1]] );
		my_del_lcbs.resize(100);	// buffer for use during lcb removal score computation
		for( size_t i = 0; i < 3; ++i )
		{
			internal_lcb_score_diff[i].resize( boost::extents[pairwise_adjacencies.shape()[0]][pairwise_adjacencies.shape()[1]] );
			internal_lcb_removed_count[i].resize( boost::extents[pairwise_adjacencies.shape()[0]][pairwise_adjacencies.shape()[1]] );
		}
		lsd_zeros.resize( internal_lcb_score_diff[0].num_elements(), 0 );
		lrc_zeros.resize( internal_lcb_removed_count[0].num_elements(), 0 );
		using_lsd = -1;
		size_t max_pair_adj_size = 0;
		for( size_t i = 0; i < seqI_count; ++i )
		{
			for( size_t j = 0; j < seqJ_count; ++j )
			{
				pairwise_lcb_count[i][j] = pairwise_adjacencies[i][j].size();
				pairwise_lcb_score[i][j] = 0;
				max_pair_adj_size = (std::max)(max_pair_adj_size, pairwise_adjacencies[i][j].size());
				for( size_t lcbI = 0; lcbI < pairwise_adjacencies[i][j].size(); ++lcbI )
					pairwise_lcb_score[i][j] += pairwise_adjacencies[i][j][lcbI].weight;
			}
		}
		bogus_scores.resize(max_pair_adj_size+10);
	  };

	/**
	 * Returns the number of possible moves a search algorithm may make from the current 
	 * location in LCB search space.  In this case it's simply the total number of pairwise LCBs
	 */
	size_t getMoveCount()
	{
		size_t move_count = 0;
		for( size_t i = seqI_first; i < seqI_last; ++i )
			for( size_t j = seqJ_first; j < seqJ_last; ++j )
				move_count += pairwise_adjacencies[i][j].size();
		return move_count;
	}

	/** returns the score of the current state */
	double score()
	{
		// score is the sum of all pairwise LCB scores,
		// minus the sum of all pairwise breakpoint penalties
		double score = 0;
		for( size_t seqI = seqI_first; seqI < seqI_last; ++seqI )
		{
			for( size_t seqJ = seqJ_first; seqJ < seqJ_last; ++seqJ )
			{
				const double pw_lcb_score = pairwise_lcb_score[seqI][seqJ];
				// add LCB scores
				score += pairwise_lcb_score[seqI][seqJ];
				// subtract breakpoint penalty
				// subtract 1 from number of LCBs so that a single circular LCB doesn't get penalized
				score -= (bp_penalty * (1-bp_weights[seqI][seqJ]) * (1-conservation_weights[seqI][seqJ]) * (pairwise_lcb_count[seqI][seqJ]-1));
				if( !(score > -1e200 && score < 1e200) )
				{
					genome::breakHere();
					cerr << "bp_weights[seqI][seqJ] " << bp_weights[seqI][seqJ] << endl;
					cerr << "conservation_weights[seqI][seqJ] " << conservation_weights[seqI][seqJ] << endl;
					cerr << "pairwise_lcb_count[seqI][seqJ] " << pairwise_lcb_count[seqI][seqJ] << endl;
					cerr << "pairwise_lcb_score[seqI][seqJ] " << pw_lcb_score << endl;
					cerr << "Invalid score!!\n";
				}
			}
		}
		return score;
	}

	/** scores a move */
	double operator()( pair< double, size_t >& the_move  )
	{
		size_t new_move_count;
		vector< pair< double, size_t > > new_move_list;
		using_lsd++;
		std::copy(lsd_zeros.begin(),lsd_zeros.end(),internal_lcb_score_diff[using_lsd].data());
		std::copy(lrc_zeros.begin(),lrc_zeros.end(),internal_lcb_removed_count[using_lsd].data());
		remove( the_move, false, internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd], false, new_move_list, new_move_count );
		applyScoreDifference( internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd] );
		double m_score = score();
		undoScoreDifference( internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd] );
		using_lsd--;
		return m_score;
	}

	bool isValid( pair< double, size_t >& the_move )
	{
		using_lsd++;
		std::copy(lsd_zeros.begin(),lsd_zeros.end(),internal_lcb_score_diff[using_lsd].data());
		std::copy(lrc_zeros.begin(),lrc_zeros.end(),internal_lcb_removed_count[using_lsd].data());
		vector< pair< double, size_t > > new_move_list;
		size_t new_move_count;
		bool success = remove( the_move, false, internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd], false, new_move_list, new_move_count );
		using_lsd--;
		return success;
	}

	bool remove( pair< double, size_t >& the_move, vector< pair< double, size_t > >& new_move_list, size_t& new_move_count )
	{
		using_lsd++;
		std::copy(lsd_zeros.begin(),lsd_zeros.end(),internal_lcb_score_diff[using_lsd].data());
		std::copy(lrc_zeros.begin(),lrc_zeros.end(),internal_lcb_removed_count[using_lsd].data());
		bool success = remove( the_move, true, internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd], true, new_move_list, new_move_count );
		if( success )
			applyScoreDifference( internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd] );
		using_lsd--;
/*		cout << "score_diff:\n";
		print2d_matrix( internal_lcb_score_diff, std::cout );
		cout << "\nlcb_removed_count:\n";
		print2d_matrix( internal_lcb_removed_count, std::cout );
		cout << endl;
*/
		return success;
	}

	void applyScoreDifference( boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count )
	{
		size_t nelems = pairwise_lcb_count.num_elements();
		for( size_t elemI = 0; elemI < nelems; elemI++ )
		{
			if( !(lcb_score_diff.data()[elemI] > -1e200 && lcb_score_diff.data()[elemI] < 1e200) )
			{
				genome::breakHere();
				cerr << "Invalid score!!\n";
			}
			pairwise_lcb_count.data()[elemI] -= lcb_removed_count.data()[elemI];
			pairwise_lcb_score.data()[elemI] -= lcb_score_diff.data()[elemI];
			if( !(pairwise_lcb_score.data()[elemI] > -1e200 && pairwise_lcb_score.data()[elemI] < 1e200) )
			{
				genome::breakHere();
				cerr << "Invalid score!!\n";
			}
		}
	}

	void undoScoreDifference( boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count )
	{
		size_t nelems = pairwise_lcb_count.num_elements();
		for( size_t elemI = 0; elemI < nelems; elemI++ )
		{
			if( !(lcb_score_diff.data()[elemI] > -1e200 && lcb_score_diff.data()[elemI] < 1e200) )
			{
				genome::breakHere();
				cerr << "Invalid score!!\n";
			}
			pairwise_lcb_count.data()[elemI] += lcb_removed_count.data()[elemI];
			pairwise_lcb_score.data()[elemI] += lcb_score_diff.data()[elemI];
			if( !(pairwise_lcb_score.data()[elemI] > -1e200 && pairwise_lcb_score.data()[elemI] < 1e200) )
			{
				genome::breakHere();
				cerr << "Invalid score!!\n";
			}
		}
	}

	size_t getMaxNewMoveCount()
	{
		return 20 * seqI_count * seqJ_count;
	}

	/** call to indicate that the given LCB has been removed 
	  * returns false if the move was invalid
	  */
	bool remove( pair< double, size_t >& the_move, bool really_remove, boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count, bool score_new_moves, vector< pair< double, size_t > >& new_move_list, size_t& new_move_count )
	{
		if( score_new_moves && !really_remove )
		{
			cerr << "Error: Incompatible options in the breakpoint scorer!!!\n";
			throw "oh shit!";
		}
		new_move_count = 0;
		// figure out which lcb we're being asked to delete
		size_t moveI = the_move.second;
		size_t move_count = 0;
		size_t move_base = 0;
		size_t seqI = 0;
		size_t seqJ = 0;
		for( seqI = seqI_first; seqI < seqI_last; ++seqI )
		{
			for( seqJ = seqJ_first; seqJ < seqJ_last; ++seqJ )
			{
				all_id_remaps[seqI][seqJ].clear();
				full_impact_list[seqI][seqJ].clear();
			}
		}

		for( seqI = seqI_first; seqI < seqI_last; ++seqI )
		{
			for( seqJ = seqJ_first; seqJ < seqJ_last; ++seqJ )
			{
				move_count += pairwise_adjacencies[seqI][seqJ].size();
				if( move_count > moveI )
					break;
				move_base = move_count;
			}
			if( move_count > moveI )
				break;
		}
		// score deletion of the LCB at (moveI - move_base) from the pairwise alignment of seqI and seqJ
		size_t del_lcb = moveI - move_base;
		if( pairwise_adjacencies[seqI][seqJ][del_lcb].lcb_id != del_lcb && really_remove )
		{
			if( pairwise_adjacencies[seqI][seqJ][del_lcb].lcb_id == LCB_UNASSIGNED )
				cerr << "bad movement, dirty dancing\n";
			return false;	// this is an invalid move -- already deleted or coalesced with another
		}
		if( pairwise_adjacencies[seqI][seqJ][del_lcb].lcb_id != del_lcb )
		{
			return false;	// this is an invalid move -- already deleted
		}
		
		vector< TrackingMatch* > matches(pairwise_adjacencies[seqI][seqJ][del_lcb].matches);
		double cur_score = score();

		if( really_remove )
		{
			deleted_tracking_matches.insert( deleted_tracking_matches.end(), matches.begin(), matches.end() );
		}

		for( size_t i = seqI_first; i < seqI_last; ++i )
		{
			for( size_t j = seqJ_first; j < seqJ_last; ++j )
			{
				lcb_score_diff[i][j] = 0;
				vector< TrackingLCB< TrackingMatch* > >& adjs = pairwise_adjacencies[i][j];
				// create a list of LCBs affected by deletion of this match
				// check whether any of them will have all of their matches removed
				if( lcb_ids.size() < matches.size() )
					lcb_ids.resize( matches.size() + 100 );
				for( size_t mI = 0; mI < matches.size(); ++mI )
					lcb_ids[mI] = tm_lcb_id_array[matches[mI]->match_id][i][j];
				size_t lcb_id_count = matches.size();
				std::sort(lcb_ids.begin(), lcb_ids.begin()+lcb_id_count);
				vector< size_t >::iterator last = std::unique(lcb_ids.begin(), lcb_ids.begin()+lcb_id_count);
				lcb_id_count = last - lcb_ids.begin();
				// delete the last one if its unassigned
				if( lcb_ids[lcb_id_count-1] == LCB_UNASSIGNED )
					lcb_id_count--;

				vector< pair< size_t, vector< TrackingMatch* > > > aff_lcbs(lcb_id_count);
				for( size_t lI = 0; lI < lcb_id_count; ++lI )
					aff_lcbs[lI].first = lcb_ids[lI];

				// organize the deleted matches
				for( size_t mI = 0; mI < matches.size(); ++mI )
				{
					size_t id = tm_lcb_id_array[matches[mI]->match_id][i][j];
					if( id == LCB_UNASSIGNED )
						continue;
					vector< pair< size_t, vector< TrackingMatch* > > >::iterator iter = std::lower_bound( aff_lcbs.begin(), aff_lcbs.end(), make_pair(id,vector< TrackingMatch* >() ) );
					iter->second.push_back( matches[mI] );
				}

				// actually delete the matches and keep a list of LCBs that get completely deleted
				size_t my_del_count = 0;
				for( size_t lI = 0; lI < aff_lcbs.size(); ++lI )
				{
					vector< TrackingMatch* >& cur_matches = adjs[lcb_ids[lI]].matches;
					size_t diff = cur_matches.size() - aff_lcbs[lI].second.size();
					if( diff == 0 )
					{
						if( my_del_count + 1 >= my_del_lcbs.size() )
							my_del_lcbs.resize(2*my_del_lcbs.size());
						my_del_lcbs[my_del_count++] = lcb_ids[lI];
						adjs[lcb_ids[lI]].to_be_deleted = true;
						lcb_score_diff[i][j] += adjs[lcb_ids[lI]].weight;
						if( really_remove )
						{
							adjs[lcb_ids[lI]].weight = 0;
							cur_matches.clear();
						}
						continue;
					}

					// update the LCB score
					double del_score_sum = 0;
					for( size_t mI = 0; mI < aff_lcbs[lI].second.size(); ++mI )
						del_score_sum += tm_score_array[aff_lcbs[lI].second[mI]->match_id][i][j];
					lcb_score_diff[i][j] += del_score_sum;
					full_impact_list[i][j].push_back( aff_lcbs[lI].first );

					if( really_remove )
					{
						adjs[lcb_ids[lI]].weight -= del_score_sum;
					
						// remove the deleted matches
						vector< TrackingMatch* > dest( diff );
						std::set_difference( cur_matches.begin(), cur_matches.end(), 
							aff_lcbs[lI].second.begin(), aff_lcbs[lI].second.end(), dest.begin() );
						swap( dest, cur_matches );
					}
				}

				lcb_removed_count[i][j] = 0;

				// now remove each LCB that needs to be deleted
				std::vector< std::pair< uint, uint > >& fid_remaps = all_id_remaps[i][j];
				std::vector< uint >& fimp_list = full_impact_list[i][j];
				for( size_t delI = 0; delI < my_del_count; ++delI )
				{
					if( adjs[my_del_lcbs[delI]].lcb_id != my_del_lcbs[delI] )
						continue;	// skip this one if it's already been deleted

					std::vector< std::pair< uint, uint > > id_remaps;
					std::vector< uint > impact_list;
					uint removed_count = RemoveLCBandCoalesce( my_del_lcbs[delI], adjs, bogus_scores, id_remaps, impact_list );
					fid_remaps.insert( fid_remaps.end(), id_remaps.begin(), id_remaps.end() );
					fimp_list.insert( fimp_list.end(), impact_list.begin(), impact_list.end() );

					lcb_removed_count[i][j] += removed_count;
					// only do this part if we're really deleting
					if( really_remove )
					{
						// move all matches to the new LCB
						for( size_t rI = 0; rI < id_remaps.size(); ++rI )
						{
							if( id_remaps[rI].second == -1 )
								continue;	// deletion
							vector< TrackingMatch* >& src_matches = adjs[id_remaps[rI].first].matches;
							vector< TrackingMatch* >& dest_matches = adjs[id_remaps[rI].second].matches;
							for( size_t mI = 0; mI < src_matches.size(); ++mI )
								tm_lcb_id_array[src_matches[mI]->match_id][i][j] = id_remaps[rI].second;
							dest_matches.insert( dest_matches.end(), src_matches.begin(), src_matches.end() );
							std::sort( dest_matches.begin(), dest_matches.end() );
							src_matches.clear();
						}
					}
				}
			}
		}

		// will be undone later
		applyScoreDifference( lcb_score_diff, lcb_removed_count );
		double new_score = score();

		if( score_new_moves )
		{
			size_t mbase = 0;
			for( size_t i = seqI_first; i < seqI_last; ++i )
			{
				for( size_t j = seqJ_first; j < seqJ_last; ++j )
				{
					vector< TrackingLCB< TrackingMatch* > >& adjs = pairwise_adjacencies[i][j];
					std::vector< uint >& fimp_list = full_impact_list[i][j];
					sort( fimp_list.begin(), fimp_list.end() );
					vector< uint >::iterator iter = std::unique( fimp_list.begin(), fimp_list.end() );
					fimp_list.erase( iter, fimp_list.end() );
					for( size_t fI = 0; fI < fimp_list.size(); fI++ )
					{
						if( adjs[fimp_list[fI]].lcb_id != fimp_list[fI] )
						{
							new_move_list[new_move_count++] = make_pair( -(std::numeric_limits<double>::max)(), mbase + fimp_list[fI] );
							continue;	// this one got trashed
						}
						// score removal of this block
						pair< double, size_t > p( 0, mbase + fimp_list[fI] );
						double scorediff = (*this)(p) - new_score;
						p.first = scorediff;
						new_move_list[new_move_count++] = p;
					}
					mbase += adjs.size();
				}
			}
		}


		// if we're not really removing, undo all the removals
		if( !really_remove )
			for( size_t i = seqI_first; i < seqI_last; ++i )
				for( size_t j = seqJ_first; j < seqJ_last; ++j )
					undoLcbRemoval( pairwise_adjacencies[i][j], all_id_remaps[i][j] );

		undoScoreDifference( lcb_score_diff, lcb_removed_count );

		// if the change in score doesn't match then this is an invalid move!!
		// allow for some numerical instability
		bool valid = true;
		if( new_score - cur_score < the_move.first - 0.00001 ||
			new_score - cur_score  > the_move.first + 0.00001 )
			valid = false;

		return valid;
	}

	vector< TrackingMatch* > getResults() 
	{
		std::sort(deleted_tracking_matches.begin(), deleted_tracking_matches.end());
		vector< TrackingMatch* > result_matches(tracking_matches.size()-deleted_tracking_matches.size());
		std::set_difference( tracking_matches.begin(), tracking_matches.end(), deleted_tracking_matches.begin(), deleted_tracking_matches.end(), result_matches.begin() );
		return result_matches;
//		return tracking_matches;
	}

	bool validate()
	{
		vector< TrackingMatch* > trams = getResults();	// need to apply any deletions...
		bool success = true;	// be optimistic!
		// make sure all the tracking matches point to the right LCBs
		for( size_t tmI = 0; tmI < trams.size(); tmI++ )
		{
			TrackingMatch* tm = trams[tmI];
			for( size_t i = 0; i < tm_lcb_id_array.shape()[1]; ++i )
				for( size_t j = 0; j < tm_lcb_id_array.shape()[2]; ++j )
				{
					// skip this match if it's not defined
					if( tm->node_match->LeftEnd(n1_des[i]) == NO_MATCH ||
						tm->node_match->LeftEnd(n2_des[j]) == NO_MATCH ||
						tm_lcb_id_array[tm->match_id][i][j] == LCB_UNASSIGNED)
						continue;
					// find the tracking match in this LCB
					size_t id = tm_lcb_id_array[tm->match_id][i][j];
					vector< TrackingMatch* >& matches = pairwise_adjacencies[i][j][id].matches;
					vector< TrackingMatch* >::iterator iter = std::lower_bound( matches.begin(), matches.end(), tm );
					if( iter == matches.end() || *iter != tm )
					{
						cerr << "Missing match!!\n";
						cerr << "lcb_id: " << id << endl;
						cerr << "match: " << tm << endl;
						genome::breakHere();
						success = false;
					}
				}
		}
		// make sure all the LCBs point to valid tracking matches
		for( size_t i = 0; i < pairwise_adjacencies.shape()[0]; ++i )
			for( size_t j = 0; j < pairwise_adjacencies.shape()[1]; ++j )
			{
				vector< TrackingLCB< TrackingMatch* > >& adjs = pairwise_adjacencies[i][j];
				for( size_t lcbI = 0; lcbI < adjs.size(); lcbI++ )
				{
					for( size_t mI = 0; mI < adjs[lcbI].matches.size(); ++mI )
					{
						vector< TrackingMatch* >::iterator iter = std::lower_bound( trams.begin(), trams.end(), adjs[lcbI].matches[mI] );
						if( *iter != adjs[lcbI].matches[mI] )
						{
							cerr << "Missing match:  in adjacencies but not tracking_matches!!\n";
							cerr << "lcb_id: " << tm_lcb_id_array[adjs[lcbI].matches[mI]->match_id][i][j] << endl;
							genome::breakHere();
							success = false;
						}
					}
				}
			}

		// make sure that the number of breakpoints matches up with what tracking_matches suggests
		vector< TrackingMatch* > final = trams;
		// convert back to an LCB list
		vector< AbstractMatch* > new_matches(final.size());
		for( size_t mI = 0; mI < final.size(); ++mI )
			new_matches[mI] = final[mI]->original_match;

		vector< gnSeqI > breakpoints;
		IdentifyBreakpoints( new_matches, breakpoints );
		vector< vector< AbstractMatch* > > LCB_list;
		IdentifyBreakpoints( new_matches, breakpoints );
		ComputeLCBs_v2( new_matches, breakpoints, LCB_list );
		cout << "breakpoints.size(): " << breakpoints.size() << "\tpairwise_lcb_count[0][0]: " << pairwise_lcb_count[0][0] << endl;
		if( breakpoints.size() != pairwise_lcb_count[0][0] )
			success = false;
		size_t adjI = 0;
		vector< TrackingLCB< TrackingMatch* > >& adjs = pairwise_adjacencies[0][0];
		for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
		{
			// make sure each LCB exists...
			while( adjI != -1 && adjI != adjs[adjI].lcb_id )
				adjI++;

			// compare matches...
			vector< AbstractMatch* > ms(adjs[adjI].matches.size()+LCB_list[lcbI].size(), (AbstractMatch*)NULL);
			std::sort( LCB_list[lcbI].begin(), LCB_list[lcbI].end() );
			vector< AbstractMatch* > asdf(adjs[adjI].matches.size());
			for( size_t mI = 0; mI < adjs[adjI].matches.size(); ++mI )
				asdf[mI] = adjs[adjI].matches[mI]->original_match;
			std::sort( asdf.begin(), asdf.end() );
			std::set_symmetric_difference( LCB_list[lcbI].begin(), LCB_list[lcbI].end(), asdf.begin(), asdf.end(), ms.begin() );
			// this should throw a fit if the sets aren't equal.
			if( ms[0] != NULL )
			{
				cerr << "In adjacencies:\n";
				for( size_t asdfI = 0; asdfI < asdf.size(); asdfI++ )
				{
					printMatch(asdf[asdfI], cerr);
					cerr << endl;
				}
				cerr << "\nIn LCB_list:\n";
				for( size_t mI = 0; mI < LCB_list[lcbI].size(); mI++ )
				{
					printMatch(LCB_list[lcbI][mI], cerr);
					cerr << endl;
				}
				cerr << "\nAll matches ssc1\n";
				SingleStartComparator<AbstractMatch> ssc1(1);
				std::sort(new_matches.begin(), new_matches.end(), ssc1);
				for( size_t mI = 0; mI < new_matches.size(); mI++ )
				{
					printMatch(new_matches[mI], cerr);
					cerr << endl;
				}

				cerr << "\nAll matches ssc0\n";
				SingleStartComparator<AbstractMatch> ssc0(0);
				std::sort(new_matches.begin(), new_matches.end(), ssc0);
				for( size_t mI = 0; mI < new_matches.size(); mI++ )
				{
					printMatch(new_matches[mI], cerr);
					cerr << endl;
				}
				genome::breakHere();
			}
			adjI++;
		}

		return success;
	}

protected:
	double bp_penalty;
	boost::multi_array<double,2> bp_weights;
	boost::multi_array<double,2> conservation_weights;
	vector< TrackingMatch* > tracking_matches;
	PairwiseLCBMatrix pairwise_adjacencies;
	std::vector<node_id_t> n1_des;
	std::vector<node_id_t> n2_des;

	boost::multi_array< size_t, 2 > pairwise_lcb_count;
	boost::multi_array< double, 2 > pairwise_lcb_score;

	vector< TrackingMatch* > deleted_tracking_matches;

private:
	// avoid continuous size lookup
	const size_t seqI_count;
	const size_t seqJ_count;

	// variables used during score computation
	boost::multi_array< std::vector< std::pair< uint, uint > >, 2 > all_id_remaps;
	boost::multi_array< std::vector< uint >, 2 > full_impact_list;
	boost::multi_array< double, 2 > internal_lcb_score_diff[3];
	boost::multi_array< size_t, 2 > internal_lcb_removed_count[3];
	int using_lsd;
	std::vector< double > lsd_zeros;
	std::vector< size_t > lrc_zeros;
	vector< double > bogus_scores;
	vector< size_t > my_del_lcbs;
	vector< size_t > lcb_ids;

	boost::multi_array< double, 3 >& tm_score_array;
	boost::multi_array< size_t, 3 >& tm_lcb_id_array;

	// limit to a range of sequences
	const size_t seqI_first;
	const size_t seqJ_first;
	const size_t seqI_last;
	const size_t seqJ_last;
};


class MoveScoreHeapComparator
{
public:
	bool operator()( const pair< double, size_t >& a, const pair< double, size_t >& b )
	{
		return a.first < b.first;	// want to order by > instead of <
	}
};

/** finds the best anchoring, returns the anchoring score */
template< class SearchScorer >
double greedySearch( SearchScorer& spbs )
{
	double prev_score = spbs.score();
	uint report_frequency = 10;
	uint moves_made = 0;
	if( debug_aligner )
		spbs.validate();
	size_t move_count = spbs.getMoveCount();
	vector< double > current_moves( spbs.getMoveCount() );
	// use double the size for the move heap to avoid an almost instant reallocation
	// when a new move gets pushed onto the heap
	size_t heap_end = spbs.getMoveCount();
	vector< pair< double, size_t > > move_heap( spbs.getMoveCount() * 2 );
	vector< pair< double, size_t > > new_moves( spbs.getMaxNewMoveCount() + 10 );
	for( size_t moveI = 0; moveI < move_count; ++moveI )
	{
		pair< double, size_t > p( 0, moveI );
		double scorediff = spbs(p) - prev_score;
		p.first = scorediff;
		move_heap[moveI] = p;
		current_moves[moveI] = p.first;
	}

	if( debug_aligner )
		spbs.validate();
	// make a heap of moves ordered by score
	// repeatedly:
	// 1) pop the highest scoring move off the heap
	// 2) attempt to apply the move
	// 3) add any new moves to the heap
	// 4) stop when the highest scoring move no longer increases the score
	MoveScoreHeapComparator mshc;
	std::make_heap( move_heap.begin(), move_heap.begin() + heap_end, mshc );
	double successful = 0;
	double invalids = 0;
	int progress = 0;
	int prev_progress = -1;
	while( heap_end > 0 )
	{
		std::pop_heap( move_heap.begin(), move_heap.begin()+heap_end, mshc );
		pair< double, size_t > best_move = move_heap[--heap_end];
		if( best_move.first < 0 )
			break;	// can't improve score

		if( best_move.first != current_moves[best_move.second] )
			continue;

		if( !spbs.isValid(best_move) )
		{
			invalids++;
			continue;
		}

		size_t new_move_count = 0;
		bool success = spbs.remove(best_move, new_moves, new_move_count);
		if( !success )
		{
			cerr << "numerical instability?  need to investigate this...\n";
//			genome::breakHere();
			invalids++;
			continue;
		}

		successful++;
		if( debug_aligner )
			spbs.validate();

		current_moves[ best_move.second ] = -(std::numeric_limits<double>::max)();
		for( size_t newI = 0; newI < new_move_count; newI++ )
			current_moves[ new_moves[newI].second ] = new_moves[newI].first;

		for( size_t newI = 0; newI < new_move_count; newI++ )
		{
			if( heap_end < move_heap.size() )
			{
				heap_end++;
				move_heap[heap_end-1] = new_moves[newI];
				std::push_heap( move_heap.begin(), move_heap.begin()+heap_end, mshc );
			}else{
				// just push the rest on all at once
				move_heap.resize( (std::min)((size_t)(heap_end * 1.6), heap_end + new_move_count) );
				std::copy( new_moves.begin() + newI, new_moves.begin() + new_move_count, move_heap.begin()+heap_end );
				for( size_t newdI = 0; newdI < new_move_count-newI; newdI++ )
					std::push_heap( move_heap.begin(), move_heap.begin()+heap_end+newdI+1, mshc );
				heap_end = move_heap.size();
				break;
			}
		}

		moves_made++;
		prev_progress = progress;
		progress = (100 * moves_made) / move_count;
		printProgress( prev_progress, progress, cout );
//		if( moves_made % report_frequency == 0 )
//			cout << "move: " << moves_made << " alignment score " << cur_score << " success ratio " << successful / invalids << endl;
	}

	return spbs.score();
}


/**
 * A breakpoint scorer that applies a fixed penalty for each breakpoint that exists in a set of
 * two or more sequences 
 */
class SimpleBreakpointScorer
{
public:
//	SimpleBreakpointScorer(){};
	SimpleBreakpointScorer( std::vector< LCB >& adjacencies, double breakpoint_penalty ) : 
	  adjs( adjacencies ),
	  bp_penalty( breakpoint_penalty )
	  {
		scores = std::vector< double >(adjs.size(), 0);
		total_weight = 0;
		bp_count = adjs.size();
		for( size_t lcbI = 0; lcbI < adjs.size(); lcbI++ )
			total_weight += adjs[lcbI].weight;
	  }

	size_t getMoveCount() 
	{
		return adjs.size();
	}

	double score()
	{
		double bp_score = (double)bp_count * bp_penalty;
		return total_weight - bp_score;
	}

	bool isValid( size_t lcbI, double move_score )
	{
		if( adjs[lcbI].lcb_id != lcbI )
			return false;
		return (*this)(lcbI) == move_score;
	}

	/** return the relative change in score if lcbI were to be removed */
	double operator()( size_t lcbI )
	{
		double cur_score = score();
		std::vector< std::pair< uint, uint > > id_remaps;
		std::vector< uint > impact_list;
		uint bp_removed = RemoveLCBandCoalesce( lcbI, adjs, scores, id_remaps, impact_list );
		undoLcbRemoval( adjs, id_remaps );
		double bp_score = (double)(bp_count - bp_removed) * bp_penalty;
		double move_score = total_weight - adjs[lcbI].weight - bp_score;
		return move_score - cur_score;
	}

	/** call to indicate that the given LCB has been removed */
	void remove( uint lcbI, vector< pair< double, size_t > >& new_moves )
	{
		std::vector< std::pair< uint, uint > > id_remaps;
		std::vector< uint > impact_list;
		uint bp_removed = RemoveLCBandCoalesce( lcbI, adjs, scores, id_remaps, impact_list );
		total_weight -= adjs[lcbI].weight;
		bp_count -= bp_removed;
		for( size_t impI = 0; impI < impact_list.size(); impI++ )
		{
			if( adjs[impact_list[impI]].lcb_id != impact_list[impI] )
				continue;
			double scorediff = (*this)(impact_list[impI]);
			new_moves.push_back(make_pair(scorediff, impact_list[impI]));
		}
	}

private:
	std::vector< LCB > adjs;
	double bp_penalty;
	std::vector< double > scores;
	double total_weight;
	size_t bp_count;
};


class GreedyRemovalScorer
{
public:
//	GreedyRemovalScorer(){};
	GreedyRemovalScorer( std::vector< LCB >& adjacencies, double minimum_weight ) : 
	  adjs( adjacencies ),
	  min_weight( minimum_weight )
	  {
		scores = std::vector< double >(adjs.size(), 0);
		total_weight = 0;
		for( size_t lcbI = 0; lcbI < adjs.size(); lcbI++ )
			total_weight += adjs[lcbI].weight - min_weight;
	  }

	size_t getMoveCount() 
	{
		return adjs.size();
	}

	double score()
	{
		return total_weight;
	}

	bool isValid( size_t lcbI, double move_score )
	{
		if( adjs[lcbI].lcb_id != lcbI )
			return false;
		return (*this)(lcbI) == move_score;
	}

	/** return the relative change in score if lcbI were to be removed */
	double operator()( size_t lcbI )
	{
		return -(adjs[lcbI].weight-min_weight);
	}

	/** call to indicate that the given LCB has been removed */
	void remove( uint lcbI, vector< pair< double, size_t > >& new_moves )
	{
		std::vector< std::pair< uint, uint > > id_remaps;
		std::vector< uint > impact_list;
		uint bp_removed = RemoveLCBandCoalesce( lcbI, adjs, scores, id_remaps, impact_list );
		total_weight -= (adjs[lcbI].weight-min_weight);
		for( size_t impI = 0; impI < impact_list.size(); impI++ )
		{
			if( adjs[impact_list[impI]].lcb_id != impact_list[impI] )
				continue;
			double scorediff = (*this)(impact_list[impI]);
			new_moves.push_back(make_pair(scorediff, impact_list[impI]));
		}
	}

private:
	std::vector< LCB > adjs;
	double min_weight;
	std::vector< double > scores;
	double total_weight;
};

template< class BreakpointScorerType >
int64 greedyBreakpointElimination_v4( vector< LCB >& adjacencies, vector< double >& scores, BreakpointScorerType& bp_scorer, ostream* status_out, size_t g1_tag = 0, size_t g2_tag = 0 ){
	// repeatedly remove the low weight LCBs until the minimum weight criteria is satisfied
	uint lcb_count = adjacencies.size();
	double total_initial_lcb_weight = 0;
	for( size_t wI = 0; wI < scores.size(); wI++ )
		total_initial_lcb_weight += scores[wI];
	double total_current_lcb_weight = total_initial_lcb_weight;

	if( adjacencies.size() == 0 )
		return 0;	// nothing can be done
	uint seq_count = adjacencies[0].left_end.size();
	
	double prev_score = bp_scorer.score();
	uint report_frequency = 10;
	uint moves_made = 0;

	size_t move_count = bp_scorer.getMoveCount();
	vector< pair< double, size_t > > move_heap( move_count * 2 );
	size_t heap_end = move_count;
	for( size_t moveI = 0; moveI < move_count; ++moveI )
	{
		move_heap[moveI].first = bp_scorer(moveI);
		move_heap[moveI].second = moveI;
	}

#ifdef LCB_WEIGHT_LOSS_PLOT
	vector< double >::iterator min_iter = std::min_element(scores.begin(), scores.end());
	double mins = *min_iter;
	if( status_out != NULL )
	{
		(*status_out) << g1_tag << '\t' << g2_tag << '\t' << lcb_count << '\t' << 1 - (total_current_lcb_weight / total_initial_lcb_weight) << '\t' << mins << endl;
	}
#endif

	// make a heap of moves ordered by score
	// repeatedly:
	// 1) pop the highest scoring move off the heap
	// 2) attempt to apply the move
	// 3) add any new moves to the heap
	// 4) stop when the highest scoring move no longer increases the score
	MoveScoreHeapComparator mshc;
	std::make_heap( move_heap.begin(), move_heap.end(), mshc );
	while( heap_end > 0 )
	{
		std::pop_heap( move_heap.begin(), move_heap.begin()+heap_end, mshc );
		heap_end--;
		pair< double, size_t > best_move = move_heap[ heap_end ];
#ifdef LCB_WEIGHT_LOSS_PLOT
		if( total_current_lcb_weight == scores[best_move.second] )
			break;	// don't remove the last LCB
#else
		if( (best_move.first < 0 ) ||
			total_current_lcb_weight == scores[best_move.second] )
			break;	// can't improve score
#endif

		vector< pair< double, size_t > > new_moves;
		bool success = bp_scorer.isValid(best_move.second, best_move.first);
		if( !success )
			continue;
		bp_scorer.remove(best_move.second, new_moves);

		
		for( size_t newI = 0; newI < new_moves.size(); newI++ )
		{
			if( heap_end < move_heap.size() )
			{
				heap_end++;
				move_heap[heap_end-1] = new_moves[newI];
				std::push_heap( move_heap.begin(), move_heap.begin()+heap_end, mshc );
			}else{
				// just push the rest on all at once
				size_t prev_size = move_heap.size();
				move_heap.insert( move_heap.end(), new_moves.begin()+newI, new_moves.end() );
				for( size_t newdI = 0; newdI < new_moves.size()-newI; newdI++ )
					std::push_heap( move_heap.begin(), move_heap.begin()+prev_size+newdI+1, mshc );
				heap_end = move_heap.size();
				break;
			}
		}

		total_current_lcb_weight -= scores[best_move.second];
		std::vector< std::pair< uint, uint > > id_remaps;
		std::vector< uint > impact_list;
		lcb_count -= RemoveLCBandCoalesce( best_move.second, adjacencies, scores, id_remaps, impact_list );
#ifdef LCB_WEIGHT_LOSS_PLOT
		mins = scores[best_move.second];
		if( status_out != NULL )
		{
			(*status_out) << g1_tag << '\t' << g2_tag << '\t' << lcb_count << '\t' << 1 - (total_current_lcb_weight / total_initial_lcb_weight) << '\t' << mins << endl;
		}
#endif
		double cur_score = bp_scorer.score();
		prev_score = cur_score;
		moves_made++;
#ifndef LCB_WEIGHT_LOSS_PLOT
		if( status_out != NULL && moves_made % report_frequency == 0 )
			(*status_out) << "move: " << moves_made << " alignment score " << cur_score << endl;
#endif
	}

	return 0;
}

/*
 * A progressive alignment algorithm for genomes with rearrangements.
 * Start simple, add complexity later.
 * (1) Compute gene content distance matrix and phylogeny
 * (2) Align the two most similar genomes
 *     - Do ordinary Mauve alignment but after the first part of LCB
 *       extension, do an anchored extension off the ends of each LCB
 *       - If two anchored extensions overlap, give the overlapping region
 *         to the higher scoring alignment -- calculate up to 100 bp of overlapping
 *         alignment...
 *     - recursively align and complete an alignment
 * (3) Add a third genome, calculate avg. distance between each genome
 *     and the current pair, select the genome with min avg. distance.
 *     - Do pairwise anchoring between the new genome and the aligned genomes
 *       - discard any overlapping anchors that are /inconsistent/ with
 *         the existing alignment
 * (4) Extend LCBs in the third genome -- find matches outside existing LCBs
 *     Then do the windowed extension -- find matches within a 10k (or some
 *     other fixed size) window off the end of each LCB
 * (5) Repeat steps (3) and (4) for each additional genome.
 * (6) Feed windows of the alignment to MUSCLE for iterative refinement
 */

ProgressiveAligner::ProgressiveAligner( uint seq_count ) :
Aligner( seq_count ),
breakpoint_penalty( -1 ),
debug(false),
refine(true),
scoring_scheme(ExtantSumOfPairsScoring),
use_weight_scaling(true),
conservation_dist_scale(.5),
bp_dist_scale(.5),
max_gapped_alignment_length(20000),
bp_dist_estimate_score(-1)
{
	gapped_alignment = true;
	max_window_size = max_gapped_alignment_length;
}

void ProgressiveAligner::SetMaxGappedAlignmentLength( size_t len )
{ 
	max_gapped_alignment_length = len; 
	max_window_size = max_gapped_alignment_length;
}

/** determine which extant sequences have been aligned at a given node */
void ProgressiveAligner::getAlignedChildren( node_id_t node, vector< node_id_t >& descendants )
{
	// do a depth first search along edges that have been aligned
	stack< node_id_t > node_stack;
	node_stack.push( node );
	vector< bool > visited( alignment_tree.size(), false );
	descendants.clear();
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		if(progress_msgs) cout << "Evaluating aligned nodes linked to node " << cur_node << endl;
		node_stack.pop();
		visited[cur_node] = true;
		for( uint childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			node_id_t child_id = alignment_tree[cur_node].children[childI];
			if( alignment_tree[cur_node].children_aligned[childI] && !visited[child_id])
				node_stack.push( child_id );
		}
		if( alignment_tree[ cur_node ].sequence != NULL )
			descendants.push_back( cur_node );
	}
}


/** determine which extant sequences have been aligned at a given node */
void ProgressiveAligner::getPath( node_id_t first_n, node_id_t last_n, vector< node_id_t >& path )
{
	// do a depth first search along edges that have been aligned
	stack< node_id_t > node_stack;
	node_stack.push( last_n );
	vector< bool > visited( alignment_tree.size(), false );
	while( node_stack.top() != first_n )
	{
		node_id_t cur_node = node_stack.top();
		size_t pre_size = node_stack.size();
		visited[cur_node] = true;
		for( uint childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			node_id_t child_id = alignment_tree[cur_node].children[childI];
			if(!visited[child_id])
			{
				node_stack.push( child_id );
				break;
			}
		}
		if( pre_size != node_stack.size() )
			continue;
		for( uint parentI = 0; parentI < alignment_tree[cur_node].parents.size(); parentI++ )
		{
			node_id_t parent_id = alignment_tree[cur_node].parents[parentI];
			if(!visited[parent_id])
			{
				node_stack.push( parent_id );
				break;
			}
		}
		if( pre_size != node_stack.size() )
			continue;
		node_stack.pop();	// didn't make any progress
	}
	path = vector< node_id_t >( node_stack.size() );
	for( size_t pI = 0; pI < path.size(); pI++ )
	{
		path[pI] = node_stack.top();
		node_stack.pop();
	}
}



void SuperInterval::CropLeft( gnSeqI amount )
{
	reference_iv.CropStart(amount);

	left_end += amount;
	length -= amount;

	if(debug_aligner)
		ValidateSelf();
}

void SuperInterval::CropRight( gnSeqI amount )
{
	reference_iv.CropEnd(amount);
	length -= amount;

	if(debug_aligner)
		ValidateSelf();
}

void SuperInterval::ValidateSelf() const
{
	vector< bitset_t > aln_mat;
	reference_iv.GetAlignment(aln_mat);
	if( aln_mat[0].size() != reference_iv.AlignmentLength() )
	{
		breakHere();
		cerr << "trouble! aln_mat[0].size() is: " << aln_mat[0].size() << " while reference_iv.AlignmentLength() is: " << reference_iv.AlignmentLength() << endl;
		cerr << "mult: " << reference_iv.Multiplicity() << endl;
		cerr << "matches.size(): " << reference_iv.GetMatches().size() << endl;
	}
	for( size_t i = 0; i < aln_mat.size(); i++ )
	{
		gnSeqI lenny = 0;
		for( size_t j = 0; j < aln_mat[i].size(); j++ )
			if( aln_mat[i][j] )
				lenny++;
		if( lenny != reference_iv.Length(i) )
		{
			cerr << "krudunkle, ref_iv.Length(" << i << "): " << reference_iv.Length(i) << "\n";
			cerr << "should be: " << lenny << endl;
			breakHere();
		}
	}
	if( reference_iv.LeftEnd(0) != NO_MATCH && reference_iv.Length(0) == 0 )
	{
		cerr << "brozooka\n";
		breakHere();
	}
	if( reference_iv.LeftEnd(1) != NO_MATCH && reference_iv.Length(1) == 0 )
	{
		cerr << "brokazooka\n";
		breakHere();
	}

	if( Length() != reference_iv.AlignmentLength() )
	{
		breakHere();
		cerr << "crapola\n";
	}
}

template< class MatchType = mems::AbstractMatch >
class GenericMatchSeqManipulator
{
public:
	GenericMatchSeqManipulator( uint seq ) : m_seq(seq) {}
	gnSeqI LeftEnd(MatchType*& m) const{ return m->LeftEnd(m_seq); }
	gnSeqI Length(MatchType*& m) const{ return m->Length(m_seq); }
	void CropLeft(MatchType*& m, gnSeqI amount ) const{ m->CropLeft(amount, m_seq); }
	void CropRight(MatchType*& m, gnSeqI amount ) const{ m->CropRight(amount, m_seq); }
	template< typename ContainerType >
	void AddCopy(ContainerType& c, MatchType*& m) const{ c.push_back( m->Copy() ); }
private:
	uint m_seq;
};

typedef GenericMatchSeqManipulator<> AbstractMatchSeqManipulator;

class SuperIntervalManipulator
{
public:
	gnSeqI LeftEnd(const SuperInterval& siv) const{ return siv.LeftEnd(); }
	gnSeqI Length(const SuperInterval& siv) const{ return siv.Length(); }
	void CropLeft( SuperInterval& siv, gnSeqI amount ) const{ siv.CropLeft( amount );}
	void CropRight( SuperInterval& siv, gnSeqI amount ) const{ siv.CropRight( amount );}
	template< typename ContainerType >
	void AddCopy(ContainerType& c, const SuperInterval& siv) const{ c.push_back( siv ); }
};

// iv_list is a container class that contains pointers to intervals or 
// matches of some sort
// precondition: both bp_list and intervals *must* be sorted
template< class T, class Maniplator >
void applyBreakpoints( vector< gnSeqI >& bp_list, vector<T>& iv_list, Maniplator& manip )
{

	size_t iv_count = iv_list.size();
	size_t bpI = 0;
	size_t ivI = 0;
	while( ivI < iv_count && bpI < bp_list.size() )
	{
		if( manip.LeftEnd(iv_list[ivI]) == NO_MATCH )
		{
			++ivI;
			continue;	// undefined in seqI, so no breakpoint here
		}
		//  -(ivI)----
		//  -------|--
		if( manip.LeftEnd(iv_list[ivI]) + manip.Length(iv_list[ivI]) <= bp_list[bpI] )
		{
			++ivI;
			continue;
		}
		//  -----(ivI)-
		//  --|--------
		if( bp_list[bpI] <= manip.LeftEnd(iv_list[ivI]) )
		{
			++bpI;
			continue;
		}

		// if split_at isn't 0 then we need to split cur_iv
		// put the left side in the new list and crop cur_iv
		gnSeqI crop_amt = bp_list[bpI] - manip.LeftEnd(iv_list[ivI]);
		manip.AddCopy( iv_list, iv_list[ivI] );
		T& left_iv = iv_list.back();
		if( progress_msgs ) cout << "crop_amt is: " << crop_amt << endl;

		manip.CropLeft( iv_list[ivI], crop_amt );
		manip.CropRight( left_iv, manip.Length(left_iv)-crop_amt );
		// restore ordering
		size_t nextI = ivI + 1;
		while( nextI < iv_count && manip.LeftEnd( iv_list[nextI-1] ) > manip.LeftEnd( iv_list[nextI] ) )
		{
			swap( iv_list[nextI-1], iv_list[nextI] );
			nextI++;
		}

// assume that crop works correctly and that it's okay to pass matches with NO_MATCH		
/**/
		if( manip.Length( iv_list[ivI] ) == 0 )
		{
			cerr << "Big fat generic zero 1\n";
			breakHere();
		}
		if( manip.Length( left_iv ) == 0 )
		{
			cerr << "Big fat generic zero 2\n";
			breakHere();
		}
		if( manip.LeftEnd( iv_list[ivI] ) == 0 )
		{
			cerr << "uh oh\n";
			breakHere();
		}
		if( manip.LeftEnd( left_iv ) == 0 )
		{
			cerr << "uh oh 2\n";
			breakHere();
		}
/**/
	}
}

template<class MatchType>
void ProgressiveAligner::propagateDescendantBreakpoints( node_id_t node1, uint seqI, std::vector<MatchType*>& iv_list )
{
	SSC<MatchType> ilc(seqI);
	sort( iv_list.begin(), iv_list.end(), ilc );
	vector< SuperInterval >& ord = alignment_tree[ node1 ].ordering;
	vector<gnSeqI> bp_list;
	for( size_t sI = 0; sI < ord.size(); sI++ )
		bp_list.push_back( ord[sI].LeftEnd() );

	GenericMatchSeqManipulator<MatchType> ism( seqI );
	applyBreakpoints( bp_list, iv_list, ism );
}

// T should be a pointer type
template<class T, class Manipulator>
void applyAncestralBreakpoints( const vector< SuperInterval >& siv_list, vector<T>& ord, uint seqI, Manipulator& m )
{
	// make bp list
	vector<gnSeqI> bp_list(siv_list.size()*2, 0);
	size_t cur = 0;
	for( size_t i = 0; i < siv_list.size(); i++ )
	{
		if( siv_list[i].reference_iv.Start(seqI) == NO_MATCH )
			continue;
		bp_list[cur++] = siv_list[i].reference_iv.LeftEnd(seqI);
		bp_list[cur++] = siv_list[i].reference_iv.LeftEnd(seqI) + siv_list[i].reference_iv.Length(seqI);
	}
	bp_list.resize(cur);
	// sort the breakpoints and apply...
	sort( bp_list.begin(), bp_list.end() );
	applyBreakpoints( bp_list, ord, m );
}


// assuming breakpoints have been propagated in both directions
// there should now be a 1-to-1 correspondence between superintervals
// in the ancestor and descendants.
void ProgressiveAligner::linkSuperIntervals( node_id_t node1, uint seqI, node_id_t ancestor )
{
	// TODO: speed this up by implementing O(N) instead of O(N^2)
	vector<SuperInterval>& a_ord = alignment_tree[ancestor].ordering;
	vector<SuperInterval>& c_ord = alignment_tree[node1].ordering;
	// initialize all linkages to nothing
	for( size_t aI = 0; aI < a_ord.size(); aI++ )
		if( seqI == 0 )
			a_ord[aI].c1_siv = (std::numeric_limits<size_t>::max)();
		else
			a_ord[aI].c2_siv = (std::numeric_limits<size_t>::max)();
	for( size_t cI = 0; cI < c_ord.size(); cI++ )
		c_ord[cI].parent_siv = (std::numeric_limits<size_t>::max)();

	for( size_t aI = 0; aI < a_ord.size(); aI++ )
	{
		if( a_ord[aI].reference_iv.LeftEnd(seqI) == NO_MATCH )
			continue;
		size_t cI = 0;
		for( ; cI < c_ord.size(); cI++ )
		{
			if( absolut(a_ord[aI].reference_iv.Start(seqI)) != c_ord[cI].LeftEnd() )
				continue;
			if( a_ord[aI].reference_iv.Length(seqI) != c_ord[cI].Length() )
			{
				breakHere();
				cerr << "mapping length mismatch\n";
				cerr << "ancestor: " << ancestor << "\t node1: " << node1 << endl;
				cerr << "a_ord[" << aI << "].reference_iv.Length(" << seqI << "): " << a_ord[aI].reference_iv.Length(seqI) << endl;
				cerr << "a_ord[" << aI << "].reference_iv.LeftEnd(" << seqI << "): " << a_ord[aI].reference_iv.LeftEnd(seqI) << endl;
				cerr << "c_ord[" << cI << "].Length(): " << c_ord[cI].Length() << endl;
				cerr << "c_ord[" << cI << "].LeftEnd(): " << c_ord[cI].LeftEnd() << endl;
				cerr << "";
				cerr << "";
			}
			// link these
			if( seqI == 0 )
				a_ord[aI].c1_siv = cI;
			else
				a_ord[aI].c2_siv = cI;
			c_ord[cI].parent_siv = aI;
			break;
		}
		if( cI == c_ord.size() )
		{
			breakHere();
			cerr << "error no mapping\n";
		}
	}
}


void ProgressiveAligner::translateGappedCoordinates( vector<AbstractMatch*>& ml, uint seqI, node_id_t extant, node_id_t ancestor )
{
	// determine the path that must be traversed
	vector< node_id_t > trans_path;
	getPath( extant, ancestor, trans_path );

	// set seqI to forward orientation 
	for( size_t mI = 0; mI < ml.size(); mI++ )
		if( ml[mI]->Orientation(seqI) == AbstractMatch::reverse )
			ml[mI]->Invert();

	// for each node on the path, construct a complete coordinate translation
	for( size_t nI = 1; nI < trans_path.size(); nI++ )
	{
		// first sort matches on start pos and make them all forward oriented
		// then split them on superinterval boundaries and assign each to a superinterval
		// then convert each match's coordinates to be superinterval-local
		// then apply the coordinate translation with transposeCoordinates
		// then shift each match's coordinates to the global ancestral coordinate space
		SSC<AbstractMatch> ssc(seqI);
		sort(ml.begin(), ml.end(), ssc);

		// split on superinterval boundaries
		vector< SuperInterval >& siv_list = alignment_tree[trans_path[nI]].ordering;
		vector< vector< AbstractMatch* > > siv_matches = vector< vector< AbstractMatch* > >(siv_list.size());
		size_t cur_child = 0;
		if( alignment_tree[trans_path[nI]].children[0] == trans_path[nI-1] )
			cur_child = 0;
		else if( alignment_tree[trans_path[nI]].children[1] == trans_path[nI-1] )
			cur_child = 1;
		else 
		{
			breakHere();
			cerr << "forest fire\n";
		}

		AbstractMatchSeqManipulator amsm( seqI );
		applyAncestralBreakpoints(siv_list, ml, cur_child, amsm );
		
		// sort matches again because new ones were added at the end
		sort(ml.begin(), ml.end(), ssc);

		// assign each match to a siv, and convert coords to siv-local
		for( size_t mI = 0; mI < ml.size(); mI++ )
		{
			if( ml[mI]->LeftEnd(seqI) == 0 )
			{
				breakHere();
				cerr << "fefefe";
			}
			size_t sivI = 0;
			for( ; sivI < siv_list.size(); sivI++ )
			{
				if( siv_list[sivI].reference_iv.LeftEnd(cur_child) == NO_MATCH )
					continue;
				if( ml[mI]->LeftEnd(seqI) >= siv_list[sivI].reference_iv.LeftEnd(cur_child) &&
					ml[mI]->LeftEnd(seqI) < siv_list[sivI].reference_iv.LeftEnd(cur_child) + siv_list[sivI].reference_iv.Length(cur_child) )
					break;
			}
			if( sivI == siv_list.size() )
			{
				cerr << "nI is: "<< nI << endl;
				cerr << "trans_path: ";
				for( size_t ttI = 0; ttI < trans_path.size(); ttI++ )
					cerr << "  " << trans_path[ttI];
				cerr << endl;
				cerr << "problem seq: " << seqI << std::endl;
				cerr << "ml[" << mI << "]->Start(0) == " << ml[mI]->Start(0) << endl;
				cerr << "ml[" << mI << "]->Length(0) == " << ml[mI]->Length(1) << endl;
				cerr << "ml[" << mI << "]->Start(1) == " << ml[mI]->Start(0) << endl;
				cerr << "ml[" << mI << "]->Length(1) == " << ml[mI]->Length(1) << endl;
				cerr << "ml.size(): " << ml.size() << endl;
				for( sivI = 0; sivI < siv_list.size(); sivI++ )
				{
					cerr << "siv_list[" << sivI << "] left end 0: " << siv_list[sivI].reference_iv.LeftEnd(0)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(0) != 0 )
						cerr << "siv_list[" << sivI << "] right end 0: " << siv_list[sivI].reference_iv.LeftEnd(0) + siv_list[sivI].reference_iv.Length(0) << endl;
					cerr << "siv_list[" << sivI << "] left end 1: " << siv_list[sivI].reference_iv.LeftEnd(1)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(1) != 0 )
						cerr << "siv_list[" << sivI << "] right end 1: " << siv_list[sivI].reference_iv.LeftEnd(1) + siv_list[sivI].reference_iv.Length(1) << endl;
				}
				breakHere();
			}
			if( ml[mI]->LeftEnd(seqI) + ml[mI]->Length(seqI) > 
				siv_list[sivI].reference_iv.LeftEnd(cur_child) + siv_list[sivI].reference_iv.Length(cur_child) )
			{
				cerr << "doesn't fit\n";
				cerr << "ml[" << mI << "]->LeftEnd(" << seqI << "): " << ml[mI]->LeftEnd(seqI) << endl;
				cerr << "ml[" << mI << "]->RightEnd(" << seqI << "): " << ml[mI]->RightEnd(seqI) << endl;
				cerr << "siv_list[" << sivI << "] left end 0: " << siv_list[sivI].reference_iv.LeftEnd(0)  << endl;
				if( siv_list[sivI].reference_iv.LeftEnd(0) != 0 )
					cerr << "siv_list[" << sivI << "] right end 0: " << siv_list[sivI].reference_iv.LeftEnd(0) + siv_list[sivI].reference_iv.Length(0) << endl;
				cerr << "siv_list[" << sivI << "] left end 1: " << siv_list[sivI].reference_iv.LeftEnd(1)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(1) != 0 )
						cerr << "siv_list[" << sivI << "] right end 1: " << siv_list[sivI].reference_iv.LeftEnd(1) + siv_list[sivI].reference_iv.Length(1) << endl;
				cerr << "ml.size(): " << ml.size() << endl;
				cerr << "siv_list.size(): " << siv_list.size() << endl;
				cerr << "trans_path:";
				for( size_t tI = 0; tI < trans_path.size(); tI++ )
					cerr << " " << trans_path[tI];
				cerr << endl;
				cerr << "trans_path[" << nI << "]: " << trans_path[nI] << endl;
				breakHere();
			}

			ml[mI]->SetLeftEnd( seqI, ml[mI]->LeftEnd(seqI) - siv_list[sivI].reference_iv.LeftEnd(cur_child) + 1 );
			// if this interval matches the reverse strand then we should effectively invert all matches
			if( siv_list[sivI].reference_iv.Start(cur_child) < 0 )
			{
				int64 new_lend = siv_list[sivI].reference_iv.Length(cur_child) - ml[mI]->LeftEnd(seqI);
				new_lend -= ml[mI]->Length( seqI ) - 2;
				new_lend *= ml[mI]->Orientation(seqI) == AbstractMatch::forward ? 1 : -1;
				ml[mI]->Invert();
				ml[mI]->SetStart( seqI, new_lend ); 
			}
			siv_matches[sivI].push_back( ml[mI] );
		}

		// apply the coordinate translation
		ml.clear();
		for( size_t sivI = 0; sivI < siv_matches.size(); sivI++ )
		{
			if( siv_matches[sivI].size() == 0 )
				continue;
			
			// get a CompactGappedAlignment<> for this interval
			CompactGappedAlignment<>* siv_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_list[sivI].reference_iv.GetMatches()[0]);
			if( siv_list[sivI].reference_iv.GetMatches().size() > 1 )
				siv_cga = NULL;
			bool alloc_new_siv = false;
			CompactGappedAlignment<> tmp_cga;
			if( siv_cga == NULL )
			{
				alloc_new_siv = true;
				siv_cga = tmp_cga.Copy();
				CompactGappedAlignment<> dorkas(siv_list[sivI].reference_iv);
				*siv_cga = dorkas;
			}

			// now translate each match...
			for( size_t mI = 0; mI < siv_matches[sivI].size(); mI++ )
			{
				CompactGappedAlignment<>* match_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_matches[sivI][mI]);
				bool alloc_new = false;
				if( match_cga == NULL )
				{
					match_cga = tmp_cga.Copy();
					*match_cga = CompactGappedAlignment<>(*(siv_matches[sivI][mI]));
					alloc_new = true;
				}
				siv_cga->translate( *match_cga, seqI, cur_child );

				if( alloc_new )
				{
					siv_matches[sivI][mI]->Free();
					siv_matches[sivI][mI] = match_cga;
				}
			}

			// shift coordinates back to global space
			for( size_t mI = 0; mI < siv_matches[sivI].size(); mI++ )
			{
				int64 cur_start = siv_matches[sivI][mI]->Start(seqI);
				if( cur_start > 0 )
					siv_matches[sivI][mI]->SetStart( seqI, cur_start + siv_list[sivI].LeftEnd() - 1 );
				else
					siv_matches[sivI][mI]->SetStart( seqI, cur_start - siv_list[sivI].LeftEnd() + 1);
				if( (siv_matches[sivI][mI]->LeftEnd(seqI) + siv_matches[sivI][mI]->Length(seqI) > siv_list.back().LeftEnd() + siv_list.back().Length() )
					 )
				{
					// is there something wrong with the translation table?
					cerr << "siv left is: " << siv_list[sivI].LeftEnd() << endl;
					cerr << "siv right is: " << siv_list[sivI].LeftEnd() + siv_list[sivI].Length() << endl;
					cerr << "match right is: " << siv_matches[sivI][mI]->LeftEnd(seqI) + siv_matches[sivI][mI]->Length(seqI) << endl;
					cerr << "superseq right is: " << siv_list.back().LeftEnd() + siv_list.back().Length() << endl;
					cerr << "";
					breakHere();
				}
				if( debug_aligner && siv_matches[sivI][mI]->Start(seqI) == 0 )
				{
					breakHere();
				}
			}
			if(alloc_new_siv)
				siv_cga->Free();
			ml.insert( ml.end(), siv_matches[sivI].begin(), siv_matches[sivI].end() );
		}
	}
	// restore forward orientation seqI
	for( size_t mI = 0; mI < ml.size(); mI++ )
		if( ml[mI]->Orientation(seqI) == AbstractMatch::reverse )
			ml[mI]->Invert();
}

class SuperIntervalPtrComp
{
public:
	bool operator()( const SuperInterval* a, const SuperInterval* b )
	{
		return (*a) < (*b);
	}
};

void ProgressiveAligner::recursiveApplyAncestralBreakpoints( node_id_t ancestor )
{
	stack<node_id_t> node_stack;
	node_stack.push(ancestor);
	while( node_stack.size() > 0 )
	{
		// pop the current node, apply ancestral breakpoints, recurse on children
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		SuperIntervalManipulator sim;
		if( progress_msgs ) cout << "cur node: " << cur_node << endl;
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			AlignmentTreeNode& atn = alignment_tree[alignment_tree[cur_node].children[childI]];
			if( progress_msgs ) cout << "childI " << childI << " aab\n";
			applyAncestralBreakpoints( alignment_tree[cur_node].ordering, atn.ordering, childI, sim );
			if( progress_msgs ) cout << "sort childI " << childI << "\n";
			vector<SuperInterval*> siv_ptr_list(atn.ordering.size());
			for( size_t sivI = 0; sivI < atn.ordering.size(); ++sivI )
				siv_ptr_list[sivI] = &(atn.ordering[sivI]);
			SuperIntervalPtrComp sipc;
			sort( siv_ptr_list.begin(), siv_ptr_list.end(), sipc );
			vector< SuperInterval > siv_list;
			for( size_t sivI = 0; sivI < siv_ptr_list.size(); ++sivI )
				siv_list.push_back(*siv_ptr_list[sivI]);
			swap(siv_list, atn.ordering);
			node_stack.push( alignment_tree[cur_node].children[childI] );
		}
		if( debug_aligner && alignment_tree[cur_node].children.size() > 0 )
			validateSuperIntervals(alignment_tree[cur_node].children[0], alignment_tree[cur_node].children[1], cur_node);
		if( progress_msgs ) cout << "linking node " << cur_node << "'s" << alignment_tree[cur_node].ordering.size() << " superintervals\n"; 
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
			linkSuperIntervals( alignment_tree[cur_node].children[childI], childI, cur_node );
	}
}


void ProgressiveAligner::pairwiseAnchorSearch( MatchList& r_list, Match* r_begin, Match* r_end )
{
	try
	{
		uint seqI = 0;
		MatchList gap_list;
		vector< int64 > starts;
	// 
	//	Get the sequence in the intervening gaps between these two matches
	//
		for( seqI = 0; seqI < 2; seqI++ )
		{
			int64 gap_end = 0;
			int64 gap_start = 0;
			getInterveningCoordinates( r_list.seq_table, r_begin, r_end, seqI, gap_start, gap_end );
			int64 diff = gap_end - gap_start;
			diff = 0 < diff ? diff : 0;

			starts.push_back( gap_start );
			gnSequence* new_seq = new gnSequence( r_list.seq_table[ seqI ]->subseq( gap_start, diff ) );
			gap_list.seq_table.push_back( new_seq );
			gap_list.sml_table.push_back( new DNAMemorySML() );
		}

		gnSeqI avg_len = (gap_list.seq_table[0]->length() + gap_list.seq_table[1]->length())/2;
		uint search_seed_size = getDefaultSeedWeight( avg_len );
		gap_mh.Clear();

		//
		//	Create sorted mer lists for the intervening gap region
		//
		uint64 default_seed = getSeed( search_seed_size );
		if( search_seed_size < MIN_DNA_SEED_WEIGHT )
		{
			for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ )
				delete gap_list.seq_table[ seqI ];
			for( uint seqI = 0; seqI < gap_list.sml_table.size(); seqI++ )
				delete gap_list.sml_table[ seqI ];
			return;
		}
		for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ ){
			gap_list.sml_table[ seqI ]->Clear();
			gap_list.sml_table[ seqI ]->Create( *(gap_list.seq_table[ seqI ]), default_seed );
		}

		//
		//	Find all matches in the gap region
		//
		gap_mh.ClearSequences();
		gap_mh.FindMatches( gap_list );
		
		EliminateOverlaps_v2( gap_list );
		
		// for anchor accuracy, throw out any anchors that are shorter than the minimum
		// anchor length after EliminateOverlaps()
		gap_list.LengthFilter( MIN_ANCHOR_LENGTH );

		if( gap_list.size() > 0 )
		{
			// shift all the matches that were found
			vector< Match* >::iterator mum_iter = gap_list.begin();
			for( ; mum_iter != gap_list.end(); ){
				boolean add_ok = true;
				for( uint seqI = 0; seqI < (*mum_iter)->SeqCount(); seqI++ ){
					int64 gap_start;
					if( (*mum_iter)->Start( seqI ) < 0 ){
						gap_start = r_begin != NULL ? -r_begin->End( seqI ) : 0;
						if( gap_start > 0 )
							gap_start = r_end != NULL ? r_end->Start( seqI ) - r_end->Length() + 1 : 0;
						else if( r_begin )
							add_ok = false;
						(*mum_iter)->SetStart( seqI, (*mum_iter)->Start( seqI ) + gap_start );
					}else{
						// insert them all before mem_iter
						gap_start = r_begin != NULL ? r_begin->End( seqI ) : 0;
						if( gap_start < 0 ){
							gap_start = r_end != NULL ? r_end->Start( seqI ) - r_end->Length() + 1 : 0;
							add_ok = false;
						}
						(*mum_iter)->SetStart( seqI, (*mum_iter)->Start( seqI ) + gap_start );
					}
				}
				if( add_ok )
					r_list.push_back( *mum_iter );
				else{
					(*mum_iter)->UnlinkSelf();
					(*mum_iter)->Free();
					(*mum_iter) = NULL;
				}
				++mum_iter;
			}
		}

		// delete sequences and smls
		for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ )
			delete gap_list.seq_table[ seqI ];
		for( uint seqI = 0; seqI < gap_list.sml_table.size(); seqI++ )
			delete gap_list.sml_table[ seqI ];
			
	}catch( gnException& gne ){
		cerr << gne << endl;
	}catch( exception& e ){
		cerr << e.what() << endl;
	}catch(...){
		cerr << "When I say 'ohhh' you say 'shit'!\n";
	}
}


template<class GappedAlignmentType>
void ProgressiveAligner::recurseOnPairs( const vector<node_id_t>& node1_seqs, const vector<node_id_t>& node2_seqs, const GappedAlignmentType& iv, Matrix<MatchList>& matches, Matrix< std::vector< search_cache_t > >& search_cache_db, Matrix< std::vector< search_cache_t > >& new_cache_db )
{
	matches = Matrix<MatchList>(node1_seqs.size(),node2_seqs.size());

	vector<node_id_t>::const_iterator n1_iter = node1_seqs.begin();
	std::vector< bitset_t > aln_matrix;
	iv.GetAlignment(aln_matrix);
	Match tmp(2);
	for( size_t n1 = 0; n1 < node1_seqs.size(); n1++ )
	{
		if( n1 > 0 )
			++n1_iter;
		vector<node_id_t>::const_iterator n2_iter = node2_seqs.begin();
		for( size_t n2 = 0; n2 < node2_seqs.size(); n2++ )
		{
			if( n2 > 0 )
				++n2_iter;
			uint seqI = node_sequence_map[*n1_iter];
			uint seqJ = node_sequence_map[*n2_iter];
			MatchList& mlist = matches(n1, n2);
			std::vector< search_cache_t >& cache = search_cache_db(n1, n2);
			std::vector< search_cache_t >& new_cache = new_cache_db(n1, n2);
			mlist.seq_table.push_back( alignment_tree[*n1_iter].sequence );
			mlist.seq_table.push_back( alignment_tree[*n2_iter].sequence );

			gnSeqI charI = 0;
			gnSeqI charJ = 0;
			gnSeqI prev_charI = 0;
			gnSeqI prev_charJ = 0;
			bool in_gap = false;
			for( uint colI = 0; colI <= iv.AlignmentLength(); colI++ )
			{
				if( colI == iv.AlignmentLength() || 
					(aln_matrix[seqI].test(colI) && aln_matrix[seqJ].test(colI)) )
				{
					if( in_gap && 
						charI - prev_charI > min_recursive_gap_length &&
						charJ - prev_charJ > min_recursive_gap_length )
					{

						Match* l_match = NULL;
						l_match = tmp.Copy();
						if(iv.Orientation(seqI) == AbstractMatch::forward)
							l_match->SetLeftEnd(0, iv.LeftEnd(seqI)+prev_charI);
						else
						{
							l_match->SetLeftEnd(0, iv.RightEnd(seqI)-prev_charI);
							l_match->SetOrientation(0, AbstractMatch::reverse );
						}
						if(iv.Orientation(seqJ) == AbstractMatch::forward)
							l_match->SetLeftEnd(1, iv.LeftEnd(seqJ)+prev_charJ);
						else
						{
							l_match->SetLeftEnd(1, iv.RightEnd(seqJ)-prev_charJ);
							l_match->SetOrientation(1, AbstractMatch::reverse );
						}
						l_match->SetLength(0);

						Match* r_match = NULL;
						if( charJ != iv.RightEnd(seqJ) && charI != iv.RightEnd(seqI) )
						{
							r_match = tmp.Copy();
							if(iv.Orientation(seqI) == AbstractMatch::forward)
								r_match->SetLeftEnd(0, iv.LeftEnd(seqI)+charI);
							else
							{
								r_match->SetLeftEnd(0, iv.RightEnd(seqI)-charI);
								r_match->SetOrientation(0, AbstractMatch::reverse );
							}
							if(iv.Orientation(seqJ) == AbstractMatch::forward)
								r_match->SetLeftEnd(1, iv.LeftEnd(seqJ)+charJ);
							else
							{
								r_match->SetLeftEnd(1, iv.RightEnd(seqJ)-charJ);
								r_match->SetOrientation(1, AbstractMatch::reverse );
							}
							r_match->SetLength(0);
						}

						if( iv.Orientation(seqI) == AbstractMatch::reverse )
						{
							swap(l_match,r_match);
							if( l_match != NULL ) l_match->Invert();
							if( r_match != NULL ) r_match->Invert();
						}

						// check whether the current cache already has the searched region
						search_cache_t cacheval = make_pair( l_match, r_match );
						std::vector< search_cache_t >::iterator cache_entry = std::upper_bound( cache.begin(), cache.end(), cacheval, mems::cache_comparator );
						if( cache_entry == cache.end() || 
							(mems::cache_comparator( cacheval, *cache_entry ) || mems::cache_comparator( *cache_entry, cacheval )) )
						{
							// search this region
							pairwiseAnchorSearch(mlist, l_match, r_match);
						}
						new_cache.push_back( cacheval );
					}
					prev_charI = charI;
					prev_charJ = charJ;
					in_gap = false;
				}
				else
					in_gap = true;
				if( colI < iv.AlignmentLength() )
				{
					if( aln_matrix[seqI].test(colI) )
						++charI;
					if( aln_matrix[seqJ].test(colI) )
						++charJ;
				}
			}
		}
	}
}

void ProgressiveAligner::getAncestralMatches( const vector< node_id_t > node1_seqs, const vector< node_id_t > node2_seqs, node_id_t node1, node_id_t node2, node_id_t ancestor, std::vector< AbstractMatch* >& ancestral_matches )
{
	// to save memory, always make node1_seqs the bigger vector
//	if( node1_seqs.size() < node2_seqs.size() )
//		swap( node1_seqs, node2_seqs );

	// for each pair of genomes, extract pairwise matches and translate up
	// eliminate overlaps
	for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
	{
		uint ii = this->node_sequence_map[node1_seqs[seqI]];
		vector< AbstractMatch* > seqI_matches;

		for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
		{
			uint jj = this->node_sequence_map[node2_seqs[seqJ]];
			vector< AbstractMatch* > cur_matches;
			for( size_t mI = 0; mI < original_ml.size(); mI++ )
			{
				if( original_ml[mI]->LeftEnd(ii) == NO_MATCH )
					continue;
				if( original_ml[mI]->LeftEnd(jj) == NO_MATCH )
					continue;
				Match mm( 2 );
				Match* new_m = mm.Copy();
				new_m->SetStart( 0, original_ml[mI]->Start(ii));
				new_m->SetStart( 1, original_ml[mI]->Start(jj));
				new_m->SetLength(original_ml[mI]->Length());
				if( new_m->Start(0) < 0 )
					new_m->Invert();	// assign reference orientation to seq 0
				cur_matches.push_back( new_m );
			}
			// now translate cur_matches
			translateGappedCoordinates( cur_matches, 1, node2_seqs[seqJ], node2 );
			seqI_matches.insert( seqI_matches.end(), cur_matches.begin(), cur_matches.end() );
		}
		EliminateOverlaps_v2( seqI_matches );
		translateGappedCoordinates( seqI_matches, 0, node1_seqs[seqI], node1 );
		ancestral_matches.insert( ancestral_matches.end(), seqI_matches.begin(), seqI_matches.end() );
	}
	EliminateOverlaps_v2( ancestral_matches );
}


void ProgressiveAligner::getPairwiseMatches( const vector< node_id_t >& node1_seqs, const vector< node_id_t >& node2_seqs, Matrix<MatchList>& pairwise_matches )
{
	pairwise_matches = Matrix< MatchList >( node1_seqs.size(), node2_seqs.size() );

	// copy sequence tables
	for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
	{
		for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
		{
			uint ii = this->node_sequence_map[node1_seqs[seqI]];
			uint jj = this->node_sequence_map[node2_seqs[seqJ]];
			pairwise_matches(seqI, seqJ).seq_table.push_back(original_ml.seq_table[ii]);
			pairwise_matches(seqI, seqJ).seq_table.push_back(original_ml.seq_table[jj]);
			pairwise_matches(seqI, seqJ).seq_filename.push_back(original_ml.seq_filename[ii]);
			pairwise_matches(seqI, seqJ).seq_filename.push_back(original_ml.seq_filename[jj]);
		}
	}

	// now copy pairwise matches
	for( size_t mI = 0; mI < original_ml.size(); mI++ )
	{
		for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
		{
			uint ii = this->node_sequence_map[node1_seqs[seqI]];
			if( original_ml[mI]->LeftEnd(ii) == NO_MATCH )
				continue;
			for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			{
				uint jj = this->node_sequence_map[node2_seqs[seqJ]];
				if( original_ml[mI]->LeftEnd(jj) == NO_MATCH )
					continue;
				Match mm( 2 );
				Match* new_m = mm.Copy();
				new_m->SetStart( 0, original_ml[mI]->Start(ii));
				new_m->SetStart( 1, original_ml[mI]->Start(jj));
				new_m->SetLength(original_ml[mI]->Length());
				if( new_m->Start(0) < 0 )
					new_m->Invert();	// assign reference orientation to seq 0
				pairwise_matches(seqI,seqJ).push_back( new_m );
			}
		}
	}
}


int IsDenseEnough( GappedAlignment* gal_iter )
{
	double total_len = 0;
	gnSeqI seqs = 0;
	for( uint seqI = 0; seqI < gal_iter->SeqCount(); seqI++ )
	{
		if( gal_iter->LeftEnd(seqI) == NO_MATCH )
			continue;
		total_len += gal_iter->Length(seqI);
	}
	double density = total_len / (gal_iter->AlignmentLength() * (double)gal_iter->Multiplicity());
	// density of 1 is ideal
	// the shorter the alignment, the closer we should be to 1 to allow splitting
	// use a linear threshold with (min_window_size,1) and (max_window_size,min_gappiness)
	// as endpoints of the threshold line
	
	// determine the density threshold for the given alignment length
	double threshold = ((max_density - min_density)/(min_window_size - max_window_size)) * ( (double)gal_iter->AlignmentLength() - max_window_size ) + min_density;
	if( density > max_density )	// don't bother aligning this, it's so dense we'll wait until iterative refinement.
		return 2;
	if( density > threshold )
		return 1;
	return 0;
}

void splitGappedAlignment( const GappedAlignment& ga, GappedAlignment& ga1, GappedAlignment& ga2, std::vector<uint>& seqs1, std::vector<uint>& seqs2 )
{
	const vector< string >& aln = GetAlignment( ga, std::vector<gnSequence*>(ga.SeqCount()) );
	ga1 = ga;
	ga2 = ga;
	for( uint seqI = 0; seqI < seqs1.size(); seqI++ )
		ga2.SetLeftEnd(seqs1[seqI], NO_MATCH);
	for( uint seqI = 0; seqI < seqs2.size(); seqI++ )
		ga1.SetLeftEnd(seqs2[seqI], NO_MATCH);
}

void removeLargeGapsPP( GappedAlignment& gal, list< GappedAlignment* >& gal_list, vector<bool>& gap_iv, const vector< size_t >& group1, const vector< size_t >& group2 )
{
	// scan through and remove any section where members of group1 aren't aligned to members of group2
	// for more than some number of nucleotides
	gap_iv.clear();
	gal_list.clear();
	const vector< string >& aln_matrix = GetAlignment(gal, vector<gnSequence*>(gal.SeqCount(),NULL));
	size_t gap_cols = 0;
	size_t last_aln_col = (std::numeric_limits<size_t>::max)();
	for( size_t colI = 0; colI < gal.AlignmentLength(); colI++ )
	{
		 size_t g1 = 0;
		 size_t g2 = 0;
		 for( ; g1 < group1.size(); ++g1 )
		 {
			 if( aln_matrix[group1[g1]][colI] != '-' )
				 break;
		 }
		 for( ; g2 < group2.size(); ++g2 )
		 {
			 if( aln_matrix[group2[g2]][colI] != '-' )
				 break;
		 }
		 if( g1 < group1.size() && g2 < group2.size() )
		 {
			 // it's an aligned col
			 if( gap_cols > max_gap_length )
			 {
				// crop out the middle gapped section
				gnSeqI split_point = 0;
				if( last_aln_col != (std::numeric_limits<size_t>::max)() )
				{
					split_point = last_aln_col + lcb_hangover;
					gal_list.push_back( gal.Copy() );
					gap_iv.push_back(false);
					gal_list.back()->CropEnd( gal.AlignmentLength()-split_point );
					gal.CropStart( split_point );
				}
				split_point = colI - lcb_hangover - split_point;
				gal_list.push_back( gal.Copy() );
				gap_iv.push_back(true);
				gal_list.back()->CropEnd( gal.AlignmentLength()-split_point );
				gal.CropStart( split_point );
			 }
			 gap_cols = 0;
			 last_aln_col = colI;
		 }else
			 ++gap_cols;
	}

	if( gap_cols > max_gap_length )
	{
		gnSeqI split_point = 0;
		if( last_aln_col != (std::numeric_limits<size_t>::max)() )
		{
			split_point = last_aln_col + lcb_hangover;
			gal_list.push_back( gal.Copy() );
			gap_iv.push_back(false);
			gal_list.back()->CropEnd( gal.AlignmentLength()-split_point );
			gal.CropStart( split_point );
		}
		gap_iv.push_back(true);
	}else
		gap_iv.push_back(false);
	gal_list.push_back( gal.Copy() );
}

void ProgressiveAligner::refineAlignment( GappedAlignment& gal, node_id_t ancestor, bool profile_aln, AlnProgressTracker& apt )
{
	// divide the gapped alignment up into windows of a given size and have
	// muscle refine the alignments
	// when anchors are dense use smaller windows to improve speed efficiency
	list< GappedAlignment* > gal_list;
	vector<bool> gap_iv;
	std::vector<node_id_t> nodes1;
	std::vector<node_id_t> nodes2;
	getAlignedChildren( alignment_tree[ancestor].children[0], nodes1 );
	getAlignedChildren( alignment_tree[ancestor].children[1], nodes2 );
	std::vector<uint> seqs1( nodes1.size() );
	std::vector<uint> seqs2( nodes2.size() );
	for( size_t nI = 0; nI < nodes1.size(); nI++ )
		seqs1[nI] = node_sequence_map[nodes1[nI]];
	for( size_t nI = 0; nI < nodes2.size(); nI++ )
		seqs2[nI] = node_sequence_map[nodes2[nI]];
//	if( profile_aln )
//	{
//		removeLargeGapsPP( gal, gal_list, gap_iv, seqs1, seqs2 );
//	}else{
		gal_list.push_back( gal.Copy() );
		gap_iv.push_back(false);
//	}
	list< GappedAlignment* >::iterator gal_iter = gal_list.begin();
	vector<bool>::iterator gap_iter = gap_iv.begin();
	while(gal_iter != gal_list.end())
	{
		int density = IsDenseEnough( *gal_iter );
		if( (density == 0 && (*gal_iter)->AlignmentLength() > max_window_size / 3) ||
			(density == 1 && (*gal_iter)->AlignmentLength() > max_window_size ) ||
			(density == 2 && (*gal_iter)->AlignmentLength() > max_window_size * 3 )

//			  || ( (*gal_iter)->AlignmentLength() > min_window_size && density == 1 && profile_aln == true ) 
			  )
		{
			// split in half
			gnSeqI split_point = (*gal_iter)->AlignmentLength() / 2;
			list< GappedAlignment* >::iterator ins_iter = gal_iter;
			++ins_iter;
			ins_iter = gal_list.insert(ins_iter, (*gal_iter)->Copy());
			vector<bool>::iterator gap_ins_iter = gap_iter;
			size_t gap_off = gap_iter - gap_iv.begin();
			++gap_ins_iter;
			gap_iv.insert( gap_ins_iter, *gap_iter );
			gap_iter = gap_iv.begin() + gap_off;
			(*gal_iter)->CropEnd( split_point );
			(*ins_iter)->CropStart( (*ins_iter)->AlignmentLength() - split_point );
			continue;
		}

		++gal_iter;
		++gap_iter;
	}
	MuscleInterface& mi = MuscleInterface::getMuscleInterface();
	// now that the alignment is all split up use muscle to refine it
	gnSeqI new_len = 0;

	gap_iter = gap_iv.begin();
	for( gal_iter = gal_list.begin(); gal_iter != gal_list.end(); ++gal_iter )
	{
		apt.cur_leftend += (*gal_iter)->AlignmentLength();
		if( profile_aln && !(*gap_iter) )
		{
			GappedAlignment ga1;
			GappedAlignment ga2;
			splitGappedAlignment( **gal_iter, ga1, ga2, seqs1, seqs2 );
			if( ga1.Multiplicity() > 0 && ga2.Multiplicity() > 0 )
				mi.ProfileAlign( ga1, ga2, **gal_iter, true );
		}else
		{
			int density = IsDenseEnough( *gal_iter );
			if( density == 0 )
				mi.Refine( **gal_iter );
			else if( density == 1 )
				mi.Refine( **gal_iter, 500 );
			else
				mi.Refine( **gal_iter, 200 );
		}

		new_len += (*gal_iter)->AlignmentLength();
		++gap_iter;

		// print a progress message
		double cur_progress = ((double)apt.cur_leftend / (double)apt.total_len)*100.0;
		printProgress(apt.prev_progress, cur_progress, cout);
		apt.prev_progress = cur_progress;
	}

	// put humpty dumpty back together
	vector< string > aln_matrix( gal.SeqCount(), string( new_len, '-' ) );
	vector< string::size_type > pos( gal.SeqCount(), 0 );
	for( gal_iter = gal_list.begin(); gal_iter != gal_list.end(); ++gal_iter )
	{
		const vector< string >& tmp_mat = GetAlignment(**gal_iter, vector<gnSequence*>( gal.SeqCount() ) );
		for( uint seqI = 0; seqI < tmp_mat.size(); seqI++ )
		{
			if( gal.LeftEnd(seqI) == 0 )
				continue;
			aln_matrix[seqI].replace(pos[seqI], tmp_mat[seqI].size(), tmp_mat[seqI]);
			pos[seqI] += tmp_mat[seqI].size();
		}
		(*gal_iter)->Free();
	}
	gal.SetAlignment(aln_matrix);
}

void ProgressiveAligner::doGappedAlignment( node_id_t ancestor, bool profile_aln )
{
	AlnProgressTracker apt;
	gnSeqI total_len = 0;
	for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
		total_len += alignment_tree[ancestor].ordering[aI].Length();
	apt.total_len = total_len;
	apt.prev_progress = 0;

	printProgress(-1, 0, cout);
	apt.cur_leftend = 1;

	for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
	{
		if( alignment_tree[ancestor].ordering[aI].reference_iv.Multiplicity() == 1 )
		{
			apt.cur_leftend += alignment_tree[ancestor].ordering[aI].reference_iv.AlignmentLength();
			continue;	// don't bother re-refining intervals that didn't get aligned here
		}

		GappedAlignment gal;
		extractAlignment(ancestor, aI, gal);
		if( gal.Multiplicity() > 1 )	// no point in refining intervals that are unaligned anyways
			refineAlignment( gal, ancestor, profile_aln, apt );
		else
			apt.cur_leftend += gal.AlignmentLength();
		ConstructSuperIntervalFromMSA(ancestor, aI, gal);

		// print a progress message
		double cur_progress = ((double)apt.cur_leftend / (double)apt.total_len)*100.0;
		printProgress(apt.prev_progress, cur_progress, cout);
		apt.prev_progress = cur_progress;
	}
	FixLeftEnds(ancestor);

	if( debug_aligner )
		validateSuperIntervals(alignment_tree[ancestor].children[0], alignment_tree[ancestor].children[1], ancestor);
	cout << "\ndone.\n";
}

void ProgressiveAligner::FixLeftEnds( node_id_t ancestor )
{
	// fixes all SuperInterval left-end coordinates for nodes below ancestor
	stack< node_id_t > node_stack;
	node_stack.push( ancestor );
	vector<bool> visited( alignment_tree.size(), false );
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		// visit post-order
		if( !visited[cur_node] )
		{
			for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
				node_stack.push( alignment_tree[cur_node].children[childI] );
			visited[cur_node] = true;
			continue;
		}
		node_stack.pop();
		if( alignment_tree[cur_node].sequence != NULL )
			continue;	// don't do anything on leaf nodes

		vector< SuperInterval >& siv_list = alignment_tree[cur_node].ordering;
		gnSeqI left_end = 1;
		for( size_t sivI = 0; sivI < siv_list.size(); sivI++ )
		{
			siv_list[sivI].SetLeftEnd(left_end);
			siv_list[sivI].SetLength(siv_list[sivI].reference_iv.AlignmentLength());
			left_end += siv_list[sivI].reference_iv.AlignmentLength();
			CompactGappedAlignment<>* m_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_list[sivI].reference_iv.GetMatches()[0]);
			
			// this one wasn't refined, just move it appropriately
			if( m_cga == NULL || siv_list[sivI].reference_iv.GetMatches().size() > 1 )
			{
				for( uint childI = 0; childI <= 1; childI++ )
				{
					size_t cur_siv = childI == 0 ? alignment_tree[cur_node].ordering[sivI].c1_siv : alignment_tree[cur_node].ordering[sivI].c2_siv;
					if( cur_siv == (std::numeric_limits<size_t>::max)() )
						continue;
					const SuperInterval& c_siv = alignment_tree[ alignment_tree[cur_node].children[childI] ].ordering[ cur_siv ];
					int64 diff = c_siv.LeftEnd() - siv_list[sivI].reference_iv.LeftEnd(childI);
					siv_list[sivI].reference_iv.SetLeftEnd(childI, c_siv.LeftEnd());
					const vector< AbstractMatch* >& matches = siv_list[sivI].reference_iv.GetMatches();
					for( size_t mI = 0; mI < matches.size(); mI++ )
					{
						if( matches[mI]->LeftEnd(childI) != NO_MATCH )
							matches[mI]->SetLeftEnd(childI, matches[mI]->LeftEnd(childI) + diff);
					}
				}

			}else{

				size_t c1_siv = alignment_tree[cur_node].ordering[sivI].c1_siv;
				if( c1_siv != (std::numeric_limits<size_t>::max)() )
				{
					const SuperInterval& c_siv = alignment_tree[ alignment_tree[cur_node].children[0] ].ordering[ c1_siv ];
					m_cga->SetLeftEnd(0, c_siv.LeftEnd());
					siv_list[sivI].reference_iv.SetLeftEnd(0, c_siv.LeftEnd());
					m_cga->SetLength(c_siv.Length(), 0);
					siv_list[sivI].reference_iv.SetLength(c_siv.Length(), 0);
					siv_list[sivI].reference_iv.SetOrientation(0, m_cga->Orientation(0));
				}
				size_t c2_siv = alignment_tree[cur_node].ordering[sivI].c2_siv;
				if( c2_siv != (std::numeric_limits<size_t>::max)() )
				{
					const SuperInterval& c_siv = alignment_tree[ alignment_tree[cur_node].children[1] ].ordering[ c2_siv ];
					m_cga->SetLeftEnd(1, c_siv.LeftEnd());
					siv_list[sivI].reference_iv.SetLeftEnd(1, c_siv.LeftEnd());
					m_cga->SetLength(c_siv.Length(), 1);
					siv_list[sivI].reference_iv.SetLength(c_siv.Length(), 1);
					siv_list[sivI].reference_iv.SetOrientation(1, m_cga->Orientation(1));
				}
			}
			if( debug_cga && m_cga && !m_cga->validate() )
//			if( m_cga && !m_cga->validate() )
				cerr << "oh junkedy\n";

			if( siv_list[sivI].reference_iv.Length(0) > 20000000 )
			{
				breakHere();
				cerr << "explosion 1\n";
			}
			if( siv_list[sivI].reference_iv.Length(1) > 20000000 )
			{
				breakHere();
				cerr << "explosion 2\n";
			}
		}
	}
}

void propagateInvert( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t ancestor, size_t ans_siv )
{
	stack< pair< node_id_t, size_t > > node_siv_stack;
	node_siv_stack.push( make_pair(ancestor, ans_siv) );
	while( node_siv_stack.size() > 0 )
	{
		pair< node_id_t, size_t > cur = node_siv_stack.top();
		node_siv_stack.pop();
		node_id_t cur_node = cur.first;
		if( alignment_tree[cur_node].ordering[cur.second].c1_siv != (std::numeric_limits<size_t>::max)() )
			node_siv_stack.push( make_pair( alignment_tree[cur_node].children[0], alignment_tree[cur_node].ordering[cur.second].c1_siv ) );
		if( alignment_tree[cur_node].ordering[cur.second].c2_siv != (std::numeric_limits<size_t>::max)() )
			node_siv_stack.push( make_pair( alignment_tree[cur_node].children[1], alignment_tree[cur_node].ordering[cur.second].c2_siv ) );
		if( cur_node == ancestor )
			continue;	// don't do anything at the ancestor
		if( alignment_tree[cur_node].sequence != NULL )
			continue;	// don't do anything on leaf nodes

		// reverse the homology structure at this node
		Interval& ref_iv = alignment_tree[cur_node].ordering[cur.second].reference_iv;
		vector< AbstractMatch* > matches;
		ref_iv.StealMatches( matches );
		AbstractMatch::orientation o0 = matches[0]->Orientation(0);
		AbstractMatch::orientation o1 = matches[0]->Orientation(1);
		matches[0]->Invert();
		if( o0 != AbstractMatch::undefined )
			matches[0]->SetOrientation(0,o0);
		if( o1 != AbstractMatch::undefined )
			matches[0]->SetOrientation(1,o1);
		ref_iv.SetMatches( matches );
		if( o0 != AbstractMatch::undefined )
		{
			ref_iv.SetOrientation(0,o0);
			ref_iv.SetLeftEnd(0,0);
		}
		if( o1 != AbstractMatch::undefined )
		{
			ref_iv.SetOrientation(1,o1);
			ref_iv.SetLeftEnd(1,0);
		}
	}
}


void ProgressiveAligner::ConstructSuperIntervalFromMSA( node_id_t ancestor, size_t ans_siv, GappedAlignment& gal )
{
	const vector< string >& aln_matrix = GetAlignment( gal, vector< gnSequence* >() );
	stack< pair< node_id_t, size_t > > node_siv_stack;
	node_siv_stack.push( make_pair(ancestor, ans_siv) );
	vector<bool> visited( alignment_tree.size(), false );
	while( node_siv_stack.size() > 0 )
	{
		pair< node_id_t, size_t > cur = node_siv_stack.top();
		node_id_t cur_node = cur.first;
		// visit post-order
		if( !visited[cur_node] )
		{
			if( alignment_tree[cur_node].ordering[cur.second].c1_siv != (std::numeric_limits<size_t>::max)() )
				node_siv_stack.push( make_pair( alignment_tree[cur_node].children[0], alignment_tree[cur_node].ordering[cur.second].c1_siv ) );
			if( alignment_tree[cur_node].ordering[cur.second].c2_siv != (std::numeric_limits<size_t>::max)() )
				node_siv_stack.push( make_pair( alignment_tree[cur_node].children[1], alignment_tree[cur_node].ordering[cur.second].c2_siv ) );
			visited[cur_node] = true;
			continue;
		}
		node_siv_stack.pop();
		if( alignment_tree[cur_node].sequence != NULL )
			continue;	// don't do anything on leaf nodes

		// build a super-interval
		vector< node_id_t > node1_seqs;	/**< the node id's of extant sequences below node 1 */
		vector< node_id_t > node2_seqs;	/**< the node id's of extant sequences below node 2 */
		getAlignedChildren( alignment_tree[cur_node].children[0], node1_seqs );
		getAlignedChildren( alignment_tree[cur_node].children[1], node2_seqs );
		vector< bitset_t > m_aln(2, bitset_t( aln_matrix[0].size(), false ) );
		gnSeqI seqI_len = 0;
		gnSeqI seqJ_len = 0;
		gnSeqI cur_col = 0;
		for( size_t colI = 0; colI < aln_matrix[0].size(); colI++ )
		{
			uint seqI = 0;
			uint seqJ = 0;
			for( ; seqI < node1_seqs.size(); ++seqI )
				if( aln_matrix[node_sequence_map[node1_seqs[seqI]]][colI] != '-' )
					break;
			for( ; seqJ < node2_seqs.size(); ++seqJ )
				if( aln_matrix[node_sequence_map[node2_seqs[seqJ]]][colI] != '-' )
					break;

			if( seqI == node1_seqs.size() && seqJ == node2_seqs.size() )
				continue;	// nothing in this column
			if( seqI != node1_seqs.size() )
			{
				seqI_len++;
				m_aln[0].set(cur_col);
			}
			if( seqJ != node2_seqs.size() )
			{
				seqJ_len++;
				m_aln[1].set(cur_col);
			}
			cur_col++;
		}
		m_aln[0].resize(cur_col);
		m_aln[1].resize(cur_col);
		CompactGappedAlignment<> tmp_cga(m_aln.size(), cur_col);
		CompactGappedAlignment<>* cga = tmp_cga.Copy();
		cga->SetLeftEnd(0, seqI_len > 0 ? 1 : 0);	// at this point we have no idea where the left end should really be
		cga->SetLeftEnd(1, seqJ_len > 0 ? 1 : 0);
		if( cga->LeftEnd(0) != NO_MATCH )
			cga->SetOrientation(0, alignment_tree[cur_node].ordering[cur.second].reference_iv.Orientation(0));
		if( cga->LeftEnd(1) != NO_MATCH )
			cga->SetOrientation(1, alignment_tree[cur_node].ordering[cur.second].reference_iv.Orientation(1));
		cga->SetLength(seqI_len,0);
		cga->SetLength(seqJ_len,1);
		cga->SetAlignment(m_aln);	// do this afterwords so that it can create the bitcount

		// the alignment may need to be reversed if the aligned parent is reverse
		size_t p_siv = alignment_tree[cur_node].ordering[cur.second].parent_siv;
		bool reverse_me = false;
		if( p_siv != (std::numeric_limits<size_t>::max)() )
		{
			size_t p_node = alignment_tree[cur_node].parents[0];
			int p_child = alignment_tree[p_node].children[0] == cur_node ? 0 : 1;
			if( alignment_tree[p_node].ordering[p_siv].reference_iv.Orientation(p_child) == AbstractMatch::reverse )
				reverse_me = true;
		}
		if( reverse_me )
		{
			cga->Invert();
			if( cga->LeftEnd(0) != NO_MATCH )
				cga->SetOrientation(0, alignment_tree[cur_node].ordering[cur.second].reference_iv.Orientation(0));
			if( cga->LeftEnd(1) != NO_MATCH )
				cga->SetOrientation(1, alignment_tree[cur_node].ordering[cur.second].reference_iv.Orientation(1));
			propagateInvert( alignment_tree, cur_node, cur.second );
		}

		alignment_tree[cur_node].ordering[cur.second].reference_iv = Interval();
		vector< AbstractMatch* > am_list(1, cga);
		alignment_tree[cur_node].ordering[cur.second].reference_iv.SetMatches( am_list );
		// set these to zero so they don't interfere with coordinate translation
		alignment_tree[cur_node].ordering[cur.second].reference_iv.SetLeftEnd(0, 0);
		alignment_tree[cur_node].ordering[cur.second].reference_iv.SetLeftEnd(1, 0);
	}
}

typedef boost::tuple<CompactGappedAlignment<>*, vector< bitset_t >*, AbstractMatch* > _sort_tracker_type;

template< class CompType >
class CgaBsComp
{
public:
	CgaBsComp( CompType& c ) : comp(c) {};
	bool operator()( const _sort_tracker_type& a, const _sort_tracker_type& b )
	{
		return comp( a.get<0>(), b.get<0>() );
	}
protected:
	CompType& comp;
};

template< typename MatchVector >
void multFilter( MatchVector& matches, uint mult = 2 )
{
	// apply a multiplicity filter
	size_t cur = 0;
	for( size_t mI = 0; mI < matches.size(); ++mI )
	{
		if( matches[mI]->Multiplicity() == mult )
			matches[cur++] = matches[mI];
		else
			matches[mI]->Free();
	}
	matches.erase(matches.begin()+cur, matches.end());
}

bool debugging_cltm = false;
void ProgressiveAligner::constructLcbTrackingMatches( 
	node_id_t ancestral_node, 
	vector< AbstractMatch* >& ancestral_matches, 
	vector< LcbTrackingMatch< AbstractMatch* > >& tracking_matches 
	)
{
	node_id_t child_0 = alignment_tree[ancestral_node].children[0];
	node_id_t child_1 = alignment_tree[ancestral_node].children[1];
	// split up matches at descendant's breakpoints
	propagateDescendantBreakpoints( child_0, 0, ancestral_matches );
	propagateDescendantBreakpoints( child_1, 1, ancestral_matches );

	// store alignment bitvectors for each match...
	vector< bitset_t > bs_tmp(alignment_tree.size());
	vector< vector< bitset_t > > bs(ancestral_matches.size(), bs_tmp);
	vector< _sort_tracker_type > cga_list;
	// initialize alignment bitvectors
	for( size_t mI = 0; mI < ancestral_matches.size(); mI++ )
	{
		vector< bitset_t > aln( alignment_tree.size(), bitset_t(ancestral_matches[mI]->AlignmentLength() ) );
		swap( bs[mI], aln );
		ancestral_matches[mI]->GetAlignment(aln);
		swap( bs[mI][child_0], aln[0] );
		swap( bs[mI][child_1], aln[1] );
		CompactGappedAlignment<> c(alignment_tree.size(),0);
		c.SetLeftEnd(child_0, ancestral_matches[mI]->LeftEnd(0));
		c.SetOrientation(child_0, ancestral_matches[mI]->Orientation(0));
		c.SetLength(ancestral_matches[mI]->Length(0), child_0);
		c.SetLeftEnd(child_1, ancestral_matches[mI]->LeftEnd(1));
		c.SetOrientation(child_1, ancestral_matches[mI]->Orientation(1));
		c.SetLength(ancestral_matches[mI]->Length(1), child_1);
		cga_list.push_back(make_tuple(c.Copy(), &bs[mI], ancestral_matches[mI]));
	}

	stack<node_id_t> node_stack;
	node_stack.push(child_0);
	node_stack.push(child_1);
	while(node_stack.size() > 0)
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() == 0 )
			continue;
		node_stack.push(alignment_tree[cur_node].children[0]);
		node_stack.push(alignment_tree[cur_node].children[1]);

		// do processing for cur_node...
		// 1. determine which interval in the current node each match falls into
		// 2. determine the offset of this match in that interval
		// 3. translate with that interval

		vector< SuperInterval >& siv_list = alignment_tree[cur_node].ordering;
		SingleStartComparator< CompactGappedAlignment<> > ssc(cur_node);
		CgaBsComp< SingleStartComparator< CompactGappedAlignment<> > > comp( ssc );
		sort(cga_list.begin(), cga_list.end(), comp);
		size_t mI = 0;
		size_t sivI = 0;
		while( mI < cga_list.size() && sivI < siv_list.size() )
		{
			CompactGappedAlignment<>* cur_match = cga_list[mI].get<0>();
			if( cur_match->Start(cur_node) == 0 )
			{
				mI++;
				continue;	// this one doesn't match in this lineage!!
			}
			if( cur_match->LeftEnd(cur_node) >= siv_list[sivI].LeftEnd() + siv_list[sivI].Length() )
			{
				sivI++;
				continue;
			}

			if( cur_match->LeftEnd(cur_node) + cur_match->Length(cur_node) > 
				siv_list[sivI].LeftEnd() + siv_list[sivI].Length() )
			{
				cerr << "doesn't fit\n";
				cerr << "cga_list[" << mI << "]->LeftEnd(" << cur_node << "): " << cur_match->LeftEnd(cur_node) << endl;
				cerr << "cga_list[" << mI << "]->RightEnd(" << cur_node << "): " << cur_match->RightEnd(cur_node) << endl;
				breakHere();
			}

			// extract the region of the siv matched by the current match
			CompactGappedAlignment<>* siv_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_list[sivI].reference_iv.GetMatches()[0]);
			if( siv_list[sivI].reference_iv.GetMatches().size() > 1 )
				siv_cga = NULL;
			if( siv_cga == NULL )
			{
				CompactGappedAlignment<> tmp_cga;
				siv_cga = tmp_cga.Copy();
				*siv_cga = CompactGappedAlignment<>(siv_list[sivI].reference_iv);
				vector<AbstractMatch*> tmp_matches(1,siv_cga);
				siv_list[sivI].reference_iv.SetMatches(tmp_matches);
			}
			CompactGappedAlignment<> new_cga;
			siv_cga->copyRange(new_cga, cur_match->LeftEnd(cur_node) - siv_list[sivI].LeftEnd(), cur_match->Length(cur_node));
			if( cur_match->Orientation(cur_node) == AbstractMatch::reverse )
				new_cga.Invert();
			if( new_cga.Multiplicity() == 0 )
			{
				cerr << "impossible!  there's no match!\n";
				genome::breakHere();
			}
			// set the leftend in cga_list
			for( uint cur_child = 0; cur_child < 2; cur_child++ )
			{
				node_id_t sweet_child = alignment_tree[cur_node].children[cur_child];
				cur_match->SetLeftEnd(sweet_child, new_cga.LeftEnd(cur_child));
				if( new_cga.LeftEnd(cur_child) != NO_MATCH )
				{
					cur_match->SetOrientation(sweet_child, new_cga.Orientation(cur_child));
					cur_match->SetLength(new_cga.Length(cur_child), sweet_child);
				}
			}

			// prepare a cga for translation
			CompactGappedAlignment<> c(1,(*cga_list[mI].get<1>())[cur_node].size());
			c.SetLeftEnd(0,1);
			c.SetLength((*cga_list[mI].get<1>())[cur_node].count(),0);
			vector<bitset_t> bivouac(1, (*cga_list[mI].get<1>())[cur_node]);
			c.SetAlignment(bivouac);

			// now translate each child
			for( uint cur_child = 0; cur_child < 2; cur_child++ )
			{
				if( new_cga.Orientation(cur_child) == AbstractMatch::undefined )
					continue;
				CompactGappedAlignment<> cga_tmp = new_cga;
				cga_tmp.SetStart(cur_child, 1);
				c.translate(cga_tmp, cur_child, 0, false);
				// adjust for end-gaps
				bitset_t bs = (cga_tmp.GetAlignment())[cur_child];
				bs.resize(c.GetAlignment()[0].size(), false);
				bs <<= c.GetAlignment()[0].find_first();
				node_id_t sweet_child = alignment_tree[cur_node].children[cur_child];
				swap( (*cga_list[mI].get<1>())[sweet_child], bs );
				for( size_t testI = 0; testI < cga_tmp.SeqCount(); ++testI )
				{
					if( ((*cga_list[mI].get<1>())[testI].size() != 0 && (*cga_list[mI].get<1>())[testI].size() != (*cga_list[mI].get<1>())[sweet_child].size() ) )
					{
						cerr << "bj0rk3l\n";
						genome::breakHere();
					}
				}
			}

			debugging_cltm = false;
			mI++;	// advance to the next match
		}
	}
	tracking_matches.resize( cga_list.size() );
	// finally, construct CompactGappedAlignments out of the bitsets
	for( size_t bsI = 0; bsI < cga_list.size(); ++bsI )
	{
		cga_list[bsI].get<0>()->SetAlignment(*cga_list[bsI].get<1>());
		cga_list[bsI].get<0>()->validate();
		TrackingMatch& ltm = tracking_matches[bsI];
		ltm.node_match = cga_list[bsI].get<0>();
		ltm.original_match = cga_list[bsI].get<2>();
		ltm.match_id = bsI;

		bool found_extant = false;
		for( size_t i = 0; i < alignment_tree.size()-1; ++i )
		{
			size_t im = node_sequence_map[i];
			if( im == (std::numeric_limits<size_t>::max)() )
				continue;
			if( ltm.node_match->LeftEnd(i) != NO_MATCH )
				found_extant = true;
		}
		if( !found_extant )
		{
			cout << "orig aln len: " << ltm.original_match->AlignmentLength() << endl;
			cout << "orig lend 0: " << ltm.original_match->Start(0) << endl;
			cout << "orig lend 1: " << ltm.original_match->Start(1) << endl;
			cout << "orig length 0: " << ltm.original_match->Length(0) << endl;
			cout << "orig length 1: " << ltm.original_match->Length(1) << endl;

			cerr << "this is an ungrounded match!!!\n";
			genome::breakHere();
		}
	}
}

size_t countUnrefined( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t ancestor )
{
	stack< node_id_t > node_stack;
	node_stack.push(ancestor);
	size_t unrefined_count = 0;
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() > 0 )
			for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); ++childI )
				node_stack.push( alignment_tree[cur_node].children[childI] );
		if( !alignment_tree[cur_node].refined )
			unrefined_count++;
	}
	return unrefined_count;
}

void markAsRefined( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t ancestor )
{
	stack< node_id_t > node_stack;
	node_stack.push(ancestor);
	size_t refined_count = 0;
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() > 0 )
			for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); ++childI )
				node_stack.push( alignment_tree[cur_node].children[childI] );
		alignment_tree[cur_node].refined = true;
	}
	alignment_tree[ancestor].refined = false;
}

template< class MatchVector >
double GetPairwiseLcbScore( MatchVector& lcb, vector< gnSequence* >& seq_table, SeedOccurrenceList& sol_1, SeedOccurrenceList& sol_2 ){
	double lcb_score = 0;
	typename MatchVector::iterator match_iter = lcb.begin();
	for( ; match_iter != lcb.end(); ++match_iter )
	{
		typedef typename MatchVector::value_type MatchPtrType;
		MatchPtrType m = *match_iter;
		vector< score_t > scores(m->AlignmentLength(), 0);
		vector< string > et;
		GetAlignment(*m, seq_table, et);

		// get substitution/gap score
		computeMatchScores( et[0], et[1], scores );
		computeGapScores( et[0], et[1], scores );
		double m_score = 0;
		for( size_t i = 0; i < scores.size(); ++i )
			if( scores[i] != INVALID_SCORE )
				m_score += scores[i];

		if( !( m_score > -1000000000 && m_score < 1000000000 ) )
		{
			cerr << "scoring error\n";
			genome::breakHere();
		}
		size_t merI = 0;
		size_t merJ = 0;
		double uni_count = 0;
		double uni_score = 0;
		for( size_t colI = 0; colI < m->AlignmentLength(); ++colI )
		{
			if(et[0][colI] != '-' && et[1][colI] != '-' )
			{
				double uni1 = (double)sol_1.getFrequency(m->LeftEnd(0) + merI - 1);
				double uni2 = (double)sol_2.getFrequency(m->LeftEnd(1) + merJ - 1);
				uni_score += 1 / uni1;
				uni_score += 1 / uni2;
				uni_count += 2;
			}
			if(et[0][colI] != '-')
				merI++;
			if(et[1][colI] != '-')
				merJ++;
		}

		double avg_uni = uni_count > 0 ? uni_score / uni_count : 1;	// if there were no aligned columns then m_score should already be very negative
		lcb_score = m_score * avg_uni;
	}
	return lcb_score;
}


void ProgressiveAligner::pairwiseScoreTrackingMatches( 
						std::vector< TrackingMatch >& tracking_matches, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants,
						boost::multi_array< double, 3 >& tm_score_array)
{
	tm_score_array.resize( boost::extents[tracking_matches.size()][node1_descendants.size()][node2_descendants.size()] );
	for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
	{
		TrackingMatch* cur_match = &tracking_matches[mI];
		AbstractMatch* node_match = cur_match->node_match;

		vector<bitset_t> aln_mat;
		node_match->GetAlignment(aln_mat);

		// build a match among extant seqs only
		bitset_t tmp_bs( aln_mat[0].size(), false );
		vector< bitset_t > extant_mat( aln_mat.size(), tmp_bs );
		vector< gnSequence* > extant_seqs( aln_mat.size(), (gnSequence*)NULL );
		CompactGappedAlignment<> extant_cga( aln_mat.size(), aln_mat[0].size() );
		for( size_t nI = 0; nI < node1_descendants.size(); ++nI )
		{
			if( node_sequence_map[node1_descendants[nI]] == (std::numeric_limits<uint>::max)()  ||
				node_match->LeftEnd(node1_descendants[nI]) == NO_MATCH )
				continue;
			extant_mat[node1_descendants[nI]] = aln_mat[ node1_descendants[nI] ];
			extant_seqs[node1_descendants[nI]] = alignment_tree[ node1_descendants[nI] ].sequence;
			extant_cga.SetStart( node1_descendants[nI], node_match->Start( node1_descendants[nI] ) );
			extant_cga.SetLength( node_match->Length( node1_descendants[nI] ), node1_descendants[nI] );
		}
		for( size_t nJ = 0; nJ < node2_descendants.size(); ++nJ )
		{
			if( node_sequence_map[node2_descendants[nJ]] == (std::numeric_limits<uint>::max)()  ||
				node_match->LeftEnd(node2_descendants[nJ]) == NO_MATCH )
				continue;
			extant_mat[ node2_descendants[nJ] ] = aln_mat[ node2_descendants[nJ] ];
			extant_seqs[node2_descendants[nJ]] = alignment_tree[ node2_descendants[nJ] ].sequence;
			extant_cga.SetStart( node2_descendants[nJ], node_match->Start( node2_descendants[nJ] ) );
			extant_cga.SetLength( node_match->Length( node2_descendants[nJ] ), node2_descendants[nJ] );
		}
		extant_cga.SetAlignment( extant_mat );

		vector< std::string > et;
		GetAlignment( extant_cga, extant_seqs, et );
		const vector< bitset_t >& cga_mat = extant_cga.GetAlignment();

		double match_sum = 0;
		vector< score_t > scores( et[0].size() );
		for( size_t nI = 0; nI < node1_descendants.size(); ++nI )
		{
			if( node_sequence_map[node1_descendants[nI]] == (std::numeric_limits<uint>::max)()  ||
				node_match->LeftEnd(node1_descendants[nI]) == NO_MATCH )
				continue;
			for( size_t nJ = 0; nJ < node2_descendants.size(); ++nJ )
			{
				if( node_sequence_map[node2_descendants[nJ]] == (std::numeric_limits<uint>::max)() ||
					node_match->LeftEnd(node2_descendants[nJ]) == NO_MATCH )
					continue;	// not extant or no match between this pair

				double uni_score = 0;
				double uni_count = 0;
				node_id_t cur_n1 = node1_descendants[nI];
				node_id_t cur_n2 = node2_descendants[nJ];

				// get substitution/gap score
				std::fill( scores.begin(), scores.end(), 0 );
				computeMatchScores( et[cur_n1], et[cur_n2], scores );
				computeGapScores( et[cur_n1], et[cur_n2], scores );
				double m_score = 0;
				for( size_t i = 0; i < scores.size(); ++i )
					if( scores[i] != INVALID_SCORE )
						m_score += scores[i];

				if( !( m_score > -1000000000 && m_score < 1000000000 ) )
				{
					cerr << "scoring error\n";
					genome::breakHere();
				}
				size_t merI = 0;
				size_t merJ = 0;
				for( size_t colI = 0; colI < node_match->AlignmentLength(); ++colI )
				{
					if(cga_mat[cur_n1].test(colI) && cga_mat[cur_n2].test(colI) )
					{
						double uni1 = (double)sol_list[node_sequence_map[cur_n1]].getFrequency(extant_cga.LeftEnd(cur_n1) + merI - 1);
						double uni2 = (double)sol_list[node_sequence_map[cur_n2]].getFrequency(extant_cga.LeftEnd(cur_n2) + merJ - 1);
						uni_score += 1 / uni1;
						uni_score += 1 / uni2;
						uni_count += 2;
					}
					if(cga_mat[cur_n1].test(colI))
						merI++;
					if(cga_mat[cur_n2].test(colI))
						merJ++;
				}

				double avg_uni = uni_count > 0 ? uni_score / uni_count : 1;	// if there were no aligned columns then m_score should already be very negative
				tm_score_array[mI][nI][nJ] = m_score * avg_uni;
				if( !( m_score > -1000000000 && m_score < 1000000000 ) )
				{
					cerr << "scoring error\n";
					genome::breakHere();
				}
			}
		}
	}
	computeAvgAncestralMatchScores(tracking_matches, node1_descendants, node2_descendants, tm_score_array);
}

void ProgressiveAligner::computeAvgAncestralMatchScores( 
						std::vector< TrackingMatch >& tracking_matches, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants,
						boost::multi_array< double, 3 >& tm_score_array)
{
	// now build up the consensus (ancestral) match scores and bp distances
	for( uint nodeI = 0; nodeI < node1_descendants.size(); nodeI++ )
	{
		for( uint nodeJ = 0; nodeJ < node2_descendants.size(); nodeJ++ )
		{
			node_id_t n1 = node1_descendants[nodeI];
			node_id_t n2 = node2_descendants[nodeJ];

			vector<node_id_t> n1_ext;
			vector<node_id_t> n2_ext;
			getAlignedChildren(n1, n1_ext);
			getAlignedChildren(n2, n2_ext);
			if( n1_ext.size() == 1 && n2_ext.size() == 1 )
				continue;	// this node has two extant nodes below it and was already scored

			// map the nodes in n1_ext to their indices in n1_descendants
			vector< node_id_t > n1_ext_map(n1_ext.size());
			for( size_t i = 0; i < n1_ext.size(); ++i )
			{
				vector< node_id_t >::iterator iter = std::find( node1_descendants.begin(), node1_descendants.end(), n1_ext[i] );
				n1_ext_map[i] = iter - node1_descendants.begin();
			}
			vector< node_id_t > n2_ext_map(n2_ext.size());
			for( size_t i = 0; i < n2_ext.size(); ++i )
			{
				vector< node_id_t >::iterator iter = std::find( node2_descendants.begin(), node2_descendants.end(), n2_ext[i] );
				n2_ext_map[i] = iter - node2_descendants.begin();
			}

			// compute scores for all matches at this node
			for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
			{
				uint tally = 0;
				double score_sum = 0;
				for( size_t i = 0; i < n1_ext.size(); ++i )
				{
					if( tracking_matches[mI].node_match->LeftEnd(n1_ext[i]) == NO_MATCH )
						continue;
					for( size_t j = 0; j < n2_ext.size(); ++j )
					{
						if( tracking_matches[mI].node_match->LeftEnd(n2_ext[j]) == NO_MATCH )
							continue;
						++tally;
						score_sum += tm_score_array[mI][n1_ext_map[i]][n2_ext_map[j]];
					}
				}
				if( tally > 0 )
					tm_score_array[mI][nodeI][nodeJ] = score_sum / (double)tally;
			}
		}
	}
}

void initTrackingMatchLCBTracking( 
	const std::vector< TrackingMatch >& tracking_matches, 
	size_t n1_count, 
	size_t n2_count, 
	boost::multi_array< size_t, 3 >& tm_lcb_id_array )
{
	tm_lcb_id_array.resize( boost::extents[tracking_matches.size()][n1_count][n2_count] );
	for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
	{
		for( size_t nI = 0; nI < n1_count; ++nI )
			for( size_t nJ = 0; nJ < n2_count; ++nJ )
				tm_lcb_id_array[mI][nI][nJ] = LCB_UNASSIGNED;
	}
}

void ProgressiveAligner::computeInternalNodeDistances( 
						boost::multi_array<double, 2>& bp_dist_mat, 
						boost::multi_array<double, 2>& cons_dist_mat, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants)
{
	// bp distances for the current node.
	bp_dist_mat.resize(boost::extents[node1_descendants.size()][node2_descendants.size()]);
	cons_dist_mat.resize(boost::extents[node1_descendants.size()][node2_descendants.size()]);
	for( size_t nI = 0; nI < node1_descendants.size(); ++nI )
	{
		if( node_sequence_map[node1_descendants[nI]] == (std::numeric_limits<uint>::max)() )
			continue;
		for( size_t nJ = 0; nJ < node2_descendants.size(); ++nJ )
		{
			if( node_sequence_map[node2_descendants[nJ]] == (std::numeric_limits<uint>::max)() )
				continue;
			size_t i = node_sequence_map[node1_descendants[nI]];
			size_t j = node_sequence_map[node2_descendants[nJ]];
			bp_dist_mat[nI][nJ] = this->bp_distance[i][j];
			cons_dist_mat[nI][nJ] = this->conservation_distance[i][j];
		}
	}

	// now build up the consensus (ancestral) bp distances
	for( uint nodeI = 0; nodeI < node1_descendants.size(); nodeI++ )
	{
		for( uint nodeJ = 0; nodeJ < node2_descendants.size(); nodeJ++ )
		{
			node_id_t n1 = node1_descendants[nodeI];
			node_id_t n2 = node2_descendants[nodeJ];

			vector<node_id_t> n1_ext;
			vector<node_id_t> n2_ext;
			getAlignedChildren(n1, n1_ext);
			getAlignedChildren(n2, n2_ext);
			if( n1_ext.size() == 1 && n2_ext.size() == 1 )
				continue;	// this node has two extant nodes below it, so already has a dist

			// map the nodes in n1_ext to their indices in n1_descendants
			vector< node_id_t > n1_ext_map(n1_ext.size());
			for( size_t i = 0; i < n1_ext.size(); ++i )
			{
				vector< node_id_t >::iterator iter = std::find( node1_descendants.begin(), node1_descendants.end(), n1_ext[i] );
				n1_ext_map[i] = iter - node1_descendants.begin();
			}
			vector< node_id_t > n2_ext_map(n2_ext.size());
			for( size_t i = 0; i < n2_ext.size(); ++i )
			{
				vector< node_id_t >::iterator iter = std::find( node2_descendants.begin(), node2_descendants.end(), n2_ext[i] );
				n2_ext_map[i] = iter - node2_descendants.begin();
			}

			// compute average bp distance
			for( size_t i = 0; i < n1_ext.size(); ++i )
			{
				for( size_t j = 0; j < n2_ext.size(); ++j )
				{
					bp_dist_mat[nodeI][nodeJ] += bp_dist_mat[n1_ext_map[i]][n2_ext_map[j]];
					cons_dist_mat[nodeI][nodeJ] += cons_dist_mat[n1_ext_map[i]][n2_ext_map[j]];
				}
			}
			bp_dist_mat[nodeI][nodeJ] /= (double)(n1_ext.size() * n2_ext.size());
			cons_dist_mat[nodeI][nodeJ] /= (double)(n1_ext.size() * n2_ext.size());
		}
	}

}

double computeID( GappedAlignment& gal, size_t seqI, size_t seqJ )
{
	const vector< string >& aln_mat = GetAlignment( gal, vector< gnSequence* >(gal.SeqCount(), NULL ));
	double id = 0;
	double possible = 0;
	for( size_t colI = 0; colI < gal.AlignmentLength(); colI++ )
	{
		if( aln_mat[seqI][colI] == '-' || aln_mat[seqJ][colI] == '-' )
			continue;
		if( toupper(aln_mat[seqI][colI]) == toupper(aln_mat[seqJ][colI]))
			id++;
		possible++;
	}
	return id / possible;
}


// options for reducing total number of pairwise match data structures during the 
// translate to ancestral phase:
// -- for each leaf below node A, call current x
//    - identify all internal nodes below B
//    - call the lowest node y
//    - create a pairwise match for each of x, des(y) multi-matches
//    - translate x, des(y) matches to x, y pairwise matches
//    - eliminate overlaps in x, y matches
//    - pick next descendent of B and call it y, repeat
//  - pick next leaf below A
// Analysis: if we select A and B such that A has fewer leaves then 
//
//
// different option -- just pick a representative from leaf(A) and leaf(B) to translate

/*
void translateToAncestral(  PhyloTree< AlignmentTreeNode >& t, node_id_t node1, node_id_t node2 )
{
	//
	// do a depth first? traversal of the tree starting at both node1 and node2
	// translate up to each cross-pair of internal nodes and eliminate overlaps
	// 
	stack< node_id_t > n1_stack;
	stack< node_id_t > n2_stack;
	n1_stack.push(node1);
	n2_stack.push(node2);
	bitset_t visited( t.size() );
	while( n1_stack.size() > 0 )
	{
		node_id_t cur_n1 = n1_stack.top();

		if( t[cur_n1].children.size() == 0 )
		{
			n1_stack.pop();
			visited[cur_n1].set();
			continue;
		}
		// both are internal nodes
		// do we need to visit children?
		node_id_t n1_c1 = t[cur_n1].children[0];
		node_id_t n1_c2 = t[cur_n1].children[1];
		if( !visited[n1_c1] && !visited[n1_c2] )
		{
			n1_stack.push(n1_c1);
			n1_stack.push(n1_c2);
			continue;
		}else if( !visited[n1_c1] || !visited[n1_c2] )
		{
			cerr << "bad tree topology\n";
		}

		while( n2_stack.size() > 0 )
		{
			node_id_t cur_n2 = n2_stack.top();
			if( t[cur_n2].children.size() == 0 )
			{
				n2_stack.pop();
				visited[cur_n1].set();
				continue;
			}
			node_id_t n2_c1 = t[cur_n2].children[0];
			node_id_t n2_c2 = t[cur_n2].children[1];
			if( !visited[n2_c1] && !visited[n2_c2] )
			{
				n1_stack.push(n2_c1);
				n1_stack.push(n2_c2);
				continue;
			}else if( !visited[n2_c1] || !visited[n2_c2] )
			{
				cerr << "bad tree topology n2\n";
			}

			// all children have been visited, ready to translate up the tree
			vector< AbstractMatch* > am_list( pairwise_matches(seqI, seqJ).begin(), pairwise_matches(seqI, seqJ).end() );
			pairwise_matches(seqI, seqJ).clear();
			translateGappedCoordinates( am_list, 1, node_sequence_map[seqJ], node2 );
			translateGappedCoordinates( am_list, 0, node1_seqs[seqI], node1 );
			ancestral_matches.insert( ancestral_matches.end(), am_list.begin(), am_list.end() );
		}
	}
}
*/

void ProgressiveAligner::getRepresentativeAncestralMatches( const vector< node_id_t > node1_seqs, const vector< node_id_t > node2_seqs, node_id_t node1, node_id_t node2, node_id_t ancestor, std::vector< AbstractMatch* >& ancestral_matches )
{
	// for each match, extract a representative match from any pair of genomes in node1_seqs and node2_seqs
	// translate up the resulting set of matches and eliminate overlaps
	vector< AbstractMatch* > cur_matches;
	boost::multi_array< vector< AbstractMatch* >, 2 > seq_matches( boost::extents[node1_seqs.size()][node2_seqs.size()] );
	for( size_t mI = 0; mI < original_ml.size(); mI++ )
	{
		for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
		{
			uint ii = this->node_sequence_map[node1_seqs[seqI]];
			if( original_ml[mI]->LeftEnd(ii) == NO_MATCH )
				continue;

			for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			{
				uint jj = this->node_sequence_map[node2_seqs[seqJ]];
				if( original_ml[mI]->LeftEnd(jj) == NO_MATCH )
					continue;
				Match mm( 2 );
				Match* new_m = mm.Copy();
				new_m->SetStart( 0, original_ml[mI]->Start(ii));
				new_m->SetStart( 1, original_ml[mI]->Start(jj));
				new_m->SetLength(original_ml[mI]->Length());
				if( new_m->Start(0) < 0 )
					new_m->Invert();	// assign reference orientation to seq 0
				seq_matches[seqI][seqJ].push_back( new_m );
				break;
			}
			break;
		}
	}
	for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
		for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
		{
			translateGappedCoordinates( seq_matches[seqI][seqJ], 0, node1_seqs[seqI], node1 );
			translateGappedCoordinates( seq_matches[seqI][seqJ], 1, node2_seqs[seqJ], node2 );
			ancestral_matches.insert( ancestral_matches.end(), seq_matches[seqI][seqJ].begin(), seq_matches[seqI][seqJ].end() );
		}

	EliminateOverlaps_v2( ancestral_matches );
}

#ifdef WIN32
#include <windows.h>
#include <PSAPI.h>
void printMemUsage()
{
	DWORD proclist[500];
	DWORD cbNeeded;
	BOOL rval;
	rval = EnumProcesses( proclist, sizeof(proclist), &cbNeeded );
	int p_count = cbNeeded / sizeof(DWORD);
	HANDLE phand;
	HMODULE hMod;
	char szFileName[MAX_PATH];
	for( int p = 0; p < p_count; p++ )
	{
		phand = OpenProcess( PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, 0, proclist[p] );
		DWORD dwSize2;
		if (EnumProcessModules(phand, &hMod, sizeof(hMod), &dwSize2)) 
		{

			// Get the module name
			if ( !GetModuleBaseName(phand, hMod, szFileName, sizeof(szFileName)) )
				szFileName[0] = 0;
			if( strncmp( szFileName, "progressiveMauve", 16 ) == 0 )
				break;	// found the right module
		}
		CloseHandle(phand);
	}

	PROCESS_MEMORY_COUNTERS mem_info;

	if( GetProcessMemoryInfo( phand, &mem_info, sizeof(mem_info) ) )
	{
		cerr << "Working set size: " << mem_info.WorkingSetSize / (1024 * 1024) << " Mb\n";
//		cerr << "Paged pool usage: " << mem_info.QuotaPagedPoolUsage << endl;
//		cerr << "Non-Paged pool usage: " << mem_info.QuotaNonPagedPoolUsage << endl;
		cerr << "Pagefile usage: " << mem_info.PagefileUsage / (1024 * 1024) << " Mb\n";
	}
}
#else
void printMemUsage()
{};
#endif

void ProgressiveAligner::alignProfileToProfile( node_id_t node1, node_id_t node2, node_id_t ancestor )
{
	// 1) find all pairwise matches
	// 2) convert to pairwise matches among the ancestral sequences
	//    - delete inconsistently aligned regions?
	// 3) perform greedy b.p. elimination on the pairwise matches
	// 4) extend LCBs
	// 5)  if total alignment weight hasn't changed, go to (8)
	// 6) search for additional matches between each match among extant sequences
	// 7) go back to 2
	// 8) perform a MUSCLE/Clustal alignment of each intervening region

	vector< node_id_t > node1_seqs;	/**< the node id's of extant sequences below node 1 */
	vector< node_id_t > node2_seqs;	/**< the node id's of extant sequences below node 2 */
	getAlignedChildren( node1, node1_seqs );
	getAlignedChildren( node2, node2_seqs );
	uint seqI, seqJ;
	gnSeqI prev_ancestral_seq_len = (std::numeric_limits<gnSeqI>::max)();

	printMemUsage();
	cerr << "get ancestral matches\n";

	Matrix<MatchList> pairwise_matches( node1_seqs.size(), node2_seqs.size() );
//	getPairwiseMatches( node1_seqs, node2_seqs, pairwise_matches );
	vector< AbstractMatch* > anc_pairwise_matches;
	getRepresentativeAncestralMatches( node1_seqs, node2_seqs, node1, node2, ancestor, anc_pairwise_matches );
	printMemUsage();
	
	PhyloTree< AlignmentTreeNode > aln_tree_backup;

	/** A cache of regions that were searched in the previous round of recursion */
	Matrix< std::vector< search_cache_t > > search_cache_db(node1_seqs.size(), node2_seqs.size());
	double prev_anchoring_score = -(std::numeric_limits<double>::max)();
	double cur_anchoring_score = -(std::numeric_limits<double>::max)();

	while(true)
	{
		vector<AbstractMatch*> ancestral_matches;
		if( anc_pairwise_matches.size() > 0 )
		{
			ancestral_matches.insert( ancestral_matches.begin(), anc_pairwise_matches.begin(), anc_pairwise_matches.end() );
			anc_pairwise_matches.clear();
		}

		// part 2, construct pairwise matches to the ancestral sequence
		// A)  for each pairwise match, translate its
		//     coordinates to the ancestral genome
		//	   -- try to use translateCoordinates
		//     -- build a translation table for translateCoordinates

		for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
		{
			for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			{
				cout << node_sequence_map[node1_seqs[seqI]] << "," << node_sequence_map[node2_seqs[seqJ]] << " has " << pairwise_matches(seqI,seqJ).size() << " pairwise matches\n";
				cout.flush();

				vector< AbstractMatch* > am_list( pairwise_matches(seqI, seqJ).begin(), pairwise_matches(seqI, seqJ).end() );
				pairwise_matches(seqI, seqJ).clear();
				translateGappedCoordinates( am_list, 1, node2_seqs[seqJ], node2 );
				translateGappedCoordinates( am_list, 0, node1_seqs[seqI], node1 );
				ancestral_matches.insert( ancestral_matches.end(), am_list.begin(), am_list.end() );
			}
		}

		// include any matches from a previous iteration of this loop
		for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
		{
			Interval& ref_iv = alignment_tree[ancestor].ordering[aI].reference_iv;
			if( ref_iv.Multiplicity() == 2 )
				for( size_t mI = 0; mI < ref_iv.GetMatches().size(); mI++ )
					if( ref_iv.GetMatches()[mI]->Multiplicity() > 1 )
						ancestral_matches.push_back( ref_iv.GetMatches()[mI]->Copy() );
		}

		// set seq 0 to forward ref. orientation
		for( size_t mI = 0; mI < ancestral_matches.size(); ++mI )
			if( ancestral_matches[mI]->Start(0) < 0 )
				ancestral_matches[mI]->Invert();

		// eliminate overlaps as they correspond to inconsistently or
		// multiply aligned regions
		EliminateOverlaps_v2( ancestral_matches );
		
		multFilter( ancestral_matches );

		vector< vector< AbstractMatch* > > LCB_list;
		vector< LCB > adjacencies;
		vector< gnSeqI > breakpoints;

		if( !collinear_genomes )
		{
			cout << "Performing Sum-of-pairs Greedy Breakpoint Elimination\n";
			cout.flush();
			// project the pairwise matches at this node to all-possible pairs matches at descendant nodes
			// keep a mapping of ancestral to extant matches so that when an ancestral match gets removed
			// the match among extant nodes also gets removed
			// how should candidate matches to remove be generated?
			// one possibility is to remove entire ancestral LCBs...  this may be problematic since ancestral
			// LCBs don't correspond to the pairwise LCBs thus an ancestral LCB could be removed with no useful
			// change in alignment score
			//
			//
			// translate the matches into LcbTrackingMatches
			printMemUsage();
			cerr << "construct LCB tracking matches\n";
			vector< TrackingMatch > tracking_matches;
			boost::multi_array< size_t, 3 > tm_lcb_id_array;
			boost::multi_array< double, 3 > tm_score_array;
			constructLcbTrackingMatches( ancestor, ancestral_matches, tracking_matches );

			vector<node_id_t> node1_descendants;
			vector<node_id_t> node2_descendants;
			if( scoring_scheme == ExtantSumOfPairsScoring )
			{
				node1_descendants = node1_seqs;
				node2_descendants = node2_seqs;
			}else{
				getDescendants(alignment_tree, node1, node1_descendants);
				getDescendants(alignment_tree, node2, node2_descendants);
			}

			//
			// score the matches
			//
			printMemUsage();
			cerr << "init tracking match LCB tracking\n";
			initTrackingMatchLCBTracking( tracking_matches, node1_descendants.size(), node2_descendants.size(), tm_lcb_id_array );
			printMemUsage();
			cerr << "pairwise score tracking matches\n";
			pairwiseScoreTrackingMatches( tracking_matches, node1_descendants, node2_descendants, tm_score_array );
			printMemUsage();

			// compute bp distances for the current node.
			// ancestral nodes take the average distance of extant nodes
			boost::multi_array<double, 2> bp_dist_mat;
			boost::multi_array<double, 2> cons_dist_mat;
			computeInternalNodeDistances( bp_dist_mat, cons_dist_mat, node1_descendants, node2_descendants);

			vector< TrackingMatch* > t_matches(tracking_matches.size());
			for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
				t_matches[mI] = &tracking_matches[mI];

			// now sort these out into pairwise LCBs
			cerr << "get pairwise LCBs\n";
			PairwiseLCBMatrix pairwise_adj_mat(boost::extents[node1_descendants.size()][node2_descendants.size()]);
			for( uint nodeI = 0; nodeI < node1_descendants.size(); nodeI++ )
				for( uint nodeJ = 0; nodeJ < node2_descendants.size(); nodeJ++ )
					getPairwiseLCBs( node1_descendants[nodeI], node2_descendants[nodeJ], nodeI, nodeJ, t_matches, pairwise_adj_mat[nodeI][nodeJ], tm_score_array, tm_lcb_id_array );

			printMemUsage();
			sort( t_matches.begin(), t_matches.end() );

			// other possibility, choose pairwise LCBs to remove.  a score improvement is always guaranteed
			// compute LCBs among descendant nodes
			// this is a good idea.  it factors out ancestral breakpoint decisions entirely
			// need a data structure to track all pairwise LCBs that contain a given match
			// template <class MatchType>
			// class LcbTrackingMatch <MatchType> 
			// { 
			// public:
			//  MatchType node_match;
			//	boost::multi_array< size_t, 2 > lcb_id;
			// }
			// all pairwise LCBs would be evaluated for removal and the one that provides the greatest
			// overall score improvement gets removed.
			// upon removal, matches associated with that LCB would get removed, and any LCBs in other 
			// genomes would get removed if they no longer had any matches
			// to pull this off, the LCB struct needs to store the set of matches directly
			// 
			// but what about small cycles that appear only in 3 or more-way comparisons?  are these
			// important?  umm, yeah, but only if you believe in evolution.
			// 
			// so here's the dilly-oh: score against the ancestral ordering(s) *and* all pairwise orderings
			// for an ancestor.  ancestor contributes the sum of all descendants to the score and breakpoints
			// are penalized as the sum of /participating/ descendants.  a descendant is participating
			// if it has some matching region defined within the LCB and if removal of that matching region
			// eliminates a breakpoint in the pairwise comparison
/*
			if( node1 == 0 && node2 == 1 )
			{
				for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
				{
//					if( tracking_matches[mI].original_match->LeftEnd(0) > 220000 && 
//						tracking_matches[mI].original_match->LeftEnd(0) < 310000 )

					{
						cout << "match (" << tracking_matches[mI].original_match->LeftEnd(0) << ", " << tracking_matches[mI].original_match->RightEnd(0) << ")\t";
						cout << "(" << tracking_matches[mI].original_match->LeftEnd(1) << ", " << tracking_matches[mI].original_match->RightEnd(1) << ")\n";
						cout << "match " << mI << " score: " << tracking_matches[mI].score[0][0] << endl;
					}
				}
				cerr << "debugme!\n";
				for( size_t adjI = 0; adjI < pairwise_adj_mat[0][0].size(); ++adjI )
				{
					cout << "LCB " << adjI << " weight: " << pairwise_adj_mat[0][0][adjI].weight << endl;
					cout << "Boundaries in seq 0: " << pairwise_adj_mat[0][0][adjI].left_end[0] << ", " << pairwise_adj_mat[0][0][adjI].right_end[0] << endl;
					cout << "Boundaries in seq 1: " << pairwise_adj_mat[0][0][adjI].left_end[1] << ", " << pairwise_adj_mat[0][0][adjI].right_end[1] << endl;
				}
				cerr << "debugme!\n";

				cerr << "bp_dist_mat:\n";
				print2d_matrix( bp_dist_mat, cerr );
				cerr << "cons_dist_mat:\n";
				print2d_matrix( cons_dist_mat, cerr );
//				debug_aligner = true;
			}
*/
			cerr << "Greedy BPE\n";
			vector< TrackingMatch* > final;
			if(scoring_scheme == AncestralScoring)
			{
				vector<node_id_t>::iterator d1_iter = std::find( node1_descendants.begin(), node1_descendants.end(), node1 );
				vector<node_id_t>::iterator d2_iter = std::find( node2_descendants.begin(), node2_descendants.end(), node2 );
				size_t d1_index = d1_iter - node1_descendants.begin();
				size_t d2_index = d2_iter - node2_descendants.begin();
				EvenFasterSumOfPairsBreakpointScorer spbs( breakpoint_penalty, bp_dist_mat, cons_dist_mat, 
					t_matches, pairwise_adj_mat, node1_descendants, node2_descendants, 
					tm_score_array, tm_lcb_id_array, d1_index, d1_index+1, d2_index, d2_index+1 );
				cur_anchoring_score = greedySearch( spbs );
				final = spbs.getResults();
			}else{
				EvenFasterSumOfPairsBreakpointScorer spbs( breakpoint_penalty, bp_dist_mat, cons_dist_mat, 
					t_matches, pairwise_adj_mat, node1_descendants, node2_descendants, 
					tm_score_array, tm_lcb_id_array, 0, node1_descendants.size(), 0, node2_descendants.size() );
				cur_anchoring_score = greedySearch( spbs );
				final = spbs.getResults();
			}
			cout << "done\n";

			// free memory used by pairwise projections
			for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
				tracking_matches[mI].node_match->Free();

			ancestral_matches.clear();

			// free memory from deleted matches here
			std::sort(final.begin(), final.end());
			vector< TrackingMatch* > deleted_t_matches( t_matches.size(), NULL );
			std::set_difference( t_matches.begin(), t_matches.end(), final.begin(), final.end(), deleted_t_matches.begin() );
			for( size_t delI = 0; delI < deleted_t_matches.size(); ++delI )
			{
				if( deleted_t_matches[delI] == NULL )
					break;
				deleted_t_matches[delI]->original_match->Free();
			}

			// convert back to an LCB list
			vector< AbstractMatch* > new_matches(final.size());
			for( size_t mI = 0; mI < final.size(); ++mI )
				new_matches[mI] = final[mI]->original_match;

			IdentifyBreakpoints( new_matches, breakpoints );
			ComputeLCBs_v2( new_matches, breakpoints, LCB_list );

		} // end if !collinear
		else
		{	// if we are assuming all genomes are collinear, then we don't need the 
			// sophisticated pairwise breakpoint scoring and can get by with simple breakpoint
			// penalties
			IdentifyBreakpoints( ancestral_matches, breakpoints );
			ComputeLCBs_v2( ancestral_matches, breakpoints, LCB_list );

			vector< double > lcb_scores( LCB_list.size() );
			double score_sum = 100;	// anything > 0 would work.  this will be the breakpoint penalty
			for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
			{
				lcb_scores[lcbI] = SimpleGetLCBCoverage( LCB_list[lcbI] );
				score_sum += lcb_scores[lcbI];
			}

			computeLCBAdjacencies_v3( LCB_list, lcb_scores, adjacencies );

			// want to eliminate all breakpoints
			SimpleBreakpointScorer wbs( adjacencies, score_sum );
			cur_min_coverage = greedyBreakpointElimination_v4( adjacencies, lcb_scores, wbs, NULL, false );
			vector<AbstractMatch*> deleted_matches;
			filterMatches_v2( adjacencies, LCB_list, lcb_scores, deleted_matches );
			for( size_t delI = 0; delI < deleted_matches.size(); ++delI )
				deleted_matches[delI]->Free();
		}
		printMemUsage();

		ancestral_matches.clear();

		cout << "Arrived at " << LCB_list.size() << " intervals\n";
		// create an ancestral ordering
		vector< Interval* > pairwise_intervals;
		Interval tmp_iv;
		for( size_t lcbI = 0; lcbI < LCB_list.size(); lcbI++ )
		{
			pairwise_intervals.push_back( tmp_iv.Copy() );
			pairwise_intervals.back()->SetMatches( LCB_list[lcbI] );
		}
		LCB_list.clear();

		vector<gnSeqI> seq_lengths = vector<gnSeqI>(2,0);
		for( size_t aI = 0; aI < alignment_tree[node1].ordering.size(); ++aI )
			seq_lengths[0] += alignment_tree[node1].ordering[aI].Length();
		for( size_t aI = 0; aI < alignment_tree[node2].ordering.size(); ++aI )
			seq_lengths[1] += alignment_tree[node2].ordering[aI].Length();

		cout << "Adding unaligned intervals\n";
		addUnalignedIntervals_v2(pairwise_intervals, set<uint>(), seq_lengths);

		cout << "addUnalignedIntervals yields " << pairwise_intervals.size() << " intervals\n";

		bool borked = false;
		if(debug_aligner)
			borked = validatePairwiseIntervals(node1, node2, pairwise_intervals);

		// merge unaligned intervals
		cout << "Merging unaligned intervals\n";
		cout.flush();
		vector<Interval*> new_list1;
		vector<Interval*> merged_intervals;
		mergeUnalignedIntervals( 1, pairwise_intervals, new_list1 );
		mergeUnalignedIntervals( 0, new_list1, merged_intervals );
		cout << "Marbling gaps\n";
		cout.flush();
		for( size_t ivI = 0; ivI < merged_intervals.size(); ivI++ )
			merged_intervals[ivI]->Marble(50);

		cout << "Propagating descendant breakpoints\n";

		// split up intervals at descendant's breakpoints
		propagateDescendantBreakpoints( node1, 0, merged_intervals );
		propagateDescendantBreakpoints( node2, 1, merged_intervals );

		cout << "descendant 0(" << node1 << ") has " << alignment_tree[node1].ordering.size() << " intervals\n";
		cout << "descendant 1(" << node2 << ") has " << alignment_tree[node2].ordering.size() << " intervals\n";
		cout << "propagateDescendantBreakpoints yields " << merged_intervals.size() << " intervals\n";

		if(debug_aligner)
			borked = validatePairwiseIntervals(node1, node2, merged_intervals);
		cout << "Creating ancestral ordering\n";
		alignment_tree[ancestor].ordering.clear();
		createAncestralOrdering( merged_intervals, alignment_tree[ancestor].ordering );
		for( size_t ivI = 0; ivI < merged_intervals.size(); ivI++ )
			merged_intervals[ivI]->Free();
		merged_intervals.clear();	// free up some memory

		if(debug_aligner)
			validateSuperIntervals( node1, node2, ancestor );

		// if we're not making any progress then bail out...
		gnSeqI cur_ancestral_seq_len = 0;
		for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
			cur_ancestral_seq_len += alignment_tree[ancestor].ordering[aI].Length();

		if( !collinear_genomes )
			cout << "Previous anchoring score: " << prev_anchoring_score << ", new anchor score: " << cur_anchoring_score << endl;
		else
			cout << "Prev alignment len: " << prev_ancestral_seq_len << ", new alignment length: " << cur_ancestral_seq_len << endl;
		// if cur_seq_len has decreased then we're improving
		// if not, then we're done finding matches
		if( collinear_genomes && cur_ancestral_seq_len >= prev_ancestral_seq_len )
			break;

		if( !collinear_genomes && cur_anchoring_score <= prev_anchoring_score )
			break;
		prev_anchoring_score = cur_anchoring_score;
		prev_ancestral_seq_len = cur_ancestral_seq_len;

		// accept the new alignment tree...
		cout << "Backing up alignment tree...\n";
		cout.flush();
		aln_tree_backup = alignment_tree;

		cout << "propagating ancestral breakpoints\n";
		cout.flush();
		recursiveApplyAncestralBreakpoints(ancestor);


		if( debug_me ) 	 
		{
			for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ ) 	 
			{
				GappedAlignment gal; 	 
				extractAlignment(ancestor, aI, gal); 	 

				bool check = false;
				for( size_t ii = 0; ii < gal.SeqCount(); ++ii )
				{
					if( gal.LeftEnd(ii) == 0 )
						continue;
					for( size_t jj = 0; jj < gal.SeqCount(); ++jj )
					{
						if( gal.LeftEnd(jj) == 0 )
							continue;
						check = check || computeID( gal, ii, jj ) < .5;
					}
				}
				if( check )
					cerr << "check iv " << aI << " dbg_count " << dbg_count << endl;
				else
					continue;

				const vector< string >& aln_mat = GetAlignment(gal, this->original_ml.seq_table); 	 
				gnSequence seq; 	 
				for( size_t seqI = 0; seqI < gal.SeqCount(); ++seqI ) 	 
					if( gal.LeftEnd(seqI) != NO_MATCH ) 	 
						seq += aln_mat[seqI]; 	 

				stringstream dbg_fname; 	 
				dbg_fname << "prof_dbg_iv_" << aI << ".dbg." << dbg_count++ << ".fas"; 	 
				ofstream debug_file( dbg_fname.str().c_str() ); 	 
				gnFASSource::Write( seq, debug_file, false ); 	 
				debug_file.close(); 	 
			} 	 
		}

		if(debug_aligner)
			validateSuperIntervals( node1, node2, ancestor );

		// search for additional alignment anchors
		cout << "recursive anchor search\n";
		cout.flush();
		Matrix<MatchList> matches;
		Matrix< std::vector< search_cache_t > > new_cache_db(node1_seqs.size(), node2_seqs.size());
		for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
		{
			CompactGappedAlignment<> cga;
			extractAlignment(ancestor, aI, cga);
			recurseOnPairs(node1_seqs, node2_seqs, cga, matches, search_cache_db, new_cache_db);

			// add any new matches to the pairwise_matches matrix
			for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
				for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
					pairwise_matches(seqI, seqJ).insert( pairwise_matches(seqI, seqJ).end(), matches(seqI, seqJ).begin(), matches(seqI, seqJ).end() );
		}
		

		for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
		{
			for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			{
				// delete the previous search cache
				swap( search_cache_db(seqI, seqJ), new_cache_db(seqI, seqJ) );
				for( size_t mI = 0; mI < new_cache_db(seqI,seqJ).size(); mI++ )
				{
					if( new_cache_db(seqI,seqJ)[mI].first != NULL )
						new_cache_db(seqI,seqJ)[mI].first->Free();
					if( new_cache_db(seqI,seqJ)[mI].second != NULL )
						new_cache_db(seqI,seqJ)[mI].second->Free();
				}
				new_cache_db(seqI,seqJ).clear();
				try{
					std::sort( search_cache_db(seqI, seqJ).begin(), search_cache_db(seqI, seqJ).end(), cache_comparator );
				}catch(...){
					cerr << "Error sorting.\n";
					cerr << "Cache has " << search_cache_db(seqI, seqJ).size() << " elements\n\n\n";
					for( size_t i = 0; i < search_cache_db(seqI, seqJ).size(); ++i )
						cerr << search_cache_db(seqI, seqJ)[i].first << ",\t" << search_cache_db(seqI, seqJ)[i].second << endl;
				}
				if( pairwise_matches(seqI,seqJ).size() > 0 )
					cout << seqI << "," << seqJ << " has an additional " << pairwise_matches(seqI,seqJ).size() << " matches\n";
			}
		}

		// restore backed up tree since we only want the final set of ancestral
		// breakpoints applied to the descendants
		cout << "Restoring backed up alignment tree...\n";
		cout.flush();
		alignment_tree = aln_tree_backup;

	}	// end while(true)

	// delete the search cache
	for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
		for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			for( size_t mI = 0; mI < search_cache_db(seqI,seqJ).size(); mI++ )
			{
				if( search_cache_db(seqI,seqJ)[mI].first != NULL )
					search_cache_db(seqI,seqJ)[mI].first->Free();
				if( search_cache_db(seqI,seqJ)[mI].second != NULL )
					search_cache_db(seqI,seqJ)[mI].second->Free();
			}

	// aln_tree_backup has the highest scoring alignment_tree
	alignment_tree = aln_tree_backup;
	cout << "propagating ancestral breakpoints\n";
	recursiveApplyAncestralBreakpoints(ancestor);

	// step 8) construct a muscle alignment in each intervening region
	if( gapped_alignment )
	{
		cout << "performing a gapped alignment\n";
		doGappedAlignment(ancestor, true);
	}else
		cerr << "skipping gapped alignment\n";
	if( refine )
	{
		size_t unrefined = countUnrefined( alignment_tree, ancestor );
		if( unrefined > 5 && ancestor != alignment_tree.root )
		{
			cout << "performing iterative refinement\n";
			doGappedAlignment(ancestor, false);
			markAsRefined( alignment_tree, ancestor );
		}
	}


	if( debug_me ) 	 
	{
		for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ ) 	 
		{ 	 

			static int dbg_count = 0; 	 
			GappedAlignment gal; 	 
			extractAlignment(ancestor, aI, gal); 	 

			bool check = false;
			for( size_t ii = 0; ii < gal.SeqCount(); ++ii )
			{
				if( gal.LeftEnd(ii) == 0 )
					continue;
				for( size_t jj = 0; jj < gal.SeqCount(); ++jj )
				{
					if( gal.LeftEnd(jj) == 0 )
						continue;
					check = check || computeID( gal, ii, jj ) < .5;
				}
			}
			if( check )
				cerr << "check iv " << aI << " dbg_count " << dbg_count << endl;
			else
				continue;

			const vector< string >& aln_mat = GetAlignment(gal, this->original_ml.seq_table); 	 
			gnSequence seq; 	 
			for( size_t seqI = 0; seqI < gal.SeqCount(); ++seqI ) 	 
				if( gal.LeftEnd(seqI) != NO_MATCH ) 	 
					seq += aln_mat[seqI]; 	 

			stringstream dbg_fname; 	 
			dbg_fname << "prof_dbg_iv_" << aI << ".dbg." << dbg_count++ << ".fas"; 	 
			ofstream debug_file( dbg_fname.str().c_str() ); 	 
			gnFASSource::Write( seq, debug_file, false ); 	 
			debug_file.close(); 	 
		} 	 
	}


}


SuperInterval::SuperInterval() :
length(0),
left_end(0),
c1_siv((std::numeric_limits<size_t>::max)()),
c2_siv((std::numeric_limits<size_t>::max)()),
parent_siv((std::numeric_limits<size_t>::max)())
{}

SuperInterval::SuperInterval( const Interval& reference_iv ) :
reference_iv(reference_iv),
length(0),
left_end(0),
c1_siv((std::numeric_limits<size_t>::max)()),
c2_siv((std::numeric_limits<size_t>::max)()),
parent_siv((std::numeric_limits<size_t>::max)())
{
}

SuperInterval::SuperInterval(const SuperInterval& siv) :
left_end(siv.left_end),
length( siv.length ),
reference_iv( siv.reference_iv ),
c1_siv(siv.c1_siv),
c2_siv(siv.c2_siv),
parent_siv(siv.parent_siv)
{
}
SuperInterval& SuperInterval::operator=(const SuperInterval& siv)
{
	left_end = siv.left_end;
	length = siv.length;
	reference_iv = siv.reference_iv;
	c1_siv = siv.c1_siv;
	c2_siv = siv.c2_siv;
	parent_siv = siv.parent_siv;
	return *this;
}



/** Sets the length of this match to @param len */
void SuperInterval::SetLength( gnSeqI len )
{
	length = len;
}

void addGuy( uint seqI, AbstractMatch::orientation orient, 
			std::vector< AbstractMatch* >& new_ivs, 
			vector<Interval*>& new_list )
{
	Interval tmp_iv;
	// set the orientation for any unaligned intervals
	if( orient == AbstractMatch::reverse )
	{
		for( size_t nI = 0; nI < new_ivs.size(); nI++ )
			if( new_ivs[nI]->LeftEnd(seqI) != NO_MATCH && new_ivs[nI]->Orientation(seqI) != orient)
				new_ivs[nI]->Invert();
	}
	// add this guy
	Interval* added_iv = tmp_iv.Copy();
	added_iv->SetMatches( new_ivs );
	new_list.push_back(added_iv);
}

void mergeUnalignedIntervals( uint seqI, vector< Interval* >& iv_list, vector< Interval* >& new_list )
{
	SSC<Interval> ivlcJ(seqI);
	sort( iv_list.begin(), iv_list.end(), ivlcJ );

	Interval tmp_iv;
	AbstractMatch::orientation orient = AbstractMatch::undefined;
	vector< AbstractMatch* > new_ivs;
	vector< Interval* > to_delete;
	for( size_t ordI = 0; ordI < iv_list.size(); ordI++ )
	{
		if( iv_list[ordI]->LeftEnd(seqI) == NO_MATCH )
		{
			new_list.push_back(iv_list[ordI]);
			iv_list[ordI] = NULL;
			continue;
		}

		if( orient == AbstractMatch::undefined && iv_list[ordI]->Multiplicity() == 2 )
		{
			orient = iv_list[ordI]->Orientation(seqI);
			vector< AbstractMatch* > matches;
			iv_list[ordI]->StealMatches( matches );
			if( orient == AbstractMatch::forward )
				new_ivs.insert( new_ivs.end(), matches.begin(), matches.end() );
			else
				new_ivs.insert( new_ivs.begin(), matches.begin(), matches.end() );

			// if it's the last one then add
			if( ordI + 1 == iv_list.size() )
				addGuy( seqI, orient, new_ivs, new_list );
			continue;
		}
		if( orient != AbstractMatch::undefined && iv_list[ordI]->Multiplicity() == 2 )
		{
			// add this guy...
			// set the orientation for any unaligned intervals
			addGuy( seqI, orient, new_ivs, new_list );

			// prepare a new one
			vector< AbstractMatch* > matches;
			orient = iv_list[ordI]->Orientation(seqI);
			iv_list[ordI]->StealMatches( matches );
			new_ivs.insert( new_ivs.end(), matches.begin(), matches.end() );
			// if it's the last one then add
			if( ordI + 1 == iv_list.size() )
				addGuy( seqI, orient, new_ivs, new_list );
			continue;
		}
		if( new_ivs.size() == 0 )
		{
			vector< AbstractMatch* > matches;
			iv_list[ordI]->StealMatches( matches );
			new_ivs.insert( new_ivs.end(), matches.begin(), matches.end() );
			continue;
		}
		// split this one in half (if its not the last one and there's something to split)...
		Interval* left_iv = iv_list[ordI]->Copy();
		to_delete.push_back( left_iv );	// make sure this gets deleted later
		bool cropped = (ordI + 1 < iv_list.size() && iv_list[ordI]->Length(seqI) > 1);
		if( cropped )
		{
			gnSeqI lendo = left_iv->AlignmentLength() / 2;
			left_iv->CropEnd( left_iv->AlignmentLength() - lendo );
			iv_list[ordI]->CropStart( lendo );
		}
		vector< AbstractMatch* > matches;
		left_iv->StealMatches( matches );
		if( orient == AbstractMatch::forward )
			new_ivs.insert( new_ivs.end(), matches.begin(), matches.end() );
		else
			new_ivs.insert( new_ivs.begin(), matches.begin(), matches.end() );

		addGuy( seqI, orient, new_ivs, new_list );
		// prepare for the next
		orient = AbstractMatch::undefined;
		if(cropped)
			ordI--;	// if we split a match, make sure we get the rest of this match on the next run through the loop
	}

	if( new_ivs.size() > 0 )
	{
		// uh-oh. there must not have been anything aligned
		addGuy( seqI, AbstractMatch::forward, new_ivs, new_list );
	}

	// free up any left_ivs that were allocated
	for( size_t delI = 0; delI < to_delete.size(); delI++ )
		to_delete[delI]->Free();

	// free up ivs left in iv_list
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
		if( iv_list[ivI] != NULL )
			iv_list[ivI]->Free();
	iv_list.clear();
}


/**
 * 
 */
void ProgressiveAligner::createAncestralOrdering( vector<Interval*>& interval_list, vector< SuperInterval >& ancestral_sequence )
{
	// construct an ancestral SuperSequence
	int64 left_end = 1;
	ancestral_sequence.resize( interval_list.size() );
	for( uint ivI = 0; ivI < interval_list.size(); ++ivI ){
		if(debug_aligner)
			interval_list[ivI]->ValidateMatches();
		vector<AbstractMatch*> matches;
		interval_list[ivI]->StealMatches(matches);
		ancestral_sequence[ivI].reference_iv.SetMatches(matches);
		ancestral_sequence[ivI].SetLeftEnd(left_end);
		ancestral_sequence[ivI].SetLength(ancestral_sequence[ivI].reference_iv.AlignmentLength());
		if(debug_aligner)
			ancestral_sequence[ivI].ValidateSelf();
		left_end += ancestral_sequence[ivI].Length();
	}
}

void markAligned( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t subject_node, node_id_t neighbor )
{
	for( uint parentI = 0; parentI < alignment_tree[subject_node].parents.size(); parentI++ )
		if( alignment_tree[subject_node].parents[parentI] == neighbor )
			alignment_tree[subject_node].parents_aligned[parentI] = true;
	for( uint childI = 0; childI < alignment_tree[subject_node].children.size(); childI++ )
		if( alignment_tree[subject_node].children[childI] == neighbor )
			alignment_tree[subject_node].children_aligned[childI] = true;
}


bool
ProgressiveAligner::validateSuperIntervals(node_id_t node1, node_id_t node2, node_id_t ancestor)
{
		// validate the ancestor
	bool borked = false;
	vector< SuperInterval >& siv_list = alignment_tree[ancestor].ordering;
	gnSeqI n1_len = 0;
	gnSeqI n2_len = 0;
	gnSeqI my_len = 0;
	gnSeqI my_iv_len = 0;
	for( size_t sivI = 0; sivI < siv_list.size(); sivI++ )
	{
		if( siv_list[sivI].reference_iv.Start(0) != 0 )
			n1_len += siv_list[sivI].reference_iv.Length(0);
		if( siv_list[sivI].reference_iv.Start(1) != 0 )
			n2_len += siv_list[sivI].reference_iv.Length(1);
		my_len += siv_list[sivI].Length();
		my_iv_len += siv_list[sivI].reference_iv.AlignmentLength();
		siv_list[sivI].ValidateSelf();
	}
	gnSeqI real_n1len = 0;
	gnSeqI real_n2len = 0;

	vector< SuperInterval >& siv1_list = alignment_tree[node1].ordering;
	for( size_t sivI = 0; sivI < siv1_list.size(); sivI++ )
	{
		if( siv1_list[sivI].Length() == 0 )
			borked = true;
		real_n1len += siv1_list[sivI].Length();
		siv1_list[sivI].ValidateSelf();
	}

	vector< SuperInterval >& siv2_list = alignment_tree[node2].ordering;
	for( size_t sivI = 0; sivI < siv2_list.size(); sivI++ )
	{
		if( siv2_list[sivI].Length() == 0 )
			borked = true;
		real_n2len += siv2_list[sivI].Length();
		siv2_list[sivI].ValidateSelf();
	}

	if( real_n1len != n1_len || real_n2len != n2_len )
			borked = true;

	// check that each picks up where the last left off
	for( size_t sivI = 1; sivI < siv1_list.size(); sivI++ )
		if( siv1_list[sivI].LeftEnd() != siv1_list[sivI-1].LeftEnd() + siv1_list[sivI-1].Length() )
		{
			borked = true;
		}
	for( size_t sivI = 1; sivI < siv2_list.size(); sivI++ )
		if( siv2_list[sivI].LeftEnd() != siv2_list[sivI-1].LeftEnd() + siv2_list[sivI-1].Length() )
		{
			borked = true;
		}

	if( my_len != my_iv_len )
		borked = true;

	if( my_len < real_n1len || my_len < real_n2len )
		borked = true;

	if( borked )
	{
		breakHere();
		cerr << "child1 has " << siv1_list.size() << " ivs totalling " << real_n1len << " nt\n";
		cerr << "child2 has " << siv2_list.size() << " ivs totalling " << real_n2len << " nt\n";
		cerr << "parent has " << siv_list.size() << " ivs, n1_len: " << n1_len << " n2_len: " << n2_len << endl;
	}
	return borked;

}

bool ProgressiveAligner::validatePairwiseIntervals(node_id_t node1, node_id_t node2, std::vector<Interval*>& pair_iv)
{
		// validate the ancestor
	bool borked = false;
	gnSeqI n1_len = 0;
	gnSeqI n2_len = 0;
	for( size_t sivI = 0; sivI < pair_iv.size(); sivI++ )
	{
		if( pair_iv[sivI]->Start(0) != 0 )
			n1_len += pair_iv[sivI]->Length(0);
		if( pair_iv[sivI]->Start(1) != 0 )
			n2_len += pair_iv[sivI]->Length(1);

		vector< bitset_t > aln_mat;
		pair_iv[sivI]->GetAlignment(aln_mat);
		if( aln_mat[0].size() != pair_iv[sivI]->AlignmentLength() )
		{
			cerr << "broked\n";
		}
		pair_iv[sivI]->ValidateMatches();
	}
	gnSeqI real_n1len = 0;
	gnSeqI real_n2len = 0;

	vector< SuperInterval >& siv1_list = alignment_tree[node1].ordering;
	for( size_t sivI = 0; sivI < siv1_list.size(); sivI++ )
	{
		if( siv1_list[sivI].Length() == 0 )
			borked = true;
		real_n1len += siv1_list[sivI].Length();
	}

	vector< SuperInterval >& siv2_list = alignment_tree[node2].ordering;
	for( size_t sivI = 0; sivI < siv2_list.size(); sivI++ )
	{
		if( siv2_list[sivI].Length() == 0 )
			borked = true;
		real_n2len += siv2_list[sivI].Length();
	}

	if( real_n1len != n1_len || real_n2len != n2_len )
			borked = true;

	// check for overlapping intervals
	vector< Interval* > tmp_iv_list = pair_iv;
	for( uint seqI = 0; seqI < 2; seqI++ )
	{
		SSC<Interval> ssc(seqI);
		sort( tmp_iv_list.begin(), tmp_iv_list.end(), ssc );
		for( size_t ivI = 1; ivI < tmp_iv_list.size(); ivI++ )
		{
			if( tmp_iv_list[ivI-1]->LeftEnd(seqI) == NO_MATCH || tmp_iv_list[ivI]->LeftEnd(seqI) == NO_MATCH )
				continue;
			if( tmp_iv_list[ivI-1]->RightEnd(seqI) >= tmp_iv_list[ivI]->LeftEnd(seqI) )
			{
				cerr << "overlap:\n";
				cerr << "tmp_iv_list[ivI-1].RightEnd(seqI): " << tmp_iv_list[ivI-1]->RightEnd(seqI) << endl;
				cerr << "tmp_iv_list[ivI].LeftEnd(seqI): " << tmp_iv_list[ivI]->LeftEnd(seqI) << endl;
				breakHere();
			}
		}
	}

	if( borked )
	{
		cerr << "child1 has " << siv1_list.size() << " ivs totalling " << real_n1len << " nt\n";
		cerr << "child2 has " << siv2_list.size() << " ivs totalling " << real_n2len << " nt\n";
		cerr << "parent has " << pair_iv.size() << " ivs, n1_len: " << n1_len << " n2_len: " << n2_len << endl;
		if( n2_len < real_n2len )
		{
			SSC<Interval> sortie(1);
			sort( pair_iv.begin(), pair_iv.end(), sortie );
			size_t prev_iv = 9999999;
			for( size_t ivI = 0; ivI < pair_iv.size(); ++ivI)
			{
				if( pair_iv[ivI]->LeftEnd(1) == NO_MATCH )
					continue;

				if( prev_iv != 9999999 )
					cerr << "diff: " << pair_iv[ivI]->LeftEnd(1) - pair_iv[prev_iv]->RightEnd(1) << endl;
				cerr << "Interval " << ivI << " LeftEnd(1): " << pair_iv[ivI]->LeftEnd(1) << " RightEnd(1): " << pair_iv[ivI]->RightEnd(1) << std::endl;
				prev_iv = ivI;
			}
		}else if( n2_len > real_n2len )
		{
			SSC<Interval> sortie(1);
			sort( pair_iv.begin(), pair_iv.end(), sortie );
			for( size_t ivI = 0; ivI < pair_iv.size(); ++ivI)
			{
				if( pair_iv[ivI]->LeftEnd(1) < real_n2len )
					continue;
				cerr << "Interval " << ivI << " LeftEnd(1): " << pair_iv[ivI]->LeftEnd(1) << " RightEnd(1): " << pair_iv[ivI]->RightEnd(1) << std::endl;
			}
		}
		breakHere();
	}
	return borked;
}

void ProgressiveAligner::alignNodes( node_id_t node1, node_id_t node2, node_id_t ancestor )
{
	cout << "Aligning node " << node1 << " to " << node2 << " via " << ancestor << "!\n";
	// if node1 and node2 are not already children of ancestor then make it so...
	if( alignment_tree[node1].parents[0] != ancestor || 
		alignment_tree[node2].parents[0] != ancestor )
	{
		breakHere();
		cerr << "rotten\n";
	}
	
	alignProfileToProfile(node1, node2, ancestor);

	// mark edges as aligned
	markAligned( alignment_tree, node1, node2 );
	markAligned( alignment_tree, node2, node1 );
	markAligned( alignment_tree, node1, ancestor );
	markAligned( alignment_tree, node2, ancestor );
	markAligned( alignment_tree, ancestor, node1 );
	markAligned( alignment_tree, ancestor, node2 );
}

void findMidpoint( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t& n1, node_id_t& n2 )
{
	// use boost's all pairs shortest path to find the longest path on the tree 
	// Then actually traverse the path to determine which edge
	// is halfway.
	double scaling_factor = 10000;
	using namespace boost;
	typedef adjacency_list<vecS, vecS, undirectedS, no_property,
	property< edge_weight_t, int, property< edge_color_t, default_color_type > > > Graph;
	const int V = alignment_tree.size();
	const std::size_t E = alignment_tree.size()-1;
	typedef std::pair < int, int >Edge;
	Edge* edge_array = new Edge[ alignment_tree.size() - 1 ];
	int* weights = new int[ alignment_tree.size() - 1 ];
	bitset_t child_found( alignment_tree.size(), false );
	size_t eI = 0;
	for( size_t vI = 0; vI < V; ++vI )
	{
		if( alignment_tree[vI].parents.size() != 0 )
		{
			edge_array[eI] = Edge( vI, alignment_tree[vI].parents[0] );
			// for some reason boost insists on using an int for weights.  need to figure that out
			weights[eI] = (int)(scaling_factor * genome::absolut(alignment_tree[vI].distance));
			eI++;
		}
	}

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	// VC++ can't handle the iterator constructor
	Graph g(V);
	for (std::size_t j = 0; j < E; ++j)
	add_edge(edge_array[j].first, edge_array[j].second, g);
#else
	Graph g(edge_array, edge_array + E, V);
#endif

	property_map < Graph, edge_weight_t >::type w = get(edge_weight, g);
	int *wp = weights;

	graph_traits < Graph >::edge_iterator e, e_end;
	for (boost::tie(e, e_end) = edges(g); e != e_end; ++e)
		w[*e] = *wp++;

	boost::multi_array<int,2> D( boost::extents[V][V] );
	bool success = johnson_all_pairs_shortest_paths(g, D);
	if( !success )
	{
		cerr << "failed, is this really a tree?\n";
		return;
	}

	// find the most distant pair of nodes
	int max_dist = 0;
	for (int i = 0; i < V; ++i) {
		for (int j = 0; j < V; ++j) {
			if( D[i][j] > max_dist )
			{
				max_dist = D[i][j];
				n1 = i;
				n2 = j;
			}
		}
	}

	typedef graph_traits<Graph>::vertex_descriptor vertex_t;
	std::vector < vertex_t > pred(num_vertices(g));
	std::vector < int > dist(num_vertices(g));
	pred[n1] = n1;

	undirected_dfs(g,
		root_vertex( vertex( n1, g ) ).
		visitor( make_dfs_visitor( make_pair(
			record_predecessors(&pred[0], on_tree_edge()),
			record_distances(&dist[0], on_tree_edge())
		))).
		edge_color_map(get(edge_color, g))
		);

	int cur_node = n2;
	int prev_node = n2;
	max_dist /= 2;
	while( cur_node != n1 && max_dist > 0 )
	{
		if( alignment_tree[cur_node].parents.size() > 0 && 
			alignment_tree[cur_node].parents[0] == pred[cur_node] )
		{
			max_dist -= (int)(scaling_factor * alignment_tree[cur_node].distance);
			prev_node = cur_node;
			cur_node = pred[cur_node];
		}else
		{
			prev_node = cur_node;
			cur_node = pred[cur_node];
			max_dist -= (int)(scaling_factor * alignment_tree[cur_node].distance);
		}
	}
	n1 = cur_node;
	n2 = prev_node;

	delete[] edge_array;
	delete[] weights;
}

void extendRootBranches( PhyloTree< AlignmentTreeNode >& alignment_tree )
{
	// find the max branch length and set the root branch lengths to twice that
	// swap children while we're at it
	node_id_t ancestor = alignment_tree.root;
	double max_blen = -(std::numeric_limits<double>::max)();
	for( size_t nI = 0; nI < alignment_tree.size(); ++nI )
	{
		if( alignment_tree[nI].distance > max_blen )
			max_blen = alignment_tree[nI].distance;
		if( alignment_tree[nI].children.size() > 0 &&
			alignment_tree[nI].children[0] > alignment_tree[nI].children[1] )
		{
			swap( alignment_tree[nI].children[0], alignment_tree[nI].children[1] );
		}
	}
	for( size_t cI = 0; cI < alignment_tree[ancestor].children.size(); ++cI )
		alignment_tree[alignment_tree[ancestor].children[cI]].distance = 2.0 * max_blen;
}

void chooseNextAlignmentPair( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t& node1, node_id_t& node2, node_id_t& ancestor )
{

	// find the nearest alignable neighbor
	node1 = 0;
	node2 = 0;
	ancestor = 0;
	double nearest_distance = (numeric_limits<double>::max)();
	for( node_id_t nodeI = 0; nodeI < alignment_tree.size(); nodeI++ )
	{
		AlignmentTreeNode& cur_node = alignment_tree[ nodeI ];

		// skip this node if it's already been completely aligned
		// or is an extant sequence
		boolean completely_aligned = true;
		for( uint alignedI = 0; alignedI < cur_node.children_aligned.size(); alignedI++ )
			completely_aligned = completely_aligned && cur_node.children_aligned[alignedI];
		for( uint alignedI = 0; alignedI < cur_node.parents_aligned.size(); alignedI++ )
			completely_aligned = completely_aligned && cur_node.parents_aligned[alignedI];
		if( cur_node.sequence != NULL || completely_aligned )
			continue;
		

		vector< node_id_t > neighbor_id;
		vector< boolean > alignable;
		vector< double > distance;
		
		for( uint parentI = 0; parentI < cur_node.parents.size(); parentI++ )
		{
			neighbor_id.push_back( cur_node.parents[parentI] );
			vector< node_id_t >::iterator cur_neighbor = neighbor_id.end() - 1;
			if( *cur_neighbor == alignment_tree.root )
			{
				// need special handling for the root since the alignment
				// tree is supposed to be unrooted
				// add all of root's children except this one
			}
			distance.push_back( cur_node.distance );
			alignable.push_back( !cur_node.parents_aligned[parentI] && (alignment_tree[*cur_neighbor].ordering.size() != 0 || alignment_tree[*cur_neighbor].sequence != NULL) );
		}

		for( uint childI = 0; childI < cur_node.children.size(); childI++ )
		{
			neighbor_id.push_back( cur_node.children[childI] );
			vector< node_id_t >::iterator cur_neighbor = neighbor_id.end() - 1;
			distance.push_back( alignment_tree[*cur_neighbor].distance );
			alignable.push_back( !cur_node.children_aligned[childI] && (alignment_tree[*cur_neighbor].ordering.size() != 0 || alignment_tree[*cur_neighbor].sequence != NULL) );
		}

		if( cur_node.ordering.size() != 0 )
		{
			// this one already has at least two sequences aligned, if another
			// is alignable then check its distance
			for( int i = 0; i < neighbor_id.size(); i++ ){
				if( !alignable[i] )
					continue;
				if( distance[i] < nearest_distance )
				{
					nearest_distance = distance[i];
					node1 = nodeI;
					node2 = neighbor_id[i];
					ancestor = nodeI;
				}
			}
		}else{
			// find the nearest alignable pair
			for( int i = 0; i < neighbor_id.size(); i++ )
			{
				if( !alignable[i] )
					continue;
				for( int j = i+1; j < neighbor_id.size(); j++ )
				{
					if( !alignable[j] )
						continue;
					if( distance[i] + distance[j] < nearest_distance )
					{
						nearest_distance = distance[i] + distance[j];
						node1 = neighbor_id[i];
						node2 = neighbor_id[j];
						ancestor = nodeI;
					}
				}
			}
		}
	}
}

/** use a list of precomputed matches instead of computing them */
void ProgressiveAligner::setPairwiseMatches( MatchList& pair_ml )
{
	original_ml = pair_ml;
	pair_ml.clear();	// ProgressiveAligner owns the matches now...
}


node_id_t createAlignmentTreeRoot( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t node1, node_id_t node2 )
{
	// create a new node and link it inline between node1 and node2
	AlignmentTreeNode atn;
	alignment_tree.push_back( atn );
	AlignmentTreeNode& old_root = alignment_tree[alignment_tree.root];
	AlignmentTreeNode& new_root = alignment_tree.back();

	if( find( alignment_tree[node1].children.begin(), alignment_tree[node1].children.end(), node2 ) !=
		alignment_tree[node1].children.end() )
	{
		new_root.children.push_back(node2);
		new_root.parents.push_back(node1);
		alignment_tree[node2].parents.push_back(alignment_tree.size()-1);
		alignment_tree[node1].children.push_back(alignment_tree.size()-1);
	}else{
		new_root.parents.push_back(node2);
		new_root.children.push_back(node1);
		alignment_tree[node2].children.push_back(alignment_tree.size()-1);
		alignment_tree[node1].parents.push_back(alignment_tree.size()-1);
	}

	// completely unlink node1 and node2 from each other
	findAndErase( alignment_tree[node1].children, node2 );
	findAndErase( alignment_tree[node2].children, node1 );
	findAndErase( alignment_tree[node1].parents, node2 );
	findAndErase( alignment_tree[node2].parents, node1 );


	// re-root the tree on the new node
	rerootTree( alignment_tree, alignment_tree.size()-1 );

	new_root.children_aligned = vector< boolean >( new_root.children.size(), false );
	old_root.children_aligned = vector< boolean >( old_root.children.size(), false );
	old_root.parents_aligned = vector< boolean >( old_root.parents.size(), false );
	new_root.sequence = NULL;

	return alignment_tree.root;
}

void ProgressiveAligner::extractAlignment( node_id_t ancestor, size_t super_iv, GappedAlignment& gal )
{
	CompactGappedAlignment<> cga;
	extractAlignment( ancestor, super_iv, cga );
	vector< string > aln;
	GetAlignment( cga, this->original_ml.seq_table, aln );
	gal = GappedAlignment(cga.SeqCount(), 0);
	for( size_t seqI = 0; seqI < cga.SeqCount(); ++seqI )
	{
		gal.SetStart(seqI, cga.Start(seqI));
		if( cga.Orientation(seqI) != AbstractMatch::undefined )
			gal.SetLength(cga.Length(seqI), seqI);
	}
	gal.SetAlignment(aln);

}

void ProgressiveAligner::extractAlignment( node_id_t ancestor, size_t super_iv, CompactGappedAlignment<>& cga )
{
	// determine the leaf node intervals below this super_iv
	vector< pair< node_id_t, size_t > > node_siv_list;
	stack< pair<node_id_t,size_t> > node_stack;
	node_stack.push(make_pair(ancestor,super_iv));
	while( node_stack.size() > 0 )
	{
		pair<node_id_t,size_t> cur = node_stack.top();
		node_id_t cur_node = cur.first;
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() == 0 )
			node_siv_list.push_back( cur );
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			if( alignment_tree[cur_node].ordering[cur.second].reference_iv.LeftEnd(childI) == NO_MATCH )
				continue;
			size_t child_siv = childI == 0 ? alignment_tree[cur_node].ordering[cur.second].c1_siv : 
				alignment_tree[cur_node].ordering[cur.second].c2_siv;
			node_stack.push(make_pair(alignment_tree[cur_node].children[childI], child_siv) );
			node_id_t n = alignment_tree[cur_node].children[childI];
			if( alignment_tree[cur_node].ordering[cur.second].reference_iv.Length(childI) != alignment_tree[n].ordering[child_siv].Length() )
			{
				breakHere();
				cerr << "alignment_tree[cur_node].ordering[cur.second].reference_iv.Length(childI): " << alignment_tree[cur_node].ordering[cur.second].reference_iv.Length(childI) << endl;
				cerr << "rotten in the state of denmark...\n";
			}
		}
	}

	// armed with the list of pairs, extract each one...

	// for each interval at the root write out the alignment
	SuperInterval& a_iv = alignment_tree[ancestor].ordering[super_iv];
	cga = CompactGappedAlignment<>(seq_count, a_iv.Length());
	vector< bitset_t > aln_mats( seq_count );

	// use translateCoordinates to map out each sequence's original coordinates
	// to the alignment coordinates
	for( size_t pairI = 0; pairI < node_siv_list.size(); pairI++ )
	{
		node_id_t nodeI = node_siv_list[pairI].first;
		size_t seq_siv = node_siv_list[pairI].second;
		
		// translate seq_siv into ancestor alignment coordinates?  
		// we can abuse translateCoordinates and the Match data structure :
		//   - add a single "match" covering the entire sequence
		//   - translate it up to alignment root coordinates
		uint seqI = node_sequence_map[nodeI];
		Match mm(2);
		mm.SetStart(0, alignment_tree[nodeI].ordering[seq_siv].LeftEnd());
		mm.SetStart(1, alignment_tree[nodeI].ordering[seq_siv].LeftEnd());
		mm.SetLength( alignment_tree[nodeI].ordering[seq_siv].Length() );
		
		vector< AbstractMatch* > aml( 1, mm.Copy() );
		translateGappedCoordinates( aml, 0, nodeI, ancestor );

		if( aml.size() > 1 )
		{
			cerr << "huh?";
			genome::breakHere();
			SingleStartComparator<AbstractMatch> ssc( 0 );
			sort( aml.begin(), aml.end(), ssc );	// huh?
		}
		CompactGappedAlignment<>* trans_cga = dynamic_cast<CompactGappedAlignment<>*>(aml[0]);
		if( trans_cga == NULL )
		{
			CompactGappedAlignment<> tmp_cga;
			trans_cga = tmp_cga.Copy();
			*trans_cga = CompactGappedAlignment<>(*aml[0]);
		}

		if( trans_cga->LeftEnd(0) + trans_cga->Length(0) > a_iv.LeftEnd() + a_iv.Length() )
		{
			cerr << "trans_cga->Start(0): " << trans_cga->Start(0) << " trans_cga->Length(0): " << trans_cga->Length(0) << endl;
			cerr << "a_iv.LeftEnd(): " << a_iv.LeftEnd() << " a_iv.Length(): " << a_iv.Length() << endl;
			breakHere();
		}
		bool parity = trans_cga->Orientation(0) == trans_cga->Orientation(1);
		cga.SetLeftEnd(seqI, trans_cga->LeftEnd(1));
		AbstractMatch::orientation o = parity ? AbstractMatch::forward : AbstractMatch::reverse;
		cga.SetOrientation(seqI, o);
		const vector< bitset_t >& tmp = trans_cga->GetAlignment();
		aln_mats[seqI] = tmp[1];

		size_t offset = trans_cga->LeftEnd(0) - a_iv.LeftEnd();
		if( aln_mats[seqI].size() < a_iv.Length() )
		{
			// need to resize and shift appropriately
			aln_mats[seqI].resize( a_iv.Length() );
			aln_mats[seqI] <<= offset;	// this is backwards in boost::dynamic_bitset for some reason...
		}
		if( trans_cga->LeftEnd(0) < a_iv.LeftEnd() )
		{
			cerr << "trans_cga->LeftEnd(0): " << trans_cga->LeftEnd(0) << endl;
			cerr << "a_iv.LeftEnd(): " << a_iv.LeftEnd() << endl;
			breakHere();
		}

		// validate match lengths
		if( trans_cga->Length(1) != alignment_tree[nodeI].ordering[seq_siv].Length() )
		{
			cerr << "b0rked\n";
			breakHere();
		}
		// set the length and alignment appropriately
		cga.SetLength(trans_cga->Length(1), seqI);

		// free storage used by trans_cga
		trans_cga->Free();
	}
	for( uint seqI = 0; seqI < aln_mats.size(); seqI++ )
		if( aln_mats[seqI].size() == 0 )
			aln_mats[seqI].resize( a_iv.Length() );
	cga.SetAlignment(aln_mats);
}

unsigned getDefaultBreakpointMax( const std::vector< genome::gnSequence* >& seq_table )
{
	double avg_len = 0;
	for( size_t seqI = 0; seqI < seq_table.size(); ++seqI )
		avg_len += seq_table[seqI]->length();
	avg_len /= (double)(seq_table.size());
	// heavily rearranged, recently diverged genomes like yersinia have up to 20 rearrangements per megabase of sequence
	avg_len /= 1000000.0;	// convert to number of megabases
	avg_len *= 20.0;	// "lots" of rearrangement
	return (unsigned)avg_len;
}

// get a pairwise bp distance
void ProgressiveAligner::CreatePairwiseBPDistance( boost::multi_array<double, 2>& bp_distmat )
{
	uint seq_count = original_ml.seq_table.size();
	bp_distmat.resize(boost::extents[seq_count][seq_count]);
	for( size_t i = 0; i < seq_count; ++i )
		for( size_t j = 0; j < seq_count; ++j )
			bp_distmat[i][j] = 1;

#ifdef LCB_WEIGHT_LOSS_PLOT
	stringstream pair_bp_ofname;
	pair_bp_ofname << "pair_bp_log.txt";
	ofstream pair_bp_out( pair_bp_ofname.str().c_str() );
#endif

	for( uint seqI = 0; seqI < seq_count; seqI++ )
	{
		for( uint seqJ = seqI + 1; seqJ < seq_count; seqJ++ )
		{
			vector<uint>::iterator n1 = find( node_sequence_map.begin(), node_sequence_map.end(), seqI );
			vector<uint>::iterator n2 = find( node_sequence_map.begin(), node_sequence_map.end(), seqJ );
			vector<node_id_t> n1_seqs( 1, n1-node_sequence_map.begin() );
			vector<node_id_t> n2_seqs( 1, n2-node_sequence_map.begin() );
			Matrix<MatchList> mml;
			getPairwiseMatches(n1_seqs, n2_seqs, mml);
			MatchList& ml = mml(0,0);

			// eliminate overlaps as they correspond to inconsistently or
			// multiply aligned regions
			EliminateOverlaps_v2( ml );
			ml.MultiplicityFilter(2);

			// do greedy b.p. elimination on the matches
			vector< MatchList > LCB_list;
			vector< LCB > adjacencies;
			vector< gnSeqI > breakpoints;
			IdentifyBreakpoints( ml, breakpoints );
			ComputeLCBs_v2( ml, breakpoints, LCB_list );
			vector< double > lcb_scores( LCB_list.size() );
			for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
				lcb_scores[lcbI] = GetPairwiseLcbScore( LCB_list[lcbI], ml.seq_table, sol_list[seqI], sol_list[seqJ] );

			computeLCBAdjacencies_v3( LCB_list, lcb_scores, adjacencies );

			// want to discard all low-weight LCBs
			// to arrive at a set of reliable LCBs
			GreedyRemovalScorer wbs( adjacencies, bp_dist_estimate_score );
#ifdef LCB_WEIGHT_LOSS_PLOT
			cur_min_coverage = greedyBreakpointElimination_v4( adjacencies, lcb_scores, wbs, &pair_bp_out, seqI, seqJ );
			pair_bp_out.flush();
#else
			cur_min_coverage = greedyBreakpointElimination_v4( adjacencies, lcb_scores, wbs, NULL );
#endif
			MatchList deleted_matches;
			filterMatches_v2( adjacencies, LCB_list, lcb_scores, deleted_matches );
			cout << "Pair (" << seqI << "," << seqJ << ") has " << LCB_list.size() << " well-supported breakpoints\n";
			
			// now set the distance entry
			bp_distmat[seqI][seqJ] = LCB_list.size();
			bp_distmat[seqJ][seqI] = LCB_list.size();

			// free the matches
			for( size_t dI = 0; dI < ml.size(); dI++ )
				ml[dI]->Free();
		}
	}
	// normalize to [0,1]
	double bp_max = 0;
	for( uint i = 0; i < bp_distmat.shape()[0]; ++i )
		for( uint j = 0; j < bp_distmat.shape()[1]; ++j )
		{
			if( bp_distmat[i][j] > bp_max )
				bp_max = bp_distmat[i][j];
		}

	double default_max = getDefaultBreakpointMax(original_ml.seq_table);
	bp_max = bp_max > default_max ? bp_max : default_max;

	for( uint i = 0; i < bp_distmat.shape()[0]; ++i )
		for( uint j = 0; j < bp_distmat.shape()[1]; ++j )
		{
			if( i != j )
				bp_distmat[i][j] /= bp_max;
			bp_distmat[i][j] *= bp_dist_scale;
		}
}

template< typename MatchListType >
void makeAlignmentTree( PhyloTree< AlignmentTreeNode >& alignment_tree, MatchListType& mlist, vector< uint >& node_sequence_map )
{
	// initialize all nodes to unaligned
	for( node_id_t nodeI = 0; nodeI < alignment_tree.size(); nodeI++ )
	{
		alignment_tree[nodeI].children_aligned = vector< boolean >( alignment_tree[nodeI].children.size(), false );
		alignment_tree[nodeI].parents_aligned = vector< boolean >( alignment_tree[nodeI].parents.size(), false );
		alignment_tree[nodeI].sequence = NULL;
		alignment_tree[nodeI].refined = false;
	}

	// set the sequence appropriately for extant sequences
	node_sequence_map = vector< uint >( alignment_tree.size(), -1 );
	for( uint seqI = 0; seqI < mlist.seq_table.size(); seqI++ )
	{
		stringstream seq_name;
		seq_name << "seq" << seqI + 1;
		node_id_t nodeI = 0;
		for( ; nodeI < alignment_tree.size(); nodeI++ )
		{
			if( seq_name.str() == alignment_tree[nodeI].name )
			{
				alignment_tree[nodeI].sequence = mlist.seq_table[seqI];
				Match mm(1);
				Match* m = mm.Copy();
				m->SetStart(0,1);
				m->SetLength(alignment_tree[nodeI].sequence->length(),0);
				vector<AbstractMatch*> tmp(1,m);
				Interval iv( tmp.begin(), tmp.end() );
				m->Free();
				SuperInterval si( iv );
				si.SetLeftEnd(1);
				si.SetLength(alignment_tree[nodeI].sequence->length());
				alignment_tree[nodeI].ordering.push_back( si );
				node_sequence_map[nodeI] = seqI;
				break;
			}
		}
		if( nodeI == alignment_tree.size() )
			throw "Phylogenetic tree names unrecognized.  Should follow seqN naming format\n";
	}
}

void DistanceMatrix( IntervalList& iv_list, NumericMatrix<double>& distmat )
{
	IdentityMatrix( iv_list, distmat );
	TransformDistanceIdentity(distmat);
}

void projectIntervalList( IntervalList& iv_list, vector< uint >& projection, vector< vector< MatchProjectionAdapter* > >& LCB_list, vector< LCB >& projected_adjs )
{
	vector< size_t > proj(projection.size());
	for( size_t i = 0; i < projection.size(); ++i )
		proj[i] = projection[i];
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
		MatchProjectionAdapter mpa_tmp( &iv_list[corI], proj );
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

void makeSuperIntervals( IntervalList& iv_list, PhyloTree< TreeNode >& alignment_tree, vector< uint >& node_sequence_map )
{
	std::stack< node_id_t > node_stack;
	node_stack.push( alignment_tree.root );
	bitset_t visited( alignment_tree.size(), false );
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		// visit post-order
		for( size_t cI = 0; cI < alignment_tree[cur_node].children.size(); ++cI )
		{
			if( !visited[alignment_tree[cur_node].children[cI]] )
				node_stack.push(alignment_tree[cur_node].children[cI]);
		}
		if( node_stack.top() != cur_node )
			continue;
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() == 0 )
			continue;	// only process internal nodes

		// process this node
		// construct pairwise LCBs

		uint seqI = node_sequence_map[alignment_tree[cur_node].children[0]];
		uint seqJ = node_sequence_map[alignment_tree[cur_node].children[0]];
		vector< uint > projection( 2 );
		projection[0] = seqI;
		projection[1] = seqJ;

		vector< vector< MatchProjectionAdapter* > > LCB_list;
		vector< LCB > projected_adjs;
		projectIntervalList( iv_list, projection, LCB_list, projected_adjs );

		// create a superinterval for each adj 
/*		alignment_tree[cur_node].ordering.resize(adjs.size());
		for( size_t adjI = 0; adjI < adjs.size(); ++adjI )
		{
			SuperInterval& siv = alignment_tree[cur_node].ordering[adjI];
			Match mleft(2);
			mleft.SetStart(0,adjI);
			mleft.SetStart(1,adjI);
			mleft.SetLength(1);
			siv.SetLeftEnd( adjI );
			siv.SetLength(1);
		}
*/
	}
}

void ProgressiveAligner::alignPP(IntervalList& prof1, IntervalList& prof2, IntervalList& interval_list )
{
	if( debug_aligner )
	{
		debug_interval = true;
		debug_cga = true;
	}

	seq_count = prof1.seq_table.size() + prof2.seq_table.size();

	if( this->breakpoint_penalty == -1 )
		this->breakpoint_penalty = getDefaultBreakpointPenalty( original_ml.seq_table );

	if( this->bp_dist_estimate_score == -1 )
		this->bp_dist_estimate_score = getDefaultBpDistEstimateMinScore( original_ml.seq_table );
	cout << "using default bp penalty: " << breakpoint_penalty << endl;
	cout << "using default bp estimate min score: " << bp_dist_estimate_score << endl;

	if( this->collinear_genomes )
		this->breakpoint_penalty = -1;

	if( collinear_genomes )
		cout << "\nAssuming collinear genomes...\n";
		
	EliminateOverlaps_v2( original_ml );
	// use existing pairwise matches
	MatchList mlist;
	mlist.clear();
	mlist = original_ml;
	cout << "Starting with " << mlist.size() << " multi-matches\n";

//
// Step 1) Compute guide trees for each profile and join them
//
	NumericMatrix< double > distance1;
	DistanceMatrix( prof1, distance1 );
	NumericMatrix< double > distance2;
	DistanceMatrix( prof2, distance2 );

	// Make a phylogenetic tree
	// use the identity matrix method and convert to a distance matrix
	ClustalInterface& ci = ClustalInterface::getClustalInterface();	
	string guide_tree_fname1 = CreateTempFileName("guide_tree");
	ci.SetDistanceMatrix( distance1, guide_tree_fname1 );
	string guide_tree_fname2 = CreateTempFileName("guide_tree");
	ci.SetDistanceMatrix( distance2, guide_tree_fname2 );

	// read the trees
	ifstream tree_file1( guide_tree_fname1.c_str() );
	if( !tree_file1.is_open() )
		throw "Error opening guide tree file";
	PhyloTree< AlignmentTreeNode > tree1;
	tree1.readTree( tree_file1 );
	tree_file1.close();
	ifstream tree_file2( guide_tree_fname2.c_str() );
	if( !tree_file2.is_open() )
		throw "Error opening guide tree file";
	PhyloTree< AlignmentTreeNode > tree2;
	tree2.readTree( tree_file2 );
	tree_file2.close();


	// compute pairwise distances among all nodes
	NumericMatrix< double > distance;
	DistanceMatrix( mlist, distance );
	conservation_distance.resize(boost::extents[seq_count][seq_count]);
	for( uint seqI = 0; seqI < seq_count; ++seqI )
		for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
			if( seqJ > seqI )
				conservation_distance[seqI][seqJ] = distance(seqI,seqJ);
			else
				conservation_distance[seqI][seqJ] = distance(seqJ,seqI);


	if( !collinear_genomes )
	{
		cout << "Calculating pairwise breakpoint distances\n";
		CreatePairwiseBPDistance(bp_distance);
		cout << "bp distance matrix:\n";
		print2d_matrix(bp_distance, cout);
		cout << endl;
	}

	// rescale the conservation distance
	double conservation_range = 2;
	double bp_range = 2;
	for( uint seqI = 0; seqI < seq_count; ++seqI )
		for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
			conservation_distance[seqI][seqJ] = distance(seqI,seqJ) / conservation_range;

	if( !(collinear_genomes && seq_count > 20 ) )
	{
		cout << "genome content distance matrix:\n";
		print2d_matrix(conservation_distance, cout);
		cout << endl;
	}

//
// construct the alignment tree by joining trees from each profile
//
	vector< uint > nsmap1;
	vector< uint > nsmap2;
	makeAlignmentTree( tree1, prof1, nsmap1 );
//	prepareAlignmentTree(tree1);
	makeAlignmentTree( tree2, prof2, nsmap2 );
//	prepareAlignmentTree(tree2);

	alignment_tree.resize( tree1.size() + tree2.size() + 1 );
	// set the sequence appropriately for extant sequences
	node_sequence_map = vector< uint >( alignment_tree.size(), -1 );

	// initialize all nodes to unaligned
	for( node_id_t nodeI = 0; nodeI < alignment_tree.size()-1; nodeI++ )
	{
		if( nodeI < tree1.size() )
		{
			alignment_tree[nodeI].sequence = tree1[nodeI].sequence;
			alignment_tree[nodeI].children = tree1[nodeI].children;
			alignment_tree[nodeI].parents = tree1[nodeI].parents;
			alignment_tree[nodeI].ordering = tree1[nodeI].ordering;
			alignment_tree[nodeI].distance = tree1[nodeI].distance;
			alignment_tree[nodeI].name = tree1[nodeI].name;
			node_sequence_map[nodeI] = nsmap1[nodeI];
		}else{
			alignment_tree[nodeI].sequence = tree2[nodeI-tree1.size()].sequence;
			alignment_tree[nodeI].children = tree2[nodeI-tree1.size()].children;
			alignment_tree[nodeI].parents = tree2[nodeI-tree1.size()].parents;
			alignment_tree[nodeI].ordering = tree2[nodeI-tree1.size()].ordering;
			alignment_tree[nodeI].distance = tree2[nodeI-tree1.size()].distance;
			alignment_tree[nodeI].name = tree2[nodeI-tree1.size()].name;
			for( size_t cI = 0; cI < alignment_tree[nodeI].children.size(); cI++ )
				alignment_tree[nodeI].children[cI] += tree1.size();
			for( size_t pI = 0; pI < alignment_tree[nodeI].parents.size(); pI++ )
				alignment_tree[nodeI].parents[pI] += tree1.size();
			node_sequence_map[nodeI] = nsmap2[nodeI-tree1.size()];
			if( node_sequence_map[nodeI] != (std::numeric_limits<uint>::max)() )
				node_sequence_map[nodeI] += prof1.seq_table.size();
		}

		alignment_tree[nodeI].children_aligned = vector< boolean >( alignment_tree[nodeI].children.size(), true );
		alignment_tree[nodeI].parents_aligned = vector< boolean >( alignment_tree[nodeI].parents.size(), true );
		alignment_tree[nodeI].refined = true;
	}

	alignment_tree.back().children.push_back( tree1.size()-1 );
	alignment_tree.back().children.push_back( alignment_tree.size()-2 );
	alignment_tree.back().distance = 100;
	alignment_tree.back().children_aligned = vector< boolean >( alignment_tree.back().children.size(), true );
	alignment_tree.back().parents_aligned = vector< boolean >( alignment_tree.back().parents.size(), true );
	alignment_tree.back().refined = false;


	getAlignment( interval_list );

}

void ProgressiveAligner::getAlignment( IntervalList& interval_list )
{
	cout << "Aligning...\n";
	// pick each pair of sequences and align until none are left
	while(true)
	{
		node_id_t node1;
		node_id_t node2;
		node_id_t ancestor;
		chooseNextAlignmentPair( alignment_tree, node1, node2, ancestor );
		if( node1 == node2 )
			break;	// all pairs have been aligned

		// this is the last alignable pair in the unrooted tree
		// create a root from which the complete alignment can be extracted
		alignNodes( node1, node2, ancestor );
		if( ancestor == alignment_tree.root )
			break;  // all done
	}

	if( refine )
	{
		// perform iterative refinement
		cout << "Performing final pass iterative refinement\n";
		doGappedAlignment(alignment_tree.root, false);
	}

	// peel off the alignment from the root node
	cout << "root alignment has " << alignment_tree[alignment_tree.root].ordering.size() << " superintervals\n";
	vector< SuperInterval >& a_ivs = alignment_tree[alignment_tree.root].ordering;
	gnSeqI len = 0;
	for( size_t ivI = 0; ivI < a_ivs.size(); ivI++ )
	{
		len += a_ivs[ivI].Length();
	}
	cout << "root alignment length: " << len << endl;


	// for each interval at the root write out the alignment
	for( size_t ivI = 0; ivI < a_ivs.size(); ivI++ )
	{
		GappedAlignment ga(seq_count, a_ivs[ivI].Length());
		extractAlignment(alignment_tree.root, ivI, ga);
		vector<AbstractMatch*> tmp(1, &ga);
		interval_list.push_back( Interval(tmp.begin(), tmp.end()) );
	}
}

template< typename MatchVector >
void getBpList( MatchVector& mvect, uint seq, vector< gnSeqI >& bp_list )
{
	bp_list.clear();
	for( size_t ivI = 0; ivI < mvect.size(); ivI++ )
	{
		if( mvect[ivI]->LeftEnd(seq) == NO_MATCH )
			continue;
		bp_list.push_back( mvect[ivI]->LeftEnd(seq) );
		bp_list.push_back( mvect[ivI]->RightEnd(seq)+1 );
	}
	std::sort( bp_list.begin(), bp_list.end() );
}

template< typename MatchVector >
void createMap( const MatchVector& mv_from, const MatchVector& mv_to, vector< size_t >& map )
{
	typedef typename MatchVector::value_type MatchPtr;
	vector< pair< MatchPtr, size_t > > m1(mv_from.size());
	vector< pair< MatchPtr, size_t > > m2(mv_to.size());
	for( size_t i = 0; i < mv_from.size(); ++i )
		m1[i] = make_pair( mv_from[i], i );
	for( size_t i = 0; i < mv_to.size(); ++i )
		m2[i] = make_pair( mv_to[i], i );
	std::sort( m1.begin(), m1.end() );
	std::sort( m2.begin(), m2.end() );
	map.resize( m1.size() );
	for( size_t i = 0; i < m1.size(); ++i )
		map[m1[i].second] = m2[i].second;
}

typedef pair< size_t, Interval* > iv_tracker_t;
class IvTrackerComp
{
public:
	IvTrackerComp( uint seq ) : ssc( seq ) {}
	bool operator()( const iv_tracker_t& a, const iv_tracker_t& b )
	{
		return ssc(a.second, b.second);
	}
private:
	SingleStartComparator<Interval> ssc;
};

const int LEFT_NEIGHBOR = -1;
const int RIGHT_NEIGHBOR = 1;
typedef vector< size_t > neighbor_t;

neighbor_t& getNeighbor( pair< neighbor_t, neighbor_t >& entry, int direction )
{
	if( direction == RIGHT_NEIGHBOR )
		return entry.first;
	else
		return entry.second;
}


void collapseCollinear( IntervalList& iv_list )
{
	const size_t seq_count = iv_list.seq_table.size();
	std::vector< Interval* > iv_ptrs(iv_list.size());
	size_t lilI = 0;
	for( size_t i = 0; i < iv_list.size(); ++i )
	{
		// ignore unaligned regions
		if( iv_list[i].Multiplicity() < 2 )
			continue;
		iv_ptrs[lilI++] = &iv_list[i];
	}
	iv_ptrs.resize(lilI);
	const size_t NEIGHBOR_UNKNOWN = (std::numeric_limits<size_t>::max)();
	neighbor_t lefties_tmp( seq_count, NEIGHBOR_UNKNOWN );
	pair< neighbor_t, neighbor_t > neighbor_pair( lefties_tmp, lefties_tmp );
	vector< pair< neighbor_t, neighbor_t > > neighbor_list( iv_ptrs.size(), neighbor_pair );
	vector< iv_tracker_t > iv_tracker( iv_ptrs.size() );
	for( size_t i = 0; i < iv_ptrs.size(); ++i )
	{
		iv_tracker[i] = make_pair( i, iv_ptrs[i] );
	}
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		IvTrackerComp ivc( seqI );
		sort( iv_tracker.begin(), iv_tracker.end(), ivc );
		size_t prev_i = NEIGHBOR_UNKNOWN;
		size_t cur_i = NEIGHBOR_UNKNOWN;
		for( size_t i = 0; i < iv_tracker.size(); ++i )
		{
			if( iv_tracker[i].second->LeftEnd(seqI) == NO_MATCH )
				continue;
			if( cur_i != NEIGHBOR_UNKNOWN )
			{
				neighbor_list[cur_i].first[seqI] = prev_i;
				neighbor_list[cur_i].second[seqI] = iv_tracker[i].first;
			}
			prev_i = cur_i;
			cur_i = iv_tracker[i].first;
		}
		// get the last one
		neighbor_list[cur_i].first[seqI] = prev_i;
		neighbor_list[cur_i].second[seqI] = NEIGHBOR_UNKNOWN;
	}

	// now look for neighbor pair entries which can be merged
	for( int d = -1; d < 2; d+= 2 )	// iterate over both directions
	{
		size_t unknown_count = 0;
		for( size_t nI = 0; nI < neighbor_list.size(); ++nI )
		{
			size_t nayb = NEIGHBOR_UNKNOWN;
			size_t seqI = 0;
			bool parity = false;
			size_t ct = 0;
			for( ; seqI < seq_count; ++seqI )
			{
				if( iv_ptrs[nI]->Orientation(seqI) == AbstractMatch::undefined )
					continue;
				int orient = iv_ptrs[nI]->Orientation(seqI) == AbstractMatch::forward ? 1 : -1;

				if( nayb == NEIGHBOR_UNKNOWN )
				{
					nayb = getNeighbor( neighbor_list[nI], d * orient * -1 )[seqI];
					if( nayb != NEIGHBOR_UNKNOWN )
						parity = iv_ptrs[nI]->Orientation(seqI) == iv_ptrs[nayb]->Orientation(seqI);
				}
				else if( nayb != getNeighbor( neighbor_list[nI], d * orient * -1 )[seqI] )
					break;
				else if( parity != (iv_ptrs[nI]->Orientation(seqI) == iv_ptrs[nayb]->Orientation(seqI)) )
					break;
				if( nayb != NEIGHBOR_UNKNOWN )
					ct++;
			}
			if( seqI < seq_count || ct < iv_ptrs[nI]->Multiplicity() )
				continue;	// not collinear
			if( nayb == NEIGHBOR_UNKNOWN )
				continue;

			// merge nI and nayb
			uint fs = iv_ptrs[nI]->FirstStart();
			gnSeqI nI_lend_fs = iv_ptrs[nI]->LeftEnd(fs);
			gnSeqI nayb_lend_fs = iv_ptrs[nayb]->LeftEnd(fs);
			AbstractMatch::orientation o = iv_ptrs[nI]->Orientation(fs);
			vector< AbstractMatch* > nI_matches;
			iv_ptrs[nI]->StealMatches( nI_matches );
			vector< AbstractMatch* > nayb_matches;
			iv_ptrs[nayb]->StealMatches( nayb_matches );
			if( !parity )
			{
				std::reverse( nI_matches.begin(), nI_matches.end() );
				for( size_t i = 0; i < nI_matches.size(); ++i )
					nI_matches[i]->Invert();
				o = o == AbstractMatch::forward ? AbstractMatch::reverse : AbstractMatch::forward;
			}
			if( (o == AbstractMatch::forward && nI_lend_fs > nayb_lend_fs) ||
				(o == AbstractMatch::reverse && nI_lend_fs < nayb_lend_fs))
				nayb_matches.insert( nayb_matches.end(), nI_matches.begin(), nI_matches.end() );
			else
				nayb_matches.insert( nayb_matches.begin(), nI_matches.begin(), nI_matches.end() );

			iv_ptrs[nayb]->SetMatches( nayb_matches );

			// update all pointers to point to nayb
			seqI = 0;
			for( ; seqI < seq_count; ++seqI )
			{
				if( getNeighbor( neighbor_list[nI], -1 )[seqI] == NEIGHBOR_UNKNOWN &&
					getNeighbor( neighbor_list[nI], 1 )[seqI] == NEIGHBOR_UNKNOWN )
					continue;
				int orient = iv_ptrs[nayb]->Orientation(seqI) == AbstractMatch::forward ? 1 : -1;
				size_t other_nayb = getNeighbor( neighbor_list[nI], d * orient * (parity ? 1 : -1) )[seqI];
				if( other_nayb != NEIGHBOR_UNKNOWN )
				{
					if( getNeighbor( neighbor_list[other_nayb], 1 )[seqI] == nI )
						getNeighbor( neighbor_list[other_nayb], 1 )[seqI] = nayb;
					else if( getNeighbor( neighbor_list[other_nayb], -1 )[seqI] == nI )
						getNeighbor( neighbor_list[other_nayb], -1 )[seqI] = nayb;
					else
					{
						cerr << "serious programmer error\n";
						genome::breakHere();
					}
				}
				if( getNeighbor( neighbor_list[nayb], 1 )[seqI] == nI )
					getNeighbor( neighbor_list[nayb], 1 )[seqI] = other_nayb;
				else if( getNeighbor( neighbor_list[nayb], -1 )[seqI] == nI )
					getNeighbor( neighbor_list[nayb], -1 )[seqI] = other_nayb;
				else
				{
					cerr << "inexcusable programmer error\n";
					genome::breakHere();
				}
				neighbor_list[nI].first[seqI] = NEIGHBOR_UNKNOWN;
				neighbor_list[nI].second[seqI] = NEIGHBOR_UNKNOWN;
			}
		}
	}

	IntervalList new_list;
	new_list.seq_filename = iv_list.seq_filename;
	new_list.seq_table = iv_list.seq_table;
	new_list.resize( iv_ptrs.size() );
	size_t newI = 0;
	for( size_t ivI = 0; ivI < iv_ptrs.size(); ++ivI )
	{
		vector< AbstractMatch* > matches;
		iv_ptrs[ivI]->StealMatches( matches );
		if( matches.size() > 0 )
			new_list[newI++].SetMatches( matches );
	}
	new_list.resize(newI);
	swap( iv_list, new_list );
	addUnalignedRegions(iv_list);
}


void applyIslands( IntervalList& iv_list, backbone_list_t& bb_list, score_t score_threshold )
{
	// collapse any intervals that are trivially collinear
	collapseCollinear( iv_list );

	uint seq_count = iv_list.seq_table.size();
	boost::multi_array< vector< CompactGappedAlignment<>* >, 2> island_array;
	// indexed by seqI, seqJ, ivI, hssI (left col, right col)
	boost::multi_array< vector< pair< size_t, size_t > >, 3 > hss_cols(boost::extents[seq_count][seq_count][iv_list.size()]);

	// ugg.  need CompactGappedAlignment for its SeqPosToColumn
	vector< CompactGappedAlignment<>* > iv_ptrs(iv_list.size());
	for( size_t i = 0; i < iv_list.size(); ++i )
	{
		CompactGappedAlignment<> tmp_cga;
		iv_ptrs[i] = tmp_cga.Copy();
		new (iv_ptrs[i])CompactGappedAlignment<>( iv_list[i] );
	}
	vector< CompactGappedAlignment<>* > iv_orig_ptrs(iv_ptrs);

	// make pairwise projections of intervals and find LCBs...
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		for( size_t seqJ = seqI+1; seqJ < seq_count; ++seqJ )
		{
			vector< uint > projection;
			projection.push_back( seqI );
			projection.push_back( seqJ );
			vector< vector< MatchProjectionAdapter* > > LCB_list;
			vector< LCB > projected_adjs;
			projectIntervalList( iv_list, projection, LCB_list, projected_adjs );
			// make intervals
			IntervalList pair_ivs;
			pair_ivs.seq_table.push_back( iv_list.seq_table[seqI] );
			pair_ivs.seq_table.push_back( iv_list.seq_table[seqJ] );
			pair_ivs.resize( LCB_list.size() );
			for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
				pair_ivs[lcbI].SetMatches( LCB_list[lcbI] );
			LCB_list.clear();

			vector< CompactGappedAlignment<>* > pair_cgas( pair_ivs.size() );
			for( size_t lcbI = 0; lcbI < pair_ivs.size(); ++lcbI )
			{
				CompactGappedAlignment<> tmp_cga;
				pair_cgas[lcbI] = tmp_cga.Copy();
				new (pair_cgas[lcbI])CompactGappedAlignment<>( pair_ivs[lcbI] );
			}

			vector< CompactGappedAlignment<>* > hss_list;
			// now find islands
			findHssRandomWalkCga( pair_cgas, pair_ivs.seq_table, score_threshold, hss_list );
			for( size_t cgaI = 0; cgaI < pair_cgas.size(); ++cgaI )
				pair_cgas[cgaI]->Free();
			pair_cgas.clear();

			// now split up on iv boundaries
			vector< gnSeqI > bp_list;
			getBpList( iv_ptrs, seqI, bp_list );
			GenericMatchSeqManipulator< CompactGappedAlignment<> > gmsm(0);
			SingleStartComparator< CompactGappedAlignment<> > ssc(0);
			std::sort(hss_list.begin(), hss_list.end(), ssc );
			applyBreakpoints( bp_list, hss_list, gmsm );
			// and again on seqJ
			getBpList( iv_ptrs, seqJ, bp_list );
			GenericMatchSeqManipulator< CompactGappedAlignment<> > gmsm1(1);
			SingleStartComparator< CompactGappedAlignment<> > ssc1(1);
			std::sort(hss_list.begin(), hss_list.end(), ssc1 );
			applyBreakpoints( bp_list, hss_list, gmsm1 );

			// now transform into interval-specific columns
			std::sort(hss_list.begin(), hss_list.end(), ssc );

			SingleStartComparator< CompactGappedAlignment<> > ivcomp(seqI);
			std::sort( iv_ptrs.begin(), iv_ptrs.end(), ivcomp );
			vector< size_t > iv_map;
			createMap( iv_ptrs, iv_orig_ptrs, iv_map );
			size_t ivI = 0;
			while( ivI < iv_ptrs.size() && iv_ptrs[ivI]->LeftEnd(0) == NO_MATCH )
				++ivI;
			for( size_t hssI = 0; hssI < hss_list.size(); ++hssI )
			{
				if( hss_list[hssI]->LeftEnd(0) == NO_MATCH || hss_list[hssI]->Length(0) == 0 )
					continue;
				if( ivI == iv_ptrs.size() )
				{
					cerr << "huh?\n";
					cerr << hss_list[hssI]->LeftEnd(0) << endl;
					cerr << hss_list[hssI]->RightEnd(0) << endl;
					cerr << iv_ptrs.back()->LeftEnd(seqI) << endl;
					cerr << iv_ptrs.back()->RightEnd(seqI) << endl;
				}
				while( ivI < iv_ptrs.size() && 
					(iv_ptrs[ivI]->LeftEnd(seqI) == NO_MATCH ||
					hss_list[hssI]->LeftEnd(0) > iv_ptrs[ivI]->RightEnd(seqI) ) )
					++ivI;
				if( ivI == iv_ptrs.size() )
				{
					cerr << "hssI fit!!\n";
					genome::breakHere();
				}
				// check for containment in seqJ
				if( iv_ptrs[ivI]->LeftEnd(seqJ) == NO_MATCH ||
					iv_ptrs[ivI]->RightEnd(seqJ) < hss_list[hssI]->LeftEnd(1) ||
					hss_list[hssI]->RightEnd(1) < iv_ptrs[ivI]->LeftEnd(seqJ) )
					continue;	// this hss falls to an invalid range in seqJ

				if( hss_list[hssI]->RightEnd(0) < iv_ptrs[ivI]->LeftEnd(seqI) )
				{
					cerr << "huh 2?\n";
					cerr << hss_list[hssI]->LeftEnd(0) << endl;
					cerr << hss_list[hssI]->RightEnd(0) << endl;
					cerr << iv_ptrs[ivI]->LeftEnd(seqI) << endl;
					cerr << iv_ptrs[ivI]->RightEnd(seqI) << endl;
					hssI++;
					continue;
				}

				vector< pair< size_t, size_t > >& cur_hss_cols = hss_cols[seqI][seqJ][iv_map[ivI]];

				gnSeqI left_col = iv_ptrs[ivI]->SeqPosToColumn( seqI, hss_list[hssI]->LeftEnd(0) );
				gnSeqI right_col = iv_ptrs[ivI]->SeqPosToColumn( seqI, hss_list[hssI]->RightEnd(0) );
				if(left_col > right_col && iv_ptrs[ivI]->Orientation(seqI) == AbstractMatch::reverse )
				{
					swap(left_col, right_col);	// must have been a revcomp seq
				}
				else if(left_col > right_col)
				{
					cerr << "bad cols\n";
					cerr << hss_list[hssI]->LeftEnd(0) << endl;
					cerr << hss_list[hssI]->RightEnd(0) << endl;
					cerr << iv_ptrs[ivI]->LeftEnd(seqI) << endl;
					cerr << iv_ptrs[ivI]->RightEnd(seqI) << endl;
					genome::breakHere();
				}

				if( left_col > 2000000000 || right_col > 2000000000 )
				{
					cerr << "huh 2?\n";
					cerr << hss_list[hssI]->LeftEnd(0) << endl;
					cerr << hss_list[hssI]->RightEnd(0) << endl;
					cerr << iv_ptrs[ivI]->LeftEnd(seqI) << endl;
					cerr << iv_ptrs[ivI]->RightEnd(seqI) << endl;
					genome::breakHere();
				}
				cur_hss_cols.push_back( make_pair( left_col, right_col ) );
			}
			for( size_t hssI = 0; hssI < hss_list.size(); ++hssI )
				hss_list[hssI]->Free();
		}
	}


	//
	// FINALLY!  ready to merge.  how to do it?
	// make an empty list of UngappedLocalAlignments
	// start with the first seq and create a ULA for every col
	// range.  Then continue to the second seq, and when
	// a col range overlaps a pre-existing ULA, create a new ULA
	// for the intersected region and a smaller ULA for the non-intersected region
	vector< vector< ULA* > > ula_list( iv_list.size() );
	for( size_t ivI = 0; ivI < iv_list.size(); ++ivI )
	{
		vector< ULA* >& iv_ulas = ula_list[ivI];
		for( size_t seqI = 0; seqI < seq_count; ++seqI )
		{
			for( size_t seqJ = seqI+1; seqJ < seq_count; ++seqJ )
			{
				vector< pair< size_t, size_t > >& cur_hss_cols = hss_cols[seqI][seqJ][ivI];
				vector< ULA* > cur_ulas( cur_hss_cols.size() );
				ULA tmp_ula(seq_count);
				for( size_t hssI = 0; hssI < cur_hss_cols.size(); ++hssI )
				{
					cur_ulas[hssI] = tmp_ula.Copy();
					cur_ulas[hssI]->SetStart(seqI, cur_hss_cols[hssI].first+1);
					cur_ulas[hssI]->SetStart(seqJ, cur_hss_cols[hssI].first+1);
					cur_ulas[hssI]->SetLength( cur_hss_cols[hssI].second - cur_hss_cols[hssI].first + 1 );
				}

				vector< gnSeqI > iv_bp_list;
				vector< gnSeqI > cur_bp_list;
				SingleStartComparator<ULA> ulacompI(seqI);
				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompI );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompI );
				getBpList( iv_ulas, seqI, iv_bp_list );
				getBpList( cur_ulas, seqI, cur_bp_list );
				GenericMatchSeqManipulator< ULA > gmsm(seqI);
				applyBreakpoints( iv_bp_list, cur_ulas, gmsm );
				applyBreakpoints( cur_bp_list, iv_ulas, gmsm );

				SingleStartComparator<ULA> ulacompJ(seqJ);
				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompJ );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompJ );
				getBpList( iv_ulas, seqJ, iv_bp_list );
				getBpList( cur_ulas, seqJ, cur_bp_list );
				GenericMatchSeqManipulator< ULA > gmsmJ(seqJ);
				applyBreakpoints( iv_bp_list, cur_ulas, gmsmJ );
				applyBreakpoints( cur_bp_list, iv_ulas, gmsmJ );

				// do seqI a second time to propagate any breakpoints introduced by seqJ
				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompI );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompI );
				getBpList( iv_ulas, seqI, iv_bp_list );
				getBpList( cur_ulas, seqI, cur_bp_list );
				applyBreakpoints( iv_bp_list, cur_ulas, gmsm );
				applyBreakpoints( cur_bp_list, iv_ulas, gmsm );

				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompI );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompI );
				// now that cur_ulas and iv_ulas are all broken up according to each other's boundaries
				// we can simply scan along and add
				size_t iv_ulas_size = iv_ulas.size();
				size_t ivuI = 0;
				size_t curuI = 0;
				vector< ULA* > added_to( cur_ulas.size(), NULL );	// this tracks which of iv_ulas a cur_ula was added to
				vector< ULA* > to_delete;
				while( ivuI < iv_ulas_size && curuI < cur_ulas.size() )
				{
					if( iv_ulas[ivuI]->LeftEnd(seqI) == cur_ulas[curuI]->LeftEnd(seqI) )
					{
						if( added_to[curuI] == iv_ulas[ivuI] )
						{
							// do nothing
						}else if( added_to[curuI] == NULL )
						{
							iv_ulas[ivuI]->SetLeftEnd(seqJ, cur_ulas[curuI]->LeftEnd(seqJ));
							added_to[curuI] = iv_ulas[ivuI];
						}else{
							ULA* merge = added_to[curuI];
							for( size_t seqK = 0; seqK < seq_count; ++seqK )
							{
								if( merge->Start(seqK) == NO_MATCH )
									continue;
								iv_ulas[ivuI]->SetStart( seqK, merge->Start(seqK) );
							}
							to_delete.push_back( merge );
						}
						ivuI++;
					}else if( iv_ulas[ivuI]->LeftEnd(seqI) < cur_ulas[curuI]->LeftEnd(seqI) )
					{
						ivuI++;
					}else
						curuI++;
				}

				// delete to_delete...
				std::sort( to_delete.begin(), to_delete.end() );
				vector< ULA* >::iterator last = std::unique( to_delete.begin(), to_delete.end() );
				to_delete.erase( last, to_delete.end() );
				vector< ULA* > new_iv_ulas( iv_ulas.size() - to_delete.size() );
				std::sort( iv_ulas.begin(), iv_ulas.end() );
				std::set_difference( iv_ulas.begin(), iv_ulas.end(), to_delete.begin(), to_delete.end(), new_iv_ulas.begin() );
				swap( iv_ulas, new_iv_ulas );
				for( size_t delI = 0; delI < to_delete.size(); ++delI )
					to_delete[delI]->Free();

				vector< ULA* > orig_ula_order = cur_ulas;
				// now do something similar for seqJ
				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompJ );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompJ );

				vector< size_t > added_map;
				createMap( cur_ulas, orig_ula_order, added_map );

				ivuI = 0;
				curuI = 0;
				to_delete.clear();
				while( ivuI < iv_ulas_size && curuI < cur_ulas.size() )
				{
					if( iv_ulas[ivuI]->LeftEnd(seqJ) == cur_ulas[curuI]->LeftEnd(seqJ) )
					{
						if( added_to[added_map[curuI]] == iv_ulas[ivuI] )
						{
							// do nothing
						}else if( added_to[added_map[curuI]] == NULL )
						{
							iv_ulas[ivuI]->SetLeftEnd(seqI, cur_ulas[curuI]->LeftEnd(seqI));
							added_to[added_map[curuI]] = iv_ulas[ivuI];
						}else{
							ULA* merge = added_to[added_map[curuI]];
							for( size_t seqK = 0; seqK < seq_count; ++seqK )
							{
								if( merge->Start(seqK) == NO_MATCH )
									continue;
								iv_ulas[ivuI]->SetStart( seqK, merge->Start(seqK) );
							}
							to_delete.push_back( merge );
						}
						ivuI++;
					}else if( iv_ulas[ivuI]->LeftEnd(seqJ) < cur_ulas[curuI]->LeftEnd(seqJ) )
					{
						ivuI++;
					}else
					{
						curuI++;
					}
				}

				// anything with a null added_to entry needs to be added to iv_ulas
				// everything else needs to get freed
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompI );
				for( curuI = 0; curuI < cur_ulas.size(); ++curuI )
				{
					if( added_to[curuI] == NULL )
						iv_ulas.push_back( cur_ulas[curuI] );
					else
						cur_ulas[curuI]->Free();
				}
				// delete to_delete...
				std::sort( to_delete.begin(), to_delete.end() );
				last = std::unique( to_delete.begin(), to_delete.end() );
				to_delete.erase( last, to_delete.end() );
				new_iv_ulas = vector< ULA* >( iv_ulas.size() - to_delete.size() );
				std::sort( iv_ulas.begin(), iv_ulas.end() );
				std::set_difference( iv_ulas.begin(), iv_ulas.end(), to_delete.begin(), to_delete.end(), new_iv_ulas.begin() );
				swap( iv_ulas, new_iv_ulas );
				for( size_t delI = 0; delI < to_delete.size(); ++delI )
					to_delete[delI]->Free();
			}
		}
	}

	// Eliminate segments that have no representation in a genome
	for( size_t ivI = 0; ivI < ula_list.size(); ++ivI )
	{
		for( size_t mI = 0; mI < ula_list[ivI].size(); ++mI )
		{
			size_t seqI = ula_list[ivI][mI]->FirstStart();
			std::vector<gnSeqI> l_pos;
			std::vector<bool> l_column;
			std::vector<gnSeqI> r_pos;
			std::vector<bool> r_column;
			gnSeqI left_col = ula_list[ivI][mI]->LeftEnd(seqI)-1;
			gnSeqI right_col = ula_list[ivI][mI]->RightEnd(seqI)-1;
			iv_orig_ptrs[ivI]->GetColumn(left_col, l_pos, l_column);
			iv_orig_ptrs[ivI]->GetColumn(right_col, r_pos, r_column);
			for( ; seqI < ula_list[ivI][mI]->SeqCount(); ++seqI )
			{
				if( ula_list[ivI][mI]->LeftEnd(seqI) == NO_MATCH )
					continue;
				if( l_pos[seqI] == r_pos[seqI] && !l_column[seqI] && !r_column[seqI] )
					ula_list[ivI][mI]->SetStart(seqI, NO_MATCH);	// no match in this col
			}
			if( ula_list[ivI][mI]->Multiplicity() < 2 )
			{
				ula_list[ivI][mI]->Free();
				ula_list[ivI][mI] = NULL;
			}
		}
		// clean out any NULL ptrs
		std::vector< ULA* >::iterator last = std::remove( ula_list[ivI].begin(), ula_list[ivI].end(), (ULA*)NULL );
		ula_list[ivI].erase( last, ula_list[ivI].end() );
	}

	// unalign regions in the iv list that aren't contained in backbone
	for( size_t ivI = 0; ivI < ula_list.size(); ++ivI )
	{
		vector< AbstractMatch* > new_matches(ula_list[ivI].size());
		for( size_t mI = 0; mI < ula_list[ivI].size(); ++mI )
		{
			size_t seqI = ula_list[ivI][mI]->FirstStart();
			gnSeqI left_col = ula_list[ivI][mI]->LeftEnd(seqI)-1;
			CompactGappedAlignment<> tmp_cga;
			CompactGappedAlignment<>* new_cga = tmp_cga.Copy();
			iv_orig_ptrs[ivI]->copyRange(*new_cga, left_col, ula_list[ivI][mI]->Length(seqI));
			for( seqI = 0; seqI < ula_list[ivI][mI]->SeqCount(); ++seqI )
			{
				if( ula_list[ivI][mI]->LeftEnd(seqI) == NO_MATCH )
					new_cga->SetLeftEnd(seqI, NO_MATCH);
			}
			new_matches[mI] = new_cga;
		}
		if( new_matches.size() > 0 )
		{

			vector< vector< AbstractMatch* > > disjoint_subsets;
			{
				// split into multiple intervals if some sequences are completely unaligned
				// use a union-find structure to quickly figure out how many subgroups there are
				vector< uint > seq_map( seq_count );
				for( size_t sI = 0; sI < seq_map.size(); ++sI )
					seq_map[sI] = sI;
				for( size_t mI = 0; mI < new_matches.size(); ++mI )
				{
					uint sI = new_matches[mI]->FirstStart();
					uint map_to = seq_map[sI];
					while( map_to != seq_map[map_to] )
						map_to = seq_map[map_to];
					seq_map[sI] = map_to;
					for( ++sI; sI < seq_count; ++sI )
					{
						if( new_matches[mI]->LeftEnd(sI) == NO_MATCH )
							continue;
						uint map_from = seq_map[sI];
						while( map_from != seq_map[map_from] )
							map_from = seq_map[map_from];
						seq_map[map_from] = map_to;
					}
				}
				vector< vector< AbstractMatch* > > mapped_lists( seq_count );
				for( size_t mI = 0; mI < new_matches.size(); ++mI )
				{
					uint sI = new_matches[mI]->FirstStart();
					uint map_to = seq_map[sI];
					while( map_to != seq_map[map_to] )
						map_to = seq_map[map_to];
					mapped_lists[map_to].push_back( new_matches[mI] );
				}
				for( uint sI = 0; sI < seq_count; ++sI )
				{
					if( mapped_lists[sI].size() > 0 )
						disjoint_subsets.push_back( mapped_lists[sI] );
				}
			}

			for( size_t dI = 0; dI < disjoint_subsets.size(); ++dI )
			{
				vector< AbstractMatch* >& cur_d_matches = disjoint_subsets[dI];
				vector< AbstractMatch* > orig_order = cur_d_matches;
				// need to sort these.  use boost's topological sort.
				vector< size_t > id_map;
				typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, boost::property<boost::vertex_color_t, boost::default_color_type> > Graph;
				typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
				typedef std::pair< int, int > Pair;
				vector< Pair > edges;
				for( size_t seqI = 0; seqI < seq_count; ++seqI )
				{
					SingleStartComparator<AbstractMatch> ssc(seqI);
					std::sort( cur_d_matches.begin(), cur_d_matches.end(), ssc );
					createMap( cur_d_matches, orig_order, id_map );
					int prev = -1;
					int first = -1;
					bool reverse = false;
					for( int mI = 0; mI < cur_d_matches.size(); ++mI )
					{
						if( cur_d_matches[mI]->LeftEnd(seqI) == NO_MATCH )
							continue;
						if( prev != -1 )
						{
							Pair edge( id_map[prev], id_map[mI] );
							if( reverse )
								swap( edge.first, edge.second );
							edges.push_back(edge);
						}else
						{
							reverse = cur_d_matches[mI]->Start(seqI) < 0;
							first = mI;
						}
						prev = mI;
					}
					if( prev != -1 && !reverse )
						edges.push_back( Pair( id_map[prev], cur_d_matches.size() ) );
					else if( prev != -1 && reverse )
						edges.push_back( Pair( id_map[first], cur_d_matches.size() ) );
				}
				std::sort( edges.begin(), edges.end() );
				vector< Pair >::iterator ee_iter = std::unique( edges.begin(), edges.end() );
				edges.erase( ee_iter, edges.end() );
				Pair* edge_array = new Pair[edges.size()];
				for( size_t eI = 0; eI < edges.size(); ++eI )
					edge_array[eI] = edges[eI];
				typedef boost::graph_traits<Graph>::vertices_size_type v_size_t;
				Graph G(edge_array, edge_array + edges.size(), v_size_t(edges.size()));
				typedef std::vector< Vertex > container;
				container c;
				topological_sort(G, std::back_inserter(c));
				cur_d_matches.clear();
				for ( container::reverse_iterator ii=c.rbegin(); ii!=c.rend(); ++ii)
				{
					if( *ii < orig_order.size() )
						cur_d_matches.push_back( orig_order[ *ii ] );
				}
				if( dI == 0 )
					iv_list[ivI].SetMatches(cur_d_matches);
				else
				{
					Interval new_iv( cur_d_matches.begin(), cur_d_matches.end() );
					iv_list.push_back(new_iv);
				}
				delete[] edge_array;
			}
		}
		else
		{
			iv_orig_ptrs[ivI]->Free();
			iv_orig_ptrs[ivI] = NULL;
		}
	}


	// update iv_list to match the filtered iv_orig_ptrs
	size_t givI = 0;
	for( size_t iI = 0; iI < iv_orig_ptrs.size(); ++iI )
	{
		if( iv_orig_ptrs[iI] != NULL )
		{
			swap( iv_list[givI], iv_list[iI] );
			iv_orig_ptrs[iI]->Free();	// done with the CompactGappedAlignments
			iv_orig_ptrs[iI] = NULL;
			givI++;
		}
	}
	// pick up any intervals that were split in half
	for( size_t iI = iv_orig_ptrs.size(); iI < iv_list.size(); ++iI )
		swap( iv_list[givI++], iv_list[iI] );
	iv_list.erase( iv_list.begin()+givI, iv_list.end() );

	// collapse any intervals that are trivially collinear
	collapseCollinear( iv_list );

	// need to add in all the unaligned regions so the viewer doesn't throw a fit
	addUnalignedRegions( iv_list );

	// free all ULAs and reconstruct them from the new alignment column coordinates
	for( size_t ulaI = 0; ulaI < ula_list.size(); ++ulaI )
		for( size_t i = 0; i < ula_list[ulaI].size(); ++i )
			ula_list[ulaI][i]->Free();
	ula_list.clear();


	ula_list.resize( iv_list.size() );
	for( size_t ivI = 0; ivI < iv_list.size(); ++ivI )
	{
		if( iv_list[ivI].Multiplicity() < 2 )
			continue;
		const vector< AbstractMatch* >& matches = iv_list[ivI].GetMatches();
		int64 right_col = 0;
		int64 left_col = 0;
		for( size_t mI = 0; mI < matches.size(); ++mI )
		{
			left_col = right_col;
			right_col += matches[mI]->AlignmentLength();
			if( matches[mI]->Multiplicity() < 2 )
				continue;
			ULA tmp_ula(matches[mI]->SeqCount());
			ULA* mula = tmp_ula.Copy();
			for( size_t seqI = 0; seqI < matches[mI]->SeqCount(); ++seqI )
				if( matches[mI]->LeftEnd(seqI) != NO_MATCH )
					mula->SetLeftEnd( seqI, left_col+1 );
			mula->SetLength( right_col - left_col );
			ula_list[ivI].push_back(mula);
		}
	}

	iv_orig_ptrs.clear();

	bb_list.clear();
	bb_list = ula_list;

	// debug: sanity check whether there are all gap columns
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
	{
		vector< string > aln;
		mems::GetAlignment( iv_list[ivI], iv_list.seq_table, aln );
		for( size_t colI = 0; colI < aln[0].size(); ++colI )
		{
			size_t rowI = 0;
			for( ; rowI < aln.size(); ++rowI )
				if( aln[rowI][colI] != '-' )
					break;
			if( rowI == aln.size() )
			{
				cerr << "ERROR!  IV " << ivI << " COLUMN " << colI << " IS ALL GAPS!\n";
			}
		}
	}
}

void writeBackboneColumns( ostream& bb_out, backbone_list_t& bb_list )
{
	//
	// At last! write out the backbone list
	//
	for( size_t ivI = 0; ivI < bb_list.size(); ++ivI )
	{
		for( size_t mI = 0; mI < bb_list[ivI].size(); ++mI )
		{
			size_t seqI = bb_list[ivI][mI]->FirstStart();
			bb_out << ivI << '\t' << bb_list[ivI][mI]->LeftEnd(seqI) << '\t' << bb_list[ivI][mI]->Length();
			for( ; seqI < bb_list[ivI][mI]->SeqCount(); ++seqI )
			{
				if( bb_list[ivI][mI]->LeftEnd(seqI) == NO_MATCH )
					continue;
				bb_out << '\t' << seqI;
			}
			bb_out << endl;
		}
	}
}

void writeBackboneSeqCoordinates( backbone_list_t& bb_list, IntervalList& iv_list, ostream& bb_out )
{
	if( bb_list.size() == 0 )
		return;
	// find seq_count
	uint seq_count = 0;
	for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
		if( bb_list[bbI].size() > 0 )
		{
			seq_count = bb_list[bbI].front()->SeqCount();
			break;
		}

	// different format -- use real sequence coordinates...
	// print a header line first
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		if( seqI > 0 )
			bb_out << '\t';
		bb_out << "seq_" << seqI << "_leftend\t";
		bb_out << "seq_" << seqI << "_rightend";
	}
	bb_out << endl;
	for( size_t ivI = 0; ivI < bb_list.size(); ++ivI )
	{
		// there seems to be a bug in the backbone creation code that causes the CGA that gets
		// stuffed into the interval to have the wrong coordinates internally, while the interval
		// maintains the correct coordinates.  work around it by converting the whole interval to a cga
		CompactGappedAlignment<> iv_cga( iv_list[ivI] );
		for( size_t mI = 0; mI < bb_list[ivI].size(); ++mI )
		{
			uint fs = bb_list[ivI][mI]->FirstStart();
			// get the sequence positions out of the alignment
			vector< gnSeqI > left_pos;
			vector< bool > left_cols;
			iv_cga.GetColumn( bb_list[ivI][mI]->LeftEnd(fs)-1, left_pos, left_cols );
			vector< gnSeqI > right_pos;
			vector< bool > right_cols;
			iv_cga.GetColumn( bb_list[ivI][mI]->RightEnd(fs)-1, right_pos, right_cols );
			for( size_t seqI = 0; seqI < bb_list[ivI][mI]->SeqCount(); ++seqI )
			{
				if( seqI > 0 )
					bb_out << '\t';
				if( bb_list[ivI][mI]->LeftEnd(seqI) == NO_MATCH )
				{
					bb_out << "0\t0";
					continue;
				}else{
					int64 leftI = left_pos[seqI];
					int64 rightI = right_pos[seqI];
					if( iv_cga.Orientation(seqI) == AbstractMatch::forward && leftI != 0 && !left_cols[seqI] )
						leftI++;
					if( iv_cga.Orientation(seqI) == AbstractMatch::reverse && rightI != 0 && !right_cols[seqI] )
						rightI++;
					if( iv_cga.Orientation(seqI) == AbstractMatch::reverse )
					{
						swap( leftI, rightI );	// must be reverse complement
					}
					if( rightI + 1 == leftI )
					{
						bb_out << "0\t0";
						continue;
					}
					if( leftI > rightI )
					{
						cerr << "oh crahpey!\n";
						cerr << "leftI: " << leftI << endl;
						cerr << "rightI: " << rightI << endl;
						cerr << "seqI: " << seqI << endl;
						cerr << "ivI: " << ivI << endl;
					}
					if( leftI == 0 )
						leftI = iv_cga.LeftEnd(seqI);
					if( rightI == iv_cga.RightEnd(seqI)+1 )
						rightI--;
					if( iv_cga.Orientation(seqI) == AbstractMatch::reverse )
					{
						leftI *= -1;
						rightI *= -1;
					}
					bb_out << leftI << '\t' << rightI;
				}
			}
			bb_out << endl;
		}
	}
}
/**
 * 
 */

void ProgressiveAligner::align( vector< gnSequence* >& seq_table, IntervalList& interval_list ){
	if( debug_aligner )
	{
		debug_interval = true;
		debug_cga = true;
	}

	seq_count = seq_table.size();
	this->currently_recursing = false;
	interval_list.seq_table = seq_table;

	// find pairwise matches
	MatchList mlist;
	mlist.seq_table = seq_table;

	if( this->breakpoint_penalty == -1 )
		this->breakpoint_penalty = getDefaultBreakpointPenalty( seq_table );
	if( this->bp_dist_estimate_score == -1 )
		this->bp_dist_estimate_score = getDefaultBpDistEstimateMinScore( original_ml.seq_table );
	cout << "using default bp penalty: " << breakpoint_penalty << endl;
	cout << "using default bp estimate min score: " << bp_dist_estimate_score << endl;

	if( this->collinear_genomes )
		this->breakpoint_penalty = -1;

	if( collinear_genomes )
		cout << "\nAssuming collinear genomes...\n";

	mlist.clear();
	mlist = original_ml;
	cout << "Starting with " << mlist.size() << " multi-matches\n";

//
// Step 2) Compute a phylogenetic guide tree using the pairwise matches
//
	NumericMatrix< double > distance;
	SingleCopyDistanceMatrix( mlist, mlist.seq_table, distance );
	cout << "\n\nGenome conservation distance matrix: " << endl;
	distance.print(cout);
	cout << endl;

	bool input_tree_specified = input_guide_tree_fname != "";
	bool output_tree_specified = output_guide_tree_fname != "";
	if( !input_tree_specified )
	{
		ClustalInterface& ci = ClustalInterface::getClustalInterface();	
		// Make a phylogenetic guide tree
		if( !output_tree_specified )
			output_guide_tree_fname = CreateTempFileName("guide_tree");
		input_guide_tree_fname = output_guide_tree_fname;
		cout << "Writing guide tree to " << output_guide_tree_fname << endl;
		ci.SetDistanceMatrix( distance, output_guide_tree_fname );
	}

	conservation_distance.resize(boost::extents[seq_count][seq_count]);
	for( uint seqI = 0; seqI < seq_count; ++seqI )
		for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
			if( seqJ > seqI )
			{
				conservation_distance[seqI][seqJ] = distance(seqI,seqJ);
				conservation_distance[seqJ][seqI] = distance(seqI,seqJ);
			}
			else
			{
				conservation_distance[seqI][seqJ] = distance(seqJ,seqI);
				conservation_distance[seqJ][seqI] = distance(seqJ,seqI);
			}

	cout << "reading tree...\n";
	// load the guide tree
	ifstream tree_file( input_guide_tree_fname.c_str() );
	if( !tree_file.is_open() )
		throw "Error opening guide tree file";
	alignment_tree.readTree( tree_file );
	tree_file.close();

	cout << "initializing alignment tree...\n";

	makeAlignmentTree( alignment_tree, mlist, node_sequence_map );
	node_id_t node1;
	node_id_t node2;
	// midpoint root the tree
	findMidpoint( alignment_tree, node1, node2 );
	node_id_t ancestor = 0;
	if( seq_count > 2 )	// if only two sequences then the tree already has a root
		ancestor = createAlignmentTreeRoot( alignment_tree, node1, node2 );

	// write out the rooted guide tree, but don't clobber the user's input tree
	if( !input_tree_specified || output_tree_specified )
	{
		ofstream out_tree_file( output_guide_tree_fname.c_str() );
		if( !out_tree_file.is_open() )
			throw "Error opening guide tree file for write";
		alignment_tree.writeTree( out_tree_file );
		out_tree_file.close();
	}

	// ensure the root is the last to get aligned and swap children to canonical order
	extendRootBranches(alignment_tree);


	if( !collinear_genomes )
	{
		// need sol lists for scoring
		cout << "Constructing seed occurrence lists\n";
		sol_list.resize(seq_count);
		for( uint seqI = 0; seqI < seq_count; seqI++ )
			sol_list[seqI].construct(*(original_ml.sml_table[seqI]));
	}
	if( !collinear_genomes && use_weight_scaling )
	{
		cout << "Calculating pairwise breakpoint distances\n";
		CreatePairwiseBPDistance(bp_distance);
	}

	// rescale the conservation distance
	if( use_weight_scaling )
	{
		for( uint seqI = 0; seqI < seq_count; ++seqI )
			for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
				conservation_distance[seqI][seqJ] = distance(seqI,seqJ) * conservation_dist_scale;
	}else{
		bp_distance.resize(boost::extents[seq_count][seq_count]);
		for( uint seqI = 0; seqI < seq_count; ++seqI )
			for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
			{
				conservation_distance[seqI][seqJ] = 0;
				bp_distance[seqI][seqJ] = 0;
			}
	}

	if( !collinear_genomes )
	{
		cout << "genome content distance matrix:\n";
		print2d_matrix(conservation_distance, cout);
		cout << endl;
		cout << "bp distance matrix:\n";
		print2d_matrix(bp_distance, cout);
		cout << endl;
	}

	getAlignment( interval_list );
}


// broken and unused function graveyard

