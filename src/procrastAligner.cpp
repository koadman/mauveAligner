#include "libGenome/gnSequence.h"
#include "libMems/Interval.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/Islands.h"
#include "libMems/Aligner.h"
#include "libMems/MuscleInterface.h"
#include "libGenome/gnFASSource.h"
#include "libMems/Backbone.h"
#include "libMems/ProgressiveAligner.h"

#include <iostream>
#include <algorithm>
#include <cctype>

#include "MatchRecord.h"
#include "SeedMatchEnumerator.h"
//#include "procrastUtilities.h"

#include <boost/tuple/tuple.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace genome;
using namespace mems;

enum rvalue { OK=0, FAILED=1, FIXME=110}; 
int scoredropoff_matrix[10] = {0,0,4756,9144,13471,17981,25302,30945,38361,40754};
int ccount = 0;
/** A Match Position Entry stores a pointer to a match and its component for a given sequence coordinate */
typedef std::pair< MatchRecord*, size_t >	MatchPositionEntry;
/** the Match Position Lookup Table should be sized match the length of the sequence */
typedef vector< MatchPositionEntry > MatchPositionLookupTable;

/** This class stores a single entry in the neighborhood list */
class NeighborhoodListEntry
{
public:
	MatchRecord* match;
	bool relative_orientation;	/** true for identical (1) and false for opposite (-1) */
	size_t Mi_component;	/** the x value in the paper (Matching component of M_i)*/
	size_t distance;		/** the d value from the paper */
	size_t Mj_component;	/** the y value in the paper (Matching component of M_j) */
};

/** Used to sort the neighborhood list using std::sort */
class NeighborhoodListComparator
{
public:
	bool operator()( const NeighborhoodListEntry& a, const NeighborhoodListEntry& b )
	{
		if( a.match != b.match )
			return a.match < b.match;
		if( a.relative_orientation != b.relative_orientation )
			return a.relative_orientation == false;
		if( a.Mi_component != b.Mi_component )
			return a.Mi_component < b.Mi_component;
		return a.distance < b.distance; 
	}
};


bool scorecmp( pair< double, GappedMatchRecord* > a, pair< double, GappedMatchRecord* > b ) 
{
   return a.first > b.first;
 }
/** The NeighborhoodGroup contains the MatchRecord pointer, the component map to the match being extended (M_i), and a vector of distances to M_i*/
typedef boost::tuple< MatchRecord*, std::vector< size_t >, std::vector< size_t > > NeighborhoodGroup;

class NeighborhoodGroupComponentCompare
{
public:
	bool operator()( const NeighborhoodGroup& a, const NeighborhoodGroup& b )
	{
		return compare(a,b) < 0;
	}
	int compare( const NeighborhoodGroup& a, const NeighborhoodGroup& b )
	{
	// compare component map vectors
		vector< size_t > ac(a.get<1>());
		vector< size_t > bc(b.get<1>());
		std::sort(ac.begin(), ac.end());
		std::sort(bc.begin(), bc.end());
		size_t i = 0;
		for( ; i < ac.size() && i < bc.size(); ++i )
		{
			if( ac[i] != bc[i] )
				return ac[i] - bc[i];
		}
		if( i < ac.size() && ac[i] != (std::numeric_limits<size_t>::max)())
			return 1;
		else if( i < bc.size() && bc[i] != (std::numeric_limits<size_t>::max)())
			return -1;

		return 0;
	}
};

class NeighborhoodGroupCompare
{
public:
	bool operator()( const NeighborhoodGroup& a, const NeighborhoodGroup& b )
	{
		int cval = srcc.compare(a,b);
		if( cval != 0 )
			return cval < 0;

	// compare distance vectors
		vector< size_t > ad(a.get<2>());
		vector< size_t > bd(b.get<2>());
		std::sort(ad.begin(), ad.end());
		std::sort(bd.begin(), bd.end());
		size_t i = 0;
		for( ; i < ad.size() && i < bd.size(); ++i )
		{
			if( ad[i] != bd[i] )
				return ad[i] < bd[i];
		}
		if( i < ad.size() )
			return false;
		else if( i < bd.size() )
			return true;

		return false;
	}
protected:
	NeighborhoodGroupComponentCompare srcc;
};


bool extendRange( MatchRecord* M_i, MatchRecord* M_m, const vector< size_t >& component_map )
{
	bool changed = false;
		// set range to cover M_m
	for( size_t x = 0; x < M_i->Multiplicity(); ++x )
	{
		size_t z = component_map[x];
		if( M_i->LeftEnd(x) == NO_MATCH || M_m->LeftEnd(z) == NO_MATCH )
			genome::breakHere();
		int64 lend_diff = M_i->LeftEnd(x) - M_m->LeftEnd(z);
		if( lend_diff > 0 )
		{
			M_i->SetLeftEnd(x, M_i->LeftEnd(x) - lend_diff);
			M_i->SetLength(M_i->Length(x)+lend_diff, x);
			changed = true;
		}
		int64 rend_diff = M_m->RightEnd(z) - M_i->RightEnd(x);
		if( rend_diff > 0 )
		{
			M_i->SetLength( M_i->Length(x)+rend_diff, x );
			changed = true;
		}
	}
	return changed;
}

void remapComponents(const vector< size_t >& srcmap, size_t mid_multiplicity, const vector< size_t >& destmap, vector< size_t >& newmap )
{
	vector< size_t > super_map( mid_multiplicity, (std::numeric_limits<size_t>::max)() );
	for( size_t mapI = 0; mapI < destmap.size(); ++mapI )
		super_map[destmap[mapI]] = mapI;
	for( size_t mapI = 0; mapI < srcmap.size(); ++mapI )
		newmap[mapI] = super_map[srcmap[mapI]];
}

void classifyMatch( AbstractMatch* M_i, AbstractMatch* M_j, vector< size_t >& ji_component_map, bool& subsumed, bool& partial )
{
	subsumed = true;
	partial = false;
	for( size_t i = 0; i < ji_component_map.size(); ++i )
	{
		size_t x = ji_component_map[i];
		size_t y = i;
		int64 lend_diff = M_i->LeftEnd(x) - M_j->LeftEnd(y);
		int64 rend_diff = M_j->RightEnd(y) - M_i->RightEnd(x);
		if( lend_diff > 0 || rend_diff > 0 )
			subsumed = false;
		if( lend_diff <= 0 && rend_diff <= 0 )
			partial = true;
	}
}

void classifySubset( MatchRecord* M_i, NeighborhoodGroup& sr, bool& subsumed, bool& partial )
{
	classifyMatch( M_i, sr.get<0>(), sr.get<1>(), subsumed, partial );
}

void checkLink( MatchRecord*& mr )
{
	while( mr->subsuming_match != NULL )
		mr = mr->subsuming_match;
}


void checkLink( MatchLink& mlink )
{
	while( mlink.subset->subsuming_match != NULL )
	{
		vector< size_t > new_map( mlink.sub_to_super_map.size() );
		for( size_t i = 0; i < mlink.sub_to_super_map.size(); ++i )
			new_map[i] = mlink.sub_to_super_map[ mlink.subset->subsumption_component_map[i] ];
		swap( new_map, mlink.sub_to_super_map );
		mlink.subset = mlink.subset->subsuming_match;
	}
}

void checkLinkAndComponent( MatchRecord*& mr, size_t& component )
{
	while( mr->subsuming_match != NULL )
	{
		component = mr->subsumption_component_map[component];
		mr = mr->subsuming_match;
	}
}

/** returns one of the superset links f
rom a match.  direction is 1 for left, -1 for right */
MatchLink& getSuperset( MatchRecord* mr, int direction )
{
	if( direction == 1 )
		return mr->left_superset;
	return mr->right_superset;
}

/** returns the subset links for a given direction.  direction is 1 for left, -1 for right */
vector<MatchLink>& getSubsets( MatchRecord* mr, int direction )
{
	if( direction == 1 )
		return mr->left_subset_links;
	return mr->right_subset_links;
}

/** returns the extra subsets for a given direction.  direction is 1 for left, -1 for right */
vector<MatchLink>& getExtraSubsets( MatchRecord* mr, int direction )
{
	if( direction == 1 )
		return mr->extra_left_subsets;
	return mr->extra_right_subsets;
}

void unlinkSuperset( MatchRecord* mr, int direction )
{
	MatchLink& superlink = getSuperset( mr, direction );
	MatchRecord* super = superlink.superset;
	if( super != NULL )
	{
		int parity = mr->Orientation(0) == super->Orientation(superlink.sub_to_super_map[0]) ? 1 : -1;
		vector< MatchLink >& subs = getSubsets( super, -direction*parity );
		for( size_t subI = 0; subI < subs.size(); ++subI )
		{
			if( subs[subI].subset == mr )
			{
				subs.erase( subs.begin() + subI, subs.begin() + subI + 1 );
				subI--;
			}
		}
		superlink.clear();
	}
}

void unlinkSupersets( MatchRecord* mr )
{
	unlinkSuperset( mr, 1 );
	unlinkSuperset( mr, -1 );
}

template< class MatchRecordPtrType >
void validate( vector< MatchRecordPtrType >& records )
{
	// make sure all matches have non-zero components
	for( size_t recI = 0; recI < records.size(); ++recI )
	{
		size_t seqI = 0;
		for( ; seqI < records[recI]->SeqCount(); ++seqI )
			if( records[recI]->LeftEnd(seqI) == NO_MATCH )
				break;
		if( seqI < records[recI]->SeqCount() )
		{
			cerr << "missing component\n";
			genome::breakHere();
		}
	}

	// make sure all links are consistent
	for( size_t recI = 0; recI < records.size(); ++recI )
	{
		MatchRecord* mr = records[recI];
		for( int direction = 1; direction > -2; direction -= 2 )
		{
			for( size_t subI = 0; subI < getSubsets(mr, direction).size(); subI++ )
			{
				// follow any stale links
				MatchRecord* sub = getSubsets(mr, direction)[subI].subset;
				size_t sub_mult = sub->Multiplicity();
				while( sub->subsuming_match != NULL )
					sub = sub->subsuming_match;
				size_t parity_seq = getSubsets(mr, direction)[subI].sub_to_super_map[0];
				int parity = mr->Orientation(parity_seq) == sub->Orientation(0) ? 1 : -1;
				// make sure that each of the subsets in these points back to this superset in its own link
				if( getSuperset(sub, -direction*parity).superset != mr )
				{
					cerr << "ohno\n";
					genome::breakHere();
				}
				if( sub_mult != sub->Multiplicity() )
				{
					cerr << "unequal mult\n";
					genome::breakHere();
				}
				if( getSubsets(mr,direction)[subI].super_component_list.count() != getSubsets(mr,direction)[subI].sub_to_super_map.size())
				{
					cerr << "broke\n";
					genome::breakHere();
				}
			}

			// make sure the supersets have this subset
			if( getSuperset(mr,direction).superset != NULL )
			{
				MatchRecord* sup = getSuperset(mr,direction).superset;
				int parity = mr->Orientation(0) == sup->Orientation(getSuperset(mr,direction).sub_to_super_map[0]) ? 1 : -1;
				size_t subI = 0;
				for( ; subI < getSubsets(sup,-direction*parity).size(); subI++ )
				{
					if( getSubsets(sup,-direction*parity)[subI].subset == mr )
						break;
				}
				if( subI == getSubsets(sup,-direction*parity).size() )
				{
					cerr << "oh crap!\n";
					genome::breakHere();
				}
				if( getSuperset(mr,direction).super_component_list.count() != getSuperset(mr,direction).sub_to_super_map.size())
				{
					cerr << "broke 3\n";
					genome::breakHere();
				}
			}
		}
	}
}

void createNeighborhoodGroupList( vector< NeighborhoodGroup >& group_list, vector< vector< size_t > >& group_members, vector< NeighborhoodListEntry >& neighborhood_list )
{
	group_list.resize( group_members.size() );
	for( size_t gI = 0; gI < group_members.size(); gI++ )
	{
		// is this subset completely contained--is it subsumed?
		MatchRecord* M_j = neighborhood_list[group_members[gI][0]].match;

		vector< size_t > component_map(M_j->Multiplicity(), (std::numeric_limits<size_t>::max)());
		vector< size_t > distances(M_j->Multiplicity(), (std::numeric_limits<size_t>::max)());
		for( vector< size_t >::iterator rec_iter = group_members[gI].begin(); rec_iter != group_members[gI].end(); ++rec_iter )
		{
			component_map[neighborhood_list[*rec_iter].Mj_component] = neighborhood_list[*rec_iter].Mi_component;
			distances[neighborhood_list[*rec_iter].Mj_component] = neighborhood_list[*rec_iter].distance;
		}
		group_list[gI].get<0>() = M_j;
		swap( group_list[gI].get<1>(), component_map );
		swap( group_list[gI].get<2>(), distances );
	}

	static NeighborhoodGroupCompare src;
	std::sort( group_list.begin(), group_list.end(), src );
}

void inheritSuperset( MatchRecord* M_i, MatchRecord* M_j, int direction, int parity )
{
	// remap superset components
	vector< size_t > comp_map( M_i->Multiplicity() );
	for( size_t ci = 0; ci < comp_map.size(); ci++ )
		comp_map[ci] = getSuperset( M_j, direction*parity ).sub_to_super_map[ M_j->subsumption_component_map[ci] ];
	// rebuild the superset component list
	boost::dynamic_bitset<> comp_list(getSuperset( M_j, direction*parity ).superset->Multiplicity(), false);
	for( size_t compI = 0; compI < comp_map.size(); ++compI )
		comp_list.set(comp_map[compI]);
	MatchLink& slink = getSuperset(M_i, direction);
	slink = MatchLink( getSuperset( M_j, direction*parity ).superset, M_i, comp_list, comp_map );
	unlinkSuperset(M_j,direction*parity);
	int slink_parity = M_i->Orientation(0) == slink.superset->Orientation(slink.sub_to_super_map[0]) ? 1 : -1;
	getSubsets(slink.superset,-direction*slink_parity).push_back(slink);

}

vector< NeighborhoodGroup >& selectList( vector< NeighborhoodGroup >& left_list, vector< NeighborhoodGroup >& right_list, int direction )
{
	return direction == 1 ? left_list : right_list;
}

class ToUPPER
{
public:
	char operator()( char a ){ return toupper(a); }
};



/**
 * Performs a superset link extension on M_i
 */
void supersetLinkExtension( GappedMatchRecord*& M_i, int direction, int& last_linked, 
						   vector< NeighborhoodGroup >& left_deferred_subsets, 
						   vector< NeighborhoodGroup >& right_deferred_subsets )
{
	// update the left end and look for another superset to chain with
	// then extend all the way to that match
	MatchRecord* M_j = getSuperset(M_i, direction).superset;
	MatchLink ij_link = getSuperset(M_i, direction);	// make a copy for safekeeping
	int ij_parity = M_i->Orientation(0) == M_j->Orientation(ij_link.sub_to_super_map[0]) ? 1 : -1;

	//
	// Link extension part 1: 
	// extend M_i to include M_j, add M_j to the chained matches
	bool changed = extendRange( M_i, M_j, ij_link.sub_to_super_map );
	M_i->chained_matches.push_back(M_j);
	M_i->chained_component_maps.push_back(ij_link.sub_to_super_map);



	// Link extension part 2:
	// figure out whether any subsets between M_j and M_i got subsumed
	for( size_t subtypeI = 0; subtypeI < 2; subtypeI++ )
	{
		vector< MatchLink >* mjsubs;
		if( subtypeI == 0 )
			mjsubs = &getSubsets(M_j, -direction*ij_parity);
		else
			mjsubs = &getExtraSubsets(M_j, -direction*ij_parity);
		vector< MatchLink >& mj_otherside_subsets = *mjsubs;

		for( size_t leftI = 0; leftI < mj_otherside_subsets.size(); ++leftI )
		{
			if( subtypeI == 0 )
				checkLink( mj_otherside_subsets[leftI] );
			MatchLink& jk_link = mj_otherside_subsets[leftI];
			boost::dynamic_bitset<> intersect = ij_link.super_component_list & jk_link.super_component_list;
			MatchRecord* M_k = jk_link.subset;
			if( M_k == M_i )
				continue;	// been there, chained that.
			size_t inter_size = intersect.count();
			if( inter_size < 2 )
				continue;	// no match
			if( inter_size >= M_i->Multiplicity() || M_k->Multiplicity() != inter_size )
				continue;

			// has this guy already been subsumed?  if so then just skip him
			if( M_k->subsuming_match != NULL )
			{
				if( subtypeI != 1 )
					breakHere(); // this should only happen with extra subsets
				mj_otherside_subsets.erase(mj_otherside_subsets.begin()+leftI, mj_otherside_subsets.begin()+leftI+1 );
				leftI--;
				continue;
			}

			// M_k is a subset relative to M_i
			int jk_parity = M_k->Orientation(0) == M_j->Orientation(jk_link.sub_to_super_map[0]) ? 1 : -1;
			int ik_parity = ij_parity * jk_parity;

			vector< size_t > component_map( M_k->Multiplicity() );
			remapComponents(jk_link.sub_to_super_map, M_j->Multiplicity(), ij_link.sub_to_super_map, component_map );

			NeighborhoodGroup sr = boost::make_tuple( M_k, component_map, vector<size_t>( M_k->Multiplicity(), 0 ) );
			// defer it until we're done extending
			selectList( left_deferred_subsets, right_deferred_subsets, -direction ).push_back( sr );
		}
	}

	//
	// Link extension part 3:
	// classify outgoing links that share components with M_i
	unlinkSuperset(M_i,direction);
	vector< size_t > supersets;
	vector< size_t > chainable;
	vector< size_t > subsets;
	vector< size_t > novel_subsets;
	vector< MatchLink >& mj_subsets = getSubsets(M_j, direction*ij_parity);
	for( size_t leftI = 0; leftI < mj_subsets.size(); ++leftI )
	{
		checkLink( mj_subsets[leftI] );
		boost::dynamic_bitset<> intersect = ij_link.super_component_list & mj_subsets[leftI].super_component_list;
		MatchRecord* M_k = mj_subsets[leftI].subset;
		if( M_k == M_i )
			continue;	// been there, chained that.
		size_t inter_size = intersect.count();
		if( inter_size < 2 )
			continue;	// no match
			// M_k is a superset relative to M_i
		if( inter_size == M_i->Multiplicity() && M_k->Multiplicity() > inter_size )
			supersets.push_back(leftI);
		else if( inter_size == M_i->Multiplicity() && M_k->Multiplicity() == inter_size )
			chainable.push_back(leftI);
		else if( inter_size < M_i->Multiplicity() && M_k->Multiplicity() == inter_size )
			subsets.push_back(leftI);
		else
			novel_subsets.push_back(leftI);
	}

	if( supersets.size() > 0 )
	{
//#4018
		cerr << "something is wrong, we should never have supersets during link extension!\n";
		breakHere();
	}

	for( size_t cI = 0; cI < chainable.size(); ++cI )
	{
		if( chainable.size() > 1 )
		{
			cerr << "bad news bruthah\n";
			genome::breakHere();
		}
		// chain with this guy
		MatchLink& jk_link = mj_subsets[chainable[cI]];
		MatchRecord* M_k = jk_link.subset;
		if( M_k->extended )
		{
			cerr << "extensor crap\n";
			breakHere();
		}
		if( M_k == M_i )
		{
			cerr << "crap\n";
			breakHere();
		}

		// update boundary coordinates
		vector< size_t > component_map( M_i->Multiplicity() );
		remapComponents(ij_link.sub_to_super_map, M_j->Multiplicity(), jk_link.sub_to_super_map, component_map );
		bool changed = extendRange( M_i, M_k, component_map );
		if( changed )
			last_linked = 2;

		// unlink from superset
		int jk_parity = M_k->Orientation(0) == M_j->Orientation(jk_link.sub_to_super_map[0]) ? 1 : -1;
		unlinkSuperset(M_k,-direction*ij_parity*jk_parity);
		// set subsuming match ptrs
		M_k->subsuming_match = M_i;
		M_k->subsumption_component_map = component_map;
		M_i->chained_matches.push_back( M_k );
		M_i->chained_component_maps.push_back( component_map );

		// compensate for the deletion in subsets
		for( size_t subI = 0; subI < chainable.size(); subI++ )
			if( chainable[subI] > chainable[cI] )
				chainable[subI]--;
		for( size_t subI = 0; subI < subsets.size(); subI++ )
			if( subsets[subI] > chainable[cI] )
				subsets[subI]--;

		// inherit M_k's outward superset and stop chaining here
		if( getSuperset( M_k, direction*ij_parity*jk_parity ).superset != NULL )
		{
			inheritSuperset( M_i, M_k, direction, ij_parity*jk_parity );
			last_linked = 2;
			break;
		}
	}

	// process subsets
	for( size_t sI = 0; sI < subsets.size(); ++sI )
	{
		// change M_k to point at M_i
		MatchLink& jk_link = mj_subsets[subsets[sI]];
		MatchRecord* M_k = jk_link.subset;
		int jk_parity = M_k->Orientation(0) == M_j->Orientation(jk_link.sub_to_super_map[0]) ? 1 : -1;
		int ik_parity = ij_parity * jk_parity;

		vector< size_t > component_map( M_k->Multiplicity() );
		remapComponents(jk_link.sub_to_super_map, M_j->Multiplicity(), ij_link.sub_to_super_map, component_map );
		// rebuild the superset component list
		boost::dynamic_bitset<> comp_list(M_i->Multiplicity(), false);
		for( size_t compI = 0; compI < component_map.size(); ++compI )
			if(component_map[compI] != (std::numeric_limits<size_t>::max)())
				comp_list.set(component_map[compI]);
		unlinkSuperset(M_k,-1*direction*ik_parity);

		// add to the deferred subsets list
		NeighborhoodGroup sr = boost::make_tuple( M_k, component_map, vector<size_t>( M_k->Multiplicity(), 0 ) );
		vector< NeighborhoodGroup >& subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
		subset_list.push_back( sr );

		// compensate for the deletion in subsets
		for( size_t subI = 0; subI < subsets.size(); subI++ )
			if( subsets[subI] > subsets[sI] )
				subsets[subI]--;
	}
}

/**
 * Performs a neighborhood list extension
 */
void neighborhoodListLookup( GappedMatchRecord*& M_i, 
						   MatchPositionLookupTable& match_pos_lookup_table,
						   vector< NeighborhoodGroup >& superset_list, 
						   vector< NeighborhoodGroup >& chainable_list,
						   vector< NeighborhoodGroup >& subset_list, 
						   vector< NeighborhoodGroup >& novel_subset_list,
						   int direction,
						   uint seed_size,
						   uint w,
						   bitset_t& left_lookups,
						   bitset_t& right_lookups
						   )
{
	//
	// construct a neighborhood list and process the neighborhood groups
	//
	vector< NeighborhoodListEntry > neighborhood_list;
	for( size_t x = 0; x < M_i->Multiplicity(); ++x )
	{
		int o_x = M_i->Orientation(x) == AbstractMatch::forward ? 1 : -1;
		int parity = o_x * direction;
		int64 match_end = parity == 1 ? M_i->LeftEnd(x) : M_i->RightEnd(x) - seed_size + 1;

		if( match_end > 0 )
			if( (direction == 1 && left_lookups.test(match_end)) ||
				(direction == -1 && right_lookups.test(match_end)) )
			{
				cerr << "looking twice in the same place\n";
//							genome::breakHere();
			}else{
				if( direction == 1 )
					left_lookups.set(match_end);
				if( direction == -1 )
					right_lookups.set(match_end);
			}

		for( int d = 1; d <= w; ++d )
		{
			if( match_end <= parity * d )
				continue;	// we're too close to the beginning
			size_t mplt_index = match_end - parity * d;
			if( mplt_index >= match_pos_lookup_table.size() )
				continue;	// we're too close to the end!

			MatchRecord* M_j = match_pos_lookup_table[ mplt_index ].first;
			size_t y = match_pos_lookup_table[ mplt_index ].second;
			if( M_j == NULL )
				continue;	// no match at this position

			NeighborhoodListEntry nle;
			nle.match = M_j;
			nle.Mi_component = x;
			nle.Mj_component = y;
			// update the link if this one was subsumed
			checkLinkAndComponent( M_j, y );
			int o_y = ((AbstractMatch*)M_j)->Orientation(y) == AbstractMatch::forward ? 1 : -1;
			nle.relative_orientation = o_x * o_y == 1 ? true : false;
			nle.distance = d;
			neighborhood_list.push_back( nle );
			
			if( M_j == M_i )
			{
				M_i->tandem = true;
				break;	// not so fast there cowboy!  can't chain beyond ourself!
			}
		}
	}

	//
	// now classify each group of the neighborhood list and act appropriately
	// group types are superset, chainable, subset, novel subset
	//
	NeighborhoodListComparator nlc;
	std::sort( neighborhood_list.begin(), neighborhood_list.end(), nlc );

	std::vector< std::vector< size_t > > superset_groups;
	std::vector< std::vector< size_t > > chainable_groups;
	std::vector< std::vector< size_t > > subset_groups;
	std::vector< std::vector< size_t > > novel_subset_groups;

	size_t group_end = 0;
	for( size_t prev = 0; prev < neighborhood_list.size(); prev = group_end )
	{
		group_end = prev + 1;
		while( group_end < neighborhood_list.size() && 
			neighborhood_list[prev].match == neighborhood_list[group_end].match && 
			neighborhood_list[prev].relative_orientation == neighborhood_list[group_end].relative_orientation )
		{
			++group_end;
		}
		// the group is everything in the range of prev to end-1
		if( prev + 1 == group_end )
			continue;	// can't do anything with groups of size 1 -- there's no match

		// do something about ties here...???
		// this code selects the *furthest* away match (e.g. that with the largest d)
		// because that's what got sorted in last in the comparator
		// it eliminates both duplicate M_i and duplicate M_j components...
		vector< pair< size_t, size_t > > j_comp_sort_list;
		for( size_t i = prev + 1; i < group_end; ++i )
		{
			if( neighborhood_list[i-1].Mi_component == neighborhood_list[i].Mi_component )
				continue;
			j_comp_sort_list.push_back(make_pair(neighborhood_list[i-1].Mj_component, i-1));
		}
		j_comp_sort_list.push_back(make_pair(neighborhood_list[group_end-1].Mj_component, group_end-1));
		std::sort(j_comp_sort_list.begin(), j_comp_sort_list.end());
		vector<size_t> group_entries;
		for( size_t i = 1; i < j_comp_sort_list.size(); ++i )
		{
			if( j_comp_sort_list[i-1].first == j_comp_sort_list[i].first )
				continue;
			group_entries.push_back(j_comp_sort_list[i-1].second);
		}
		group_entries.push_back(j_comp_sort_list.back().second);

		// update the links in case something is subsumed
		for( size_t gI = 0; gI < group_entries.size(); ++gI )
			checkLinkAndComponent( neighborhood_list[group_entries[gI]].match, neighborhood_list[group_entries[gI]].Mj_component );

		// finally, classify the match as one of superset, subset, 
		// chainable, novel subset
		MatchRecord* M_j = neighborhood_list[prev].match;

		if( group_entries.size() == M_i->Multiplicity() && 
			M_j->Multiplicity() > M_i->Multiplicity() )
		{
			// superset
			superset_groups.push_back( group_entries );
		}else
		if( group_entries.size() == M_i->Multiplicity() && 
			M_j->Multiplicity() == M_i->Multiplicity() )
		{
			// chainable
			chainable_groups.push_back( group_entries );
		}else
		if( group_entries.size() < M_i->Multiplicity() && 
			group_entries.size() == M_j->Multiplicity() )
		{
			// subset
			subset_groups.push_back( group_entries );
		}else
		{
			// novel subset
			novel_subset_groups.push_back( group_entries );
		}

	}	// end loop that splits the neighborhood into groups

	createNeighborhoodGroupList( superset_list, superset_groups, neighborhood_list );
	createNeighborhoodGroupList( chainable_list, chainable_groups, neighborhood_list );
	createNeighborhoodGroupList( subset_list, subset_groups, neighborhood_list );
	createNeighborhoodGroupList( novel_subset_list, novel_subset_groups, neighborhood_list );
}

/**
 * Chains matches onto M_i or subsumes them as appropriate
 */
void processChainableMatches( GappedMatchRecord*& M_i, vector< NeighborhoodGroup >& chainable_list,
				  int direction, int& last_linked )
{
	// link the closest possible chainable first.
	for( size_t gI = 0; gI < chainable_list.size(); gI++ )
	{
		MatchRecord* M_j = chainable_list[gI].get<0>();

		vector< size_t >& component_map = chainable_list[gI].get<1>();

		if( M_j == M_i )
		{
			// this is an inverted overlapping repeat, skip it.
			continue;
		}
		if( M_j->extended )
		{
			// oh no!  M_i should have been swallowed up already!
//						cerr << "extensor crap 2\n";
//						breakHere();
		}

		bool subsumed;
		bool partial;
		classifySubset( M_i, chainable_list[gI], subsumed, partial );

		M_j->subsuming_match = M_i;
		M_j->subsumption_component_map = component_map;
		vector< size_t >& yx_map = chainable_list[gI].get<1>();
		vector< size_t > xy_map(yx_map.size());
		for( size_t i = 0; i < yx_map.size(); ++i )
			xy_map[ yx_map[i] ] = i;
//		for( vector< size_t >::iterator rec_iter = chainable_groups[gI].begin(); rec_iter !=  chainable_groups[gI].end(); ++rec_iter )
//			xy_map[ neighborhood_list[*rec_iter].Mi_component ] = neighborhood_list[*rec_iter].Mj_component;

		// if M_j isn't extending the boundaries of every component of M_i then
		// it may be inconsistent with already chained matches.  just subsume it without
		// chaining in that case.
		if( !subsumed && !partial )
		{
			M_i->chained_matches.push_back( M_j );
			M_i->chained_component_maps.push_back( component_map );
			// update the left-end and right-end coords
			bool changed = extendRange(M_i, M_j, xy_map);
			if( changed )
				last_linked = 2;
		}

		int parity = M_i->Orientation(0) == M_j->Orientation(xy_map[0]) ? 1 : -1;
		if( getSuperset( M_j, -direction*parity ).superset != NULL )
			unlinkSuperset(M_j,-direction*parity);	// won't be needing this anymore...

		// if M_j has a superset then inherit it and stop chaining here
		if( getSuperset( M_j, direction*parity ).superset != NULL )
		{
			inheritSuperset( M_i, M_j, direction, parity );
			last_linked = 2;	// we may do a link extension!
			break;
		}
	}
}
//void ExtendMatch(vector< GappedMatchRecord* >& final, int fI, vector< gnSequence* >& seq_table, PairwiseScoringScheme& pss)
int ExtendMatch(GappedMatchRecord*& mte, vector< gnSequence* >& seq_table, PairwiseScoringScheme& pss, unsigned w, int direction = 0)
{
	ccount +=1;
	//todo: remove alignment parameter
	//      dont pass vector? 
	static bool debug_extension = false;
	//punt on this for now..
	bool novel_hss_regions_support = false;
	bool danger_zone_active = true;
	vector< string > alignment;
	mems::GetAlignment(*mte,seq_table, alignment);	// expects one seq_table entry per matching component
	int multi = mte->Multiplicity();
	
//	int extend_length_1 = ceil(-0.78*multi+150); 
	int extend_length = max(w, sqrt(max(pow(150,2.0)-pow(2*multi,2.0),0)));
	// aced 06/25/07: should re-evaluate the situation now that the Homology HMM has been implemented
	//
//  tjt: looks like with hmm w is ok!
	//int extend_length = min(int(2*w), 130);	
	
	vector<int> left_extend_vector;
	vector<int> right_extend_vector;
	int left_extend_length = extend_length;	
	int right_extend_length = extend_length;
	if ( mte->tandem )
	{		
		cerr << "Sorry, no extension for tandem repeats.." << endl << endl;	
		return FIXME;
	}

	//careful, if mte->LeftEnd(j) < extend_length, ToString() will be disappointed...
	for( gnSeqI j = 0; j < alignment.size(); j++)
	{
		//now put check for curpos+extend_length<startpos of next match component..
		 
		if( mte->Orientation(j) == AbstractMatch::reverse )
		{
			//if leftend <= 0 set right extension to 0
			if( mte->LeftEnd(j) <= 0 || mte->LeftEnd(j) > 4000000000u )
				right_extend_vector.push_back(0);
			//if extend_length goes to far, set to maximum possible
			else if ( mte->LeftEnd(j) <= extend_length )
				right_extend_vector.push_back(mte->LeftEnd(j)-1);

			//if we run into another match, don't extend into it
			else if ( j > 0 && mte->LeftEnd(j) - extend_length <= mte->RightEnd(j-1) )
				right_extend_vector.push_back(mte->LeftEnd(j)-mte->RightEnd(j-1)-1);
			
			//else everything ok to set to preset extend_length
			else
				right_extend_vector.push_back(extend_length-1);
			
			if(mte->RightEnd(j) <= 0 || mte->RightEnd(j) > 4000000000u)
				left_extend_vector.push_back(0);
			else if ( mte->RightEnd(j) + extend_length > seq_table[0]->length() )
				left_extend_vector.push_back(seq_table[0]->length()-mte->RightEnd(j));
			else if ( j+1 < alignment.size() && mte->RightEnd(j) + extend_length >= mte->LeftEnd(j+1) )
				left_extend_vector.push_back(mte->LeftEnd(j+1)-mte->RightEnd(j)-1);
			else
				left_extend_vector.push_back(extend_length-1);
	
		}
		else
		{
			if( mte->LeftEnd(j) <= 0 || mte->LeftEnd(j) > 4000000000u )
				left_extend_vector.push_back(0);
			else if ( mte->LeftEnd(j) <= extend_length )
				left_extend_vector.push_back(mte->LeftEnd(j)-1);
			else if ( j > 0 && mte->LeftEnd(j) - extend_length <= mte->RightEnd(j-1) )
				left_extend_vector.push_back(mte->LeftEnd(j)-mte->RightEnd(j-1)-1);
			else
				left_extend_vector.push_back(extend_length-1);

			if(mte->RightEnd(j) <= 0 || mte->RightEnd(j) > 4000000000u)
				right_extend_vector.push_back(0);
			else if ( mte->RightEnd(j) + extend_length > seq_table[0]->length() )
				right_extend_vector.push_back(seq_table[0]->length()-mte->RightEnd(j));
			else if ( j+1 < alignment.size() && mte->RightEnd(j) + extend_length >= mte->LeftEnd(j+1) )
				right_extend_vector.push_back(mte->LeftEnd(j+1)-mte->RightEnd(j)-1);
			else
				right_extend_vector.push_back(extend_length-1);	
		}
	}
	left_extend_length = *(std::min_element( left_extend_vector.begin(), left_extend_vector.end() ));
	right_extend_length = *(std::min_element( right_extend_vector.begin(), right_extend_vector.end() ));

	const gnFilter* rc_filter = gnFilter::DNAComplementFilter();
	std::vector<std::string> leftExtension(multi);
	GappedAlignment leftside(multi,left_extend_length);
	std::vector<std::string> rightExtension(multi);
	GappedAlignment rightside(multi,right_extend_length);
	vector< string > leftExtension_aln;
	vector< string > rightExtension_aln;

	if ( left_extend_length > 10 && direction == 1  )
	{
		// extract sequence data
		for( gnSeqI j = 0; j < alignment.size(); j++)
		{			
			if( mte->Orientation(j) == AbstractMatch::reverse )
			{			
				seq_table[0]->ToString( leftExtension[j], left_extend_length, mte->RightEnd(j)+1 );
				leftside.SetLeftEnd(j,mte->RightEnd(j)+1);
				rc_filter->ReverseFilter(leftExtension[j]);
			}else{
				seq_table[0]->ToString( leftExtension[j], left_extend_length, mte->LeftEnd(j) - left_extend_length );
				leftside.SetLeftEnd(j,mte->LeftEnd(j) - left_extend_length);
			}
			leftside.SetOrientation(j,mte->Orientation(j));
			leftside.SetLength(left_extend_length,j);
		}
		bool align_success = false;
		align_success = mems::MuscleInterface::getMuscleInterface().CallMuscleFast( leftExtension_aln, leftExtension );
		if ( align_success ){		
			leftside.SetAlignment(leftExtension_aln);
			leftside.SetAlignmentLength(leftExtension_aln.at(0).size());
		}else{
			cerr << "Muscle failed" << endl;
			return FAILED;
		}
	}
	else if ( right_extend_length > 10 && direction == -1 )
	{
		for( gnSeqI j = 0; j < alignment.size(); j++)
		{
			if( mte->Orientation(j) == AbstractMatch::reverse )
			{			
				rightside.SetLeftEnd(j,mte->LeftEnd(j) - right_extend_length-1);
				seq_table[0]->ToString( rightExtension[j], right_extend_length, mte->LeftEnd(j) - right_extend_length-1);
				rc_filter->ReverseFilter(rightExtension[j]);
			}else{
				rightside.SetLeftEnd(j,mte->RightEnd(j)+1 );
				seq_table[0]->ToString( rightExtension[j], right_extend_length, mte->RightEnd(j)+1 );
			}
			rightside.SetOrientation(j,mte->Orientation(j));
			rightside.SetLength(right_extend_length,j);
		}
		bool align_success = false;		
		align_success = mems::MuscleInterface::getMuscleInterface().CallMuscleFast( rightExtension_aln, rightExtension );
		if ( align_success ){
			rightside.SetAlignment(rightExtension_aln);
			rightside.SetAlignmentLength(rightExtension_aln.at(0).size());
		}else{
			cerr << "Muscle failed" << endl;
			return FAILED;
		}
	}else{
		//what are you even doing here?!?
		cerr << "Extension failed" << endl;
		return FAILED;
	}

	vector< AbstractMatch* > mlist;
	if( direction == 1 )
		mlist.push_back(leftside.Copy());
	//tjt: don't use original match, only regions to the left/right
	//     since for now we won't modify mte, even if the homology detection method suggests otherwise
	//mlist.push_back(mte->Copy());
	if( direction == -1 )
		mlist.push_back(rightside.Copy());
 
	//createInterval
	Interval iv;
	iv.SetMatches(mlist);
	CompactGappedAlignment<>* cga = new CompactGappedAlignment<>( iv );
	vector< CompactGappedAlignment<>* > cga_list;
	CompactGappedAlignment<>* result;
	//detectAndApplyBackbone
	backbone_list_t bb_list;
	detectAndApplyBackbone( cga, seq_table,result,bb_list,pss, DEFAULT_ISLAND_SCORE_THRESHOLD, direction != 1, direction == 1 );
	cga->Free();

	bool boundaries_improved = false;
	if( bb_list.size() == 0 ||
		bb_list.at(0).size() == 0)
	{
		//no backbone segment found
		cerr << "Crikey!! no backbone found during extension..." << endl;
		return FAILED;
	}

	if(debug_extension)
	{
	// aced: debug printing to get bb segments right
		for( size_t bbI = 0; bb_list.size() > 0 && bbI < bb_list[0].size(); bbI++ )
			cerr << "bbI: " << bbI << '\t' << *(bb_list[0][bbI]) << endl;
	}
	
	int backbone_i = -1;
	if (bb_list.at(0).size() == 1)
		backbone_i = 0;
	int invalid_matches = 0;

	// tjt: only want to check the first/last backbone (for now)
	//      process novel homolgous regions in the near future
	CompactGappedAlignment<> tmp_cga;
	cga_list.push_back( tmp_cga.Copy() );
	if ( direction > 0 )
		result->copyRange(*(cga_list.back()),bb_list.at(0).back()->LeftEnd(0),bb_list.at(0).back()->AlignmentLength()-1);
	else
		result->copyRange(*(cga_list.back()),bb_list.at(0).front()->LeftEnd(0),bb_list.at(0).front()->AlignmentLength()-1);
	if( cga_list.back()->Multiplicity() == 0 )
	{
		// this one must have been covering an invalid region (gaps aligned to gaps)
		cga_list.back()->Free();
		cga_list.erase( cga_list.end()-1 );
		invalid_matches++;
		return FAILED;
	}
	if( (cga_list.back()->Multiplicity() == mte->Multiplicity() ) && (cga_list.back()->Length() > 1 ))
	{
		// successful extension!!
		//boundaries were improved, current match is extended original match
		cerr << "Extension worked!! Improved boundaries! Multiplicity: " << mte->Multiplicity() << endl;
		cerr << "Old boundaries: " << mte->LeftEnd(0) << " " << mte->RightEnd(0) << endl;

		vector< AbstractMatch* > matches( 1, cga_list.back());
		GappedMatchRecord* M_new = mte->Copy();
		//tjt: clobber mte's GappedMatchRecord data and set boundaries
		M_new->SetMatches(matches);
		//build the component map for processChainableMatches
		vector< size_t > component_map( mte->Multiplicity() );
		for( size_t i = 0; i < component_map.size(); ++i )
			component_map[i] = i;
		vector< NeighborhoodGroup > chainable_list;
		NeighborhoodGroup sr = boost::make_tuple( M_new, component_map, vector<size_t>( mte->Multiplicity(), 0 ) );
		chainable_list.push_back(sr);
		//tjt: does this need to be updated somewhere after call to processChainableMatches?
		int last_linked = 0;
		//use this function to correctly chain the left/right gapped extension to mte(M_i)
		processChainableMatches( mte, chainable_list, direction, last_linked );
		cerr << "New boundaries: " << mte->LeftEnd(0) << " " << mte->RightEnd(0) << endl << endl;

		if(debug_extension)
		{
			//tjt: write alignment to file
			//need to finalize before calling GetAlignment!
			mte->finalize(seq_table);
			vector<string> new_alignment;
			mems::GetAlignment(*mte, seq_table, new_alignment);	// expects one seq_table entry per matching component
			gnSequence myseq;
			for( uint j = 0; j < new_alignment.size(); j++)
				myseq += new_alignment.at(j);
			gnFASSource::Write(myseq, "coolio2.txt" );
		}
		
		return OK;

	}
	else
	{
		cerr << "Extension failed.. " << endl << endl;
		return FAILED;
	}
	

	
	return FAILED;
	
}//tjt: match should be extended!

class ProcrastinationQueue
{
public:
	template< typename MrPtrType >
	ProcrastinationQueue( vector< MrPtrType >& match_record_list )
	{
		q.resize( match_record_list.size() );
		std::copy(match_record_list.begin(), match_record_list.end(), q.begin() );
		std::make_heap( q.begin(), q.end(), mhc );
		q_end = q.size();
		q_size = q.size();
	}

	/** pops an element from the queue, maintaining heap order */
	MatchRecord* pop()
	{
		std::pop_heap( q.begin(), q.begin()+q_end, mhc );
		q_end--;
		return *(q.begin() + q_end);
	}

	/** Adds an element to the queue and restores heap order */
	void push( MatchRecord* M_n )
	{
		if( q_end < q.size() )
			q[q_end] = M_n;
		else
		{
			q.push_back(M_n);
		}
		q_size++;
		q_end++;
		std::push_heap(q.begin(), q.begin()+q_end, mhc);
	}
	/** gets the total number of elements that have been placed in the queue */
	size_t size() const{ return q_size; }
	/** returns the current number of elements in the queue */
	size_t end() const{ return q_end; }


	/** defines a multiplicity heap ordering */
	class MultiplicityHeapCompare
	{
	public:
		bool operator()( const MatchRecord* a, const MatchRecord* b )
		{
			return a->Multiplicity() < b->Multiplicity();
		}
	};

private:
	const MultiplicityHeapCompare mhc;
	std::vector< MatchRecord* > q;
	size_t q_end;
	size_t q_size;
};

/**
 * Creates novel subset matches where appropriate and adds them to the procrastination queue
 */
void processNovelSubsetMatches( GappedMatchRecord*& M_i, vector< NeighborhoodGroup >& novel_subset_list,
				bool find_novel_subsets, ProcrastinationQueue& procrastination_queue,
				vector< gnSequence* >& seq_table, int direction, uint w, int& last_linked,
				size_t& novel_subset_count )
{
	// finally process novel subset
	bool prev_linked = false;	// we only want to link the closest of a group with the same components
	int created_thisround = 0;
	static NeighborhoodGroupComponentCompare srcc;
	for( size_t gI = 0; gI < novel_subset_list.size(); gI++ )
	{
		if( !find_novel_subsets )
			continue;	// only find novel subsets if we're supposed to
		if( last_linked != 0 )
			continue;	// only generate subsets when last_linked == 0

		// be sure to handle case where:
		// --M_i-->   --M_j--   <--M_i--
		// that may cause an identical novel subset to get spawned but with
		// M_i and M_j swapped as left and right supersets

		bool same_components = false;
		if( gI > 0 )
			same_components = srcc.compare(novel_subset_list[gI], novel_subset_list[gI-1]) == 0;
		prev_linked = same_components? prev_linked : false;

		if( prev_linked )
			continue;	// don't link a subset with the same components...

		// TODO: handle the tandem repeat case
		if( M_i->tandem )
			continue;

		MatchRecord* M_j = novel_subset_list[gI].get<0>();
		// if M_j hasn't been extended then we don't do anything yet.
		// we may find this novel subset again when M_j gets extended
		if( M_j->extended == false )
			continue;

		size_t mult = 0;	// multiplicity of novel subset
		for( size_t i = 0; i < novel_subset_list[gI].get<1>().size(); ++i )
			if( novel_subset_list[gI].get<1>()[i] != (std::numeric_limits<size_t>::max)() )
				mult++;

		if( mult < 2 )
			continue;	// can't do anything if there's no match!

		UngappedMatchRecord tmper1(mult,0);
		GappedMatchRecord tmper2(tmper1);  // this is lame
		GappedMatchRecord* M_n = tmper2.Copy();

		size_t mnewi = 0;
		vector< size_t > new_to_i_map(mult);
		vector< size_t > new_to_j_map(mult);
		boost::dynamic_bitset<> ni_list(M_i->Multiplicity());
		boost::dynamic_bitset<> nj_list(M_j->Multiplicity());
		for( size_t i = 0; i < novel_subset_list[gI].get<1>().size(); ++i )
		{
			if( novel_subset_list[gI].get<1>()[i] != (std::numeric_limits<size_t>::max)() )
			{
				new_to_i_map[mnewi] = novel_subset_list[gI].get<1>()[i];
				new_to_j_map[mnewi] = i;
				ni_list.set(new_to_i_map[mnewi]);
				nj_list.set(i);
				M_n->SetStart(mnewi, M_j->Start(i));	// sets left-end and orientation
				M_n->SetLength(M_j->Length(i),mnewi);
				mnewi++;
			}
		}
		if( M_n->Orientation(0) == AbstractMatch::reverse )
			M_n->Invert();
		// before we go any further, make sure that the relevant portion of M_i is not 
		// either completely or partially subsumed by the relevant portion of M_j!
		MatchProjectionAdapter mpaa( M_i, new_to_i_map );
		vector< size_t > mpaa_to_Mn_map( new_to_i_map.size() );
		for( size_t i = 0; i < mpaa_to_Mn_map.size(); ++i )
			mpaa_to_Mn_map[i] = i;
		bool subsumed;
		bool partial;
		classifyMatch( M_n, &mpaa, mpaa_to_Mn_map, subsumed, partial );
		if( subsumed )
		{
			M_n->Free();
			continue;	// there's nothing novel about this subset...
		}
		if( partial )
		{
			// FIXME: we should really spawn a novel subset on the non-subsumed components
			M_n->Free();
			continue;
		}
		created_thisround+= M_n->Multiplicity();

		M_n->chained_matches.push_back(M_j);
		M_n->chained_component_maps.push_back(new_to_j_map);
		//tjt: need to send finalize seq_table for muscle alignment
		M_n->finalize(seq_table);	// make this one a legitimate match...
		M_n->chained_matches.clear();
		M_n->chained_component_maps.clear();

		// create links from M_n to M_i and M_j
		int ni_parity = M_n->Orientation(0) == M_i->Orientation(new_to_i_map[0]) ? 1 : -1;
		int nj_parity = M_n->Orientation(0) == M_j->Orientation(new_to_j_map[0]) ? 1 : -1;
		MatchLink& ni_link = getSuperset(M_n,-direction*ni_parity);
		ni_link = MatchLink(M_i,M_n,ni_list,new_to_i_map);
		getSubsets(M_i,direction).push_back(ni_link);
		MatchLink& nj_link = getSuperset(M_n,direction*ni_parity);
		nj_link = MatchLink(M_j,M_n,nj_list,new_to_j_map);
		getSubsets(M_j,-direction*ni_parity*nj_parity).push_back(nj_link);

		// push M_n onto the heap
		novel_subset_list.push_back(M_n);
		procrastination_queue.push(M_n);
		novel_subset_count++;
	}

	if( created_thisround > w * M_i->Multiplicity() )
	{
		cerr << "made too many!\n";
		genome::breakHere();
	}
}

void writeXmfa( MatchList& seedml, std::vector< GappedMatchRecord* >& extended_matches, const std::string& xmfa_file )
{
	GenericIntervalList<GappedMatchRecord> gmr_list;
	for( size_t gmrI = 0; gmrI < extended_matches.size(); ++gmrI )
		gmr_list.push_back(*extended_matches[gmrI]);

	if( xmfa_file.length() > 0  && xmfa_file != "-")
	{
		gmr_list.seq_filename.push_back( seedml.seq_filename[0] );
		gmr_list.seq_table.push_back( seedml.seq_table[0] );
		if( xmfa_file == "-" )
			gmr_list.WriteStandardAlignment(cout);
		else
		{
			ofstream xmfa_out(xmfa_file.c_str());
			gmr_list.WriteStandardAlignment(xmfa_out);
			xmfa_out.close();
		}
	}
}

int main( int argc, char* argv[] )
{
//	debug_interval = true;
	// Declare the supported options.

	string sequence_file = "";
	
	unsigned w = 0;
	int kmersize =0;
	uint seed_weight = 0;
	string outputfile = "";
	string output2file = "";
	string xmfa_file = "";
	string stat_file = "";
	bool find_novel_subsets = false;
	bool solid_seed = false;
	bool extend_chains = true;
	bool chain = true;
	bool two_hits = false;
	bool unalign = true;

	po::variables_map vm;
	try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "get help message")
            ("sequence", po::value<string>(&sequence_file), "FastA sequence file")
			("w", po::value<unsigned>(&w)->default_value(0), "max gap width ")
			("z", po::value<unsigned>(&seed_weight)->default_value(0), "seed weight")
			("solid", po::value<bool>(&solid_seed)->default_value(0), "use solid seed")
			("two-hits", po::value<bool>(&two_hits)->default_value(false), "require two hits within w to trigger gapped extension")
			("unalign", po::value<bool>(&unalign)->default_value(true), "unalign non-homologous sequence")
			("chain", po::value<bool>(&chain)->default_value(true), "chain matches")
			("extend", po::value<bool>(&extend_chains)->default_value(true), "extend chains")
			("novel-subsets", po::value<bool>(&find_novel_subsets)->default_value(false), "find novel subset matches ")
			("output", po::value<string>(&outputfile)->default_value(""), "output ")
			("score-out", po::value<string>(&output2file)->default_value(""), "output with corresponding score and alignment info ")
			("highest", po::value<string>(&stat_file)->default_value("procrast.highest"), "file containing highest scoring alignment for each multiplicity ")
			("xmfa", po::value<string>(&xmfa_file)->default_value(""), "XMFA format output")
        ;

		if( argc < 2 )
		{
            cout << desc << "\n";
            return 1;
		}
                
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

		
        if (vm.count("w")) {
            cout << "max gap width (w) was set to " 
                 << w << ".\n";
        } else {
            cout << "Max gap width (w) not specified\n, using default value of 10\n";
        }

		if (vm.count("z")) {
            cout << "seed weight set to " 
                 << seed_weight << ".\n";
        } else {
            cout << "Using default seed weight.\n";
        }
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
	


	// Roadmap: 
	// 1. identify seed matches using a Sorted Mer List
	// 2. create a "UngappedMatchRecord" for each seed match and put in the match record list
	// 3. create a Match Position Lookup Table
	// 4. create a multiplicity priority queue
	// 5. create an (empty) Novel Subset Match Record list
	// 6. extend all matches!
	// 7. create a procrastination queue for novel subset matches
	// 8. extend all novel subset matches!
	// 9. create a final list of matches
	// 10. score matches
	// 11. report matches


	//
	// part 1, load sequence and find seed matches using SML and a repeat class...
	//
	MatchList seedml;
	seedml.seq_filename = vector< string >( 1, sequence_file );
	seedml.sml_filename = vector< string >( 1, seedml.seq_filename[0] + ".sml" );
	//seedml.LoadSequences( &cout );
	LoadSequences( seedml, &cout );
	if( seed_weight == 0 )
	{
		seed_weight = (int)((double)getDefaultSeedWeight( seedml.seq_table[0]->length() ) * .9);
	}
	
	int seed_rank = 0;
	if ( solid_seed )
	{
		seed_rank = INT_MAX;
		std::cout << "Using solid seed" << std::endl;
	}
	seedml.LoadSMLs( seed_weight, &cout, seed_rank );
	int64 seed = getSeed( seed_weight, seed_rank);
	uint seed_size = getSeedLength( seed );
	if( w == 0 )
		w = seed_weight * 3;	// default value
	
	cout << "Using seed weight: " << seed_weight << " and w: " << w << endl;
	SeedMatchEnumerator sme;
	sme.FindMatches( seedml );
	
	//
	// part 2, convert to match records
	//
	vector< UngappedMatchRecord* > match_record_list( seedml.size() );
	size_t component_count = 0;
	for( size_t mI = 0; mI < seedml.size(); ++mI )
	{
		UngappedMatchRecord tmp( seedml[mI]->SeqCount(), seedml[mI]->AlignmentLength() );
		match_record_list[mI] = tmp.Copy();
		
		for( size_t seqI = 0; seqI < seedml[mI]->SeqCount(); seqI++ )
		{
			match_record_list[mI]->SetStart( seqI, seedml[mI]->Start( seqI ) );
			match_record_list[mI]->SetLength( seedml[mI]->Length( seqI ), seqI );
		}
		component_count += seedml[mI]->SeqCount();
		//match_record_list[mI]->
		seedml[mI]->Free();
	}
	
	//
	// part 3, create a match position lookup table
	//
	vector< pair< gnSeqI, MatchPositionEntry > > mplt_sort_list( component_count );
	size_t compI = 0;
	for( size_t mI = 0; mI < match_record_list.size(); ++mI )
	{
		UngappedMatchRecord* mr = match_record_list[mI];
		for( size_t seqI = 0; seqI < mr->SeqCount(); ++seqI )
			mplt_sort_list[compI++] = make_pair( mr->LeftEnd( seqI ), make_pair( mr, seqI ) );
	}
	// pairs get ordered on the first element by default 
	std::sort( mplt_sort_list.begin(), mplt_sort_list.end() );
	gnSeqI seq_length = seedml.seq_table[0]->length();
	MatchPositionLookupTable match_pos_lookup_table( seq_length+1, make_pair( (UngappedMatchRecord*)NULL, 0 ) );
	for( size_t i = 0; i < mplt_sort_list.size(); ++i )
		match_pos_lookup_table[ mplt_sort_list[i].first ] = mplt_sort_list[i].second;


	//
	// part 4, create a procrastination queue
	//
	ProcrastinationQueue procrastination_queue( match_record_list );

	//
	// part 5, create an (empty) Novel Subset Match Record list
	//
	vector< GappedMatchRecord* > novel_subset_list;

	size_t superset_count = 0;
	size_t chainable_count = 0;
	size_t subset_count = 0;
	size_t novel_subset_count = 0;

	boost::dynamic_bitset<> left_lookups(seedml.seq_table[0]->length(), false);
	boost::dynamic_bitset<> right_lookups(seedml.seq_table[0]->length(), false);

	//
	// part 6, extend all matches!
	//
	vector< GappedMatchRecord* > extended_matches;	/**< The extended matches will be chains of UngappedMatchRecords */

	//for extension
	PairwiseScoringScheme pss = PairwiseScoringScheme(hoxd_matrix,-100,-20);
	
	int curI = 0;
	uint curr_extensions = 0;
	uint max_extensions = 2;
	while(  procrastination_queue.end() > 0 )
	{
		
		int prevI = curI;
		curI +=1;
		
		if( (curI * 100) / procrastination_queue.size() != (prevI * 100) / procrastination_queue.size() )
		{
			cout << (curI * 100) / procrastination_queue.size() << "%..";
			cout.flush();
		}
		
		// pop the next match off the heap
		MatchRecord* umr = procrastination_queue.pop(); 
		// if the match has been subsumed then skip it
		if( umr->subsuming_match != NULL )
			continue;
		if( umr->dont_extend == true )
			continue;

//		if( umr == (MatchRecord*)0x01335878 )
//			cout << "umr:\n" << *(UngappedMatchRecord*)umr << endl;

		GappedMatchRecord* M_i = dynamic_cast<GappedMatchRecord*>(umr);
		if( M_i == NULL )
		{
			// create a new gapped match record for M_i
			GappedMatchRecord gmr( *(UngappedMatchRecord*)umr );
			M_i = gmr.Copy();
			umr->subsuming_match = M_i;
			M_i->chained_matches.push_back( umr );
			vector< size_t > component_map( M_i->Multiplicity() );
			for( size_t i = 0; i < component_map.size(); ++i )
				component_map[i] = i;
			M_i->chained_component_maps.push_back(component_map);
			swap(umr->subsumption_component_map, component_map);	// swap avoids reallocation
			// update superset and subset links
			for( int dI = 1; dI > -2; dI -= 2 )
			{
				MatchLink& ij_link = getSuperset(M_i,dI);
				if( ij_link.superset != NULL )
				{
					ij_link.subset = M_i;
					unlinkSuperset(umr,dI);
					int parity = M_i->Orientation(0) == ij_link.superset->Orientation(ij_link.sub_to_super_map[0]) ? 1 : -1;
					getSubsets(ij_link.superset,-dI*parity).push_back(ij_link);
				}
				vector< MatchLink >& subsets = getSubsets(M_i,dI);
				for( size_t subI = 0; subI < subsets.size(); ++subI )
				{
					subsets[subI].superset = M_i;
					int parity = M_i->Orientation(subsets[subI].sub_to_super_map[0]) == subsets[subI].subset->Orientation(0) ? 1 : -1;
					getSuperset(subsets[subI].subset, -dI*parity).superset = M_i;
				}
				getSubsets(umr,dI).clear();	// so that validate works...
			}
		}

		if( M_i == (MatchRecord*)0x01839da4)
			cerr << "debugme!!\n";
		if( M_i == (MatchRecord*)0x017aaf24)
			cerr << "supersetdebugme!!\n";


		M_i->extended = true;
		extended_matches.push_back( M_i );
		
		// extend the match in each direction 
		// if a superset exists use that first
		// otherwise create a neighborhood list
		int direction = 1;	// leftward == 1, rightward == -1, done == -3
		int last_linked = 0;	// stores the group type that was chained.  1 == superset, 2 == chainable, 0 == none
		vector< NeighborhoodGroup > left_deferred_subsets;
		vector< NeighborhoodGroup > right_deferred_subsets;

		score_t score = 0;
		vector< gnSequence* > seqtable( M_i->SeqCount(), seedml.seq_table[0] );
		vector< string > alignment;
		vector<score_t> scores;
		
		while( direction > -2 )
		{
			last_linked = 0;
			
			// check for superset
			if( getSuperset(M_i, direction).superset != NULL )
			{
				supersetLinkExtension( M_i, direction, last_linked, left_deferred_subsets, right_deferred_subsets );
			}	// end if there was a superset
			else
			{
				//
				// perform a neighborhood list extension, 
				// looks for neighboring matches in the match position lookup table
				// 
				vector< NeighborhoodGroup > superset_list;
				vector< NeighborhoodGroup > chainable_list;
				vector< NeighborhoodGroup > subset_list;
				vector< NeighborhoodGroup > novel_subset_list;
				neighborhoodListLookup( M_i, match_pos_lookup_table,
								superset_list, chainable_list, subset_list, novel_subset_list,
								direction, seed_size, w, left_lookups, right_lookups);

				// tallies for debugging
				superset_count += superset_list.size();
				chainable_count += chainable_list.size();
				subset_count += subset_list.size();
				
				// now process each type of neighborhood group
				// supersets are already done.  happy chrismakwanzuhkkah

				// then process chainable
				processChainableMatches( M_i, chainable_list, direction, last_linked );

				// defer subset processing
				for( size_t gI = 0; gI < subset_list.size(); gI++ )
				{
					vector< NeighborhoodGroup >& cur_subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
					cur_subset_list.push_back( subset_list[gI] );
				}

				// finally process novel subset
				processNovelSubsetMatches(M_i, novel_subset_list, find_novel_subsets, procrastination_queue, 
					seedml.seq_table, direction, w, last_linked, novel_subset_count );

			} // end if no superset was found then do neighborhood list lookup


			// if we didn't do a chaining or superset extension, try a gapped extension
			if( last_linked == 0 )
			{
				if ( !extend_chains )
				{
					direction -= 2;
					continue;
				}
				bool failed = true;
				// only extend if two matches are chained if two-hits == true
				if( !two_hits || (two_hits && M_i->chained_matches.size() > 1 ))
				{	
					cerr << "Preparing to extend Chain #" << curI << endl;
					cerr << "Direction = " << direction << endl;
					failed = ExtendMatch(M_i, seqtable, pss, w, direction);
				}
				if (failed )
				{
					//end gapped extension  whenever extension fails.
					direction -=2;
					continue;
				}

				//update links appropriately, and we can take another round
				//through the evil megaloop, possibly discovering additional chainable
				//seeds or superset links.
				
				// need to update links by looking for matches in the region that was just extended over
				vector< NeighborhoodGroup > superset_list;
				vector< NeighborhoodGroup > chainable_list;
				vector< NeighborhoodGroup > subset_list;
				vector< NeighborhoodGroup > novel_subset_list;
				neighborhoodListLookup( M_i, match_pos_lookup_table,
								superset_list, chainable_list, subset_list, novel_subset_list,
								direction, seed_size, w, left_lookups, right_lookups);
					
				// now process each type of neighborhood group
				// FIXME: need to process supersets

				// how to process supersets?
				// if we have completely extended through a superset
				//   then we want to replace that part of the alignment with the superset
				// if the superset continues beyond the end of at least one component, then 
				// we want to create a superset link for it, and process it during a link extension
				// what to do when the superset doesn't match very well with the 
				MatchRecord* M_j = getSuperset(M_i, direction).superset;
				if ( M_j != NULL )
				{
					MatchLink ij_link = getSuperset(M_i, direction);	
					int ij_parity = M_i->Orientation(0) == M_j->Orientation(ij_link.sub_to_super_map[0]) ? 1 : -1;
					mems::MatchProjectionAdapter mpaa( M_j, M_i->chained_component_maps[0] );
					int overextended_components = 0;
					for( size_t seqI = 0; seqI < mpaa.Multiplicity(); seqI++ )
					{
						if ( direction < 0 )
						{
							if( M_i->RightEnd(seqI) > mpaa.RightEnd(seqI) )
							{
								//uh oh, extended past superset!
								overextended_components+=1;
							}
						}
						else
						{
							if( M_i->LeftEnd(seqI) < mpaa.LeftEnd(seqI) )
							{
								//uh oh, extended past superset!	
								overextended_components+=1;
							}			
						}
					}
					if ( overextended_components == mpaa.Multiplicity() )
					{
						//   then we want to replace that part of the alignment with the superset
						//reset to superset boundaries
						bool changed = extendRange( M_i, M_j, ij_link.sub_to_super_map );
					}
					else if ( overextended_components > 0 )
					{
						// we want to create a superset link for it, and process it during a link extension
						inheritSuperset( M_i, M_j, direction, ij_parity );
						last_linked = 2;	// we may do a link extension!
					}
				
				}
				// then process chainable
				processChainableMatches( M_i, chainable_list, direction, last_linked );

				// defer subset processing
				for( size_t gI = 0; gI < subset_list.size(); gI++ )
				{
					vector< NeighborhoodGroup >& cur_subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
					cur_subset_list.push_back( subset_list[gI] );
				}

				// finally process novel subset
				processNovelSubsetMatches(M_i, novel_subset_list, find_novel_subsets, procrastination_queue, 
					seedml.seq_table, direction, w, last_linked, novel_subset_count );
				
			}
		}	// end loop over leftward and rightward extension

		//
		// finalize the alignment -- this resolves overlapping components into a single gapped alignment
		//
		//tjt: need to send finalize seq_table for muscle alignment
		if( M_i == (GappedMatchRecord*)0x00d37364 )
			cerr << "debugmult\n";
		//tjt: make sure finalize only gets called once!
		M_i->finalize(seedml.seq_table);
		if(extended_matches.size() % 25 == 0 && extended_matches.size() / 25 > 0 )
			writeXmfa( seedml, extended_matches, xmfa_file );
	
		//
		// process deferred subsets
		//
		for( int direction = 1; direction > -2; direction -= 2 )
		{
			vector< NeighborhoodGroup >& subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
			NeighborhoodGroupCompare ngc;
			NeighborhoodGroupComponentCompare ngcc;
			std::sort( subset_list.begin(), subset_list.end(), ngc );
			bool prev_linked = false;
			for( size_t sI = 0; sI < subset_list.size(); ++sI )
			{
				bool same_components = false;
				if( sI > 0 )
					same_components = ngcc.compare(subset_list[sI], subset_list[sI-1]) == 0;
				prev_linked = same_components? prev_linked : false;

				// check whether each of these ended up getting subsumed
				bool subsumed;
				bool partial;
				classifySubset( M_i, subset_list[sI], subsumed, partial );
				MatchRecord* M_j = subset_list[sI].get<0>();
				if( subsumed )
				{
					M_j->subsuming_match = M_i;
					M_j->subsumption_component_map = subset_list[sI].get<1>();
					unlinkSupersets(M_j);
					continue;
				}
				if( partial )
				{
					// create a novel subset record, mark this one as subsumed
					// just destroy it for now...
					M_j->dont_extend = true;
					unlinkSupersets(M_j);
					for( size_t mjI = 0; mjI < M_j->Multiplicity(); ++mjI )
						if( match_pos_lookup_table[M_j->LeftEnd(mjI)].first == M_j )
							match_pos_lookup_table[M_j->LeftEnd(mjI)] = make_pair((MatchRecord*)NULL,0);
					continue;
				}

				if( prev_linked )
				{
					// the previous subset has the same components as this one and was linked.
					// we may consider this one an 'extra' if all components are further away
					NeighborhoodGroup cur_group = subset_list[sI];
					subset_list.erase(subset_list.begin() + sI, subset_list.begin() + sI + 1);
					sI--;
					size_t dI = 0;
					for( ; dI < cur_group.get<2>().size(); ++dI )
						if( cur_group.get<2>()[dI] <= subset_list[sI].get<2>()[dI] )
							break;

					if( dI == cur_group.get<2>().size() )
					{
						// include this in a list of extra subsets
						boost::dynamic_bitset<> tmp_bs(M_i->Multiplicity());
						getExtraSubsets( M_i, direction ).push_back( MatchLink( (MatchRecord*)M_i, M_j, tmp_bs, cur_group.get<1>() ) );
						continue;
					}
					// we've got a subset tie.
					cerr << "Subset tie, erasing M_j\n";
					M_j->dont_extend = true;
					unlinkSupersets(M_j);
					for( size_t mjI = 0; mjI < M_j->Multiplicity(); ++mjI )
						if( match_pos_lookup_table[M_j->LeftEnd(mjI)].first == M_j )
							match_pos_lookup_table[M_j->LeftEnd(mjI)] = make_pair((MatchRecord*)NULL,0);
					continue;
				}

				int parity = M_i->Orientation( subset_list[sI].get<1>()[0] ) == M_j->Orientation(0) ? 1 : -1;
				// if we have the following case:
				// --M_i-->   --M_j--   <--M_i-- ... ... --M_i-->   --M_j--   <--M_i--
				// then M_j may already be linked to M_i but on the other side
				if( getSuperset(M_j,direction*parity).superset == M_i )
					continue;
				unlinkSuperset( M_j, -direction*parity );
				// it's outside, just link it in
				// rebuild the superset component list
				boost::dynamic_bitset<> comp_list(M_i->Multiplicity(), false);
				for( size_t compI = 0; compI < subset_list[sI].get<1>().size(); ++compI )
					comp_list.set(subset_list[sI].get<1>()[compI]);
				getSuperset(M_j,-direction*parity) = MatchLink( M_i, M_j, comp_list, subset_list[sI].get<1>() );
				getSubsets(M_i,direction).push_back( getSuperset(M_j,-direction*parity));
				prev_linked = true;
			}
			subset_list.clear();
		}

		
//		validate( extended_matches );
//		validate( match_record_list );
	}

	cout << "superset count: " << superset_count << endl;
	cout << "chainable count: " << chainable_count << endl;
	cout << "subset count: " << subset_count << endl;
	cout << "novel subset count: " << novel_subset_count << endl;
	
	//write output!!!!!

	// 
	// part 9, create a final list of local multiple alignments (already done in extended_matches)
	//
	vector< GappedMatchRecord* >& final = extended_matches;
	writeXmfa( seedml, extended_matches, xmfa_file );

	// 
	// part 10, score matches
	
	//create output stream
	ostream* output;
	ostream* output2;
 	ofstream score_out_file;
	ofstream aln_out_file;
	ofstream stats_out_file;
	if(stat_file != "" && stat_file != "-")
		stats_out_file.open( stat_file.c_str() );

	if(outputfile == "" || outputfile == "-")
		output = &cout;
	else
	{
		aln_out_file.open( outputfile.c_str() );
		output = &aln_out_file;
	}
	if(output2file == "" || output2file == "-")
		output2 = &cout;
	else
	{
		score_out_file.open( output2file.c_str() );
		output2 = &score_out_file;
	}
	vector< pair< double, GappedMatchRecord* > > scored( final.size() );
	vector<score_t> scores_final;
	score_t score_final = 0;
	for( size_t fI = 0; fI < final.size(); fI++ )
	{
	    vector<string> alignment;
		vector< gnSequence* > seq_table( final[fI]->SeqCount(), seedml.seq_table[0] );
		mems::GetAlignment(*final[fI], seq_table, alignment);	// expects one seq_table entry per matching component
		//send temporary output format to file if requested 
		*output << "#procrastAlignment " << fI+1 << endl << *final.at(fI) << endl;
		computeSPScore( alignment, pss, scores_final, score_final);
		scored[fI] = make_pair( score_final, final[fI] );
	}
	std::sort( scored.begin(), scored.end() );
	std::reverse( scored.begin(), scored.end() );

	// 
	// part 11, report matches in scored order
	//
	output2->setf(ios::fixed);
	output2->precision(0);
	for( size_t sI = 0; sI < scored.size(); ++sI )
	{
		*output2 << "#procrastAlignment " << sI+1 << endl << *scored[sI].second << endl;
		*output2 << "Alignment length: " << scored[sI].second->AlignmentLength() << endl;
		*output2 << "Score: " << scored[sI].first << endl;
	}
	
	///report highest scoring lma for each multiplicity
	vector<int> multiplicity_list;
	vector< pair< int,  pair<double, GappedMatchRecord*> > > ordered;
	for( size_t tI = 0; tI < scored.size(); ++tI )
	{ 
		if ( find( multiplicity_list.begin(), multiplicity_list.end(), scored[tI].second->Multiplicity() )  == multiplicity_list.end() )
		{
			multiplicity_list.push_back( scored[tI].second->Multiplicity() );
			ordered.push_back(make_pair( scored[tI].second->Multiplicity(), scored[tI] ) );
		}
	}

	
	std::sort(ordered.begin(), ordered.end());
	stats_out_file.setf(ios::fixed);
	stats_out_file.precision(0);
	for( size_t tI = 0; tI < ordered.size(); ++tI )
	{
		stats_out_file << "#" << tI+1 << ": r= " << ordered[tI].first << " l= " << ordered[tI].second.second->AlignmentLength() << " s= " << ordered[tI].second.first << endl;
	}
	
	// punt: part 12, finally add any unaligned regions to the interval list	
	//if( gapped_alignment )
	//addUnalignedIntervals( interval_list );

	//score matches again!
	//call calculateSPScore/calculateConsensusScore on each chain
	//std::sort( scored.begin(), scored.end() );
 
	
	// clean up
	for( size_t eI = 0; eI < match_record_list.size(); ++eI )
		match_record_list[eI]->Free();
	for( size_t eI = 0; eI < novel_subset_list.size(); ++eI )
		if( novel_subset_list[eI]->subsuming_match != NULL || novel_subset_list[eI]->dont_extend )
			novel_subset_list[eI]->Free();
	for( size_t eI = 0; eI < extended_matches.size(); ++eI )
		if( extended_matches[eI]->subsuming_match == NULL && !extended_matches[eI]->dont_extend )
			extended_matches[eI]->Free();

	for( size_t seqI = 0; seqI < seedml.seq_table.size(); ++seqI )
		delete seedml.seq_table[seqI];
	for( size_t seqI = 0; seqI < seedml.sml_table.size(); ++seqI )
		delete seedml.sml_table[seqI];

	return 0;
}

