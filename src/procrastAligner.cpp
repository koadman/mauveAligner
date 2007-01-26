#include "libGenome/gnSequence.h"
#include "libMems/Interval.h"
#include "ProgressiveAligner.h"

#include <iostream>
#include <algorithm>
#include <cctype>

#include "MatchRecord.h"
#include "SeedMatchEnumerator.h"
#include "procrastUtilities.h"

#include <boost/tuple/tuple.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace genome;
using namespace mems;

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

/** defines a multiplicity heap ordering */
class MultiplicityHeapCompare
{
public:
	bool operator()( const MatchRecord* a, const MatchRecord* b )
	{
		return a->Multiplicity() < b->Multiplicity();
	}
};

/** The subset record contains the subset pointer, the component map to the superset, and a vector of distances */
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

/** returns one of the superset links from a match.  direction is 1 for left, -1 for right */
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

	NeighborhoodGroupCompare src;
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

int main( int argc, char* argv[] )
{
//	debug_interval = true;
	// Declare the supported options.

	string sequence_file = "";
	unsigned w = 0;
	int kmersize =0;
	uint seed_weight = 0;
	string outputfile = "";
	string xmfa_file = "";
	bool find_novel_subsets = false;

	po::variables_map vm;
	try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "get help message")
            ("sequence", po::value<string>(&sequence_file), "FastA sequence file")
			("w", po::value<unsigned>(&w)->default_value(0), "max gap width ")
			("z", po::value<unsigned>(&seed_weight)->default_value(0), "seed weight")
			("b", po::value<int>(&kmersize)->default_value(1), "kmer background freq size ")
			("novel-subsets", po::value<bool>(&find_novel_subsets)->default_value(false), "find novel subset matches ")
			("output", po::value<string>(&outputfile)->default_value(""), "output ")
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

		if (vm.count("b")) {
            cout << "kmer background freq size (b) was set to " 
                 << kmersize << ".\n";
        } else {
            cout << "kmer background freq size (b) was not specified\n, using default value of 1\n";
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
	LoadSequences( seedml, &cout );
	if( seed_weight == 0 )
		seed_weight = (int)((double)getDefaultSeedWeight( seedml.seq_table[0]->length() ) * .9);
	seedml.LoadSMLs( seed_weight, &cout );
	int64 seed = getSeed( seed_weight );
	uint seed_size = getSeedLength( seed );
	if( w == 0 )
		w = seed_weight * 4;	// default value
	cout << "Using seed weight: " << seed_weight << " and w: " << w << endl;
	SeedMatchEnumerator sme;
	sme.FindMatches( seedml );

	// need single nuc & kmer frequency
	
	string sequence = seedml.seq_table.at(0)->ToString();
	string uppercase = sequence;
	ToUPPER tupperware;
	std::transform(sequence.begin(),sequence.end(), uppercase.begin(), tupperware);

	map<string,gnSeqI> polyfreq;
	map<string,gnSeqI> monofreq;
	map<string, gnSeqI>::iterator it;
	for (gnSeqI i = 0; i <= uppercase.size()-kmersize; i++)
	{
	   string kmer = uppercase.substr(i,kmersize);
	   string nucleotide = uppercase.substr(i,1);
	   if( nucleotide[0] != 'A' &&
			nucleotide[0] != 'C' &&
			nucleotide[0] != 'G' &&
			nucleotide[0] != 'T' )
			nucleotide[0] = 'A';
	   for( size_t kI = 0; kI < kmer.size(); kI++ )
		   if( kmer[kI] != 'A' &&
				kmer[kI] != 'C' &&
				kmer[kI] != 'G' &&
				kmer[kI] != 'T' )
				kmer[kI] = 'A';

	   polyfreq[kmer] +=1;	
	   monofreq[nucleotide] +=1;	
	   //insert( const string& val );
	   //it = find( const string& mer );
       //it->second+=1;	   
	}
	
	int K = (int)pow(double(4.0),double(kmersize));
	genome::Array<double> pu_arr(K);
	double* pu = pu_arr.data;
	for( size_t kI = 0; kI < K; kI++ )
		pu[kI] = 1;	// add pseudocounts (laplace estimates) of 1 so that we don't end up with a 0 if some k-mer doesn't exist
	int kcount = 0;
    for( map<string, gnSeqI>::iterator iter = polyfreq.begin(); iter != polyfreq.end(); iter++ ) 
	{
		cout << (*iter).first << " occurs " << (*iter).second << " times" << endl;
		//Background frequencies, calculate once from sequence string
		//if doesnt appear in original sequenced, will not be stored
		if((*iter).second > 0)
			pu[kcount] = (float)(*iter).second/(float)sequence.size();
		else
			pu[kcount] = 0.0;
		kcount++;

    }

	vector<float> nucFrequency;
	nucFrequency.push_back(float(monofreq["A"])/float(sequence.size()));
	nucFrequency.push_back(float(monofreq["T"])/float(sequence.size()));
	nucFrequency.push_back(float(monofreq["G"])/float(sequence.size()));
	nucFrequency.push_back(float(monofreq["C"])/float(sequence.size()));

	//
	// part 2, convert to match records
	//
	vector< UngappedMatchRecord* > match_record_list( seedml.size() );
	for( size_t mI = 0; mI < seedml.size(); ++mI )
	{
		UngappedMatchRecord tmp( seedml[mI]->SeqCount(), seedml[mI]->AlignmentLength() );
		match_record_list[mI] = tmp.Copy();

		for( size_t seqI = 0; seqI < seedml[mI]->SeqCount(); seqI++ )
		{
			match_record_list[mI]->SetStart( seqI, seedml[mI]->Start( seqI ) );
			match_record_list[mI]->SetLength( seedml[mI]->Length( seqI ), seqI );
		}
		seedml[mI]->Free();
	}
	
	//
	// part 3, create a match position lookup table
	//
	vector< pair< gnSeqI, MatchPositionEntry > > mplt_sort_list;
	for( size_t mI = 0; mI < match_record_list.size(); ++mI )
	{
		UngappedMatchRecord* mr = match_record_list[mI];
		for( size_t seqI = 0; seqI < mr->SeqCount(); ++seqI )
			mplt_sort_list.push_back( make_pair( mr->LeftEnd( seqI ), make_pair( mr, seqI ) ) );
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
	std::vector< MatchRecord* > procrastination_queue( match_record_list.size() );
	std::copy(match_record_list.begin(), match_record_list.end(), procrastination_queue.begin() );
	MultiplicityHeapCompare mhc;
	std::make_heap( procrastination_queue.begin(), procrastination_queue.end(), mhc );

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
	size_t queue_end = procrastination_queue.size();
	while( queue_end > 0 )
	{
		// pop the next match off the heap
		std::pop_heap( procrastination_queue.begin(), procrastination_queue.begin()+queue_end, mhc );
		queue_end--;
		MatchRecord* umr = *(procrastination_queue.begin() + queue_end);

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

		extended_matches.push_back( M_i );
		M_i->extended = true;

//		if( M_i == (GappedMatchRecord*)0x012b7658 )
//			cout << "M_i:\n" << *(GappedMatchRecord*)M_i << endl;
//		if( M_i == (GappedMatchRecord*)0x0133ab78 )
//			cout << "M_i:\n" << *(GappedMatchRecord*)M_i << endl;
//		if( M_i == (GappedMatchRecord*)0x00f57e88 )
//			cout << "M_i:\n" << *(GappedMatchRecord*)M_i << endl;

		// extend the match in each direction 
		// if a superset exists use that first
		// otherwise create a neighborhood list
		int direction = 1;	// leftward == 1, rightward == -1, done == -3
		int last_linked = 0;	// stores the group type that was chained.  1 == superset, 2 == chainable, 0 == none
		vector< NeighborhoodGroup > left_deferred_subsets;
		vector< NeighborhoodGroup > right_deferred_subsets;
		while( direction > -2 )
		{
			last_linked = 0;

			// check for superset
			if( getSuperset(M_i, direction).superset != NULL )
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

				//
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
					cerr << "something is wrong, we should never have supersets during link extension!\n";
					breakHere();
				}

				for( size_t cI = 0; cI < chainable.size(); ++cI )
				{
					if( chainable.size() > 1 )
					{
//						cerr << "bad news bruthah\n";
//						genome::breakHere();
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

			}	// end if there was a superset
			else
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

				superset_count += superset_groups.size();
				chainable_count += chainable_groups.size();
				subset_count += subset_groups.size();

				//
				// process each group
				//
				
				// process supersets first
				for( size_t gI = 0; gI < superset_groups.size(); gI++ )
				{
					// supersets are already done.  happy chrismakwanzuhkkah
				}

				// then process chainable
				// link the closest possible chainable first.
				vector< NeighborhoodGroup > cur_chainable;
				createNeighborhoodGroupList( cur_chainable, chainable_groups, neighborhood_list );
				for( size_t gI = 0; gI < cur_chainable.size(); gI++ )
				{
					MatchRecord* M_j = cur_chainable[gI].get<0>();

					vector< size_t >& component_map = cur_chainable[gI].get<1>();

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
					classifySubset( M_i, cur_chainable[gI], subsumed, partial );

					M_j->subsuming_match = M_i;
					M_j->subsumption_component_map = component_map;
					vector< size_t > xy_map(chainable_groups[gI].size());
					for( vector< size_t >::iterator rec_iter = chainable_groups[gI].begin(); rec_iter !=  chainable_groups[gI].end(); ++rec_iter )
						xy_map[ neighborhood_list[*rec_iter].Mi_component ] = neighborhood_list[*rec_iter].Mj_component;

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

				// then process subset
				// link the closest possible subset of a given type.
				vector< NeighborhoodGroup > cur_subsets;
				NeighborhoodGroupComponentCompare srcc;
				createNeighborhoodGroupList( cur_subsets, subset_groups, neighborhood_list );
				for( size_t gI = 0; gI < cur_subsets.size(); gI++ )
				{
					vector< NeighborhoodGroup >& subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
					subset_list.push_back( cur_subsets[gI] );
				}

				// finally process novel subset
				vector< NeighborhoodGroup > cur_novel;
				createNeighborhoodGroupList( cur_novel, novel_subset_groups, neighborhood_list );
				bool prev_linked = false;	// we only want to link the closest of a group with the same components
				int created_thisround = 0;
				for( size_t gI = 0; gI < cur_novel.size(); gI++ )
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
						same_components = srcc.compare(cur_novel[gI], cur_novel[gI-1]) == 0;
					prev_linked = same_components? prev_linked : false;

					if( prev_linked )
						continue;	// don't link a subset with the same components...

					// TODO: handle the tandem repeat case
					if( M_i->tandem )
						continue;

					MatchRecord* M_j = cur_novel[gI].get<0>();
					// if M_j hasn't been extended then we don't do anything yet.
					// we may find this novel subset again when M_j gets extended
					if( M_j->extended == false )
						continue;

					size_t mult = 0;
					for( size_t i = 0; i < cur_novel[gI].get<1>().size(); ++i )
						if( cur_novel[gI].get<1>()[i] != (std::numeric_limits<size_t>::max)() )
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
					for( size_t i = 0; i < cur_novel[gI].get<1>().size(); ++i )
					{
						if( cur_novel[gI].get<1>()[i] != (std::numeric_limits<size_t>::max)() )
						{
							new_to_i_map[mnewi] = cur_novel[gI].get<1>()[i];
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
					M_n->finalize();	// make this one a legitimate match...
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
					if( queue_end < procrastination_queue.size() )
						procrastination_queue[queue_end] = M_n;
					else
						procrastination_queue.push_back(M_n);
					queue_end++;
					std::push_heap(procrastination_queue.begin(), procrastination_queue.begin()+queue_end, mhc);
					novel_subset_count++;
				}
				if( created_thisround > w * M_i->Multiplicity() )
				{
					cerr << "made too many!\n";
					genome::breakHere();
				}

			} // end if no superset was found then build neighborhood list

			if( last_linked == 0 )
			{
				direction -= 2;	// if we didn't extend with a superset or chainable then change directions
			}

		}	// end loop over leftward and rightward extension


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
				// if we have he following case:
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

		//
		// finalize the alignment -- this resolves overlapping components into a single gapped alignment
		//

		M_i->finalize();

//		if( M_i == (GappedMatchRecord*)0x012b7658 )
//			cout << "M_i:\n" << *(GappedMatchRecord*)M_i << endl;
//		if( M_i == (GappedMatchRecord*)0x00f57e88 )
//			cout << "M_i:\n" << *(GappedMatchRecord*)M_i << endl;
//		if( M_i == (GappedMatchRecord*)0x00c9f218 )
//			cout << "M_i:\n" << *(GappedMatchRecord*)M_i << endl;

//		validate( extended_matches );
//		validate( match_record_list );
	}

	cout << "superset count: " << superset_count << endl;
	cout << "chainable count: " << chainable_count << endl;
	cout << "subset count: " << subset_count << endl;
	cout << "novel subset count: " << novel_subset_count << endl;
	
	// 
	// part 9, create a final list of local multiple alignments (already done in extended_matches)
	//
	vector< GappedMatchRecord* >& final = extended_matches;

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
	// 
	// part 10, score matches
	
	//create output stream
	ostream* output;
 	ofstream aln_out_file;
	if(outputfile == "" || outputfile == "-")
		output = &cout;
	else
	{
		aln_out_file.open( outputfile.c_str() );
		output = &aln_out_file;
	}
	vector< pair< double, GappedMatchRecord* > > scored( final.size() );
	string alphabet = "ATGC";
    vector<double> entropyscore;
	for( size_t fI = 0; fI < final.size(); fI++ )
	{
		//entropy score here we go!
	    vector<string> alignment;
		vector< gnSequence* > seq_table( final[fI]->SeqCount(), seedml.seq_table[0] );
		mems::GetAlignment(*final[fI], seq_table, alignment);	// expects one seq_table entry per matching component
		//send temporary output format to file if requested 
		*output << "#procrastAlignment " << fI+1 << endl << *final.at(fI) << endl;
		double entropy = 0.0;
		for( gnSeqI i = 0; i < alignment.at(0).size(); i++)
		{
			for ( int k = 0; k < alphabet.size(); k++)
			{
				double occurrences = 0;
				string letter = alphabet.substr(k,1);
				for( int s = 0; s < alignment.size(); s++)
				{
					char aln_s_i = toupper(alignment[s][i]);
					if( aln_s_i == alphabet[k] )
						occurrences +=1;
					else if( aln_s_i != 'A' &&  
							aln_s_i != 'C' && 
							aln_s_i != 'G' && 
							aln_s_i != 'T'  )
						occurrences += .25;	// unknown
				}
				//tjt: equation given in Uri Keich Bioinformatics Paper for calculating
				//	   the entropy score for a given column i of the alignment
				//double(((occurences*alignment.size())/nucFrequency(k)))
				if( occurrences != 0 )	// log of 0 is undefined
					entropy += occurrences * log( (double) (occurrences/alignment.size())/nucFrequency.at(k) );
			}
		}
		entropyscore.push_back(entropy);

		//Log-likelihood score, entropy score from alignment? int?
		double s = entropy;

		// alignment row count
		gnSeqI N = alignment[0].length();

		//Lattice size
		int Q = 16384;

		double score = s;
		string pvalue="1";
		/*
		if (!CallbagFFT("bagFFT",N, K, s, pu, Q,pvalue))
			cout << "pvalue not calculated, error" << endl;
		else
			cout << "bagFFT result: " << pvalue << endl;
		double score = atof(pvalue.c_str());
		*/
		/**FIX MEMORY LEAK BEFORE USING THIS CODE **/
/*		bagFFT_return_value result = bagFFT(N, K, s, pu, Q);
		cout << "Entropy: " << entropy << "\tp-Value = " << result.log_pvalue << " (" << exp(result.log_pvalue) << ")" << endl;
		double score = result.log_pvalue;
		**/
		
		scored[fI] = make_pair( score, final[fI] );
	}
	std::sort( scored.begin(), scored.end() );


	// 
	// part 11, report matches in scored order
	//
	for( size_t sI = 0; sI < scored.size(); ++sI )
	{
		cout << "Alignment length: " << scored[sI].second->AlignmentLength() << endl;
		cout << "Score: " << scored[sI].first << endl;
		cout << *scored[sI].second << endl;
	}
	
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

