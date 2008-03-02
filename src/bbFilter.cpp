#include "libMems/Backbone.h"
using namespace mems;
using namespace std;
using namespace genome;

typedef pair< bb_seqentry_t, size_t > labeled_bb_t;

class BbSorter
{
public:
	BbSorter( size_t seqI ){ m_seq = seqI; }
	bool operator()( const labeled_bb_t& a, const labeled_bb_t& b )
	{
		return abs(a.first[m_seq].first) < abs(b.first[m_seq].first);
	}
	size_t m_seq;
};


class BbSeqEntrySorter
{
public:
	BbSeqEntrySorter( size_t seqI ){ m_seq = seqI; }
	bool operator()( const bb_seqentry_t& a, const bb_seqentry_t& b )
	{
		return abs(a[m_seq].first) < abs(b[m_seq].first);
	}
	size_t m_seq;
};

// add unique segments of some minimum length
// FIXME: does not add begin and end segments!
void addUniqueSegments( std::vector< bb_seqentry_t >& bb_seq_list, size_t min_length = 20 )
{
	if( bb_seq_list.size() == 0 )
		return;
	vector< bb_seqentry_t > new_segs;
	uint seq_count = bb_seq_list[0].size();
	// now mark segs that are too close to each other to be considered independent
	for( size_t sI = 0; sI < seq_count; sI++ )
	{
		BbSeqEntrySorter bbs(sI);
		std::sort( bb_seq_list.begin(), bb_seq_list.end(), bbs );
		for( size_t bbI = 1; bbI < bb_seq_list.size(); bbI++ )
		{
			if( bb_seq_list[bbI][sI].first == 0 )
				continue;
			int64 diff = abs(bb_seq_list[bbI][sI].first) - abs(bb_seq_list[bbI-1][sI].second); 
			if( abs(diff) > min_length )
			{
				bb_seqentry_t newb( seq_count, make_pair( 0,0 ) );
				newb[sI].first = abs(bb_seq_list[bbI-1][sI].second) + 1;
				newb[sI].second = abs(bb_seq_list[bbI][sI].first) - 1;
				new_segs.push_back( newb );
			}
		}
	}
	cout << "Adding " << new_segs.size() << " genome-specific features\n";
	bb_seq_list.insert( bb_seq_list.end(), new_segs.begin(), new_segs.end() );
}

int main( int argc, char* argv[] )
{
#if	WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	if( argc < 4 )
	{
		cerr << "bbFilter <backbone file> <independent dist> <output file> <beast|gp> [<seq1> <seq2>...<seqN>]\n";
		cerr << "seq index starts at 0.\n";
		cerr << "\nExample:\n";
		cerr << "bbFilter my_alignment.backbone 50 my_feats.bin gp\n";
		cerr << "the above command extracts binary features from \"my_alignment.backbone\" which are separated by a minimum of 50nt sequence conserved among all taxa in the alignment.  The output is written to my_feats.bin in genoplast format\n";
		cerr << "\n\nExample 2:\nbbFilter aln.backbone 100 feats.xml beast 0 1 2 5 6\n";
		cerr << "the above command extracts binary features from \"aln.backbone\" which are separated by a minimum of 100nt sequence conserved among genomes 0,1,2,5, and 6 from the alignment.  The output is written to feats.xml in beast format\n";
		return -1;
	}
	string bbseq_fname( argv[1] );
	int indie_dist = atoi( argv[2] );
	string output_fname( argv[3] );
	string target_format( argv[4] );

	ifstream bbseq_input( bbseq_fname.c_str() );
	if( !bbseq_input.is_open() ){
		cerr << "Error opening \"" << bbseq_fname << "\"" << endl;
		return -4;
	}
	ofstream anal_output( output_fname.c_str() );
	if( !anal_output.is_open() ){
		cerr << "Error opening \"" << output_fname << "\" for writing" << endl;
		return -6;
	}
	
	// read the backbone column file	
	vector< bb_seqentry_t > bb_seq_list;
	readBackboneSeqFile( bbseq_input, bb_seq_list );

	// read the list of seqs of interest
	vector< int > seqs;
	for( int i = 5; i < argc; i++ )
		seqs.push_back(atoi(argv[i]));

	// assume all seqs are of interest
	if( seqs.size() == 0 && bb_seq_list.size() > 0 )
	{
		for( int i = 0; i < bb_seq_list[0].size(); i++ )
			seqs.push_back(i);
	}
	// now assign tracking IDs to the backbone segments
	vector< labeled_bb_t > bb_segs;
	for( size_t i = 0; i < bb_seq_list.size(); i++ )
	{
		bb_segs.push_back( make_pair( bb_seq_list[i], i ) );
	}

	bitset_t good_bb( bb_seq_list.size() );
	bitset_t nway( bb_seq_list.size() );
	bitset_t nunya( bb_seq_list.size() );

	// mark anything that has all of the seqs or none of seqs as not useful
	for( size_t bbI = 0; bbI < bb_seq_list.size(); bbI++ )
	{
		bool all = true;
		bool none = true;
		for( size_t sI = 0; sI < seqs.size(); sI++ )
		{
			if( bb_seq_list[bbI][seqs[sI]].first == 0 )
				all = false;
			else
				none = false;
		}
		if(all)
			nway.set(bbI);
		if(none)
			nunya.set(bbI);
	}
	good_bb = nway | nunya;
	good_bb.flip();
	
	// now mark segs that are too close to each other to be considered independent
	for( size_t sI = 0; sI < seqs.size(); sI++ )
	{
		BbSorter bbs(seqs[sI]);
		std::sort( bb_segs.begin(), bb_segs.end(), bbs );
		for( size_t bbI = 1; bbI < bb_segs.size()-1; bbI++ )
		{
			if( nway[bb_segs[bbI].second] )
				continue;
			if( bb_segs[bbI].first[seqs[sI]].first == 0 )
				continue;
			// ensure that it has n-way on both sides and that they are at least "indie_dist" long
			if( nway.test(bb_segs[bbI-1].second) && 
				nway.test(bb_segs[bbI+1].second) &&
				absolut(bb_segs[bbI-1].first[seqs[sI]].second - bb_segs[bbI-1].first[seqs[sI]].first) >= indie_dist &&
				absolut(bb_segs[bbI+1].first[seqs[sI]].second - bb_segs[bbI+1].first[seqs[sI]].first) >= indie_dist )
			{
			}else
				good_bb.set(bb_segs[bbI].second, false);
		}
	}

	// create site patterns, then write out the good ones
	bitset_t empty( bb_seq_list.size() );
	vector< bitset_t > spa_seqs( seqs.size(), empty );
	for( size_t bbI = 0; bbI < bb_seq_list.size(); bbI++ )
		for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
			spa_seqs[seqI].set(bbI, bb_seq_list[bbI][seqs[seqI]].first != 0);

	vector< string > binseqs( seqs.size(), string( good_bb.count(), '0' ) );
	for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
	{
		size_t goodI = 0;
		for( size_t bbI = 0; bbI < good_bb.size(); bbI++ )
			if(good_bb.test(bbI))
			{
				if(spa_seqs[seqI].test(bbI))
					binseqs[seqI][goodI] = '1';
				goodI++;
			}
	}
	// write out the seqs!!
	if( target_format == "beast" )
	{
		anal_output << "\t<taxa id=\"taxa\">\n";
		for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
		{
			anal_output << "\t\t<taxon id=\"seq" << seqI << "\"/>\n";
	
		}
		anal_output << "\t</taxa>\n";
		anal_output << "\t<alignment id=\"alignment\" dataType=\"binary\">\n";
	
		for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
		{
			anal_output << "\t\t<sequence>\n";
			anal_output << "\t\t\t<taxon idref=\"seq" << seqI << "\"/>\n";
			anal_output << "\t\t\t" << binseqs[seqI] << endl;
			anal_output << "\t\t</sequence>\n";
//			anal_output << "> seq" << seqI << endl;
//			for( size_t i = 0; i < binseqs[seqI].size(); i+=80 )
//				anal_output << binseqs[seqI].substr(i, 80) << endl;
		}
		anal_output << "\t</alignment>\n";
	}else{
		// write genoplast format
		for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
		{
			for( size_t cI = 0; cI < binseqs[seqI].size(); cI++ )
			{
				if( cI > 0 )
					anal_output << ' ';
				anal_output << binseqs[seqI][cI];
			}
			anal_output << endl;
		}
	}

	anal_output.close();
}

