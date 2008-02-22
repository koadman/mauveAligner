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


int main( int argc, char* argv[] )
{
#if	WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	if( argc < 4 )
	{
		cerr << "bbFilter <backbone file> <independent dist> <output file> <seq1> <seq2>...<seqN>\n";
		cerr << "seq index starts at 0.\n";
		return -1;
	}
	string bbseq_fname( argv[1] );
	int indie_dist = atoi( argv[2] );
	string output_fname( argv[3] );

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
	for( int i = 4; i < argc; i++ )
		seqs.push_back(atoi(argv[i]));

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
//		anal_output << "> seq" << seqI << endl;
//		for( size_t i = 0; i < binseqs[seqI].size(); i+=80 )
//			anal_output << binseqs[seqI].substr(i, 80) << endl;
	}
	anal_output << "\t</alignment>\n";
	anal_output.close();
}