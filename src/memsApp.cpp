#include "getopt.h"
#include "gn/gnSequence.h"
#include "DNAFileSML.h"
#include "MemHash.h"
#include "aligner.h"
#include "MatchList.h"
#include "MemSubsets.h"

class lcsMatch{
public:
	gnSeqI baseI;
	Match* match;

	lcsMatch* back_ptr;
	gnSeqI path_length;
};

class LisComparator {
public:
	LisComparator( unsigned seq = 0, boolean rest_decreasing = false, boolean full_comparison = false){
		m_seq = seq;
		this->rest_decreasing = rest_decreasing;
		this->full_comparison = full_comparison;
	}

	LisComparator( const LisComparator& msc ){
		*this = msc;
	}
	
	LisComparator& operator=( const LisComparator& msc ){
		m_seq = msc.m_seq;
		rest_decreasing = msc.rest_decreasing;
		full_comparison = msc.full_comparison;
		return *this;
	}

	bool operator()(const lcsMatch* a, const lcsMatch* b) const{
		const lcsMatch* a_match = a;
		const lcsMatch* b_match = b;		
		
//		int32 start_diff = max( a_match->match->FirstStart(), m_seq ) - max( b_match->match->FirstStart(), m_seq );
//		if(start_diff == 0){
			uint32 m_count = a_match->match->SeqCount();
//			m_count = m_count <= b_match->match->SeqCount() ? m_count : b_match->match->SeqCount();
			for(uint32 seqI = m_seq; seqI < m_count; seqI++){
				int64 a_start = a_match->match->Start( seqI ), b_start = b_match->match->Start( seqI );
				
				if( a_start == MEM_NO_MATCH || b_start == MEM_NO_MATCH ){
					if( rest_decreasing ){
						// set a and b so the comparison goes reverse next time thru
						a_match = b;
						b_match = a;
					}
					continue;
				}

				a_start += a_match->baseI;
				b_start += b_match->baseI;
				
				if( rest_decreasing ){
					// set a and b so the comparison goes reverse next time thru
					a_match = b;
					b_match = a;
				}
				
//				if(a_start < 0)
//					a_start = -a_start + a_match->match->Length();
//				if(b_start < 0)
//					b_start = -b_start + b_match->match->Length();

				int64 diff = a_start - b_start;
				if( diff == 0 )
					continue;
				else if( !full_comparison )
					return diff < 0;
				else if( diff >= 0 ){
					return false;
				}else if( seqI + 1 == m_count )
					return true;
			}
//		}
		return false;
//		return start_diff < 0;
	}
	
	/**
	 * Ensures that both matches have the same directional parity (forward or reverse complement )
	 */
	bool checkParity( const lcsMatch* a, const lcsMatch* b ){
		int32 start_diff = max( a->match->FirstStart(), m_seq ) - max( b->match->FirstStart(), m_seq );
		if(start_diff == 0){
			uint32 m_count = a->match->SeqCount();
			if( m_count != b->match->SeqCount() )
				return false;
			for(uint32 seqI = 0; seqI < m_count; seqI++){
				if( a->match->Start( seqI ) == MEM_NO_MATCH ||
					b->match->Start( seqI ) == MEM_NO_MATCH )
					continue;
				if( a->match->Start( seqI ) > 0 &&
					b->match->Start( seqI ) > 0 )
					continue;
				if( a->match->Start( seqI ) < 0 &&
					b->match->Start( seqI ) < 0 )
					continue;
				return false;
			}
			return true;
		}
		return false;
	}

private:
	unsigned m_seq;
	boolean rest_decreasing;
	boolean full_comparison;
};


void lcsAlign( MatchList& ml, vector<MatchList>& lcb_list ){
	if( ml.match_list.size() == 0 )
		return;

// Eliminate linked inclusions for now, come up with a better solution later
	SubsetInclusionRemover ms;
//	ms.EliminateLinkedInclusions( ml.match_list );

//
// separate all unlinked subsets without a match in the first genome
// from the rest of the matches
//
//	map<MatchID_t, Match*> match_map;
//	uint seq_count = (*match_iter)->SeqCount();
//	for(; match_iter != ml.match_list.end(); match_iter++ ){
//		match_map.insert( map<MatchID_t, Match*>::value_type( (*match_iter)->MatchID(), *match_iter ) );
//	}

	list<Match*>::iterator match_iter = ml.match_list.begin();
	MatchList filtered_list = ml;
	filtered_list.UnlinkedFirstStartFilter( 0 );

// sum up the total number of lcsMatches that will be needed and allocate them all at once
	match_iter = filtered_list.match_list.begin();
	gnSeqI total_matches = 0;
	for(; match_iter != filtered_list.match_list.end(); match_iter++ ){
		total_matches += (*match_iter)->Length();
	}
	
	Array<lcsMatch> match_buffer( total_matches );

	list< lcsMatch* > lcs_matches;
	match_iter = filtered_list.match_list.begin();
	gnSeqI matchI = 0;
	for(; match_iter != filtered_list.match_list.end(); match_iter++ ){
		for( gnSeqI baseI = 0; baseI < (*match_iter)->Length(); baseI++ ){
			lcsMatch* lcs_match = &( match_buffer.data[matchI] ); 
			lcs_match->baseI = baseI;
			lcs_match->match = *match_iter;
			lcs_match->back_ptr = NULL;
			lcs_match->path_length = (*match_iter)->Multiplicity();
			lcs_matches.push_back( lcs_match );
			matchI++;
		}
	}

//
// Do LIS on these guys
//

	// sort them in decreasing order on first sequence, increasing in every other seq
	LisComparator lsc( 0, true, false );
	lcs_matches.sort( lsc );
	
	lsc = LisComparator( 1, false, false );

	// add them to decreasing lists
	// vector of lists?
	vector< list<lcsMatch*> > lis_matches;
	map< lcsMatch*, gnSeqI, LisComparator > lis_map( lsc );

	LisComparator back_comp( 1, false, true );
	gnSeqI max_col = 0, max_entry = 0, max_len = 0;
	
	list< lcsMatch* >::iterator lcs_iter = lcs_matches.begin();
	for( ; lcs_iter != lcs_matches.end(); lcs_iter++ ){
//		if( (*lcs_iter)->match->Start( 1 ) == 491524 ){
//			cout << "shite goin down\n";
//		}
		map< lcsMatch*, gnSeqI, LisComparator >::iterator insert_loc = lis_map.lower_bound( *lcs_iter );
		gnSeqI column = 0;
		if( insert_loc != lis_map.end() ){
			// add to the end of the column
			column = insert_loc->second;
			lis_map.erase( insert_loc );
			lis_map.insert( map< lcsMatch*, gnSeqI >::value_type( *lcs_iter, column ) );
			lis_matches[ column ].push_back( *lcs_iter );
		}else{
			// start a new column
			list<lcsMatch*> new_list;
			column = lis_matches.size();
			lis_matches.push_back( new_list );
			lis_matches[ column ].push_back( *lcs_iter );
			lis_map.insert( map< lcsMatch*, gnSeqI >::value_type( *lcs_iter, column ) );
		}

		// search for something to point back to
		for( gnSeqI colI = column; colI > 0; colI-- ){
			list< lcsMatch* >::iterator back_iter = lis_matches[ colI - 1 ].end();
			do{
				back_iter--;
				// this must be smaller in _every_ sequence
				bool bc = back_comp( *back_iter, *lcs_iter );
				bool parity = back_comp.checkParity( *back_iter, *lcs_iter );
				if( bc && parity ){
					(*lcs_iter)->back_ptr = *back_iter;
					(*lcs_iter)->path_length += (*back_iter)->path_length;
					if( (*lcs_iter)->path_length > max_len ){
						max_len = (*lcs_iter)->path_length;
						max_entry = lis_matches[ column ].size() - 1;
						max_col = column;
					}
					break;
				}
			}while( back_iter != lis_matches[ colI - 1  ].begin() );

			if( (*lcs_iter)->back_ptr != NULL )
				break;
		}
		
	}
	
	// lcs is found.  do something.
	// aggregate all linked lcsMatches that point to the same Match,  if the match isn't complete
	// then break it up into pieces.
	// find the next lcs with the remaining pieces.  
	
	lcs_iter = lis_matches[ max_col ].begin();
	for( uint matchI = 0; matchI < max_entry; matchI++ ){
		lcs_iter++;
	}
	
	// lcs_iter points to the end of the longest common subsequence.  trace it backwards to
	// get the complete lcs
	lcsMatch* prev_lcs = *lcs_iter;
	lcsMatch* cur_lcs = prev_lcs->back_ptr;
	Match* trimmed_match = prev_lcs->match->Clone();
	trimmed_match->CropEnd( trimmed_match->Length() - prev_lcs->baseI - 1 );
	list<Match*> lcs_match_list;
	while( cur_lcs != NULL ){
		if( cur_lcs->match != prev_lcs->match ){
			trimmed_match->CropStart( prev_lcs->baseI );
			lcs_match_list.push_back( trimmed_match );
			trimmed_match = cur_lcs->match->Clone();
			trimmed_match->CropEnd( trimmed_match->Length() - cur_lcs->baseI - 1 );
		}
		prev_lcs = cur_lcs;
		cur_lcs = prev_lcs->back_ptr;
		
	}
	trimmed_match->CropStart( prev_lcs->baseI );
	lcs_match_list.push_back( trimmed_match );
	
	lcb_list.push_back( ml );
	lcb_list[0].match_list.clear();
	lcb_list[0].match_list.insert( lcb_list[0].match_list.begin(), lcs_match_list.rbegin(), lcs_match_list.rend() );

	// cover as much range as possible subject to minimum density criteria
	// lcbs compete for area, the lcb with the better diagonal distance ratio wins.
	
}



void print_usage( const char* pname ){
	cerr << pname << " [-r] [-m mer size] [-o diff] [-d density] [-s size] "
		<< "<-f output file> <seq1 filename> <sml1 filename> ... "
		<< " <seqN filename> <smlN filename>" << endl;
	cerr << "Options:" << endl;
	cerr << "\t-l Create locally collinear blocks (LCBs)" << endl;
	cerr << "\t-r Recursive alignment of gaps (Implies -l)" << endl;
	cerr << "\t-m <number> Initial mer size" << endl;
	cerr << "\t-o <number> Maximum LCB offset difference between neighboring matches" << endl;
	cerr << "\t-d <number> Minimum LCB density" << endl;
	cerr << "\t-s <number> Minimum LCB size in base pairs" << endl;
	cerr << "\t-f <file> Output file name" << endl;
	cerr << "\t-i <file> Use specified input match file instead of searching for matches\n";
	cerr << endl;
}

/**
 * This basic application demonstrates how to find MUMs in two sequences using libMems.
 * First, each sequence and its corresponding sorted mer list is loaded.  If the sorted
 * mer list fails to load a new one is created.  Then each sequence and SML are added to
 * a MemHash which searches for exact matches.  Finally, the MemList is extracted from the
 * MemHash and is written to disk.
 */
int main( int argc, char* argv[] ){
	if( argc <= 0 ){
		print_usage( "mauveAligner" );
		return -1;
	}
	if( argc == 1 ){
		print_usage( argv[0] );
		return -1;
	}
	
	vector<string> seq_files;
	vector<string> sml_files;
	vector<gnSequence*> seq_table;
	vector<DNAFileSML*> sml_table;
	uint mer_size = DNA_MER_SIZE;	// Default
	boolean recursive = false;
	boolean create_LCBs = false;
	gnSeqI LCB_offset = 0;
	double LCB_density = 0;
	gnSeqI LCB_size = 0;
	string output_file = "";
	boolean read_matches = false;
	string match_input_file = "";
	
	// parse command line with gnu getopt
	int opt;
	int ac = argc;
	char** av = argv;
	// 'm' mer size
	// 'r' recursive
	const char* short_args= "m:o:d:s:f:i:rl";

	while( (opt = getopt( ac, av, short_args )) != EOF ){
		switch( opt ){
			case 'l':
				create_LCBs = true;
				break;
			case 'r':
				create_LCBs = true;
				recursive = true;
				break;
			case 'm':
				mer_size = atoi( optarg );
				break;
			case 'o':
				LCB_offset = atol( optarg );
				break;
			case 'd':
				LCB_density = atof( optarg );
				break;
			case 's':
				LCB_size = atol( optarg );
				break;
			case 'f':
				output_file = optarg;
				break;
			case 'i':
				read_matches = true;
				match_input_file = optarg;
				break;
			default:
				print_usage( argv[0] );
				return -1;
		}
	}
	// now read in the seq and sml file names from av
	boolean seq_name_arg = true;
	for( int optI = optind; optI < argc; optI++ ){
		if( seq_name_arg )
			seq_files.push_back( av[ optI ] );
		else
			sml_files.push_back( av[ optI ] );
		seq_name_arg = !seq_name_arg;
	}
	string dump;
	cout << "in debug mode\n";
	cin >> dump;
//	if( seq_files.size() > 2 ){
//		cout << "Cannot align more than two sequences\n";
//		return -1;
//	}
//	recursive = false;
	
	// check for incorrect invocation syntax
	if( create_LCBs && ( ( LCB_offset == 0 ) ||
		( LCB_density == 0 ) || ( LCB_size == 0 ) ) ){
		
		cerr << "The maximum LCB offset, minimum LCB density, and minimum LCB size must all be specified.\n";
		return -1;
	}
	
	if( output_file == "" ){
		cerr << "An output file must be specified!\n";
		return -1;
	}
	
	if( seq_files.size() != sml_files.size() ){
		cerr << "Error: Each sequence file must have a corresponding SML file specified.\n";
	}
	
	// done parsing and checking command line options
	// Start doing the work
	
	MatchList match_list;

	match_list.seq_filename = seq_files;
	match_list.sml_filename = sml_files;
	match_list.LoadSequences( mer_size, &cout );
	
	ofstream match_out( output_file.c_str() );
	if( !match_out.is_open() ){
		cerr << "Error opening " << output_file << endl;
		return -2;
	}

	if( read_matches ){
		ifstream match_in( match_input_file.c_str() );
		if( !match_in.is_open() ){
			cerr << "Error opening " << match_input_file << endl;
			return -2;
		}
		match_list.ReadList( match_in );
		match_list.seq_filename = seq_files;
	}else{
		MemHash match_finder;
		match_finder.LogProgress( &cout );
		match_finder.FindMatches( match_list );
	}

	if( !create_LCBs ){
		match_list.WriteList( match_out );
		return 0;
	}

	Aligner aligner( LCB_offset, LCB_density, LCB_size );
	vector<MatchList> lcb_list;
	try{
		aligner.align( match_list, lcb_list, recursive );
//		lcsAlign( match_list, lcb_list );
	}catch( gnException& gne ){
		cerr << gne << endl;
	}
//	catch( ... ){
//		cerr << "unhandled exception\n";
//	}
	aligner.writeLCBlist( lcb_list, match_out );
	return 0;
}
