#include "getopt.h"
#include "gn/gnSequence.h"
#include "DNAFileSML.h"
#include "MemHash.h"
#include "Aligner.h"
#include "MatchList.h"
#include "MemSubsets.h"

void writeRearrangmentList( const vector<MemList>& LCB_list, ostream& lcb_out );
void writeRearrangmentList( const vector<MemList>& LCB_list, ostream& lcb_out ){
	if( LCB_list.size() == 0 )
		return;
	
	const list<MemHashEntry*>& mem_list = LCB_list[0].mem_list;
	
	MemHashEntry* first_mem = *(mem_list.begin());
	unsigned int seq_count = first_mem->SeqCount();
	unsigned int mer_size = first_mem->MerSize();
	unsigned int seqI;
	lcb_out << "SequenceCount" << '\t' << seq_count << "\n";
	for(seqI = 0; seqI < seq_count; seqI++){
		if( LCB_list[0].seq_filename.size() > seqI ){
			lcb_out << "Sequence" << seqI << "File" << '\t';
			lcb_out << LCB_list[0].seq_filename[seqI];
			lcb_out << "\n";
		}
	}

	for( uint lcbI = 0; lcbI < LCB_list.size(); lcbI++ ){
		//get all the mems out of the lists and write them out
	    list<MemHashEntry*>::const_iterator mem_iter;
		mem_iter = LCB_list[ lcbI ].mem_list.begin();
		lcb_out << lcbI << " Start:";
		for( seqI = 0; seqI < seq_count; seqI++ ){
			int64 start_coord = (*mem_iter)->Start( seqI );
			if( start_coord < 0 )
				start_coord -= (*mem_iter)->Length();
			lcb_out << '\t' << start_coord;
		}

		mem_iter = LCB_list[ lcbI ].mem_list.end();
		mem_iter--;
		lcb_out << lcbI << "\tEnd:";
		for( seqI = 0; seqI < seq_count; seqI++ ){
			int64 end_coord = (*mem_iter)->Start( seqI );
			if( end_coord > 0 )
				end_coord += (*mem_iter)->Length();
			lcb_out << '\t' << end_coord;
		}
		lcb_out << endl;
	}
}


void print_usage( const char* pname ){
	cerr << pname << " [-i input file] [-l] [-m mer size] "
		<< " [-s min match size] [-f max offset diff] [-d min identity] "
		<< " [-r min range] [-o output file] "
		<< " [seq1 filename] [sml1 filename] ... [seqN filename] [smlN filename]" << endl;
	cerr << "Options:" << endl;
	cerr << "\t-m, --seed-size=<number> Initial seed match size, default is " << DNA_MER_SIZE << endl;
	cerr << "\t-o, --output=<file> Output file name (Outputs to screen if not specified)" << endl;
	cerr << "\t-l, --locate-LCBs Locate locally collinear blocks (LCBs)" << endl;
	cerr << "\t-s, --min-match-size=<number> Minimum length of matches in b.p. to use for LCB determination (Implies -l)" << endl;
	cerr << "\t-f, --max-offset-diff=<number> Maximum permissible difference in generalized offset between adjacent matches (Implies -l)" << endl;
	cerr << "\t-d, --min-identity=<number> Minimum LCB identity where <number> is a real value between 0 and 1, default is 0 (Implies -l)" << endl;
	cerr << "\t-r, --min-range=<number> Minimum LCB range in base pairs (Implies -l)" << endl;
	cerr << "\t-i, --match-input=<file> Use specified match file instead of searching for matches\n";
	cerr << endl;
}

/**
 * This basic application demonstrates how to use libMems to produce full scale multiple
 * genomic alignments.  First the command line is parsed to get the names of data files
 * and user specified options.  Next each sequence and its corresponding sorted mer list
 * are loaded.  If the sorted mer list fails to load a new one is created.  
 * If it is necessary to find matches in the sequences instead of loading them, each 
 * sequence and SML are added to a MemHash which searches for exact matches.  
 * Then LCBs are found if the user requested it.  Finally, either the MatchList or the
 * LCB list is written to disk.
 */
int main( int argc, char* argv[] ){
	try{
	if( argc <= 0 ){
		print_usage( "mauveAligner" );
		return -1;
	}
	if( argc == 1 ){
		print_usage( argv[0] );
		return -1;
	}
	
//	string debuggo;
//	cin >> debuggo;
	
	vector<string> seq_files;
	vector<string> sml_files;
	vector<gnSequence*> seq_table;
	vector<DNAFileSML*> sml_table;
	uint mer_size = DNA_MER_SIZE;	// Default
	boolean create_LCBs = false;
	double LCB_density = 0;
	int64 LCB_size = -1;
	int64 min_match_size = 0;
	int64 max_offset_diff = -1;
	string output_file = "";
	boolean read_matches = false;
	boolean create_matches = false;
	string match_input_file = "";

	// parse command line with gnu getopt
	int opt;
	int ac = argc;
	char** av = argv;

	const char* short_args= "m:d:r:s:o:i:f:l";
	struct option long_opts[] = {
		{"seed-size", required_argument, NULL, 'm'},
		{"locate-LCBs", no_argument, NULL, 'l'},
		{"max-offset-diff", required_argument, NULL, 'f'},
		{"min-match-size", required_argument, NULL, 's'},
		{"min-identity", required_argument, NULL, 'd'},
		{"min-range", required_argument, NULL, 'r'},
		{"output", required_argument, NULL, 'o'},
		{"match-input", required_argument, NULL, 'i'},
	};
	int indexptr;
	while( (opt = getopt_long( ac, av, short_args, long_opts, &indexptr )) != EOF ){
		switch( opt ){
			case 'm':
				mer_size = atoi( optarg );
				break;
			case 'l':
				create_LCBs = true;
				break;
			case 'f':
				create_LCBs = true;
				max_offset_diff = atol( optarg );
				break;
			case 's':
				create_LCBs = true;
				min_match_size = atol( optarg );
				break;
			case 'd':
				create_LCBs = true;
				LCB_density = atof( optarg );
				break;
			case 'r':
				create_LCBs = true;
				LCB_size = atol( optarg );
				break;
			case 'o':
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

	// check for incorrect invocation syntax
	if( create_LCBs ){
		if( LCB_size < 0 ) {
			cerr << "A minimum LCB range greater than 0 must be specified in order to create LCBs.\n";
			return -1;
		}
		if( max_offset_diff < -1 ){
			cerr << "The maximum offset difference for adjacent matches should be non-negative.\n";
			return -1;
		}
	}
	
	
	if( seq_files.size() != sml_files.size() ){
		cerr << "Error: Each sequence file must have a corresponding SML file specified.\n";
		return -1;
	}
	// done parsing and checking command line options
	// Start doing the work
	MatchList match_list;
	match_list.seq_filename = seq_files;
	match_list.sml_filename = sml_files;
	if( !read_matches )
		match_list.LoadSequences( mer_size, &cout );
	MemList mem_list;
	mem_list.seq_filename = seq_files;
	mem_list.sml_filename = sml_files;
	mem_list.seq_table = match_list.seq_table;
	mem_list.sml_table = match_list.sml_table;
	
	ofstream match_out;
	if( output_file != "-" && output_file != "" ){
		match_out.open( output_file.c_str() );
		if( !match_out.is_open() ){
			cerr << "Error opening " << output_file << endl;
			return -2;
		}
	}

	if( read_matches ){
		ifstream match_in( match_input_file.c_str() );
		if( !match_in.is_open() ){
			cerr << "Error opening " << match_input_file << endl;
			return -2;
		}
		mem_list.ReadList( match_in );
		mem_list.seq_filename = seq_files;
	}else{
		MemHash match_finder;
		match_finder.LogProgress( &cout );
		match_finder.FindMems( mem_list );
	}

	if( !create_LCBs ){
		mem_list.MultiplicityFilter( mem_list.seq_table.size() );
		mem_list.WriteList( match_out );
		return 0;
	}

	Aligner aligner( max_offset_diff, LCB_density, LCB_size );
	vector<MemList> lcb_list;
	try{
		aligner.align( mem_list, lcb_list, min_match_size );
	}catch( gnException& gne ){
		cerr << gne << endl;
	}

	if( output_file != "-" && output_file != "" )
		writeRearrangmentList( lcb_list, match_out );
	else
		writeRearrangmentList( lcb_list, cout );

	}catch( ... ){
		cerr << "Unhandled exception, bailing out.\n";
	}
	return 0;
}
