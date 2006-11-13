// cadena.cpp : main project file.
/*******************************************************************************
 * $Id: cadena.cpp,v 1.49 2005/11/15 00:18:45 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#include "mauveAligner.h"
#include "getopt.h"
#include "cadena.h"
#include <sstream>
#include <stdexcept>
#include "libMems/Matrix.h"
#include "libMems/NumericMatrix.h"
#include "libGenome/gnSequence.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MemHash.h"
#include "libMems/MaskedMemHash.h"
#include "libMems/Aligner.h"
//#include "libMems/MatchList.h"
#include "libMems/MemSubsets.h"
#include "libMems/RepeatHashCat.h"
#include "libMems/RepeatMatchList.h"
#include "libMems/Interval.h"
#include "libMems/IntervalList.h"
#include "libMems/gnAlignedSequences.h"
#include "libMems/Islands.h"
#include "libMems/MuscleInterface.h"

using namespace std;
using namespace genome;

Cadena::Cadena( unsigned int seq_count ) :
match_allocator( SlotAllocator< Match >::GetSlotAllocator() )
{
	debug = false;
	this->seq_count = seq_count;
	collinear_genomes = false;
	this->permutation_filename = "";
	this->permutation_weight = -1;
}

/*Cadena::Cadena( const Cadena& al )  : 
match_allocator( al.match_allocator )
{
	operator=( al );
}
*/
Cadena& Cadena::operator=( const Cadena& al )
{

	return *this;
}

void print_usage( const char* pname ){
	cerr << "Cadena usage:" << endl;
	cerr << pname << " [options] <seq1 filename> <sml filename> <seq2 filename> ... <seqN filename> " << endl;
	cerr << "Options:" << endl;
	cerr << "\t    --repeat-size=<number> Min repeat size"  << endl;
	cerr << "\t    --input-alignment=<file> Read in an XMFA format alignment from the specified file" << endl;
	cerr << "\t    --output=<file> Output file name containing chained repeats.  Prints to screen by default" << endl;
	cerr << "\t    --input-tree=<file> Input guide tree file name." << endl;
	cerr << endl;
}
int main(int argc, char* argv[] )
{


	//
	// parse command line with gnu getopt
	//
	int opt;
	int config_opt;
	int ac = argc;
	char** av = argv;


	vector<string> seq_files;
	vector<string> sml_files;
	uint repeat_size = 15;	// Use default settings
	string output_fn = "output.txt";
	bool read_matches = false;	
	string input_aln_fn = "";
	string tree_fn = "";
	
	//MatchList match_list;
	//MatchList match_list_final;
	RepeatMatchList match_list;
	RepeatMatchList match_list_final;
	
	// take as parameters
	// -seq files, sml files, mer_size, distance matrix, 
	// output file(list of chains of segmental homologies found), level recursion
	// this output file can then be passed to Mauve/M-GCAT & incorporated into Display & Alignment

	// 'm' mer size
	// 'r' recursive
	const char* short_args= "r:o:t:i:vrRE";
	enum opt_names{
		opt_repeat_size,
		opt_output,
		opt_input_tree,
		opt_input_aln,
	};
	struct option long_opts[] = {
		{"repeat-size", required_argument, &config_opt, opt_repeat_size},
		{"output", required_argument, &config_opt, opt_output},
		{"input-alignment", required_argument, &config_opt, opt_input_aln},
		{"input-tree", required_argument, &config_opt, opt_input_tree},
		{0, 0, 0, 0}	// for correct termination of option list
						// getopt_long can segfault without this
	};

	int indexptr;
	while( (opt = getopt_long( ac, av, short_args, long_opts, &indexptr )) != EOF ){
		switch( opt ){
			case 0:
				switch(config_opt){
					case opt_repeat_size:
						repeat_size = atoi( optarg );
						break;
					case opt_output:
						output_fn = optarg;
						break;
					case opt_input_aln:
						input_aln_fn = optarg;
						break;
					case opt_input_tree:
						tree_fn = optarg;
						break;
					default:
						print_usage( argv[0] );
						return -1;
					
				}
				break;
			default:
				print_usage( argv[0] );
				return -1;
		}
	}

	// now read in the seq and sml file names from av
	boolean seq_name_arg = true;
	uint seqcount = 0;
	for( int optI = optind; optI < argc; optI++ )
	{
		seqcount++;
		
		if( seq_name_arg )
			seq_files.push_back( av[ optI ] );
		else
			sml_files.push_back( av[ optI ] );
		if ( seqcount == 1 )
			seq_name_arg = false;
		else
		{
			seq_name_arg = true;
			//sml_files.push_back( "tmp.sml" );
		}
	}
	if (! seqcount )
	{
		print_usage( argv[0] );
		return -1;
	}
	
	//use to find repeats in a single sequence
	//or to find segmental homology in several sequences concatenated into a single sequence
	//pass 1 or more sequences to program, if > 1 concatenate into 1
	//track borders, find matches, detain border jumpers, rinse & repeat

	
	//MLDeleter deleter( match_list );
	
	
	

	match_list.seq_filename = seq_files;
	match_list.sml_filename = sml_files;
	match_list.LoadSequences( &cout );
	
	uint contig_count = 0;
	uint32 seq_start = 0;
	//1) concatenate sequences into 1 sequence, set boundaries
	gnSequence* file_sequence = match_list.seq_table.at(0);
	contig_count = file_sequence->contigListSize();
	RepeatHashCat repeat_finder;
	for(uint i = 1; i < match_list.seq_table.size(); i++)
	{
		repeat_finder.concat_contig_start.push_back( seq_start );
		file_sequence->append((const gnSequence &)*match_list.seq_table.at(i));
		seq_start = ((const gnSequence &)*match_list.seq_table.at(i)).length();
	}
	match_list.seq_table.clear();
	match_list.seq_table.push_back(file_sequence);
	match_list_final = match_list;
	cout << contig_count << " " << file_sequence->contigListSize() << endl;
	match_list.LoadSMLs( repeat_size, &cout );

	
	
	//repeat_finder.concat_contig_start.push_back( contig_count );
	repeat_finder.LogProgress( &cout );
 
	//2) find repeats in concatenated sequence, divide matches that span borders
	cout << "\nFinding matches & deporting border jumpers... " << endl;
	repeat_finder.FindMatches( match_list );

	vector<Match*>::const_iterator match_iter;
	match_iter = match_list.begin();
	set<Match*> cur_set;
	set<Match*>::iterator set_iter;

	cout << "\nSelecting MEMs with the highest multiplicity and Cropping overlaps..." << endl;
	for(; match_iter != match_list.end(); match_iter++)
	{
		(*match_iter)->UnlinkFromSets();
		cur_set = (*match_iter)->Supersets();

		//if no supersets present, this is the highest level multiplicity
		if ( cur_set.size() == 0 )
		{
			//if overlapping, trim
			for( uint i = 0; i < (*match_iter)->SeqCount()-1; i++)
			{
				if ( absolut( (*match_iter)->Start(i) ) + (*match_iter)->Length() > absolut((*match_iter)->Start(i+1)) )
				{
					//overlapping, crop
					(*match_iter)->CropEnd( absolut( (*match_iter)->Start(i) ) + (*match_iter)->Length() - absolut((*match_iter)->Start(i+1))  ); 
				}

			}
			match_list_final.push_back(*match_iter);
		}
		//match_file << cur_set.size();
		//set_iter = cur_set.begin();
		//for(; set_iter != cur_set.end(); set_iter++ )
		//{
		//}
	}
	
	ostream* match_out;
	ofstream* match_out_file = new ofstream( output_fn.c_str()  );
	match_out = match_out_file;
	//match_list.WriteList( *match_out );
	match_list_final.WriteList( *match_out );
	match_out->flush();

	//extract_tr(match_list);
	//3) squeeze out the lowest subset in all match sets, remove overlaps in this lowest level set
	//track as Tandem Repeats for now, later change to Approximate Tandem Repeats(ATR)

    //4) chain set of repeat matches into monotonically increasing/decreasing chain
	//la_cadena = encadenar(match_list);

	//5) Output results to file

	return 0;
}
