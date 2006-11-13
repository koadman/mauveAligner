/*******************************************************************************
 * $Id: Aligner.cpp,v 1.47 2004/04/19 23:10:30 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi.h"

// For compilers that support precompilation, includes "wx/wx.h".

#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWindows headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif


#include "libMems/Aligner.h"
#include "libMems/MemSubsets.h"
#include "libMems/Islands.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MuscleInterface.h"	// it's the default gapped aligner
#include "libGenome/gnRAWSource.h"

#include <map>
#include <fstream>	// for debugging
#include <sstream>
#include <stack>
#include <algorithm>



bool debug_shite = false;


MpiAligner::MpiAligner( uint seq_count ) :
match_allocator( SlotAllocator< Match >::GetSlotAllocator() )
{
	debug = false;
	this->seq_count = seq_count;
	min_recursive_gap_length = default_min_r_gap_size;
	collinear_genomes = false;
	gal = &(MuscleInterface::getMuscleInterface());
	this->permutation_filename = "";
	this->permutation_weight = -1;
}

MpiAligner::MpiAligner( const MpiAligner& al ) : 
match_allocator( al.match_allocator )
{
	operator=( al );
}

MpiAligner& MpiAligner::operator=( const MpiAligner& al )
{
	ml = al.ml;
	gap_mh = al.gap_mh;
	nway_mh = al.nway_mh;
	seq_count = al.seq_count;
	debug = al.debug;
	
	LCB_minimum_density = al.LCB_minimum_density;
	LCB_minimum_range = al.LCB_minimum_range;
	
	cur_min_coverage = al.cur_min_coverage;
	min_recursive_gap_length = al.min_recursive_gap_length;
	collinear_genomes = al.collinear_genomes;

	gal = al.gal;

	permutation_weight = al.permutation_weight;
	permutation_filename = al.permutation_filename;

	return *this;
}

// compute the gapped alignments between anchors in an LCB
// master waits for either (1) request for work from worker
// or (2) result message
// if (1) master sends the index of the intervening gap to align
// if (2) master receives results
// 
// worker sends request for work
// worker does an alignment on the specified gap
// worker sends results message
// 
void MpiAligner::AlignLCB( MatchList& mlist, Interval& iv ){
	// check whether this function can do anything useful...
	if( !collinear_genomes && mlist.size() < 2 ){
		iv.matches.insert( iv.matches.end(), mlist.begin(), mlist.end() );
		iv.CalculateOffset();
		return;
	}

	boolean debug_recurse = false;
	int64 config_value = 138500;
	int print_interval = 50;
	try{
	list< Match* > match_list;
	match_list.insert( match_list.end(), mlist.begin(), mlist.end() );
	mlist.clear();
	MatchList r_list = mlist;

	list< Match* >::iterator recurse_iter = match_list.begin();
	list< Match* >::iterator recurse_prev = match_list.begin();
	// scan ahead to the first n-way matches
	while( recurse_prev != match_list.end() && (*recurse_prev)->Multiplicity() != seq_count )
		recurse_prev++;

	recurse_iter = recurse_prev;
	if( !collinear_genomes ){
		if( recurse_iter != match_list.end() )
			recurse_iter++;
		while( recurse_iter != match_list.end() && (*recurse_iter)->Multiplicity() != seq_count )
			recurse_iter++;
	}else
		cout << "Assuming collinear genomes...\n";

	int rank = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Barrier(MPI_COMM_WORLD);

	if( rank == 0 ){
		// master role
		// 1.  count the number of gapped alignments to perform

		int anchor_count = 0;
		while( true ){
			anchor_count++;

			// scan ahead to the next pair of n-way matches
			recurse_prev = recurse_iter;
			if( recurse_iter != match_list.end() )
				recurse_iter++;
			while( recurse_iter != match_list.end() && (*recurse_iter)->Multiplicity() != seq_count )
				recurse_iter++;

			if( ( recurse_iter == match_list.end() && !collinear_genomes ) ||
					( recurse_prev == match_list.end() && collinear_genomes ) )
					break;
		}
		int current_alignment = 0;
		//
		// 2.  wait for a message from workers
		//
		MPI_Status probe_status;
		MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_status );
		if( probe_status.MPI_TAG == WORKER_IDLE_TAG ){
			// assign the next alignment
			MPI_Send( &current_alignment, 1, MPI_INT, probe_status.MPI_SOURCE,
					ASSIGNMENT_TAG, MPI_COMM_WORLD );
			current_alignment++;
		}else{
			// receive results from the worker

		}
	}else{

	}

	uint memI = 0;
	uint matchI = 0;
	while( true ){
		if( memI >= print_interval && memI % print_interval == 0 || debug)
			cout << "Number: " << memI << " match " << **recurse_prev << endl;
		memI++;
		if( debug_recurse ){
			cout << "Recursing on " << endl;
			if( recurse_prev != match_list.end() )
				cout << **recurse_prev << " and " << endl;
			if( recurse_iter != match_list.end() )
				cout << **recurse_iter << endl;
		}
		
		if( recurse_prev != match_list.end() && (*recurse_prev)->Start( 0 ) == config_value )
			cout << "";
		
		// recurse on a pair of matches! 
		// this function should locate all matches between the two iterators
		// and add them to r_list		
		r_list.clear();
		GappedAlignment* cr = NULL;
		boolean align_success = false;
		
		Match* r_lend = NULL;
		Match* r_rend = NULL;
		if( recurse_iter != recurse_prev )
			r_lend = *recurse_prev;
		if( recurse_iter != match_list.end() )
			r_rend = *recurse_iter;

		// attempt a clustalW alignment
		cr = new GappedAlignment();
		align_success = gal->Align( *cr, r_lend, r_rend, r_list.seq_table );

		// add the gapped alignment to the Interval
		if( r_lend != NULL )
			iv.add( r_lend );
		if( align_success )
			iv.add( cr );

		// scan ahead to the next pair of n-way matches
		recurse_prev = recurse_iter;
		if( recurse_iter != match_list.end() )
			recurse_iter++;
		while( recurse_iter != match_list.end() && (*recurse_iter)->Multiplicity() != seq_count )
			recurse_iter++;

		if( ( recurse_iter == match_list.end() && !collinear_genomes ) ||
				( recurse_prev == match_list.end() && collinear_genomes ) )
				break;
	}
	// get the last little bit at the end of the LCB.
	iv.matches.insert( iv.matches.end(), recurse_prev, recurse_iter );
	mlist.insert( mlist.end(), match_list.begin(), match_list.end() );
	iv.CalculateOffset();

	}catch( gnException& gne ){
		cerr << gne << endl;
	}catch(exception& e){
		cerr << e.what();
	}catch(...){
		cerr << "matrix exception?\n";
	}
}

// just search each intervening region once for matches, no gapped alignment...
void MpiAligner::SearchWithinLCB( MatchList& mlist ){
	// check whether this function can do anything useful...
	if( !collinear_genomes && mlist.size() < 2 )
		return;

	boolean debug_recurse = false;
	int64 config_value = 138500;
	int print_interval = 50;

	try{
	list< Match* > match_list;
	match_list.insert( match_list.end(), mlist.begin(), mlist.end() );
	mlist.clear();
	MatchList r_list = mlist;

	list< Match* >::iterator recurse_iter = match_list.begin();
	list< Match* >::iterator recurse_prev = match_list.begin();
	if( !collinear_genomes && recurse_iter != match_list.end() )
		recurse_iter++;
	
	uint memI = 0;
	uint matchI = 0;
	while( recurse_prev != match_list.end() ){
		if( memI >= print_interval && memI % print_interval == 0 || debug)
			cout << "Number: " << memI << " match " << **recurse_prev << endl;
		memI++;
		if( debug_recurse ){
			cout << "Recursing on " << endl;
			if( recurse_prev != match_list.end() )
				cout << **recurse_prev << " and " << endl;
			if( recurse_iter != match_list.end() )
				cout << **recurse_iter << endl;
		}
		
		
		// recurse on a pair of matches! 
		// this function should locate all matches between the two iterators
		// and add them to r_list		
		r_list.clear();
		Match* r_left = NULL;
		Match* r_right = NULL;
		if( recurse_iter == match_list.begin() ){
			r_left = NULL;
			r_right = *recurse_iter;
		}else if( recurse_iter == match_list.end() ){
			r_left = *recurse_prev;
			r_right = NULL;
		}else{
			r_left = *recurse_prev;
			r_right = *recurse_iter;
		}
		Recursion( r_list, r_left, r_right, true );
		if( debug_recurse ){
			vector< Match* >::iterator r_iter = r_list.begin();
			cout << "Found matches " << endl;
			for(; r_iter != r_list.end(); r_iter++ )
				cout << **r_iter << endl;
		}

		// insert any n-way matches into the match list
		for( matchI = 0; matchI < r_list.size(); matchI++ ){
			if( r_list[ matchI ]->Multiplicity() == seq_count ){
				match_list.insert( recurse_iter, r_list[ matchI ] );
			}else
				match_allocator.Free( r_list[ matchI ] );
		}

		// move ahead to the next pair of n-way matches
		recurse_prev = recurse_iter;
		if( recurse_iter != match_list.end() )
			recurse_iter++;
		
		// break early if we aren't assuming genome collinearity
		if( !collinear_genomes && recurse_iter == match_list.end() )
			break;
			
	}

	mlist.insert( mlist.begin(), match_list.begin(), match_list.end() );

	}catch( gnException& gne ){
		cerr << gne << endl;
	}catch(exception& e){
		cerr << e.what();
	}catch(...){
		cerr << "matrix exception?\n";
	}
}

