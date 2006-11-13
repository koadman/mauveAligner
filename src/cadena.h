/*******************************************************************************
 * $Id: Cadena.h 
 *******************************************************************************/

#ifndef _Cadena_h_
#define _Cadena_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "mauveAligner.h"
#include "getopt.h"
#include <sstream>
#include <stdexcept>
#include "libMems/Matrix.h"
#include "libMems/NumericMatrix.h"
#include "libGenome/gnSequence.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MemHash.h"
#include "libMems/RepeatHash.h"
#include "libMems/MaskedMemHash.h"
#include "libMems/Aligner.h"
#include "libMems/MatchList.h"
#include "libMems/MemSubsets.h"
#include "libMems/RepeatHashCat.h"
#include "libMems/Interval.h"
#include "libMems/IntervalList.h"
#include "libMems/gnAlignedSequences.h"
#include "libMems/Islands.h"
#include "libMems/MuscleInterface.h"

using namespace std;
using namespace mems;

/**
 * Used to find locally colinear blocks (LCBs) and do recursive
 * alignments on the blocks
 * To create an alignment one need only use the align method.
 * LCB lists are typically stored using the IntervalList class.  They can be
 * read and written in interval format using that class.  For input and output
 * of gapped alignments in other formats, see the gnAlignedSequences class.
 * Other methods in this class are available for experimentation.
 */
class Cadena {
public:
	/** 
	 * Constructs an aligner for the specified number of sequences.
	 * @param seq_count 	The number of sequences that will be aligned with this Aligner
	 */
	Cadena( unsigned int seq_count );
	Cadena( const Aligner& al );
	Cadena& operator=( const Cadena& cd );

	/**
	 * Chains a set of Repeats, returns highest scoring chain
	**/
	void encadenar( MatchList& mlist, IntervalList& interval_list, double LCB_minimum_density, double LCB_minimum_range, boolean recursive, boolean extend_lcbs, boolean gapped_alignment, string tree_filename = "" );
	
	void Recursion( MatchList& r_list, Match* r_begin, Match* r_end, boolean nway_only = false );
	void GetBestChain( MatchList& r_list, MatchList& best_lcb );
	//Matches may have spanned sequence boundaries in concatenated sequence, trim to bounds
	void BoundMatches( void );

	/** Set output parameters for permutation matrices */
	void SetPermutationOutput( std::string& permutation_filename, int64 permutation_weight );
	void WritePermutation( vector< LCB >& adjacencies, std::string out_filename );
protected:
	//MemHash gap_mh;			/**< Used during recursive alignment */
	//MaskedMemHash nway_mh;	/**< Used during recursive alignment to find nway matches only */
	RepeatHash repeat_finder;
	uint32 seq_count;		/**< The number of sequences this aligner is working with */
	boolean debug;			/**< Flag for debugging output */
	
	
	SlotAllocator< Match >& match_allocator;

	void consistencyCheck( uint lcb_count, vector< LCB >& adjacencies, vector< MatchList >& lcb_list, vector< int64 >& weights );
	
	boolean recursive;		/**< Set to true if a recursive anchor search/gapped alignment should be performed */
	boolean extend_lcbs;	/**< Set to true if LCB extension should be attempted */
	boolean gapped_alignment;	/**< Set to true to complete a gapped alignment */
	boolean currently_recursing;	/**< True when the recursive search has begun */
	boolean collinear_genomes;	/**< Set to true if all genomes are assumed to be collinear */
	
	GappedAligner* gal;

	std::string permutation_filename;
	int64 permutation_weight;
};

/**
 * Thrown if some error occurs during alignment
 */
CREATE_EXCEPTION( CadenaError );

#endif //cadena.h