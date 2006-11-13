#include "libMems/Interval.h"
#include "libMems/Islands.h"

using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	IntervalList iv_list;
	MatchList mlist;
	if( argc != 3 )
	{
		cerr << "Usage: <input interval file> <output interval file>";
		return -1;
	}
	ifstream in_file( argv[1] );
	if( !in_file.is_open() )
	{
		cerr << "Error opening \"argv[1]\"\n";
		return -1;
	}
	iv_list.ReadList( in_file );
	mlist.seq_filename = iv_list.seq_filename;
	mlist.LoadSequences(NULL);
	iv_list.seq_table = mlist.seq_table;
	iv_list.seq_filename = mlist.seq_filename;
	addUnalignedIntervals( iv_list );
	ofstream out_file( argv[2] );
	if( !out_file.is_open() )
	{
		cerr << "Error opening \"argv[2]\"\n";
		return -2;
	}
	iv_list.WriteList( out_file );
	return 0;
}
