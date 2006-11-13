/*******************************************************************************
 * $Id: Utilities.cpp,v 1.27 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#include "libMems/MuscleInterface.h"
#include "libMems/MuscleInterface.cpp"

#include "procrastUtilities.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnStringTools.h"

#include "boost/algorithm/string/erase.hpp"
#include "boost/format.hpp"

#include <sstream>
#include <fstream>

using namespace std;
using namespace genome;
//using namespace mems;

bool CallbagFFT( string path,  int N, int K, double s, double* pu, int Q, string& pvalue )
{

	try{
		ostringstream input_seq_stream;
		genome::gnSequence seq;
	
		string bagfft_params = str( boost::format("%1% %2% %3% %4% ") % N % K % s % Q );
		for(int i = 0; i < K; i++)
		{
			bagfft_params = str( boost::format("%1% %2%") % bagfft_params % pu[i] );
		}
		
		string bagfft_cmd = path + " " + bagfft_params;
		char** bagfft_cmdline = mems::parseCommand( bagfft_cmd );
		string output;
		string error;
		
		cout << bagfft_cmd << endl;
		bool success = mems::pipeExec( bagfft_cmdline, bagfft_cmd, input_seq_stream.str(), output, error );
		if( !success || output.size() == 0 )
			return false;

		istringstream output_bagfft_stream( output );
		string cur_line;

		// parse the bagFFT output
		bool ok = false;
		while( getline( output_bagfft_stream, cur_line ) )
		{
						
			string::size_type loc = cur_line.find( "p-Value", 0 );
		    if( loc != string::npos )
			{
			  string::size_type end = cur_line.find( "(", loc );
			  string::size_type start = cur_line.find( "=", loc );
			  pvalue = cur_line.substr(start+1,end-(start+1));
			  //pvalue = pvalue.replace(pvalue.find("e",0),2,".e");
			  ok = true;
			  break;
			}
			else
			{
				continue;
			}	
		}
	    
		return ok;
	}catch( gnException& gne ){
	}catch( exception& e ){
	}catch(...){
	}
	return false;
}
