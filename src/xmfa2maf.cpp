#include "libMems/IntervalList.h"
#include <fstream>

using namespace mems;
using namespace std;
using namespace genome;

int main(int argc, char* argv[] ){
	if(argc != 3){
		cerr << "Usage: xmfa2maf <xmfa input> <maf output>\n";
		return -1;
	}
	ifstream ifile(argv[1]);
	if(!ifile.is_open()){
		cerr << "Error reading \"" << argv[1] << "\"\n";
		return -2;
	}
	ofstream ofile(argv[2]);
	if(!ofile.is_open()){
		cerr << "Error writing to \"" << argv[2] << "\"\n";
		return -2;
	}

	IntervalList xmfa;
	xmfa.ReadStandardAlignment(ifile);
	LoadSequences(xmfa, &cout);

	ofile << "##maf version=1 program=progressiveMauve\n";
	for(int ivI=0; ivI < xmfa.size(); ivI++ ){
		ofile << "a\n";
		vector<string> aln;
		GetAlignment( xmfa[ivI], xmfa.seq_table, aln );
		for( int seqI=0; seqI < xmfa.seq_filename.size(); seqI++ ){
			if(xmfa[ivI].LeftEnd(seqI)==0)
				continue;	// sequence not defined in this block
			ofile << "s " << xmfa.seq_filename[seqI];
			
			ofile << " " << xmfa[ivI].LeftEnd(seqI)-1;
			ofile << " " << xmfa[ivI].Length(seqI);
			ofile << " " << (xmfa[ivI].Orientation(seqI) == AbstractMatch::reverse ? "-" : "+");
			ofile << " " << xmfa.seq_table[seqI]->length();
			ofile << " " << aln[seqI] << endl;
		}
		ofile << endl;
	}
	ofile.close();

	return 0;
}
