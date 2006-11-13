#ifndef __SeedMatchEnumerator_h__
#define __SeedMatchEnumerator_h__

#include "libMems/MatchFinder.h"
#include "libMems/RepeatHash.h"
#include "libMems/MemHash.h"
#include "libMems/MatchList.h"
#include "libMems/SortedMerList.h"
#include "libMems/Match.h"

/**
 * Turns every seed match into a full match without extension.
 */
class SeedMatchEnumerator : public mems::MatchFinder 
{
public:
	virtual SeedMatchEnumerator* Clone() const;

	void FindMatches( mems::MatchList& match_list )
	{
		for( size_t seqI = 0; seqI < match_list.seq_table.size(); ++seqI ){
			if( !AddSequence( match_list.sml_table[ seqI ], match_list.seq_table[ seqI ] ) ){
				genome::ErrorMsg( "Error adding " + match_list.seq_filename[seqI] + "\n");
				return;
			}
		}
		CreateMatches();
		match_list.clear();
		match_list.insert( match_list.end(), mlist.begin(), mlist.end() );
	}

	virtual boolean CreateMatches();
protected:

	virtual boolean EnumerateMatches( mems::IdmerList& match_list );
	virtual boolean HashMatch(mems::IdmerList& match_list);
	virtual mems::SortedMerList* GetSar(uint32 sarI) const;
	mems::MatchList mlist;
	void SetDirection(mems::Match& mhe);
};

SeedMatchEnumerator* SeedMatchEnumerator::Clone() const{
	return new SeedMatchEnumerator(*this);
}

inline
mems::SortedMerList* SeedMatchEnumerator::GetSar(uint32 sarI) const{
	return sar_table[0];
}

boolean SeedMatchEnumerator::CreateMatches(){
	if(seq_count == 1){
		MatchFinder::FindMatchSeeds();
		return true;
	}
	return false;
}

boolean SeedMatchEnumerator::EnumerateMatches( mems::IdmerList& match_list ){
	return HashMatch(match_list);
}

boolean SeedMatchEnumerator::HashMatch(mems::IdmerList& match_list){
	//check that there is at least one forward component
	match_list.sort(&mems::idmer_position_lessthan);
	// initialize the hash entry
	mems::Match mhe = mems::Match( match_list.size() );
	mhe.SetLength( GetSar(0)->SeedLength() );
	
	//Fill in the new Match and set direction parity if needed.
	mems::IdmerList::iterator iter = match_list.begin();

	uint32 repeatI = 0;
	for(; iter != match_list.end(); iter++)
		mhe.SetStart(repeatI++, iter->position + 1);

	SetDirection( mhe );
	if(mhe.Multiplicity() < 2){
		std::cerr << "red flag " << mhe << "\n";
	}else{
		mlist.push_back(mhe.Copy());
	}
	return true;
}

// evil, evil code duplication.

void SeedMatchEnumerator::SetDirection(mems::Match& mhe){
	//get the reference direction
	boolean ref_forward = false;
	uint32 seqI=0;
	for(; seqI < mhe.SeqCount(); ++seqI)
		if(mhe[seqI] != mems::NO_MATCH){
			ref_forward = !(GetSar(seqI)->GetMer(mhe[seqI] - 1) & 0x1);
			break;
		}
	//set directional parity for the rest
	for(++seqI; seqI < mhe.SeqCount(); ++seqI)
		if(mhe[seqI] != mems::NO_MATCH)
			if(ref_forward == (GetSar(seqI)->GetMer(mhe[seqI] - 1) & 0x1))
				mhe.SetStart(seqI, -mhe[seqI]);
}


#endif	// __SeedMatchEnumerator_h__
