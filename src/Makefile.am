
if DEBUG
D_CXXFLAGS = -Wall -g -DCOMMAND_LINE -D__GNDEBUG__ 
else
P_CXXFLAGS = -DCOMMAND_LINE
endif
AM_CXXFLAGS = $(D_CXXFLAGS) $(P_CXXFLAGS) 

LIBRARY_CL = $(DEPS_LIBS) 
INCLUDES = @DEPS_CFLAGS@ 
# ./bagFFT/

bin_PROGRAMS = mauveAligner mauveStatic scoreAlignment \
extractBackbone transposeCoordinates uniqueMerCount toRawSequence\
calculateBackboneCoverage calculateCoverage mfa2xmfa addUnalignedIntervals\
toMultiFastA getAlignmentWindows extractBackbone2 uniquifyTrees \
countInPlaceInversions toGrimmFormat joinAlignmentFiles mauveToXMFA \
stripGapColumns gappiness progressiveMauve progressiveMauveStatic \
sortContigs extractBCITrees calculateBackboneCoverage2 createBackboneMFA \
pairCompare procrastAligner scoreALU evd alignmentProjector stripSubsetLCBs \
projectAndStrip makeBadgerMatrix multiEVD randomGeneSample getOrthologList

EXTRA_PROGRAMS = mauveMpatrol mauveEfence toGBKsequence multiToRawSequence unalign

mauveAligner_SOURCES = mauveAligner.cpp mauveAligner.h
mauveAligner_LDADD = $(LIBRARY_CL)
EXTRA_mauveAligner_SOURCES = getopt.c getopt.h getopt1.c
mauveAligner_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

mauveMpatrol_SOURCES = mauveAligner.cpp
mauveMpatrol_LDFLAGS = -lmpatrol -lbfd -liberty $(LIBRARY_CL) `wx-config --libs`
EXTRA_mauveMpatrol_SOURCES = getopt.c getopt.h getopt1.c
mauveMpatrol_DEPENDENCIES = @GETOPT_LONG_SYSTEM@



mauveEfence_SOURCES = mauveAligner.cpp
mauveEfence_LDADD = -lefence $(LIBRARY_CL)
EXTRA_mauveEfence_SOURCES = getopt.c getopt.h getopt1.c
mauveEfence_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

mauveStatic_SOURCES = mauveAligner.cpp
mauveStatic_LDFLAGS = @STATIC_FLAG@ 
mauveStatic_LDADD = $(LIBRARY_CL)
EXTRA_mauveStatic_SOURCES = getopt.c getopt.h getopt1.c
mauveStatic_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

#mauveAligner4_SOURCES = mauveAligner.cpp
#mauveAligner4_LDFLAGS = @STATIC_FLAG@ 
#mauveAligner4_LDADD = $(LIBRARY_CL)
#mauveAligner4_CXXFLAGS = MAX_SEQ_COUNT=4
#EXTRA_mauveAligner4_SOURCES = getopt.c getopt.h getopt1.c
#mauveAligner4_DEPENDENCIES = @GETOPT_LONG_SYSTEM@


calculateCoverage_SOURCES = calculateCoverage.cpp
calculateCoverage_LDFLAGS = @STATIC_FLAG@ 
calculateCoverage_LDADD = $(LIBRARY_CL)
EXTRA_calculateCoverage_SOURCES = getopt.c getopt.h getopt1.c
calculateCoverage_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

scoreAlignment_SOURCES = scoreAlignment.cpp
scoreAlignment_LDFLAGS = @STATIC_FLAG@
scoreAlignment_LDADD = $(LIBRARY_CL)
EXTRA_scoreAlignment_SOURCES = getopt.c getopt.h getopt1.c
scoreAlignment_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

extractBackbone_SOURCES = extractBackbone.cpp
extractBackbone_LDFLAGS = @STATIC_FLAG@
extractBackbone_LDADD = $(LIBRARY_CL)
EXTRA_extractBackbone_SOURCES = getopt.c getopt.h getopt1.c
extractBackbone_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

toRawSequence_SOURCES = toRawSequence.cpp
toRawSequence_LDFLAGS = @STATIC_FLAG@
toRawSequence_LDADD = $(LIBRARY_CL)
EXTRA_toRawSequence_SOURCES = getopt.c getopt.h getopt1.c
toRawSequence_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

transposeCoordinates_SOURCES = transposeCoordinates.cpp
transposeCoordinates_LDFLAGS = @STATIC_FLAG@
transposeCoordinates_LDADD = $(LIBRARY_CL)
EXTRA_transposeCoordinates_SOURCES = getopt.c getopt.h getopt1.c
transposeCoordinates_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

uniqueMerCount_SOURCES = uniqueMerCount.cpp
uniqueMerCount_LDFLAGS = @STATIC_FLAG@
uniqueMerCount_LDADD = $(LIBRARY_CL)
EXTRA_uniqueMerCount_SOURCES = getopt.c getopt.h getopt1.c
uniqueMerCount_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

calculateBackboneCoverage_SOURCES = calculateBackboneCoverage.cpp
calculateBackboneCoverage_LDFLAGS = @STATIC_FLAG@
calculateBackboneCoverage_LDADD = $(LIBRARY_CL)
EXTRA_calculateBackboneCoverage_SOURCES = getopt.c getopt.h getopt1.c
calculateBackboneCoverage_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

mfa2xmfa_SOURCES = mfa2xmfa.cpp
mfa2xmfa_LDFLAGS = @STATIC_FLAG@
mfa2xmfa_LDADD = $(LIBRARY_CL)
EXTRA_mfa2xmfa_SOURCES = getopt.c getopt.h getopt1.c
mfa2xmfa_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

toGBKsequence_SOURCES = toGBKsequence.cpp
toGBKsequence_LDFLAGS = @STATIC_FLAG@
toGBKsequence_LDADD = $(LIBRARY_CL)
EXTRA_toGBKsequence_SOURCES = getopt.c getopt.h getopt1.c
toGBKsequence_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

multiToRawSequence_SOURCES = multiToRawSequence.cpp
multiToRawSequence_LDFLAGS = @STATIC_FLAG@
multiToRawSequence_LDADD = $(LIBRARY_CL)
EXTRA_multiToRawSequence_SOURCES = getopt.c getopt.h getopt1.c
multiToRawSequence_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

unalign_SOURCES = unalign.cpp
unalign_LDFLAGS = @STATIC_FLAG@
unalign_LDADD = $(LIBRARY_CL)
EXTRA_unalign_SOURCES = getopt.c getopt.h getopt1.c
unalign_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

addUnalignedIntervals_SOURCES = addUnalignedIntervals.cpp
addUnalignedIntervals_LDFLAGS = @STATIC_FLAG@
addUnalignedIntervals_LDADD = $(LIBRARY_CL)

toMultiFastA_SOURCES = toMultiFastA.cpp
toMultiFastA_LDFLAGS = @STATIC_FLAG@
toMultiFastA_LDADD = $(LIBRARY_CL)

getAlignmentWindows_SOURCES = getAlignmentWindows.cpp
getAlignmentWindows_LDFLAGS = @STATIC_FLAG@
getAlignmentWindows_LDADD = $(LIBRARY_CL)

extractBackbone2_SOURCES = extractBackbone2.cpp
extractBackbone2_LDFLAGS = @STATIC_FLAG@
extractBackbone2_LDADD = $(LIBRARY_CL)

uniquifyTrees_SOURCES = uniquifyTrees.cpp
uniquifyTrees_LDFLAGS = @STATIC_FLAG@
uniquifyTrees_LDADD = $(LIBRARY_CL)


countInPlaceInversions_SOURCES = countInPlaceInversions.cpp
countInPlaceInversions_LDFLAGS = @STATIC_FLAG@
countInPlaceInversions_LDADD = $(LIBRARY_CL)

toGrimmFormat_SOURCES = toGrimmFormat.cpp
toGrimmFormat_LDFLAGS = @STATIC_FLAG@
toGrimmFormat_LDADD = $(LIBRARY_CL)

joinAlignmentFiles_SOURCES = joinAlignmentFiles.cpp
joinAlignmentFiles_LDFLAGS = @STATIC_FLAG@
joinAlignmentFiles_LDADD = $(LIBRARY_CL)

mauveToXMFA_SOURCES = mauveToXMFA.cpp
mauveToXMFA_LDFLAGS = @STATIC_FLAG@
mauveToXMFA_LDADD = $(LIBRARY_CL)

stripGapColumns_SOURCES = stripGapColumns.cpp
stripGapColumns_LDFLAGS = @STATIC_FLAG@
stripGapColumns_LDADD = $(LIBRARY_CL)

gappiness_SOURCES = gappiness.cpp
gappiness_LDFLAGS = @STATIC_FLAG@
gappiness_LDADD = $(LIBRARY_CL)


progressiveMauve_SOURCES = progressiveMauve.cpp UniqueMatchFinder.h UniqueMatchFinder.cpp 
progressiveMauve_LDADD = $(LIBRARY_CL)

progressiveMauveStatic_SOURCES = progressiveMauve.cpp UniqueMatchFinder.h UniqueMatchFinder.cpp 
progressiveMauveStatic_LDFLAGS = @STATIC_FLAG@
progressiveMauveStatic_LDADD = $(LIBRARY_CL)

sortContigs_SOURCES = sortContigs.cpp
sortContigs_LDFLAGS = @STATIC_FLAG@
sortContigs_LDADD = $(LIBRARY_CL)

extractBCITrees_SOURCES = extractBCITrees.cpp
extractBCITrees_LDFLAGS = @STATIC_FLAG@
extractBCITrees_LDADD = $(LIBRARY_CL)

calculateBackboneCoverage2_SOURCES = calculateBackboneCoverage2.cpp
calculateBackboneCoverage2_LDFLAGS = @STATIC_FLAG@
calculateBackboneCoverage2_LDADD = $(LIBRARY_CL)

createBackboneMFA_SOURCES = createBackboneMFA.cpp
createBackboneMFA_LDFLAGS = @STATIC_FLAG@
createBackboneMFA_LDADD = $(LIBRARY_CL)

pairCompare_SOURCES = pairCompare.cpp
pairCompare_LDFLAGS = @STATIC_FLAG@
pairCompare_LDADD = $(LIBRARY_CL)

procrastAligner_SOURCES = procrastAligner.cpp MatchRecord.h SeedMatchEnumerator.h
procrastAligner_LDFLAGS = @STATIC_FLAG@
procrastAligner_LDADD = $(LIBRARY_CL) 

scoreALU_SOURCES = scoreALU.cpp
scoreALU_LDFLAGS = @STATIC_FLAG@
scoreALU_LDADD = $(LIBRARY_CL) 

evd_SOURCES = evd.cpp
evd_LDFLAGS = @STATIC_FLAG@
evd_LDADD = $(LIBRARY_CL)

alignmentProjector_SOURCES = alignmentProjector.cpp
alignmentProjector_LDFLAGS = @STATIC_FLAG@
alignmentProjector_LDADD = $(LIBRARY_CL)

stripSubsetLCBs_SOURCES = stripSubsetLCBs.cpp
stripSubsetLCBs_LDFLAGS = @STATIC_FLAG@
stripSubsetLCBs_LDADD = $(LIBRARY_CL)

projectAndStrip_SOURCES = projectAndStrip.cpp
projectAndStrip_LDFLAGS = @STATIC_FLAG@
projectAndStrip_LDADD = $(LIBRARY_CL)

makeBadgerMatrix_SOURCES = makeBadgerMatrix.cpp
makeBadgerMatrix_LDFLAGS = @STATIC_FLAG@
makeBadgerMatrix_LDADD = $(LIBRARY_CL)

multiEVD_SOURCES = multiEVD.cpp
multiEVD_LDFLAGS = @STATIC_FLAG@
multiEVD_LDADD = $(LIBRARY_CL)

randomGeneSample_SOURCES = randomGeneSample.cpp
randomGeneSample_LDFLAGS = @STATIC_FLAG@
randomGeneSample_LDADD = $(LIBRARY_CL)

getOrthologList_SOURCES = getOrthologList.cpp
getOrthologList_LDFLAGS = @STATIC_FLAG@
getOrthologList_LDADD = $(LIBRARY_CL)

