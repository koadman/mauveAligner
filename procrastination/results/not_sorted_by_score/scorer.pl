#!/usr/bin/env perl

use strict;

my @data = ("AF435921-22Kb","Z15025-38Kb","AC034110-167Kb","AC010145-199kb","homo_sapiens_chr22-1MB");
my @seedweights = (9,9,15,15,15);
for( my $dataI = 0; $dataI < @data; $dataI++ )
{
	testData($data[$dataI], $seedweights[$dataI]);
}

exit 0;

sub testData
{
	my $name = shift;
	my $z = shift;
	print "Scoring $name...\n";
	my $procrast_cmd = "./procrastAligner --z $z --output $name.procrast --sequence ../../dataset/$name.fna 2> $name.procrast.err > $name.procrast.out";
	print "Executing $procrast_cmd\n";
	my $start_time = time();
	`$procrast_cmd`;
	my $end_time = time();
	open( ALN_TIME_FILE, ">$name.time" );
	print ALN_TIME_FILE ($end_time - $start_time);
	print ALN_TIME_FILE " seconds\n";
	close ALN_TIME_FILE;

	my $score_cmd = "./scoreALU --alignment $name.procrast --alus ../../dataset/$name.alu 2> $name.score.err > $name.score.out";
	print "Executing $score_cmd\n";
	`$score_cmd`;
}



