#!/usr/bin/env perl

use strict;

my @data = ("final_p1_i0_s0","final_p1_i15_s85","final_p1_i30_s170","final_p1_i60_s340","final_p1_i75_s425","final_p1_i150_s850","final_p2_i0_s0","final_p2_i15_s85","final_p2_i30_s170","final_p2_i60_s340","final_p2_i75_s425","final_p2_i150_s850","final_p3_i0_s0","final_p3_i15_s85","final_p3_i30_s170","final_p3_i60_s340","final_p3_i75_s425","final_p3_i150_s850","final_p4_i0_s0","final_p4_i15_s85","final_p4_i30_s170","final_p4_i60_s340","final_p4_i75_s425","final_p4_i150_s850","final_p5_i0_s0","final_p5_i15_s85","final_p5_i30_s170","final_p5_i60_s340","final_p5_i75_s425","final_p5_i150_s850","final_p6_i0_s0","final_p6_i15_s85","final_p6_i30_s170","final_p6_i60_s340","final_p6_i75_s425","final_p6_i150_s850","final_p7_i0_s0","final_p7_i15_s85","final_p7_i30_s170","final_p7_i60_s340","final_p7_i75_s425","final_p7_i150_s850","final_p8_i0_s0","final_p8_i15_s85","final_p8_i30_s170","final_p8_i60_s340","final_p8_i75_s425","final_p8_i150_s850","final_p9_i0_s0","final_p9_i15_s85","final_p9_i30_s170","final_p9_i60_s340","final_p9_i75_s425","final_p9_i150_s850","final_p10_i0_s0","final_p10_i15_s85","final_p10_i30_s170","final_p10_i60_s340","final_p10_i75_s425","final_p10_i150_s850");


my @seedweights = (13,13,13,13,13,13,
									  13,13,13,13,13,13,
									  13,13,13,13,13,13,
									  13,13,13,13,13,13,
									  13,13,13,13,13,13,
									  13,13,13,13,13,13,
									  13,13,13,13,13,13,
									  13,13,13,13,13,13,
									  13,13,13,13,13,13,
									  13,13,13,13,13,13);
									  
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
	my $procrast_cmd = "./procrastAligner --z $z --solid 0 --xmfa=$name.procrast.xmfa --output $name.procrast --sequence ../datasets/$name.fna 2> $name.procrast.err > $name.procrast.out";
	my $euler_cmd = "./EulerAlign -k $z -l -i 50 -v ../datasets/$name.fna $name.euler > $name.euler.out ";
	
	print "Executing $procrast_cmd\n";
	my $start_time = time();
	`$procrast_cmd`;
	my $end_time = time();
	open( ALN_TIME_FILE, ">$name.time" );
	print ALN_TIME_FILE ($end_time - $start_time);
	print ALN_TIME_FILE " seconds\n";
	close ALN_TIME_FILE;

	my $score_alu_cmd = "./scoreALU --alignment $name.procrast --alus ../datasets/$name.alu 2> $name.procrast.score.err > $name.procrast.alu.score.out";
	print "Executing $score_alu_cmd\n";
	`$score_alu_cmd`;
	my $score_cmd = "python scoreAlignments.py ../datasets/$name.tup $name.procrast.xmfa > $name.procrast.score.out";
	print "Executing $score_cmd\n";
	`$score_cmd`;
}