#!/usr/bin/perl
use warnings;
use strict;

#Name    Completeness    Contamination
my $checkm2_file = shift(@ARGV);

my $nc = 0;
my $mq = 0;
my $lq = 0;
my $t = 0;
my $comp_sum = 0;
my $cont_sum = 0;
my $q_sum = 0;
my @comp;
my @cont;
my @qual;
open(FILE,"$checkm2_file");
while(my $x = <FILE>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $completeness = $array[1];
	my $contamination = $array[2];
	if($completeness eq "Completeness" and $contamination eq "Contamination"){
		next;
	}elsif($completeness > 90 and $contamination < 5){
		my $q = $completeness - (5 * $contamination);
		$nc++;
		$t++;
		$comp_sum += $completeness;
		$cont_sum += $contamination;
		$q_sum += $q;
		push @comp, $completeness;
		push @cont, $contamination;
		push @qual, $q;
	}elsif($completeness >= 50 and $contamination < 10){
		my $q = $completeness - (5 * $contamination);
		$mq++;
		$t++;
		$comp_sum += $completeness;
                $cont_sum += $contamination;
		$q_sum += $q;
                push @comp, $completeness;
                push @cont, $contamination;
		push @qual, $q;
	}else{
		$lq++;
		$t++;
	}
}
my @sorted_comp = sort{$b<=>$a}(@comp);
my @sorted_cont = sort{$b<=>$a}(@cont);
my @sorted_qual = sort{$b<=>$a}(@qual);
my $check = (scalar(@sorted_comp)) % 2;
my $pos = int((scalar(@sorted_comp)) / 2);
my $median_comp = 0;
my $median_cont = 0;
my $median_qual = 0;
if($check != 0){
	$median_comp = $sorted_comp[$pos];
	$median_cont = $sorted_cont[$pos];
	$median_qual = $sorted_qual[$pos];
}else{
	my $pre = $pos - 1;
	$median_comp = ($sorted_comp[$pre] + $sorted_comp[$pos]) / 2;
	$median_cont = ($sorted_cont[$pre] + $sorted_cont[$pos]) / 2;
	$median_qual = ($sorted_qual[$pre] + $sorted_qual[$pos]) / 2;
}
my $mean_comp = $comp_sum/scalar(@sorted_comp);
my $mean_cont = $cont_sum/scalar(@sorted_cont);
my $mean_qual = $q_sum/scalar(@sorted_qual);
my $max_comp = shift(@sorted_comp);
my $min_comp = pop(@sorted_comp);
my $max_cont = shift(@sorted_cont);
my $min_cont = pop(@sorted_cont);
my $max_qual = shift(@sorted_qual);
my $min_qual = pop(@sorted_qual);
print "NC --> $nc, MQ --> $mq, LQ --> $lq, Total --> $t\n";
print "Completeness: mean --> $mean_comp, min --> $min_comp, median --> $median_comp, max --> $max_comp\n";
print "Contamination: mean --> $mean_cont, min --> $min_cont, median --> $median_cont, max --> $max_cont\n";
print "Quality: mean --> $mean_qual, min --> $min_qual, median --> $median_qual, max --> $max_qual\n"; 
