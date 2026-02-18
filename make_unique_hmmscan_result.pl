#!/usr/bin/perl
use warnings;
use strict;

my $list = shift(@ARGV);

my %data;
my %seen;
my %score;

open(LIST,"$list");
while(my $x = <LIST>){
	chomp($x);
	my @array = split(/\t/,$x);
	if(!exists($seen{$array[2]})){
		$seen{$array[2]} = "added";
		$score{$array[2]} = $array[5];
		$data{$array[2]} = $x;
	}elsif($array[5] > $score{$array[2]}){
        	$score{$array[2]} = $array[5];
		$data{$array[2]} = $x;
	}
}
foreach my $x(sort keys %seen){
	print "$data{$x}\n";
}
