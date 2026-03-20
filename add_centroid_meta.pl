#!/usr/bin/perl
use warnings;
use strict;

my $meta = shift(@ARGV);
my $list = shift(@ARGV);

my @fs_list;
my %centroids;
open(CENTROIDS,"$list");
while(my $x = <CENTROIDS>){
	chomp($x);
	my ($fs,$centroid) = split(/\t/,$x);
	$centroids{$centroid} = $fs;
	push @fs_list, $fs;
}
close CENTROIDS;
open(META,"$meta");
print "Protein\tReference\tSingle\tCentroid\tFS_cluster\tpident_modifier\n";

my %ref;
my %fs;
my %mod;
while(my $x = <META>){
	chomp($x);
	my ($protein,$reference,$cluster,$modifier) = split(/\t/,$x);
	$ref{$protein} = $reference;
	$fs{$protein} = $cluster;
	$mod{$protein} = $modifier;
}
close META;
my %counts;
foreach my $x(@fs_list){
	my $c = 0;
	foreach my $k(keys %fs){
		if($fs{$k} eq $x){
			$c++;
		}
	}
	$counts{$x} = $c;
}	
foreach my $x(keys %fs){
	my $single = "No";
	if($counts{$fs{$x}} == 1){
		$single = "Yes";
	}
	my $cen = "No";
	if(exists($centroids{$x})){
		$cen = "Yes";
	}
	print "$x\t$ref{$x}\t$single\t$cen\t$fs{$x}\t$mod{$x}\n";
}

