#!/usr/bin/perl
#5
use warnings;
use strict;

my $genomad_plasmid_summary = shift(@ARGV);
my $c2c = shift(@ARGV); #0 and 3
my $min_score = shift(@ARGV);
my $site_name = shift(@ARGV);

my %taxas;
my %lens;
my $check_len = 0;
open(CAT,"$c2c");
while(my $x = <CAT>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $contig = $array[0];
	if($contig =~ m/^\#/){
		next;
	}
	my @name_array = split(/_/,$contig);
	$lens{$contig} = $name_array[3];
	$check_len += $lens{$contig};
	my $tax = $array[3];
	my $check = $array[1];
	if($check ne "taxid assigned"){
		$taxas{$contig} = $check;
	}else{
		$taxas{$contig} = $tax;
	}
}
close CAT;
open(PLASMIDS,"$genomad_plasmid_summary");
my $counter = 0;
my %tax_counts;
my %tax_bases;
while(my $x = <PLASMIDS>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $contig = $array[0];
	my $score = $array[5];
	if($score eq "plasmid_score"){
		next;
	}elsif($score >= $min_score){
		$counter++;
		if(!exists($tax_counts{$taxas{$contig}})){
			$tax_counts{$taxas{$contig}} = 1;
			$tax_bases{$taxas{$contig}} = $lens{$contig}; 
		}else{
			$tax_counts{$taxas{$contig}}++;
			$tax_bases{$taxas{$contig}} += $lens{$contig};
		}
	}
}
close PLASMIDS;
print "Site\tTaxonomy\tPlasmid_count\tPlasmid_bases\n";
foreach my $z(sort keys %tax_counts){
	print "$site_name\t$z\t$tax_counts{$z}\t$tax_bases{$z}\n";
}
print STDERR "$site_name\t$counter\t$check_len\n";
