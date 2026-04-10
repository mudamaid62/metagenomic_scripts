#!/usr/bin/perl
#10
use warnings;
use strict;

my $genomad_virus_summary = shift(@ARGV);
my $c2c = shift(@ARGV);
my $min_score = shift(@ARGV);
my $tax_comparisons = shift(@ARGV); #output genomad/CAT tax comparison
my $site_name = shift(@ARGV);

my %c2c_taxas;
open(CAT,"$c2c");
while(my $x = <CAT>){
        chomp($x);
        my @array = split(/\t/,$x);
        my $contig = $array[0];
        my $tax = $array[3];
        my $check = $array[1];
        if($check ne "taxid assigned"){
                $c2c_taxas{$contig} = $check;
        }else{
                $c2c_taxas{$contig} = $tax;
        }
}
close CAT;

my %tax_counts;
my %virus_taxas;
open(VIRUSES,"$genomad_virus_summary");
my $counter = 0;
while(my $x = <VIRUSES>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $score = $array[6];
	my $tax = $array[10];
	if($score eq "virus_score"){
		next;
	}elsif($score >= $min_score){
		$counter++;
		my ($contig,$useless) = split(/\|/,$array[0]);
		$virus_taxas{$contig} = $tax;
		if(!exists($tax_counts{$tax})){
			$tax_counts{$tax} = 1;
		}else{
			$tax_counts{$tax}++;
		}
	}
}
close VIRUSES;
print "TOTAL\t$counter\t100\n;";
my $check = 0;
foreach my $z(sort keys %tax_counts){
	my $frac = ($tax_counts{$z}/$counter)*100;
	print "$z\t$tax_counts{$z}\t$frac\n";
	$check += $frac;
}
print STDERR "$check\%\n";
open(COMP,">$tax_comparisons");
print COMP "Site\tContig\tGenomad_virus_tax\tCAT_bacterial_tax\n";
foreach my $z(sort keys %virus_taxas){
	#print STDERR "$z\n";
	print COMP "$site_name\t$z\t$virus_taxas{$z}\t$c2c_taxas{$z}\n";
}
