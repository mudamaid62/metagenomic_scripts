#!/usr/bin/perl
#8
use warnings;
use strict;

my $phigaro_gff3 = shift(@ARGV);

my %tax_counts;
open(VIRUSES,"$phigaro_gff3");
my $counter = 0;
my $trans_counter = 0;
while(my $x = <VIRUSES>){
	chomp($x);
	my @array = split(/\t/,$x);
	if(exists($array[2]) and $array[2] eq "prophage"){
		my @data = split(/\;/,$array[8]);
		my @tax_array = split(/\=/,$data[1]);
		my @trans_array = split(/\=/,$data[2]);
		my $tax = $tax_array[1];
		$counter++;
		if(!exists($tax_counts{$tax})){
			$tax_counts{$tax} = 1;
		}else{
			$tax_counts{$tax}++;
		}
		if($trans_array[1] eq "True"){
			$trans_counter++;
		}
	}
}
print "TOTAL\t$counter\t100\n;";
my $trans_frac = ($trans_counter/$counter)*100;
print "Transposable\t$trans_counter\t$trans_frac\n";
my $check = 0;
foreach my $z(sort keys %tax_counts){
	my $frac = ($tax_counts{$z}/$counter)*100;
	print "$z\t$tax_counts{$z}\t$frac\n";
	$check += $frac;
}
print STDERR "$check\%\n";
