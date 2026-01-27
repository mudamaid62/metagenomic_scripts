#!/usr/bin/perl
use warnings;
use strict;

my $table = shift(@ARGV);
my $extra_hill = shift(@ARGV); #get extra Hill numbers, comma-separated

my %counts;
my %samples;
my @extra_qs;
my $extra_header = "";
open(TABLE, "$table");

while(my $x = <TABLE>){
	chomp($x);
	my($otu,$sample,$freq) = split(/\t/,$x);
	my $name = "$otu~$sample";
	$counts{$name} = $freq;
	if(!exists($samples{$sample})){
		$samples{$sample} = "added";
	}
}
if(defined($extra_hill)){
	my @q_array = split(/\,/,$extra_hill);
	my @names;
 	foreach my $h(@q_array){
		push @extra_qs, $h;
		my $q_name = "Hill_number($h)";
		push @names, $q_name;
	}
	$extra_header = join "\t",@names;
}
print "$table\tRichness R\tShannon H\tSimpson λ\tSimpson dominance 1/λ\tGini 1-λ\tPielou J\tHmax\tEfficiency η\t$extra_header\n";
foreach my $s(sort keys %samples){
	my %sample_counts;
	my %otus;
	my $total_counts;
	my $richness = 0;
	foreach my $k(keys %counts){
		my($otu,$sample) = split(/~/,$k);
		if($sample eq $s){
			$sample_counts{$otu} = $counts{$k};
			if(!exists($otus{$otu})){
				$otus{$otu} = "added";
				$richness++;
			}
			$total_counts += $counts{$k};
		}
	}
	if($richness == 1){
		print "$s\t1\n";
		next;
	}
	my @pis;
	my @logs;
	foreach my $z(keys %otus){
		my $pi = ($sample_counts{$z})/$total_counts;
		push @pis,$pi;
		my $log = log($pi);
		push @logs,$log;
	}
	my $n = $richness - 1;
	my $shannon_sum = 0;
	my $eta_sum = 0;
	for my $i(0..$n){
		my $product1 = $pis[$i] * $logs[$i];
		my $product2 = ($pis[$i] * $logs[$i])/(log($richness));
		$shannon_sum += $product1;
		$eta_sum += $product2;
	}
	my $shannon_index = $shannon_sum * -1;
	my $eta_index = $eta_sum * -1;
	my $pielou_index = $shannon_index/(log($richness));
	my $pis_reference = \@pis;
	my $d2 = get_hill(2,$n,$shannon_index,$pis_reference);
	my $simpson_index = 1/$d2;
	my $gini = 1 - $simpson_index;
	my $h_max = $shannon_index / $pielou_index;
	print "$s\t$richness\t$shannon_index\t$simpson_index\t$d2\t$gini\t$pielou_index\t$h_max\t$eta_index\t";
	foreach my $j(@extra_qs){
		my $q_index = get_hill($j,$n,$shannon_index,$pis_reference);
		print "$q_index\t";
	}
	print "\n";
}
		
sub get_hill{
	my ($q,$n,$shannon_index,$pis_reference) = @_;
	my $hill_number = 0;
	if($q == 1){
		$hill_number = exp($shannon_index);
	}else{
		my @pis = @$pis_reference;
		my $sum = 0;
		foreach my $i(0..$n){
			my $exp = ($pis[$i])**$q;
			$sum += $exp;
		}
		$hill_number = ($sum)**(1/(1-$q));
	}
	return $hill_number;
}			
