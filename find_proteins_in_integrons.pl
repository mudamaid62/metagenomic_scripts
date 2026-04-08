#!/usr/bin/perl
use warnings;
use strict;

my $proteins_positions = shift(@ARGV);
my $m8_taxed = shift(@ARGV); #single metagenome
my $FS_identification_table = shift(@ARGV);
my $mode = shift(@ARGV); #SARG or VFDB
my $integron_file = shift(@ARGV);
my $out_table = shift(@ARGV);
my $out_frequencies = shift(@ARGV);

my %contigs;
my %starts;
my %ends;
my %senses;
open(PROTEINS,"$proteins_positions");
while(my $x = <PROTEINS>){
	chomp($x);
	my ($protein,$contig,$start,$end,$sense) = split(/\t/,$x);
	$contigs{$protein} = $contig;
	$starts{$protein} = $start;
	$ends{$protein} = $end;
	$senses{$protein} = $sense;
}
close PROTEINS;
open(ID,"$FS_identification_table");
my %clusters;
while(my $x = <ID>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $cluster = shift(@array);
	my $data = join "\t",@array;
	$clusters{$cluster} = $data;
}
close ID;
open(INTEGRONS,"$integron_file");
my %in_start;
my %in_end;
my %in_type;
while(my $x = <INTEGRONS>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $contig = $array[0];
	$in_start{$contig} = $array[1];
	$in_end{$contig} = $array[2];
	$in_type{$contig} = $array[3]; 
}
close INTEGRONS;
open(M8,"$m8_taxed");
my %total_observations;
my %integron_observations;
open(TABLE,">$out_table");
if($mode eq "SARG"){
	print TABLE "protein\tcontig\tstart\tend\tsense\tFS_cluster\tObserved_genes\tType\tMechanism\tsub1\tsub2\tspecial\tin_integron\tintegron_type\ttaxonomy\n";
}elsif($mode eq "VFDB"){
	print TABLE "protein\tcontig\tstart\tend\tsense\tFS_cluster\tObserved_genes\tVFC\tVF\tVF_full\tin_integron\tintegron_type\ttaxonomy\n";
}
open(FREQ,">$out_frequencies");
print FREQ "protein\ttotal_observations\tintegron_observations\tintegron_freq\n";
while(my $x = <M8>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $protein = $array[1];
	my $hit = $array[3];
	my $FS = $array[15];
	my $taxonomy = $array[16];
	if(!exists($total_observations{$hit})){
		$total_observations{$hit} = 0;
	}
	if(!exists($integron_observations{$hit})){
                $integron_observations{$hit} = 0;
        }
	$total_observations{$hit}++;
	my $in_integron = "FALSE";
	my $integron_type = "NA";
	if(exists($in_type{$contigs{$protein}})){
		if($in_start{$contigs{$protein}} <= $starts{$protein} and $in_end{$contigs{$protein}} >= $ends{$protein}){
			$integron_observations{$hit}++;
			$in_integron = "TRUE";
			$integron_type = $in_type{$contigs{$protein}};
		}
	}
	my $data = $clusters{$FS};
	my @data_array = split(/\t/,$data);
	if($mode eq "SARG"){
		print TABLE "$hit\t$contigs{$protein}\t$starts{$protein}\t$ends{$protein}\t$senses{$protein}\t$FS\t$data_array[0]\t$data_array[2]\t$data_array[3]\t$data_array[4]\t$data_array[5]\t$data_array[6]\t$in_integron\t$integron_type\t$taxonomy\n";
	}elsif($mode eq "VFDB"){
		print TABLE "$hit\t$contigs{$protein}\t$starts{$protein}\t$ends{$protein}\t$senses{$protein}\t$FS\t$data_array[0]\t$data_array[2]\t$data_array[5]\t$data_array[4]\t$in_integron\t$integron_type\t$taxonomy\n";
	}
}
foreach my $h(sort keys %total_observations){
	my $freq = $integron_observations{$h}/$total_observations{$h};
	print FREQ "$h\t$total_observations{$h}\t$integron_observations{$h}\t$freq\n";
}
