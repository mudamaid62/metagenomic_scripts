#!/usr/bin/perl
use warnings;
use strict;

my $proteins_positions = shift(@ARGV);
my $m8_taxed = shift(@ARGV); #single metagenome
my $FS_identification_table = shift(@ARGV);
my $mode = shift(@ARGV); #SARG or VFDB
my $genomad_plasmid_summary = shift(@ARGV);
my $score_threshold = shift(@ARGV);
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
open(PLASMIDS,"$genomad_plasmid_summary");
my %plasmid_score;
my %plasmid_fdr;
while(my $x = <PLASMIDS>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $contig = $array[0];
	if($contig eq "seq_name"){
		next;
	}
	my $score = $array[5];
	my $fdr = $array[6];
	if($score >= $score_threshold){
		$plasmid_score{$contig} = $score;
		$plasmid_fdr{$contig} = $fdr;
	}
}
close PLASMIDS;
open(M8,"$m8_taxed");
my %total_observations;
my %plasmid_observations;
open(TABLE,">$out_table");
if($mode eq "SARG"){
	print TABLE "protein\tcontig\tstart\tend\tsense\tFS_cluster\tObserved_genes\tType\tMechanism\tsub1\tsub2\tspecial\tin_plasmid\tplasmid_score\tplasmid_fdr\ttaxonomy\n";
}elsif($mode eq "VFDB"){
	print TABLE "protein\tcontig\tstart\tend\tsense\tFS_cluster\tObserved_genes\tVFC\tVF\tVF_full\tin_plasmid\tplasmid_score\tplasmid_fdr\ttaxonomy\n";
}
open(FREQ,">$out_frequencies");
print FREQ "protein\ttotal_observations\tplasmid_observations\tplasmid_freq\n";
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
	if(!exists($plasmid_observations{$hit})){
                $plasmid_observations{$hit} = 0;
        }
	$total_observations{$hit}++;
	my $in_plasmid = "FALSE";
	my $score = "NA";
	my $fdr = "NA";
	if(exists($plasmid_score{$contigs{$protein}})){
		$plasmid_observations{$hit}++;
		$score = $plasmid_score{$contigs{$protein}};
		$fdr = $plasmid_fdr{$contigs{$protein}}; 
		$in_plasmid = "TRUE";
	}
	my $data = $clusters{$FS};
	my @data_array = split(/\t/,$data);
	if($mode eq "SARG"){
		print TABLE "$hit\t$contigs{$protein}\t$starts{$protein}\t$ends{$protein}\t$senses{$protein}\t$FS\t$data_array[0]\t$data_array[2]\t$data_array[3]\t$data_array[4]\t$data_array[5]\t$data_array[6]\t$in_plasmid\t$score\t$fdr\t$taxonomy\n";
	}elsif($mode eq "VFDB"){
		print TABLE "$hit\t$contigs{$protein}\t$starts{$protein}\t$ends{$protein}\t$senses{$protein}\t$FS\t$data_array[0]\t$data_array[2]\t$data_array[5]\t$data_array[4]\t$in_plasmid\t$score\t$fdr\t$taxonomy\n";	
	}
}
foreach my $h(sort keys %total_observations){
	my $freq = $plasmid_observations{$h}/$total_observations{$h};
	print FREQ "$h\t$total_observations{$h}\t$plasmid_observations{$h}\t$freq\n";
}
