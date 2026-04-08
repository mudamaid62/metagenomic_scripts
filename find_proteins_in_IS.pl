#!/usr/bin/perl
use warnings;
use strict;

my $proteins_positions = shift(@ARGV);
my $m8_taxed = shift(@ARGV); #single metagenome
my $FS_identification_table = shift(@ARGV);
my $mode = shift(@ARGV); #SARG or VFDB
my $isescan_tsv = shift(@ARGV);
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
open(ISES,"$isescan_tsv");
my %is_start;
my %is_end;
my %is_cluster;
my %is_family;
my %is_type;
my $contig_uniq = 0;
my @uniqs;
while(my $x = <ISES>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $contig = $array[0];
	$contig = "$contig~$contig_uniq";
	print "$contig\n";
	$is_start{$contig} = $array[3];
	$is_end{$contig} = $array[4];
	$is_type{$contig} = $array[21];
	$is_cluster{$contig} = $array[2];
	$is_family{$contig} = $array[1];
	push @uniqs, $contig;
	$contig_uniq++;
}
close ISES;
open(M8,"$m8_taxed");
my %total_observations;
my %ises_observations;
open(TABLE,">$out_table");
if($mode eq "SARG"){
	print TABLE "protein\tcontig\tstart\tend\tsense\tFS_cluster\tObserved_genes\tType\tMechanism\tsub1\tsub2\tspecial\tclose_to_is\tis_type\tis_family\tis_cluster\ttaxonomy\n";
}elsif($mode eq "VFDB"){
	print TABLE "protein\tcontig\tstart\tend\tsense\tFS_cluster\tObserved_genes\tVFC\tVF\tVF_full\tclose_to_is\tis_type\tis_family\tis_cluster\ttaxonomy\n";
}
open(FREQ,">$out_frequencies");
print FREQ "protein\ttotal_observations\tis_observations\tis_freq\n";
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
	if(!exists($ises_observations{$hit})){
                $ises_observations{$hit} = 0;
        }
	$total_observations{$hit}++;
	my $close_to_is = "FALSE";
	my $type = "NA";
	my $family = "NA";
	my $cluster = "NA";
	my $contig_name = $contigs{$protein};
	foreach my $c(@uniqs){
		if($close_to_is eq "TRUE"){
			next;
		}
		my($og_contig,$uniq) = split(/~/,$c);
		if($og_contig eq $contig_name){
			#$is_start{$contigs{$protein}} $starts{$protein} $is_end{$contigs{$protein}} $ends{$protein}
			my $distance = 5001;
			if($ends{$protein} <= $is_start{$c}){
				$distance = $is_start{$c} - $ends{$protein};
			}elsif($is_end{$c} <= $starts{$protein}){
				$distance = $starts{$protein} - $is_end{$c};
			}elsif($starts{$protein} >= $is_start{$c} and $ends{$protein} <= $is_end{$c}){
				$distance = 0;
			}
			if($distance <= 5000){
				$ises_observations{$hit}++;
				$close_to_is = "TRUE";
				$type = $is_type{$c};
        			$family = $is_family{$c};
        			$cluster = $is_cluster{$c};
			}
		}
	}
	my $data = $clusters{$FS};
	my @data_array = split(/\t/,$data);
	if($mode eq "SARG"){
		print TABLE "$hit\t$contigs{$protein}\t$starts{$protein}\t$ends{$protein}\t$senses{$protein}\t$FS\t$data_array[0]\t$data_array[2]\t$data_array[3]\t$data_array[4]\t$data_array[5]\t$data_array[6]\t$close_to_is\t$type\t$family\t$cluster\t$taxonomy\n";
	}elsif($mode eq "VFDB"){
		print TABLE "$hit\t$contigs{$protein}\t$starts{$protein}\t$ends{$protein}\t$senses{$protein}\t$FS\t$data_array[0]\t$data_array[2]\t$data_array[5]\t$data_array[4]\t$close_to_is\t$type\t$family\t$cluster\t$taxonomy\n";
	}
}
foreach my $h(sort keys %total_observations){
	my $freq = $ises_observations{$h}/$total_observations{$h};
	print FREQ "$h\t$total_observations{$h}\t$ises_observations{$h}\t$freq\n";
}
