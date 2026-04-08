#!/usr/bin/perl
use warnings;
use strict;

my $proteins_positions = shift(@ARGV);
my $m8_taxed = shift(@ARGV); #single metagenome
my $FS_identification_table = shift(@ARGV);
my $mode = shift(@ARGV); #SARG or VFDB
my $genomad_virus_summary = shift(@ARGV);
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
	$senses{$protein} = $sense;
	$starts{$protein} = $start;
	$ends{$protein} = $end;
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
open(VIRUSES,"$genomad_virus_summary");
my %virus_score;
my %virus_fdr;
my %virus_coord;
my %virus_tax;
while(my $x = <VIRUSES>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $contig = $array[0];
	if($contig eq "seq_name"){
		next;
	}
	my $score = $array[6];
	my $fdr = $array[7];
	my $coord = $array[3];
	my $tax = $array[10];
	if($coord ne "NA"){
		my @c_array = split(/\|/,$contig);
		$contig = $c_array[0];
		print "$contig\n";
	} 
	if($score >= $score_threshold){
		$virus_score{$contig} = $score;
		$virus_fdr{$contig} = $fdr;
		$virus_coord{$contig} = $coord;
		$virus_tax{$contig} = $tax;
	}
}
close VIRUSES;
open(M8,"$m8_taxed");
my %total_observations;
my %virus_observations;
my %provirus_observations;
open(TABLE,">$out_table");
if($mode eq "SARG"){
	print TABLE "protein\tcontig\tstart\tend\tsense\tFS_cluster\tObserved_genes\tType\tMechanism\tsub1\tsub2\tspecial\tin_virus\tvirus_score\tvirus_fdr\tis_provirus\tcontig_taxonomy\tvirus_taxonomy\n";
}elsif($mode eq "VFDB"){
	print TABLE "protein\tcontig\tstart\tend\tsense\tFS_cluster\tObserved_genes\tVFC\tVF\tVF_full\tin_virus\tvirus_score\tvirus_fdr\tis_provirus\tcontig_taxonomy\tvirus_taxonomy\n";
}
open(FREQ,">$out_frequencies");
print FREQ "protein\ttotal_observations\tvirus_observations\tprovirus_observations\tvirus_freq\tprovirus_freq\n";
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
	if(!exists($virus_observations{$hit})){
                $virus_observations{$hit} = 0;
        }
	if(!exists($provirus_observations{$hit})){
                $provirus_observations{$hit} = 0;
        }
	$total_observations{$hit}++;
	my $in_virus = "FALSE";
	my $score = "NA";
	my $fdr = "NA";
	my $virus_tax = "NA";
	my $is_provirus = "FALSE";
	if(exists($virus_score{$contigs{$protein}})){
		if($virus_coord{$contigs{$protein}} ne "NA"){
			my($v_start,$v_end) = split(/-/,$virus_coord{$contigs{$protein}});
			if($v_start <= $starts{$protein} and $v_end >= $ends{$protein}){
				$virus_observations{$hit}++;
				$provirus_observations{$hit}++;
				$score = $virus_score{$contigs{$protein}};
				$fdr = $virus_fdr{$contigs{$protein}}; 
				$in_virus = "TRUE";
				$is_provirus = "TRUE";
				$virus_tax = $virus_tax{$contigs{$protein}};
			}
		}else{
			$virus_observations{$hit}++;
			$score = $virus_score{$contigs{$protein}};
			$fdr = $virus_fdr{$contigs{$protein}};
			$in_virus = "TRUE";
			$virus_tax = $virus_tax{$contigs{$protein}};
		}
	}
	my $data = $clusters{$FS};
	my @data_array = split(/\t/,$data);
	if($mode eq "SARG"){
		print TABLE "$hit\t$contigs{$protein}\t$starts{$protein}\t$ends{$protein}\t$senses{$protein}\t$FS\t$data_array[0]\t$data_array[2]\t$data_array[3]\t$data_array[4]\t$data_array[5]\t$data_array[6]\t$in_virus\t$score\t$fdr\t$is_provirus\t$taxonomy\t$virus_tax\n";
	}elsif($mode eq "VFDB"){
		print TABLE "$hit\t$contigs{$protein}\t$starts{$protein}\t$ends{$protein}\t$senses{$protein}\t$FS\t$data_array[0]\t$data_array[2]\t$data_array[5]\t$data_array[4]\t$in_virus\t$score\t$fdr\t$is_provirus\t$taxonomy\t$virus_tax\n";
	}
}
foreach my $h(sort keys %total_observations){
	my $v_freq = $virus_observations{$h}/$total_observations{$h};
	my $p_freq = $provirus_observations{$h}/$total_observations{$h};
	print FREQ "$h\t$total_observations{$h}\t$virus_observations{$h}\t$provirus_observations{$h}\t$v_freq\t$p_freq\n";
}
