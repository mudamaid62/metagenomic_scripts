#!/usr/bin/perl
use warnings;
use strict;

my $result = shift(@ARGV);
my $similarity_matrix = shift(@ARGV);
my $go_matrix = shift(@ARGV);
my $similarity_table = shift(@ARGV);

my %proteins;
my @protein_list;
my %go_terms;
my @go_terms_list;
my %scores;
my %vectors;

open(RESULT,"$result");
print STDERR "STEP 1: Build protein and GO-terms lists\n";
while(my $x = <RESULT>){
	chomp($x);
	if($x eq "Proteins,GO-terms,Scores"){
		next;
	}
	my($seq,$go,$score) = split(/\,/,$x);
	if(!exists($proteins{$seq})){
		$proteins{$seq} = "added";
		#print STDERR "Added protein $seq\n";
		push @protein_list,$seq;
	}
	if(!exists($go_terms{$go})){
		$go_terms{$go} = "added";
		#print STDERR "Added GO-term $go\n";
		push @go_terms_list,$go;
	}
	my $name = "$seq~$go";
	$scores{$name} = $score;
}
my @sorted_proteins = sort(@protein_list);
my @sorted_go_terms = sort(@go_terms_list);
my $protein_header = join "\t",@sorted_proteins;
my $go_header = join "\t",@sorted_go_terms;
print STDERR "STEP 2: Making protein GO-terms vectors\n";
foreach my $p(@sorted_proteins){
	my @score_vector;
	foreach my $g(@sorted_go_terms){
		my $name = "$p~$g";
		if(exists($scores{$name})){
			push @score_vector, $scores{$name};
			#print STDERR "Score of $g for $p is $scores{$g}\n";
		}else{
			push @score_vector, 0;
			#print STDERR "Score of $g for $p does not exist\n";
		}
	}
	my $p_vector = join "\t",@score_vector;
	$vectors{$p} = $p_vector;
	print STDERR "$p vector done\n";
}
close RESULT;
open(GO_MATRIX,">$go_matrix");
open(SIMILARITY_MATRIX, ">$similarity_matrix");
print GO_MATRIX "\t$go_header\n";
foreach my $v(sort keys %vectors){
	print GO_MATRIX "$v\t$vectors{$v}\n";
}
close GO_MATRIX;
open(SIMILARITY_TABLE,">$similarity_table");
print STDERR "STEP 3: Getting cosine-similarity distances\n";
my %similarities;
foreach my $i(@sorted_proteins){
	print STDERR "Calculating similarities for $i\n";
	foreach my $j(@sorted_proteins){
		my $pair = "$i~$j";
		my $same = "$j~$i";
		if(exists($similarities{$pair}) or exists($similarities{$same})){
			next;
		}else{
			my $i_vector = $vectors{$i};
			my $j_vector = $vectors{$j};
			my $similarity = get_cosine($i_vector,$j_vector);
			$similarities{$pair} = $similarity;
			$similarities{$same} = $similarity;
		}
		#print STDERR "$pair cosine similarity --> $similarities{$pair}\n"; 
	}
}
print SIMILARITY_MATRIX "\t$protein_header\n";
print SIMILARITY_TABLE "Protein_i\tProtein_j\tFunctional_Similarity\n";
my %check_print;
foreach my $k(@sorted_proteins){
	my @k_vector;
	foreach my $l(@sorted_proteins){
		my $pair = "$k~$l";
		my $same = "$l~$k";
		if(!exists($check_print{$pair}) or !exists($check_print{$same})){
			print SIMILARITY_TABLE "$k\t$l\t$similarities{$pair}\n";
			$check_print{$pair} = "printed";
			$check_print{$same} = "printed";
		}
		push @k_vector, $similarities{$pair};
	}
	my $k_line = join "\t",@k_vector;
	print SIMILARITY_MATRIX "$k\t$k_line\n";
}
close SIMILARITY_MATRIX;
close SIMILARITY_TABLE;

sub get_cosine{
	my($i,$j) = @_;
	my @i_vector = split(/\t/,$i);
	my @j_vector = split(/\t/,$j);
	my $length = (scalar(@i_vector)) - 1;
	my $i_square_sum = 0;
	my $j_square_sum = 0;
	my $product_sum = 0;
	for my $k(0..$length){
		my $i_square = ($i_vector[$k])**2;
		my $j_square = ($j_vector[$k])**2;
		my $product = $i_vector[$k] * $j_vector[$k];
		$i_square_sum += $i_square;
		$j_square_sum += $j_square;
		$product_sum += $product;
	}
	my $i_norm = ($i_square_sum)**(1/2);
	my $j_norm = ($j_square_sum)**(1/2);
	my $cosine = ($product_sum/($i_norm * $j_norm));
	return $cosine;
}
