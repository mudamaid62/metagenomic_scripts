#!/usr/bin/perl
use warnings;
use strict;

if(!exists($ARGV[0])){
        die "No arguments\n";
}

my $result = shift(@ARGV);
my $cluster_file = shift(@ARGV); #foldseek clusters with metadata
my $checked_file = shift(@ARGV); #checked foldseek clusters
my $not_list = shift(@ARGV); #list with proteins less than 30 aa
my $go_matrix = shift(@ARGV);
my $similarity_table = shift(@ARGV);
my $tmp_dir = shift(@ARGV);
my $go_terms_list = shift(@ARGV);

my %not_proteins;
open(NOT,"$not_list");
while(my $n = <NOT>){
        chomp($n);
        $not_proteins{$n} = "skip";
}
close NOT;
my @go_terms;
open(GO_TERMS,"$go_terms_list");
while(my $g = <GO_TERMS>){
	chomp($g);
	push @go_terms,$g;
}
my @sorted_go_terms = sort(@go_terms);
my $sorted_go_ref = \@sorted_go_terms;
my $go_header = join "\t",@sorted_go_terms;
close GO_TERMS;
#mkdir $tmp_dir;
my %clusters;
my %metadata;
my %clusters_list;
open(CLUSTERS,"$cluster_file");
while(my $k = <CLUSTERS>){
        chomp($k);
        my @k_array = split(/\t/,$k);
        my $cluster = shift(@k_array);
        if(!exists($clusters_list{$cluster})){
                $clusters_list{$cluster} = "added";
        }
        my $protein = shift(@k_array);
        my $meta = join "\t",@k_array;
        $clusters{$protein} = $cluster;
        $metadata{$protein} = $meta;
}
close CLUSTERS;
#foreach my $r(sort keys %clusters_list){
#	print "Creating $r\_annopro_data\n";
#	open(R_TEMP,">$tmp_dir/$r\_annopro_data");
#	print R_TEMP "Proteins,GO-terms,Scores";
#	close R_TEMP;
#}
#open(RESULT,"$result");
#my $z_count = 0;
#my $total_lines = 24417676;
#while(my $z = <RESULT>){
#	chomp($z);
#      	if($z eq "Proteins,GO-terms,Scores"){
#		next;
#	}
#	$z_count++;
#	my $percent = ($z_count/$total_lines)*100;
#	my($seq,$go,$score) = split(/\,/,$z);
#	if(exists($clusters{$seq})){
#		my $cluster = $clusters{$seq};
#		system "echo $z >> $tmp_dir/$cluster\_annopro_data";
#		print "$z_count/$total_lines --> $percent\%\n";
#	}
#}
open(CHECKED,"$checked_file");
my %checked_metadata;
while(my $f = <CHECKED>){
        chomp($f);
        my @f_array = split(/\t/,$f);
        my $f_cluster = shift(@f_array);
        my $f_metadata = join "\t",@f_array;
        $checked_metadata{$f_cluster} = $f_metadata;
}
close CHECKED;
my %used;
open(TABLE,">$similarity_table");
print TABLE "Protein_i\tProtein_j\tFunctional_similarity\tStructural_cluster_i\tDB_i\tSARG_type_i\tSARG_subtype_i\tVFDB_VF_i\tVFDB_VFC_i\tVFDB_Pathogen_i\tStructural_cluster_j\tDB_j\tSARG_type_j\tSARG_subtype_j\tVFDB_VF_j\tVFDB_VFC_j\tVFDB_Pathogen_j\n";
print STDERR "Structural_cluster\tCluster_type\tKGI_number\tSARG_number\tVFDB_number\tGeneral_type\tSpecific_type\tMin_func_sim\tMedian_func_sim\tMax_func_sim\n";
open(GO_MATRIX,">$go_matrix");
print GO_MATRIX "\t$go_header\n";

foreach my $q(sort keys %clusters_list){
        print "WORKING with cluster $q\n";
        my @q_array;
        foreach my $w(keys %clusters){
                if(exists($used{$w})){
                        next;
                }elsif($clusters{$w} eq $q){
                        if(exists($not_proteins{$w})){
                                $used{$w} = "yes";
                                print TABLE "$w\t-\tNA:less than 30 aa\t$clusters{$w}\t$metadata{$w}\t-\t-\t-\t-\t-\t-\t-\n";
                        }else{
                                push @q_array,$w;
                                $used{$w} = "yes";
                        }
                }
        }
        my @sorted = sort(@q_array);
	my $q_ref = \@sorted;
	my $file = "$q\_annopro_data";
	my %vectors = get_vectors($file,$tmp_dir,$sorted_go_ref,$q_ref);
	my %norms;
	my @test_vector = split(/\t/,$vectors{$sorted[0]});
	my $vector_length = (scalar(@test_vector)) - 1;
	foreach my $v(sort keys %vectors){
		print GO_MATRIX "$v\t$vectors{$v}\n";
		$norms{$v} = get_norm($vector_length, $vectors{$v});
	} 
	my $vec_ref = \%vectors;
	my $norm_ref = \%norms;
        my %similarities = get_similarities($vector_length,$vec_ref,$norm_ref,$q_ref);
        my %printed;
        my @values;
        my $sim_number = 0;
	foreach my $s(sort keys %similarities){
                my($i,$j) = split(/~/,$s);
                my $pair = "$i~$j";
                my $same = "$j~$i";
                if(!exists($printed{$pair}) or !exists($printed{$same})){
                        my $text = "$i\t$j\t$similarities{$s}\t$clusters{$i}\t$metadata{$i}\t$clusters{$j}\t$metadata{$j}\n";
                        push @values, $similarities{$s};
                        $printed{$pair} = "yes";
                        $printed{$same} = "yes";
                        print TABLE "$text";
                        $sim_number++;
                }
        }
	if(scalar(@values) != 0){
                my @sorted_values = sort{$a<=>$b}(@values);
                my $parity = $sim_number % 2;
                my $median = 0;
		my $min = 0;
		my $max = 0;
                if($parity == 0){
                        my $center_1 = $sim_number/2;
                        my $center_2 = $center_1 - 1;
                        $median = ($sorted_values[$center_1] + $sorted_values[$center_2])/2;
                }else{
                        my $center = int($sim_number/2);
                        $median = $sorted_values[$center];
                }
		if(scalar(@values) == 1){
                	$min = $sorted_values[0];
                	$max = $min;
		}else{
			$min = $sorted_values[0];
			$max = pop(@sorted_values);
		}
		print STDERR "$q\t$checked_metadata{$q}\t$min\t$median\t$max\n";
        }else{
		print STDERR "$q\t$checked_metadata{$q}\tNA\tNA\tNA\n";
        }
}
close TABLE;
close GO_MATRIX;

sub get_vectors{
	my ($file,$tmp_dir,$sorted_go_ref,$q_ref) = @_;
	my %go_terms;
	my %scores;
	my %vectors;
	open(RESULT,"$tmp_dir/$file");
	print "STEP 1: Build protein and GO-terms lists for file $file\n";
	while(my $x = <RESULT>){
        	chomp($x);
        	if($x eq "Proteins,GO-terms,Scores"){
                	next;
        	}
        	my($seq,$go,$score) = split(/\,/,$x);
		my $name = "$seq~$go";
        	$scores{$name} = $score;
	}
	my @sorted_proteins = @$q_ref;
	my @sorted_go_terms = @$sorted_go_ref;
	my $vector_length = scalar(@sorted_go_terms) - 1;
	print "STEP 2: Making protein GO-terms vectors for file $file\n";
	foreach my $p(@sorted_proteins){
        	my $vector_check = 0;
        	my @score_vector;
        	foreach my $g(@sorted_go_terms){
                	my $name = "$p~$g";
                	if(exists($scores{$name})){
                        	push @score_vector, $scores{$name};
                        	$vector_check += $scores{$name};
                        	#print "Score of $g for $p is $scores{$g}\n";
                	}else{
                        	push @score_vector, 0;
                        	#print "Score of $g for $p does not exist\n";
                	}
        	}
       		if($vector_check == 0){
                	die "ERROR: vector for $p is a zero-vector";
        	}
        	my $p_vector = join "\t",@score_vector;
        	$vectors{$p} = $p_vector;
	}
	return %vectors;
}
sub get_cosine{
        my($length,$i_norm,$j_norm,$i,$j) = @_;
        my @i_vector = split(/\t/,$i);
        my @j_vector = split(/\t/,$j);
        my $product_sum = 0;
        for my $k(0..$length){
                my $product = $i_vector[$k] * $j_vector[$k];
                $product_sum += $product;
        }
        my $cosine = ($product_sum/($i_norm * $j_norm));
        return $cosine;
}
sub get_norm{
        my ($length,$i) = @_;
        my @i_vector = split(/\t/,$i);
        my $i_square_sum = 0;
        for my $k(0..$length){
                my $i_square = ($i_vector[$k])**2;
                $i_square_sum += $i_square;
        }
        my $i_norm = ($i_square_sum)**(1/2);
        return $i_norm;
}
sub get_similarities{
        my ($vector_length,$vector_ref,$norm_ref,$sorted_proteins_ref) = @_;
        my %vectors = %$vector_ref;
        my %norms = %$norm_ref;
        my @sorted_proteins = @$sorted_proteins_ref;
        my %similarities;
        foreach my $i(@sorted_proteins){
                my @start = localtime(time);
                print "[$start[2]\:$start[1]\:$start[0]] Calculating similarities for $i\n";
                foreach my $j(@sorted_proteins){
                        my $pair = "$i~$j";
                        my $same = "$j~$i";
                        if(exists($similarities{$pair}) or exists($similarities{$same})){
                                next;
                        }else{
                                my $i_vector = $vectors{$i};
                                my $j_vector = $vectors{$j};
                                my $i_norm = $norms{$i};
                                my $j_norm = $norms{$j};
                                my $similarity = get_cosine($vector_length,$i_norm,$j_norm,$i_vector,$j_vector);
                                $similarities{$pair} = $similarity;
                                $similarities{$same} = $similarity;
                                #print "$pair cosine --> $similarity\n";
                        }
                }
        }
        return %similarities;
}
