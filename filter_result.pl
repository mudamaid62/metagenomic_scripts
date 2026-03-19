#!/usr/bin/perl
use warnings;
use strict;

my $result = shift(@ARGV); #input FS_table
my $filtered_table = shift(@ARGV); #output FS_table
my $filtered_list = shift(@ARGV); #output list of proteins in clusters
my $filtered_log = shift(@ARGV); 
my $threshold = shift(@ARGV);
my $single_protein_list = shift(@ARGV); #output list of single proteins
my $next_iteration_list = shift(@ARGV); #output list of proteins excluded from clusters

my %clusters;
my %proteins;
my %similarities;
my %meta;
my $KGI_singles = 0;
open(RESULT,"$result");
open(SINGLE,">$single_protein_list");

while(my $x = <RESULT>){
	chomp($x);
	my @array = split(/\t/,$x);
	if($array[0] eq "Protein_i"){
		next;
	}
	if($array[2] eq "NA\:less than 30 aa"){
		if(!exists($proteins{$array[0]})){
			if($array[0] =~ m/SARG_/ or $array[0] =~ m/VFG/){	
				print SINGLE "$array[0]\n";
				$proteins{$array[0]} = "single";
			}else{
				$proteins{$array[0]} = "single";
				$KGI_singles++;
			}
		}
		next;
	}
	if(!exists($clusters{$array[3]})){
		$clusters{$array[3]} = "$array[0]";
		$proteins{$array[0]} = $array[3];
	}
	if(!exists($proteins{$array[0]})){
		$proteins{$array[0]} = $array[3];
		$clusters{$array[3]} .= "~$array[0]";
	}
	if(!exists($proteins{$array[1]})){
                $proteins{$array[1]} = $array[3];
                $clusters{$array[3]} .= "~$array[1]";
        }
	my $pair = "$array[0]~$array[1]";
	my $same = "$array[1]~$array[0]";
	$similarities{$pair} = $array[2];
	$similarities{$same} = $array[2];
	my $i_data = "$array[3]~$array[4]~$array[5]~$array[6]~$array[7]~$array[8]~$array[9]";
	my $j_data = "$array[10]~$array[11]~$array[12]~$array[13]~$array[14]~$array[15]~$array[16]";
	$meta{$array[0]} = $i_data;
	$meta{$array[1]} = $j_data;
}
close RESULT;
open(TABLE,">$filtered_table");
open(LIST,">$filtered_list");
open(LOG, ">$filtered_log");
open(NEXT,">$next_iteration_list");
print TABLE "Protein_i\tProtein_j\tFunctional_similarity\tStructural_cluster_i\tDB_i\tSARG_type_i\tSARG_subtype_i\tVFDB_VF_i\tVFDB_VFC_i\tVFDB_Pathogen_i\tStructural_cluster_j\tDB_j\tSARG_type_j\tSARG_subtype_j\tVFDB_VF_j\tVFDB_VFC_j\tVFDB_Pathogen_j\n";
print STDERR "Structural_cluster\tCluster_type\tKGI_number\tSARG_number\tVFDB_number\tTOTAL\tGeneral_type\tSpecific_type\tMin_func_sim\tMedian_func_sim\tMax_func_sim\n";

my $sep = "-" x 20;
foreach my $z(sort keys %clusters){
	print "$sep\n";
	print LOG "$sep\n";
	my @cluster = split(/~/,$clusters{$z});
	if(scalar(@cluster) == 1){
		print "Cluster $z only has 1 protein\nCluster $z will be DISCARDED\n";
		print LOG "Cluster $z only has 1 protein\nCluster $z will be DISCARDED\n";
		if($cluster[0] =~ m/SARG_/ or $cluster[0] =~ m/VFG/){
			print SINGLE "$cluster[0]\n";
		}else{
			$KGI_singles++;
		}
		next;
	}
	my @cluster_distances;
        my %distance_sum;
	my %used;
	foreach my $p(@cluster){
		my $sum = 0;
		foreach my $q(@cluster){
			if($p eq $q){
				next;
			}
			my $name = "$p~$q";
			my $same = "$q~$p";
			$sum += $similarities{$name};
			if(!exists($used{$name}) or !exists($used{$same})){
				push @cluster_distances, $similarities{$name};
				$used{$name} = "yes";
				$used{$same} = "yes";
			}
		}
		$distance_sum{$p} = $sum;
	}
	my @sorted_distances = sort{$a<=>$b}(@cluster_distances);
	my $current_min = $sorted_distances[0];
	print "For cluster $z, current min is $current_min\n";
	print LOG "For cluster $z, current min is $current_min\n";
	if($current_min < $threshold and scalar(@cluster) <= 2){
		my $elements = scalar(@cluster);
		print "For cluster $z current minimal functional similarity $current_min is less than $threshold, but $z only has $elements proteins\n";
		print LOG "For cluster $z current minimal functional similarity $current_min is less than $threshold, but $z only has $elements proteins\n";
		print "Cluster $z will be DISCARDED\n";
		print LOG "Cluster $z will be DISCARDED\n";
		if($cluster[0] =~ m/SARG_/ or $cluster[0] =~ m/VFG/){
			print SINGLE "$cluster[0]\n";
		}else{
			$KGI_singles++;
		}
		if($cluster[1] =~ m/SARG_/ or $cluster[1] =~ m/VFG/){
                        print SINGLE "$cluster[1]\n";
                }else{
                        $KGI_singles++;
		}
		next;
	}
	until($current_min >= $threshold or scalar(@cluster) <= 2){
		my @kept;
		my $min_sum = 0;
		my $current = "placeholder";
		foreach my $k(keys %distance_sum){
			#print STDERR "For $k distance sum is $distance_sum{$k}\n";
			if($min_sum == 0 and $current eq "placeholder"){
				$min_sum = $distance_sum{$k};
				$current = $k;
			}elsif($distance_sum{$k} < $min_sum){
				$min_sum = $distance_sum{$k};
                                $current = $k;		
			}
		}
		print "Will try filtering $current\n";
		print LOG "Will try filtering $current\n";
		print NEXT "$current\n";
		foreach my $c(@cluster){
			if($c eq $current){
				next;
			}else{
				push @kept,$c;
			}
		}
		@cluster = ();
        	@cluster_distances = ();
        	%distance_sum = () ;
        	%used = () ;
		@sorted_distances = ();
		foreach my $r(@kept){
			push @cluster,$r;
		}
		foreach my $p(@cluster){
                	my $sum = 0;
                	foreach my $q(@cluster){
                        	if($p eq $q){
                                	next;
                        	}
                        	my $name = "$p~$q";
                        	my $same = "$q~$p";
                        	$sum += $similarities{$name};
                        	if(!exists($used{$name}) or !exists($used{$same})){
                                	push @cluster_distances, $similarities{$name};
                                	$used{$name} = "yes";
                                	$used{$same} = "yes";
                        	}
                	}
                	$distance_sum{$p} = $sum;
        	}
        	@sorted_distances = sort{$a<=>$b}(@cluster_distances);
        	$current_min = $sorted_distances[0];
		print "For cluster $z, current min is $current_min\n";
		print LOG "For cluster $z, current min is $current_min\n";
	}
	if($current_min < $threshold and scalar(@cluster) <= 2){
        	my $elements = scalar(@cluster);
        	print "For cluster $z current minimal functional similarity $current_min is less than $threshold, but $z only has $elements proteins\n";
        	print LOG "For cluster $z current minimal functional similarity $current_min is less than $threshold, but $z only has $elements proteins\n";
        	print "Cluster $z will be DISCARDED\n";
        	print LOG "Cluster $z will be DISCARDED\n";
		if($cluster[0] =~ m/SARG_/ or $cluster[0] =~ m/VFG/){
                        print SINGLE "$cluster[0]\n";
                }else{
                        $KGI_singles++;
                }
                if($cluster[1] =~ m/SARG_/ or $cluster[1] =~ m/VFG/){
                        print SINGLE "$cluster[1]\n";
                }else{
                        $KGI_singles++;
                }
                next;
        }
	print "Current minimal functional similarity $current_min is greater or equal than $threshold\n";
	print LOG "Current minimal functional similarity $current_min is greater or equal than $threshold\n";
	print "Filtering COMPLETE\n";
	print LOG "Filtering COMPLETE\n";
	my %printed;
	my $KGI_number = 0;
	my $SARG_number = 0;
	my $VFDB_number = 0;
	my @general_types;
	my @specific_types;
	foreach my $w(@cluster){
		print LIST "$w\n";
		my @w_meta = split(/~/,$meta{$w});
		if($w_meta[1] eq "KGI"){
			$KGI_number++;
			my $dummy = "placeholder";
			push @general_types,$dummy;
                        push @specific_types,$dummy;
		}elsif($w_meta[1] eq "SARG"){
			$SARG_number++;
			push @general_types,$w_meta[2];
			push @specific_types,$w_meta[3];
		}elsif($w_meta[1] eq "VFDB"){
			$VFDB_number++;
			push @general_types,$w_meta[5];
			push @specific_types,$w_meta[4];
		}
		foreach my $v(@cluster){
			if($w eq $v){
				next;
			}
			my @v_meta = split(/~/,$meta{$v});
			my $name = "$w~$v";
			my $same = "$v~$w";
			if(!exists($printed{$name}) or !exists($printed{$same})){
				my $w_printed = join "\t",@w_meta;
				my $v_printed = join "\t",@v_meta;
				print TABLE "$w\t$v\t$similarities{$name}\t$w_printed\t$v_printed\n";
				$printed{$name} = "yes";
				$printed{$same} = "yes";
			}
		}
	}
	my $total = $KGI_number + $SARG_number + $VFDB_number;
	my $g_conflict = "no";
	my $s_conflict = "no";
	my @g_array;
	my @s_array;
	my %g_check;
	my %s_check;
	my $g_number = 0;
	my $s_number = 0;
	my $g_type = "-";
	my $s_type = "-";
	foreach my $g(@general_types){
		if(!exists($g_check{$g})){
			$g_check{$g} = "added";
			unless($g eq "placeholder"){
				push @g_array,$g;
				$g_number++;
			}
		}
	}
	if($g_number > 0){
		$g_type = join ";",@g_array;
	}
	if($g_number > 1){
		$g_conflict = "yes";
	}
	foreach my $s(@specific_types){
                if(!exists($s_check{$s})){
                        $s_check{$s} = "added";
                        unless($s eq "placeholder"){
                                push @s_array,$s;
                                $s_number++;
                        }
                }
        }
        if($s_number > 0){
                $s_type = join ";",@s_array;
        }
        if($s_number > 1){
                $s_conflict = "yes";
        }
	my $cluster_type = "none";
	if($s_conflict eq "no" and $KGI_number > 0 and $SARG_number > 0 and $VFDB_number > 0){
		$cluster_type = "Specific_KGI_SARG_VFDB";
	}elsif($s_conflict eq "no" and $KGI_number > 0 and $SARG_number > 0 and $VFDB_number == 0){ 	
		$cluster_type = "Specific_KGI_SARG";
	}elsif($s_conflict eq "no" and $KGI_number > 0 and $SARG_number == 0 and $VFDB_number > 0){
		$cluster_type = "Specific_KGI_VFDB";
	}elsif($s_conflict eq "no" and $KGI_number == 0 and $SARG_number > 0 and $VFDB_number > 0){
                $cluster_type = "Specific_SARG_VFDB";
	}elsif($s_conflict eq "no" and $KGI_number == 0 and $SARG_number > 0 and $VFDB_number == 0){
                $cluster_type = "Specific_SARG";
	}elsif($s_conflict eq "no" and $KGI_number == 0 and $SARG_number == 0 and $VFDB_number > 0){
                $cluster_type = "Specific_VFDB";
	}elsif($s_conflict eq "yes" and $g_conflict eq "no" and $KGI_number > 0 and $SARG_number > 0 and $VFDB_number > 0){
                $cluster_type = "General_KGI_SARG_VFDB";
	}elsif($s_conflict eq "yes" and $g_conflict eq "no" and $KGI_number > 0 and $SARG_number > 0 and $VFDB_number == 0){
                $cluster_type = "General_KGI_SARG";
	}elsif($s_conflict eq "yes" and $g_conflict eq "no" and $KGI_number > 0 and $SARG_number == 0 and $VFDB_number > 0){
                $cluster_type = "General_KGI_VFDB";
	}elsif($s_conflict eq "yes" and $g_conflict eq "no" and $KGI_number == 0 and $SARG_number > 0 and $VFDB_number > 0){
                $cluster_type = "General_SARG_VFDB";
	}elsif($s_conflict eq "yes" and $g_conflict eq "no" and $KGI_number == 0 and $SARG_number > 0 and $VFDB_number == 0){
                $cluster_type = "General_SARG";
	}elsif($s_conflict eq "yes" and $g_conflict eq "no" and $KGI_number == 0 and $SARG_number == 0 and $VFDB_number > 0){
                $cluster_type = "General_VFDB";
	}elsif($g_conflict eq "yes" and $KGI_number > 0 and $SARG_number > 0 and $VFDB_number > 0){
		$cluster_type = "Conflicting_KGI_SARG_VFDB";
	}elsif($g_conflict eq "yes" and $KGI_number > 0 and $SARG_number > 0 and $VFDB_number == 0){
                $cluster_type = "Conflicting_KGI_SARG";
	}elsif($g_conflict eq "yes" and $KGI_number > 0 and $SARG_number == 0 and $VFDB_number > 0){
                $cluster_type = "Conflicting_KGI_VFDB";
	}elsif($g_conflict eq "yes" and $KGI_number == 0 and $SARG_number > 0 and $VFDB_number > 0){
                $cluster_type = "Conflicting_SARG_VFDB";
	}elsif($g_conflict eq "yes" and $KGI_number == 0 and $SARG_number > 0 and $VFDB_number == 0){
                $cluster_type = "Conflicting_SARG";
	}elsif($g_conflict eq "yes" and $KGI_number == 0 and $SARG_number == 0 and $VFDB_number > 0){
                $cluster_type = "Conflicting_VFDB";
	}
	my $distances_number = scalar(@sorted_distances);
	my $parity = $distances_number % 2;
	my $half = $distances_number / 2;
	my $median = 0;
	if($parity == 0){
		my $pre = $half - 1;
		$median = ($sorted_distances[$half] + $sorted_distances[$pre])/2;
	}else{
		my $int = int($half);
		$median = $sorted_distances[$int];
	}
	my $max = pop(@sorted_distances);
	print STDERR "$z\t$cluster_type\t$KGI_number\t$SARG_number\t$VFDB_number\t$total\t$g_type\t$s_type\t$current_min\t$median\t$max\n";
}
print "KGI single proteins --> $KGI_singles\n";
print LOG "KGI single proteins --> $KGI_singles\n";
close TABLE;
close LIST;
close LOG;
close SINGLE;
close NEXT;
