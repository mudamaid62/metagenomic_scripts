#!/usr/bin/perl
use warnings;
use strict;

#STDOUT redirects centroids list, members_number, distribution of identities observed in the cluster
my $blast = shift(@ARGV); #query target pident bitscore
my $table = shift(@ARGV); #protein ref(yes/no) cluster
my $out_table = shift(@ARGV);

my %proteins;
my %cluster_list;
my %is_ref;
open(TABLE,"$table");
while(my $x = <TABLE>){
	chomp($x);
	my($protein,$ref,$cluster) = split(/\t/,$x);
	$proteins{$protein} = $cluster;
	$is_ref{$protein} = $ref;
	if(!exists($cluster_list{$cluster})){
		$cluster_list{$cluster} = "added";
	}
}
close TABLE;
open(OUT,">$out_table");
print OUT "Protein\tReference\tFS_cluster\tpident_modifier\n";
my %pairs;
my %p_bits;
open(BLAST,"$blast");
while(my $z = <BLAST>){
	chomp($z);
	my($query,$target,$pident,$bits) = split(/\t/,$z);
	my $pair = "$query~$target";
	my $frac = $pident/100;
	$pairs{$pair} = $frac;
	$p_bits{$pair} = $bits;
}
close BLAST;
print "FS_cluster\tCentroid_protein\tNumber_of_members\tMin_observed_id\tP25\tMedian_observed_id\tP75\tMax_observed_id\tAverage_observed_id\tStandard_deviation\n";
foreach my $c(sort keys %cluster_list){
	my @identities;
	my @members;
	foreach my $p(keys %proteins){
		if($proteins{$p} eq $c){
			push @members,$p;
		}
	}
	my $member_number = scalar(@members);
	if($c eq "SINGLE_KGI"){
		foreach my $p(@members){
			print OUT "$p\t$is_ref{$p}\t$c\t0\n";
		}
	}elsif($member_number == 1){
		foreach my $p(@members){
			print OUT "$p\t$is_ref{$p}\t$c\t1\n";
			print "$c\t$p\t1\t1\t1\t1\t1\t1\t1\t0\n";
		}
	}elsif($member_number <= 0){
		die "ERROR: $c cluster has 0 members\n";
	}else{
		my %sum;
		foreach my $v(@members){
			my $v_sum = 0;
			foreach my $w(@members){
				if($v eq $w){
					next;
				}
				my $name = "$v~$w";
				my $same = "$w~$v";
				my $bits = 0;
				my $id = 0;
				if(exists($pairs{$name}) and $p_bits{$name} > $bits){
					$bits = $p_bits{$name};
					$id = $pairs{$name};
				}
				if(exists($pairs{$same}) and $p_bits{$same} > $bits){
                                        $bits = $p_bits{$same};
					$id = $pairs{$same};
				}
				push @identities, $id;
				$v_sum += $bits;
			}
			$sum{$v} = $v_sum;
		}
		my $best_sum = 0;
		my $centroid = "placeholder";
		foreach my $s(keys %sum){
			if($best_sum == 0 and $centroid eq "placeholder"){
				$best_sum = $sum{$s};
				$centroid = $s;
				print STDERR "for $c best is $centroid with $best_sum\n"; 
			}elsif($sum{$s} > $best_sum){
				print STDERR "for $c, $s with sum $sum{$s} better than $best_sum\n";
				$best_sum = $sum{$s};
				$centroid = $s;
			}
		}
		foreach my $q(@members){
			my $pair = "$q~$centroid";
			my $same = "$centroid~$q";
			my $value = 0;
			if(exists($pairs{$pair}) and $pairs{$pair} > $value){
				$value = $pairs{$pair};
			}
			if(exists($pairs{$same}) and $pairs{$same} > $value){
                                $value = $pairs{$same};
			}
			print OUT "$q\t$is_ref{$q}\t$c\t$value\n";
		}
		my @sorted_ids = sort{$a<=>$b} @identities;
                my $sim_number = scalar(@sorted_ids);
                my $parity = $sim_number % 2;
                my $median = 0;
                my $min = 0;
                my $max = 0;
                if($parity == 0){
                        my $center_1 = $sim_number/2;
                        my $center_2 = $center_1 - 1;
                        $median = ($sorted_ids[$center_1] + $sorted_ids[$center_2])/2;
                }else{
                        my $center = int($sim_number/2);
                        $median = $sorted_ids[$center];
                }
		my $avg_sum = 0;
		foreach my $i(@sorted_ids){
			$avg_sum += $i;
		}
		my $mean = $avg_sum/$sim_number;
		my $squared_sum = 0;
		foreach my $i(@sorted_ids){
			my $square_dis = ($i - $mean)**2;
			$squared_sum += $square_dis;
		}
		my $std = ($squared_sum / $sim_number)**(1/2);
		if($sim_number == 1){
                        $min = $sorted_ids[0];
                        $max = $min;
                }else{
                        $min = $sorted_ids[0];
                        $max = pop(@sorted_ids);
                }
		my $cat = join "\t",@sorted_ids;
		my $first_quartile = get_percentile(25,$cat);
		my $third_quartile = get_percentile(75,$cat);
		print "$c\t$centroid\t$member_number\t$min\t$first_quartile\t$median\t$third_quartile\t$max\t$mean\t$std\n";
	}
}

sub get_percentile{
        my ($number,$values) = @_;
        my @array = split(/\t/,$values);
        my $len = scalar(@array);
        my $rank = int(($number/100)*$len);
        my $index = $rank - 1;
        my $value = $array[$index];
        return $value;
}
