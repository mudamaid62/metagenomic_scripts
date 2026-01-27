#!/usr/bin/perl
use warnings;
use strict;

my $metadata_table = shift(@ARGV); #Protein Reference FS_cluster Redundant ab_string
my $abundance_table = shift(@ARGV);
my $diamond_table = shift(@ARGV);
my $consolidation_identity = shift(@ARGV);
my $consolidation_coverage = shift(@ARGV);

my %priority;
my %abundance;
my %bitscore;
my %full_register;
my @all_proteins;

open(ALL,"$metadata_table");
print STDERR "STEP 0: Setting database\n";
while(my $x = <ALL>){
	chomp($x);
	my @array = split(/\t/,$x);
	if($array[0] eq "Protein"){
		next;
	}
	$priority{$array[0]} = $array[1];
	if(!exists($array[3])){
		push @array, "-";
	}
	my $protein = $array[0];
	my $register = "$array[1]\t$array[2]\t$array[3]";
	push @all_proteins, $protein;
	$full_register{$protein} = $register;
	if(!exists($array[4])){
		$abundance{$protein} = "placeholder";
	}else{
		my @ab_array;
		for my $i(4..14){
			push @ab_array, $array[$i];
		}
		my $ab_string = join "\t",@ab_array;
		$abundance{$protein} = $ab_string;
	}
	#print STDERR "$protein --> $full_register{$protein}\n$abundance{$protein}\n";
}
close ALL;
open(ABUNDANCE, "$abundance_table");
while(my $x = <ABUNDANCE>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $protein = shift(@array);
	my $ab_string = join "\t", @array;
	if(exists($abundance{$protein}) and $abundance{$protein} eq "placeholder"){
		$abundance{$protein} = $ab_string;
	}
}
close ABUNDANCE;
my %consolidation_hits;
my %used;
my @diamond_hits;
open(DIAMOND,"$diamond_table");
while(my $j = <DIAMOND>){
	chomp($j);
	my @j_array = split(/\t/,$j);
	my $pair = "$j_array[0]~$j_array[2]";
	my $pair_bitscore = $j_array[9];
	#$bitscore{$pair} = $pair_bitscore;
	if($j_array[10] >= $consolidation_identity and $j_array[11] >= $consolidation_coverage and $j_array[12] >= $consolidation_coverage and $j_array[0] ne $j_array[2]){
		push @diamond_hits, $j;
		$bitscore{$pair} = $pair_bitscore;
	}
}
close DIAMOND;
my $separator = "-" x 30;
my $protein_number = scalar(@all_proteins);
my $i = 0;
foreach my $z(@all_proteins){
	if($abundance{$z} eq "placeholder"){
		$abundance{$z} = "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
	}
	if(exists($used{$z})){
		next;
	}
	$i++;
	my $percent = int(($i/$protein_number)*100);
	print STDERR "$separator\nSTEP 1: checking $z hits ($i / $protein_number -> $percent \%)\n";
	my @z_array;
	foreach my $y(@diamond_hits){
		my @y_array = split(/\t/,$y);
		if(exists($used{$y_array[2]})){
			next;
		}
		if($y_array[0] eq $z){
			push @z_array, $y_array[2];
			$i++;
			$used{$y_array[2]} = "yes";
		}
	}
	my $z_value = scalar(@z_array);
	if($z_value > 0){
		my $selected = "";
		push @z_array, $z;
                $used{$z} = "yes";
		my @winners;
		my %values;
		my $winner_score = 0;
		foreach my $u(@z_array){
			my $u_value = 0;
			foreach my $v(@z_array){
				if(exists($bitscore{"$u~$v"})){
					$u_value += $bitscore{"$u~$v"};
				}			
			}
			$values{$u} = $u_value;
			print STDERR "$u --> $values{$u}\n";
			if($u_value >= $winner_score){
				$winner_score = $u_value;
			}
		}
		foreach my $w(@z_array){
			if($values{$w} == $winner_score){
				push @winners,$w;
			}
		}
		my @prioritized;
		my @sorted = sort(@winners);
		foreach my $h(@sorted){
			if($priority{$h} eq "yes"){
				push @prioritized, $h;
			}
		}
		if(scalar(@prioritized > 0)){
			my @sort_p = sort(@prioritized);
			$selected = shift(@sort_p);
		}else{
			$selected = shift(@sorted);
		}
		my $z_group = join "\t",@z_array;
		$consolidation_hits{$selected} = $z_group;
		print STDERR "For $z, selected $selected\n";
	}else{
		print STDERR "$z is unique\n";
	}
}
print "Protein\tReference\tFS_cluster\tRedundant names\tSCp1\tSCp2\tSCp3\tsemifrozen\tMp1\tMp3\tmoraine\tRp1\tRp2\tRp3\trizo\n";
my $group_number = 0;
my $consolidated_proteins = 0;

foreach my $q(sort keys %consolidation_hits){
	$group_number++;
	print STDERR "$separator\nSTEP 2: Group $q\n";
	my @q_array = split(/\t/,$consolidation_hits{$q});
	my $SCp1 = 0;
	my $SCp2 = 0;
	my $SCp3 = 0;
	my $semifrozen = 0;
	my $Mp1 = 0;
	my $Mp3 = 0;
	my $moraine = 0;
	my $Rp1 = 0;
	my $Rp2 = 0;
	my $Rp3 = 0;
	my $rizo = 0;
	my @redundant_array;
	foreach my $w(@q_array){
		my @site_ab = split(/\t/,$abundance{$w});
		print STDERR "$w --> $abundance{$w}\n";
		$SCp1 += $site_ab[0];
		$SCp2 += $site_ab[1];
		$SCp3 += $site_ab[2];
		$semifrozen += $site_ab[3];
		$Mp1 += $site_ab[4];
		$Mp3 += $site_ab[5];
		$moraine += $site_ab[6];
		$Rp1 += $site_ab[7];
		$Rp2 += $site_ab[8];
		$Rp3 += $site_ab[9];
		$rizo += $site_ab[10];
		my @w_register = split(/\t/, $full_register{$w});
		if($w_register[2] ne "-"){
			print STDERR "$w --> $w_register[2]\n";
			my @w_array = split(/\;/,$w_register[2]);
			foreach my $v(@w_array){
				push @redundant_array, $v;
			}
		}
		if($w ne $q){
			push @redundant_array, $w;
		}
	}
	my @sorted = sort(@redundant_array);
	my $final_name = $q;
	my @new_array;
	if($priority{$q} eq "no"){
		push @redundant_array, $q;
		foreach my $s(@sorted){
			if(exists($priority{$s}) and $priority{$s} eq "yes"){
				$final_name = $s;
				print STDERR "Old rep $q was replaced by $final_name\n";
				last;
			}
		}
		@sorted = ();
		@sorted = sort(@redundant_array);
		foreach my $s(@sorted){
			if($s ne $final_name){
				push @new_array, $s;
			}
		}
	}else{
		@new_array = @sorted;
	}
	$consolidated_proteins += scalar(@new_array);
	$consolidated_proteins++;
	my $extra = join "\;", @new_array;
	my ($reference,$FS_cluster,$useless) = split(/\t/,$full_register{$q});
	my $final_register = "$priority{$final_name}\t$FS_cluster\t$extra";
	print STDERR "Total $final_name group abundance --> $SCp1\t$SCp2\t$SCp3\t$semifrozen\t$Mp1\t$Mp3\t$moraine\t$Rp1\t$Rp2\t$Rp3\t$rizo\n";
	print "$final_name\t$final_register\t$SCp1\t$SCp2\t$SCp3\t$semifrozen\t$Mp1\t$Mp3\t$moraine\t$Rp1\t$Rp2\t$Rp3\t$rizo\n";
}
my $previously_consolidated = 0;
my $uniques = 0;
foreach my $n(@all_proteins){
	if(!exists($used{$n})){
		my @n_register = split(/\t/, $full_register{$n});
		if($n_register[2] ne "-"){
			print STDERR "STEP 3: $n was already consolidated in a previous iteration\n";
                        print STDERR "$n --> $n_register[2]\n";
                        my @redundant_array = split(/\;/,$n_register[2]);
			$previously_consolidated += scalar(@redundant_array);
			$previously_consolidated++;
			print "$n\t$full_register{$n}\t$abundance{$n}\n";
                }else{
			print STDERR "STEP 3: $n was never used\n";
			print "$n\t$full_register{$n}\t$abundance{$n}\n";
			$uniques++;			
		}
	}
}	
my $total_proteins = $consolidated_proteins + $previously_consolidated + $uniques;
print STDERR "$separator\nSUMMARY\nConsolidation Groups -> $group_number\nConsolidated Proteins -> $consolidated_proteins\nPreviously Consolidated Proteins --> $previously_consolidated\nUnique Proteins -> $uniques\nTotal Proteins -> $total_proteins\n$separator\n";
