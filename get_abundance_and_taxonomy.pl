#!/usr/bin/perl
use warnings;
use strict;

my $smf_file = shift(@ARGV);
my $rat_r2c = shift(@ARGV);
my $m8_file = shift(@ARGV);
my $otu_table = shift(@ARGV);

my $cell_number = 0;
open(SMF,"$smf_file");
while(my $x = <SMF>){
	chomp($x);
	my @array = split(/\t/,$x);
	if($array[1] eq "bacterial_archaeal_bases"){
		next;
	}else{
		$cell_number = $array[1]/$array[4];
	}
}
close SMF;
my %reads;
open(RAT,"$rat_r2c");
while(my $y = <RAT>){
	chomp($y);
	if($y =~ m/^\#/){
		next;
	}
	my @y_array = split(/\t/,$y);
	my $current_tax = "no taxid assigned";
	if($y_array[1] eq "no taxid assigned"){
		$current_tax = "no taxid assigned";
	}elsif($y_array[2] ne ""){
		$current_tax = $y_array[2];
	}elsif($y_array[3] ne ""){
                $current_tax = $y_array[3];
        }elsif($y_array[4] ne ""){
                $current_tax = $y_array[4];
        }else{
		my @tax_array = split(/ /,$y_array[5]);
		my $useless = shift(@tax_array);
		my $tax = join " ", @tax_array;
		$current_tax = $tax;
	}
	if(!exists($reads{$y_array[0]})){
		$reads{$y_array[0]} = $current_tax;
	}else{
		my $tax_to_solve = "$reads{$y_array[0]}~$current_tax";
		$reads{$y_array[0]} = solve_tax($tax_to_solve);
	}
}
close RAT;
open(MM,"$m8_file");
my %best_hit_bitscore;
my %best_hit;
my %valid_reads;
my %selected_targets;
my %abundances;
my %counts;
print "Protein\tCounts\tAbundance (copies/cell)\tBest_taxonomy\n";
while(my $z = <MM>){
	chomp($z);
	my($query,$target,$pident,$qcov,$tcov,$evalue,$bits,$qlen,$tlen,$alnlen) = split(/\t/,$z);
	if($evalue <= 1e-7 and $pident >= 80 and $qcov >= 0.75 and $alnlen >= 25){
		$valid_reads{$query} = $qlen;
		if(!exists($selected_targets{$target})){
			$abundances{$target} = 0;
			$counts{$target} = 0;
			$selected_targets{$target} = ($tlen*3)+3;
		}
		if(!exists($best_hit{$query})){
			$best_hit{$query} = $target;
			$best_hit_bitscore{$query} = $bits;
		}elsif($bits > $best_hit_bitscore{$query}){
			$best_hit{$query} = $target;
			$best_hit_bitscore{$query} = $bits;
		}
	}
}		
close MM;
my %tax_to_solve;
my $spaciator = "-" x 30;
my %observed_taxonomies;
foreach my $q(keys %valid_reads){
        my $value = ($valid_reads{$q})/($selected_targets{$best_hit{$q}});
        $abundances{$best_hit{$q}} += $value;
        $counts{$best_hit{$q}}++;
        my $pair = "$best_hit{$q}~$reads{$q}";
        if(!exists($observed_taxonomies{$pair})){
                $observed_taxonomies{$pair} = $value;
        }else{
                $observed_taxonomies{$pair} += $value;
        }
        if(!exists($tax_to_solve{$best_hit{$q}})){
                $tax_to_solve{$best_hit{$q}} = $reads{$q};
        }else{
                my $current = $tax_to_solve{$best_hit{$q}};
                my $new = $reads{$q};
                $tax_to_solve{$best_hit{$q}} = "$current~$new";
        }
}
my %final_taxonomy;
foreach my $i(sort keys %abundances){
        if($abundances{$i} == 0){
                next;
        }
        my $final_abundance = $abundances{$i}/$cell_number;
        print STDERR "$spaciator\nSolving tax for $i\n";
        $final_taxonomy{$i} = solve_tax($tax_to_solve{$i});
        print STDERR "Solved, best tax is $final_taxonomy{$i}\n";
        print "$i\t$counts{$i}\t$final_abundance\t$final_taxonomy{$i}\n";
}
open(OTU,">$otu_table");
my %otu_taxonomies;
foreach my $o(sort keys %observed_taxonomies){
        my $final_abundance = $observed_taxonomies{$o}/$cell_number;
        my($protein,$obs_tax) = split(/~/,$o);
	if(exists($final_taxonomy{$protein})){
		my $group_solved_tax = $final_taxonomy{$protein};
		my $this_solved_tax = solve_tax("$group_solved_tax~$obs_tax");
		if($this_solved_tax eq $group_solved_tax){
			$otu_taxonomies{"$protein~$group_solved_tax"} += $final_abundance;
		}else{
			$otu_taxonomies{"$protein~$this_solved_tax"} += $final_abundance;
		}
	}
}
foreach my $o(sort keys %otu_taxonomies){
        my($protein,$obs_tax) = split(/~/,$o);
	print OTU "$obs_tax\t$protein\t$otu_taxonomies{$o}\n";
}
close OTU;

sub solve_tax{
	my $string = shift;
	my @taxonomies = split(/~/,$string);
	my $best_tax = "no taxid assigned";
	if(scalar(@taxonomies) == 1){
		$best_tax = $taxonomies[0];
		#print STDERR "No need to solve, tax is $best_tax\n";
		return $best_tax;
	}else{
		my $conflict = "no conflict";
		my @roots;
		my @domains;
		my @phyla;
		my @classes;
		my @orders;
		my @families;
		my @genera;
		my @species;
		#print STDERR "Observed taxa:\n";
		foreach my $t(@taxonomies){
			#print STDERR "$t\n";
			if($t eq "no taxid assigned"){
				next;
			}
			my @levels = split(/\;/,$t);
			if(defined($levels[0])){
				push @roots,$levels[0];
			}
			if(defined($levels[1])){
                                push @domains,$levels[1];
                        }
			if(defined($levels[2])){
                                push @phyla,$levels[2];
                        }
			if(defined($levels[3])){
                                push @classes,$levels[3];
                        }
			if(defined($levels[4])){
                                push @orders,$levels[4];
                        }
			if(defined($levels[5])){
                                push @families,$levels[5];
                        }
			if(defined($levels[6])){
                                push @genera,$levels[6];
                        }
			if(defined($levels[7])){
                                push @species,$levels[7];
                        }
		}
		my $i = 0;
		while($conflict eq "no conflict" and $i <= 7){
			if($i == 0){
				if(scalar(@roots) > 0){
					my $current = $roots[0];
					foreach my $k(@roots){
						if($k ne $current){
							$conflict = "found conflict";
						}
					}
					if($conflict eq "no conflict"){
						$best_tax = $current;
					}
				}
			}elsif($i == 1){
                                if(scalar(@domains) > 0){
                                        my $current = $domains[0];
                                        foreach my $k(@domains){
                                                if($k ne $current){
                                                        $conflict = "found conflict";
                                                }
                                        }
					if($conflict eq "no conflict"){
                                                $best_tax .= ";$current";
                                        }
                                }
                        }elsif($i == 2){
                                if(scalar(@phyla) > 0){
                                        my $current = $phyla[0];
                                        foreach my $k(@phyla){
                                                if($k ne $current){
                                                        $conflict = "found conflict";
                                                }
                                        }
                                        if($conflict eq "no conflict"){
                                                $best_tax .= ";$current";
                                        }
                                }
                        }elsif($i == 3){
                                if(scalar(@classes) > 0){
                                        my $current = $classes[0];
                                        foreach my $k(@classes){
                                                if($k ne $current){
                                                        $conflict = "found conflict";
                                                }
                                        }
                                        if($conflict eq "no conflict"){
                                                $best_tax .= ";$current";
                                        }
                                }
                        }elsif($i == 4){
                                if(scalar(@orders) > 0){
                                        my $current = $orders[0];
                                        foreach my $k(@orders){
                                                if($k ne $current){
                                                        $conflict = "found conflict";
                                                }
                                        }
                                        if($conflict eq "no conflict"){
                                                $best_tax .= ";$current";
                                        }
                                }
                        }elsif($i == 5){
                                if(scalar(@families) > 0){
                                        my $current = $families[0];
                                        foreach my $k(@families){
                                                if($k ne $current){
                                                        $conflict = "found conflict";
                                                }
                                        }
                                        if($conflict eq "no conflict"){
                                                $best_tax .= ";$current";
                                        }
                                }
                        }elsif($i == 6){
                                if(scalar(@genera) > 0){
                                        my $current = $genera[0];
                                        foreach my $k(@genera){
                                                if($k ne $current){
                                                        $conflict = "found conflict";
                                                }
                                        }
                                        if($conflict eq "no conflict"){
                                                $best_tax .= ";$current";
                                        }
                                }
                        }elsif($i == 7){
                                if(scalar(@species) > 0){
                                        my $current = $species[0];
                                        foreach my $k(@species){
                                                if($k ne $current){
                                                        $conflict = "found conflict";
                                                }
                                        }
                                        if($conflict eq "no conflict"){
                                                $best_tax .= ";$current";
                                        }
                                }
                        }
			$i++;
		}
		return "$best_tax";
	}
}
