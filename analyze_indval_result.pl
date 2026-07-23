#!/usr/bin/perl
use warnings;
use strict;

my $indval = shift(@ARGV);
my $data_table = shift(@ARGV); #SARG or VFDB final table
my $mode = shift(@ARGV); #SARG or VFDB


#Site    Protein FS_cluster      Observed_genes  Consensus_protein       VFC     VF      VF_full Consolidated_VF Counts  Abundance       Taxonomy
#Site    Protein FS_cluster      Observed_genes  Consensus_protein       Type    Mechanism       Sub1    Sub2    Special Counts  Abundance       Taxonomy

my %site_protein_tax;
my %protein_metadata;
my %site_protein_counts;
my %site_protein_abundance;

open(DATA,"$data_table");
while(my $x = <DATA>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $site = $array[0];
	if($site eq "Site"){
		next;
	}
	my $protein = $array[1];
	my $cluster = $array[2];
	my $metadata = "placeholder";
	my $taxonomy = "no taxid assigned";
	my $counts = 0;
	my $abundance = 0;
	if($mode eq "SARG"){
		$metadata = "$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\t$array[8]\t$array[9]";
		$counts = $array[10];
		$abundance = $array[11];
		$taxonomy = $array[12];
	}elsif($mode eq "VFDB"){
		$metadata = "$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\t$array[8]";
		$counts = $array[9];
		$abundance = $array[10];
		$taxonomy = $array[11];
	}
	$protein_metadata{$protein} = $metadata;
	my $name = "$site~$protein";
	if(exists($site_protein_tax{$name})){
		die "Wrong formatted table > $name\n";
	}
	$site_protein_tax{$name} = $taxonomy;
	$site_protein_counts{$name} = $counts;
	$site_protein_abundance{$name} = $abundance;
}
close DATA;
open(INDVAL,"$indval");
while(my $x = <INDVAL>){
	chomp($x);
	my($protein,$a_value,$b_value,$stat,$p_value,$sig,$group) = split(/\t/,$x);
	my @sites = split(/\,/,$group);
	my $tax_to_solve = "placeholder";
	my @abundances;
	my @counts;
	my @names = ("SCp1","SCp2","SCp3","Mp1","Mp3","Rp1","Rp2","Rp3");
	foreach my $n(@names){
		my $name = "$n~$protein";
		my $default = 0;
		if(exists($site_protein_tax{$name})){
			if($tax_to_solve eq "placeholder"){
				$tax_to_solve = $site_protein_tax{$name};
			}else{
				my $new_tax = "$tax_to_solve~$site_protein_tax{$name}";
				$tax_to_solve = $new_tax;
			}
			push @abundances, $site_protein_abundance{$name};
			push @counts, $site_protein_counts{$name};
		}else{
			push @abundances, $default;
			push @counts, $default;
		}
	}
	my $taxonomy = solve_tax($tax_to_solve);
	my $counts_string = join "\t",@counts;
	my $abundances_string = join "\t",@abundances;
	print "$x\t$protein_metadata{$protein}\t$taxonomy\t$counts_string\t$abundances_string\n";
}

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
