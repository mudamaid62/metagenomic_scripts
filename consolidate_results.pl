#!/usr/bin/perl
use warnings;
use strict;

my $site = shift(@ARGV); #site name
my $abundance_gini = shift(@ARGV);
my $plasmids_file = shift(@ARGV); #freqs
my $viruses_file = shift(@ARGV); #freqs
my $ises_file = shift(@ARGV); #freqs
my $integrons_file = shift(@ARGV); #freqs
my $prophages_file = shift(@ARGV); #freqs
my $mode = shift(@ARGV); #SARG or VFDB
my $identification_table = shift(@ARGV);
my $concat_table = shift(@ARGV); #SARG or VFDB

my %FS_clusters;
open(CONCAT,"$concat_table");
while(my $x = <CONCAT>){
        chomp($x);
        my @array = split(/\t/,$x);
        my $name = shift(@array);
        my $fs = pop(@array);
        $FS_clusters{$name} = $fs;
}
my %cluster_data;
open(ID,"$identification_table");
while(my $x = <ID>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $fs = shift(@array);
	my $data = join "\t",@array;
	$cluster_data{$fs} = $data;
}
close ID;
my %total_obs;
my %plasmid_obs;
my %plasmid_freq;
open(PLASMIDS,"$plasmids_file");
while(my $x = <PLASMIDS>){
	chomp($x);
	my @array = split(/\t/,$x);
	$total_obs{$array[0]} = $array[1];
	$plasmid_obs{$array[0]} = $array[2];
	$plasmid_freq{$array[0]} = $array[3];
}
close PLASMIDS;
my %virus_obs;
my %virus_freq;
open(VIRUSES,"$viruses_file");
while(my $x = <VIRUSES>){
        chomp($x);
        my @array = split(/\t/,$x);
        $virus_obs{$array[0]} = $array[2];
        $virus_freq{$array[0]} = $array[4];
}
close VIRUSES;
my %ises_obs;
my %ises_freq;
open(ISES,"$ises_file");
while(my $x = <ISES>){
        chomp($x);
        my @array = split(/\t/,$x);
        $ises_obs{$array[0]} = $array[2];
        $ises_freq{$array[0]} = $array[3];
}
close ISES;
my %integrons_obs;
my %integrons_freq;
open(INTEGRONS,"$integrons_file");
while(my $x = <INTEGRONS>){
        chomp($x);
        my @array = split(/\t/,$x);
        $integrons_obs{$array[0]} = $array[2];
        $integrons_freq{$array[0]} = $array[3];
}
close INTEGRONS;
my %prophages_obs;
my %prophages_freq;
open(PROPHAGES,"$prophages_file");
while(my $x = <PROPHAGES>){
        chomp($x);
        my @array = split(/\t/,$x);
        $prophages_obs{$array[0]} = $array[2];
        $prophages_freq{$array[0]} = $array[3];
}
close PROPHAGES;
if($mode eq "SARG"){
	print "Site\tProtein\tFS_cluster\tObserved_genes\tConsensus_protein\tType\tMechanism\tSub1\tSub2\tSpecial\tCounts\tAbundance\tTaxonomy\tContig_observations\tPlasmid_obs\tVirus_obs\tIS_obs\tIntegron_obs\tProphage_obs\tGini_index\tPlasmid_freq\tVirus_freq\tIS_freq\tIntegron_freq\tProphages_freq\tMobility_score\tScaled_score\n";
}elsif($mode eq "VFDB"){
	print "Site\tProtein\tFS_cluster\tObserved_genes\tConsensus_protein\tVFC\tVF\tVF_full\tConsolidated_VF\tCounts\tAbundance\tTaxonomy\tContig_observations\tPlasmid_obs\tVirus_obs\tIS_obs\tIntegron_obs\tProphage_obs\tGini_index\tPlasmid_freq\tVirus_freq\tIS_freq\tIntegron_freq\tProphages_freq\tMobility_score\tScaled_score\n";
}else{
	die "Mode not supported\n";
}
open(GINI,"$abundance_gini");
while(my $x = <GINI>){
	chomp($x);
	my($protein,$counts,$abundance,$tax,$gini) = split(/\t/,$x);
	my $contig_obs = 0;
	my $plasmid_score = 0;
	my $virus_score = 0;
	my $is_score = 0;
	my $integron_score = 0;
	my $prophage_score = 0;
	my $plasmid_c = 0;
        my $virus_c = 0;
        my $is_c = 0;
        my $integron_c = 0;
        my $prophage_c = 0;
	if(exists($total_obs{$protein})){
		$contig_obs = $total_obs{$protein};
	}
	if(exists($plasmid_freq{$protein})){
                $plasmid_score = $plasmid_freq{$protein};
        }
	if(exists($virus_freq{$protein})){
                $virus_score = $virus_freq{$protein};
        }
	if(exists($ises_freq{$protein})){
                $is_score = $ises_freq{$protein};
        }
	if(exists($prophages_freq{$protein})){
                $prophage_score = $prophages_freq{$protein};
        }
	if(exists($integrons_freq{$protein})){
                $integron_score = $integrons_freq{$protein};
        }
	if(exists($plasmid_obs{$protein})){
                $plasmid_c = $plasmid_obs{$protein};
        }
        if(exists($virus_obs{$protein})){
                $virus_c = $virus_obs{$protein};
        }
        if(exists($ises_obs{$protein})){
                $is_c = $ises_obs{$protein};
        }
        if(exists($prophages_obs{$protein})){
                $prophage_c = $prophages_obs{$protein};
        }
        if(exists($integrons_obs{$protein})){
                $integron_c = $integrons_obs{$protein};
        }
	my $mob_score = ($gini+$plasmid_score+$virus_score+$is_score+$prophage_score+$integron_score)/6;
	my $scaled_score = $mob_score * $abundance;
	print STDERR "$protein\t$mob_score\t$scaled_score\n";
	print "$site\t$protein\t$FS_clusters{$protein}\t$cluster_data{$FS_clusters{$protein}}\t$counts\t$abundance\t$tax\t$contig_obs\t$plasmid_c\t$virus_c\t$is_c\t$integron_c\t$prophage_c\t$gini\t$plasmid_score\t$virus_score\t$is_score\t$integron_score\t$prophage_score\t$mob_score\t$scaled_score\n";
}
close GINI;




