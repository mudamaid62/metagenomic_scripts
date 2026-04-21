#!/usr/bin/perl
use warnings;
use strict;

my $VFDB_table = shift(@ARGV);
my $structure = shift(@ARGV);
my $VFDB_extra = shift(@ARGV);
my $first_clustering = shift(@ARGV);
my $second_clustering = shift(@ARGV);
my $FS_clustering = shift(@ARGV);

my %first;
my %second;
my %FS;

open(FIRST,"$first_clustering");
while(my $x = <FIRST>){
        chomp($x);
        $x =~ s/\"//g;
	my($cluster,$protein) = split(/\t/,$x);
        $first{$protein} = $cluster;
}
close FIRST;
open(SECOND,"$second_clustering");
while(my $x = <SECOND>){
	chomp($x);
	$x =~ s/\"//g;
	my($cluster,$protein) = split(/\t/,$x);
	$second{$protein} = $cluster;
}	
close SECOND;
open(FS,"$FS_clustering");
while(my $x = <FS>){
        chomp($x);
        $x =~ s/\"//g;
        my($cluster,$protein) = split(/\t/,$x);
        $FS{$protein} = $cluster;
}
close FS;
my %type;
my %subtype;
my %pathogen;
open(STRUCTURE,"$structure");
while(my $z = <STRUCTURE>){
	chomp($z);
	my($gene,$z_subtype,$z_type,$z_pathogen) = split(/\t/,$z);
	$type{$gene} = $z_type;
	$subtype{$gene} = $z_subtype;
	$pathogen{$gene} = $z_pathogen;
}
close STRUCTURE;
open(EXTRA,"$VFDB_extra");
my %full;
my %char;
my %structure;
my %function;
my %mechanism;
while(my $i = <EXTRA>){
	chomp($i);
	my @i_array = split(/\t/,$i);
	$full{$i_array[0]} = $i_array[2];
	$char{$i_array[0]} = $i_array[6];
	$structure{$i_array[0]} = $i_array[7];
	$function{$i_array[0]} = $i_array[8];
	$mechanism{$i_array[0]} = $i_array[9];
}
close EXTRA;
open(TABLE,"$VFDB_table");
#OG_header       VFG     gb      Gene    Protein VF_type VF      VF_class        VFC     Species
while(my $k = <TABLE>){
	chomp($k);
	my @k_array = split(/\t/,$k);
	if($k_array[0] eq "OG_header"){
		next;
	}
	my $inherited_second = "none";
        my $inherited_FS = "none";
        if(exists($second{$first{$k_array[1]}})){
                $inherited_second = $second{$first{$k_array[1]}};
        }
        if(exists($FS{$second{$first{$k_array[1]}}})){
                $inherited_FS = $FS{$second{$first{$k_array[1]}}};
        }
        if($inherited_FS eq "none" and exists($FS{$first{$k_array[1]}})){
                $inherited_FS = $FS{$first{$k_array[1]}};
        }
	my $k_full = "-";
	my $k_char = "-";
	my $k_structure = "-";
	my $k_function = "-";
	my $k_mechanism = "-";
	if(exists($full{$k_array[6]})){
		$k_full = $full{$k_array[6]};
	}
	if(exists($char{$k_array[6]})){
                $k_char = $char{$k_array[6]};
        }
	if(exists($structure{$k_array[6]})){
                $k_structure = $structure{$k_array[6]};
        }
	if(exists($function{$k_array[6]})){
                $k_function = $function{$k_array[6]};
        }
	if(exists($mechanism{$k_array[6]})){
                $k_mechanism = $mechanism{$k_array[6]};
        }
	print "$k_array[1]\t$k_array[3]\t$k_array[4]\t$type{$k_array[1]}\t$subtype{$k_array[1]}\t$k_full\t$k_char\t$k_structure\t$k_function\t$k_mechanism\t$pathogen{$k_array[1]}\t$first{$k_array[1]}\t$inherited_second\t$inherited_FS\n";
}
