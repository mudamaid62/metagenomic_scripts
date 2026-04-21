#!/usr/bin/perl
use warnings;
use strict;

my $SARG_table = shift(@ARGV);
my $structure = shift(@ARGV);
my $SARG_extra = shift(@ARGV);
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
open(STRUCTURE,"$structure");
while(my $z = <STRUCTURE>){
	chomp($z);
	my($gene,$z_subtype,$z_type) = split(/\t/,$z);
	$type{$gene} = $z_type;
	$subtype{$gene} = $z_subtype;
}
close STRUCTURE;
my %hmm;
my %mechanism;
my %m_sub1;
my %m_sub2;
my %special;
open(EXTRA,"$SARG_extra");
while(my $i = <EXTRA>){
	chomp($i);
	my @i_array = split(/\t/,$i);
	$hmm{$i_array[1]} = $i_array[4];
	$mechanism{$i_array[1]} = $i_array[5];
	$m_sub1{$i_array[1]} = $i_array[6];
	$m_sub2{$i_array[1]} = $i_array[7];
	$special{$i_array[1]} = $i_array[0];
}
close EXTRA;
open(TABLE,"$SARG_table");
while(my $k = <TABLE>){
	chomp($k);
	my($name,$header) = split(/\t/,$k);
	if($name eq "Sequence"){
		next;
	}
	my @h_array = split(/ /,$header);
        my $gene = $h_array[0];
	my $inherited_second = "none";
	my $inherited_FS = "none";
	if(exists($second{$first{$name}})){
                $inherited_second = $second{$first{$name}};
        }
	if(exists($FS{$second{$first{$name}}})){
		$inherited_FS = $FS{$second{$first{$name}}};
	}
	if($inherited_FS eq "none" and exists($FS{$first{$name}})){
		$inherited_FS = $FS{$first{$name}};
	}
	my $k_hmm = "-";
	my $k_mechanism = "-";
	my $k_m_sub1 = "-";
	my $k_m_sub2 = "-";
	my $k_special = "-";
	if(exists($hmm{$gene})){
		$k_hmm = $hmm{$gene};
	}
	if(exists($mechanism{$gene})){
                $k_mechanism = $mechanism{$gene};
        }
	if(exists($m_sub1{$gene})){
                $k_m_sub1 = $m_sub1{$gene};
        }
	if(exists($m_sub2{$gene})){
                $k_m_sub2 = $m_sub2{$gene};
        }
	if(exists($special{$gene})){
                $k_special = $special{$gene};
        }
	print "$name\t$gene\t$k_hmm\t$type{$gene}\t$subtype{$gene}\t$k_mechanism\t$k_m_sub1\t$k_m_sub2\t$k_special\t$first{$name}\t$inherited_second\t$inherited_FS\n";
}
