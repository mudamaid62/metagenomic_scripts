#!/usr/bin/perl
use warnings;
use strict;

my $cluster_types = shift(@ARGV); #cluster type, ex: KGI_SARG
my $protein_file = shift(@ARGV); #protein list
my $sarg_result = shift(@ARGV);
my $vfdb_result = shift(@ARGV);
my $first_clustering = shift(@ARGV); #mmseqs initial
my $second_clustering = shift(@ARGV); #mmseqs work_DB
my $FS_clustering = shift(@ARGV);

my %first;
my %second;
my %FS;

print STDERR "Reading clusterings\n";
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
my %clusters;
open(CLUSTER,"$cluster_types");
while(my $x = <CLUSTER>){
	chomp($x);
	my($cluster,$type) = split(/\t/,$x);
	$clusters{$cluster} = $type;
}
close CLUSTER;
open(PROTEINS,"$protein_file");
open(SARG,">$sarg_result");
open(VFDB,">$vfdb_result");
print STDERR "Classifying proteins by reverse inheritance\n";
while(my $z = <PROTEINS>){
	chomp($z);
	my $inherited_second = "none";
        my $inherited_FS = "none";
        if(exists($second{$first{$z}})){
                $inherited_second = $second{$first{$z}};
        }
        if(exists($FS{$second{$first{$z}}})){
                $inherited_FS = $FS{$second{$first{$z}}};
        }
        if($inherited_FS eq "none" and exists($FS{$first{$z}})){
                $inherited_FS = $FS{$first{$z}};
        }
	if($inherited_FS eq "none"){
		next;
	}else{
		if($clusters{$inherited_FS} =~ m/SARG/){	
			print SARG "$z\t-\t-\t-\t-\t-\t-\t-\t-\t$first{$z}\t$inherited_second\t$inherited_FS\n";
		}
		if($clusters{$inherited_FS} =~ m/VFDB/){
                	print VFDB "$z\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t$first{$z}\t$inherited_second\t$inherited_FS\n";
		}
	}
}
