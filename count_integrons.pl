#!/usr/bin/perl
use warnings;
use strict;

#ID_replicon     CALIN   complete        In0     topology        size
#NODE_9_length_115122_cov_10.986903      0       0       0       lin     115122
my $integron_dir = shift(@ARGV);

my $list = `ls $integron_dir/`;
my @files = split(/\n/,$list);
my $calin = 0;
my $complete = 0;
my $in0 = 0;

foreach my $x(@files){
	if($x eq "integron_finder.out"){
		next;
	}
	my ($node,$number,$extension) = split(/\./,$x);
	if($extension ne "summary"){
		next;
	}
	open(FILE,"$integron_dir/$x");
	while(my $z = <FILE>){
		chomp($z);
		my @array = split(/\t/,$z);
		if($array[0] eq "ID_replicon"){
			next;
		}else{
			$calin += $array[1];
			$complete += $array[2];
			$in0 += $array[3];
		}
	}
	close FILE;
}	
print "$integron_dir calin > $calin, complete > $complete, In0 > $in0\n";
