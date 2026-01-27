#!/usr/bin/perl
use warnings;
use strict;

my $db = shift(@ARGV);
my $candidates = shift(@ARGV);
my $temp_dir = shift(@ARGV);

mkdir "$temp_dir";
open(LIST, "$candidates") or die "$!\n";
my $candidate_number = `less $candidates | wc -l`;
chomp($candidate_number);
my $slices = int($candidate_number/10);
my $just = $candidate_number % 10;
my $prefix_len = length($slices);
print STDERR "Candidate number $candidate_number, slices $slices, just $just\n";
my @headers; 
my $s = 1;
my $k = 0;
my $m = 1;
my @final;
my @short;
while(my $x = <LIST>){
	my $pre = "0" x $prefix_len;
	until(length($pre.$m) == $prefix_len){
		chop($pre);
	}
	my $n = "$pre$m";
	chomp($x);
	print STDERR "Current s $s,k $k, m $m, pre $pre, n $n, x $x\n";
	push @headers, $x;
	if($just == 0 and $s == $slices){
		push @final, $x;
	}elsif($s > $slices){
		push @final, $x;
	}elsif($k < 10){
		push @short, $x;
		$k++;
	}elsif($k == 10){
		my $list = join "\,", @short;
		#print STDERR "CMD $n --> mmseqs view $db --id-list $list --idx-entry-type 0 --id-mode 1 > $temp_dir/mmseqs_view_tmp_$n\n";
		system "mmseqs view $db --id-list $list --idx-entry-type 0 --id-mode 1 > $temp_dir/mmseqs_view_tmp_$n 2> $temp_dir/mmseqs_errors_$n"; 
		$m++;
		for my $p(0..9){
                     	delete($short[$p]);
                }
		$s++;
		$k = 0;
		if($just == 0 and $s == $slices){
			push @final, $x;
		}elsif($s > $slices){
                	push @final, $x;
		}else{
			push @short, $x;
			$k++;
		}
	}
}
my $pre = "0" x $prefix_len;
until(length($pre.$m) == $prefix_len){
	chop($pre);
}
my $n = "$pre$m";
my $list = join "\,", @final;
print STDERR "CMD $n --> mmseqs view $db --id-list $list --idx-entry-type 0 --id-mode 1 > $temp_dir/mmseqs_view_tmp_$n\n";
system "mmseqs view $db --id-list $list --idx-entry-type 0 --id-mode 1 > $temp_dir/mmseqs_view_tmp_$n 2> $temp_dir/mmseqs_errors_$n";
system "cat $temp_dir/mmseqs_view_tmp_* > $temp_dir/complete_file";
system "cat $temp_dir/mmseqs_errors_* > get_esm_errors.log";
open (MM, "$temp_dir/complete_file");
my $i = 0;
while(my $y = <MM>){
	chomp($y);
	print ">$headers[$i]\n$y\n";
	$i++;
}
system "rm -r $temp_dir";
