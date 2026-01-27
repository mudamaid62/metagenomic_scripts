#!/usr/bin/perl
use warnings;
use strict;

my $file = shift(@ARGV);
my $out_1 = shift(@ARGV);
my $out_2 = shift(@ARGV);
my $single_out = shift(@ARGV);

my $fasta = read_fasta($file);
my %seqs = fasta_parser($fasta);
my %quality;
my @ids;
my %check;

open(OUT1, ">$out_1");
open(OUT2, ">$out_2");
foreach my $x(keys %seqs){
	my @x_array = split(/\@/,$x);
	my @sample = split(/_/,$x_array[1]);
	my $file = pop(@sample);
	if($file eq "R1" or $file eq "fwd"){
		$file = 1;
	}elsif($file eq "R2" or $file eq "rev"){
		$file = 2;
	}
	my $len = length($seqs{$x});
	$quality{$x} = "=" x $len;
	my $name = pop(@x_array);
	if(!exists($check{$name})){
		$check{$name} = $file;
	}elsif($check{$name} == 1 and $file == 2){
		$check{$name} = 3;
	}elsif($check{$name} == 2 and $file == 1){
		$check{$name} = 4;
	}
	my $id = "$name~$file~$x";
	push @ids,$id;	
}	
my @rest;
my $f = 0;
my $r = 0;
foreach my $z(sort(@ids)){
	my ($name,$file,$x) = split(/~/,$z);
	if($check{$name} == 3 or $check{$name} == 4){
		if($file == 1){
			print OUT1 "\@$name\n$seqs{$x}\n+\n$quality{$x}\n";
			$f++;
		}else{
			print OUT2 "\@$name\n$seqs{$x}\n+\n$quality{$x}\n";
			$r++;
		}
	}else{
		push @rest,$z;
	}
}
close OUT1;
close OUT2;
open(SINGLE,">$single_out");
my $s_f = 0;
my $s_r = 0;
foreach my $s(@rest){
	my ($name,$file,$x) = split(/~/,$s);
	if($file == 1){
		$s_f++;
	}else{
		$s_r++;
	}
	print SINGLE "\@$name\n$seqs{$x}\n+\n$quality{$x}\n";
}
print STDERR "FORWARD\t$f\nREVERSE\t$r\nSINGLE_F\t$s_f\nSINGLE_R\t$s_r\n";

sub read_fasta{
        my $file = shift;
        my @lines;
        open(FASTA,"$file") or die "$file not found $!";
        while(my $x = <FASTA>){
                chomp($x);
                if($x =~ m/>/){
                        my @x_array = split(/>/,$x);
                        my $white = shift(@x_array);
                        my $pre = join "_",@x_array;
                        my $y = ">$pre";
                        $x = ">$pre";
                }
                push @lines, $x;
        }
        my $fasta = join "\n",@lines;
        return $fasta;
        close FASTA;
}
sub fasta_parser{
        my $fasta = shift;
        my @seqs = split(/>/,$fasta);
        my %out;
        foreach my $x(@seqs){
                if($x eq ""){
                        next;
                }else{
                        my @f_array = split(/\n/,$x);
                        my $header = shift(@f_array);
                        #$header =~ s/[ >\:\[\/\\]+/_/g;
                        #$header =~ s/\]//g;
                        my $seq = join "",@f_array;
                        $out{$header} = $seq;
                }
        }
        return %out;
}
