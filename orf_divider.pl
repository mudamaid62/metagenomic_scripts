#!/usr/bin/perl
use warnings;
use strict;

my $protein_file = shift(@ARGV);

print STDERR "Reading fasta\n";
my $fasta = read_fasta($protein_file);
my %seqs = fasta_parser($fasta);
open(DB_1,">100_db_proteins.faa");
open(DB_2,">200_db_proteins.faa");
open(DB_3,">300_db_proteins.faa");
open(DB_4,">400_db_proteins.faa");
open(DB_5,">500_db_proteins.faa");
open(DB_6,">600_db_proteins.faa");
open(DB_7,">700_db_proteins.faa");
open(DB_8,">800+_db_proteins.faa");

foreach my $x(keys %seqs){
	$seqs{$x} =~ s/\*//g;
	my $len = length($seqs{$x});
	if($len <= 100){
		print DB_1 ">$x\n$seqs{$x}\n";
	}elsif($len > 100 and $len <= 200){
		print DB_2 ">$x\n$seqs{$x}\n";
	}elsif($len > 200 and $len <= 300){
                print DB_3 ">$x\n$seqs{$x}\n";
        }elsif($len > 300 and $len <= 400){
                print DB_4 ">$x\n$seqs{$x}\n";
        }elsif($len > 400 and $len <= 500){
                print DB_5 ">$x\n$seqs{$x}\n";
        }elsif($len > 500 and $len <= 600){
                print DB_6 ">$x\n$seqs{$x}\n";
        }elsif($len > 600 and $len <= 700){
                print DB_7 ">$x\n$seqs{$x}\n";
        }elsif($len > 700){
                print DB_8 ">$x\n$seqs{$x}\n";
	}else{
		print STDERR "$x might have an error\n";
	}
}

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
			$header =~ s/[ >\:\[\/\\]+/_/g;
			$header =~ s/\]//g;
                        my $seq = join "",@f_array;
                        $out{$header} = $seq;
                }
        }
        return %out;
}
