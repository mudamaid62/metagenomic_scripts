#!/usr/bin/perl
use warnings;
use strict;

my $bin_dir = shift(@ARGV);
my $quality_file = shift(@ARGV);
my $new_dir = shift(@ARGV);

my %quality;
open(QUALITY,"$quality_file");
while(my $x = <QUALITY>){
	chomp($x);
	my($bin,$bin_quality) = split(/\t/,$x);
	$quality{$bin} = $bin_quality;
}
close QUALITY;
my %contigs;
my %candidates;
foreach my $z(sort keys %quality){
	my $fasta = read_fasta("$bin_dir\/$z\.fasta");
	my %seqs = fasta_parser($fasta);
	foreach my $s(keys %seqs){
		if(!exists($contigs{$s})){
			$contigs{$s} = $seqs{$s};
		}
		if(!exists($candidates{$s})){
			$candidates{$s} = $z;
		}else{
			my $current_quality = $quality{$candidates{$s}};
			my $new_quality = $quality{$z};
			if($new_quality > $current_quality){
				print "For contig $s, quality $new_quality from bin $z is better than quality $current_quality from $candidates{$s}\n";
				$candidates{$s} = $z;
			}else{
				print "For contig $s, quality $new_quality from bin $z is worse than quality $current_quality from $candidates{$s}\n";
			}
		}
	}
}
mkdir($new_dir);
foreach my $i(keys %quality){
	my @i_array;
	foreach my $j(keys %candidates){
		if($candidates{$j} eq $i){
			push @i_array,$j;
		}
	}
	if(scalar(@i_array) > 0){
		open(BIN,">$new_dir/$i\.fasta");
		foreach my $k(@i_array){
			print BIN ">$k\n$contigs{$k}\n";
		}
		close BIN;
	}else{
		print STDERR "Bin $i was left without contigs\n";
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
