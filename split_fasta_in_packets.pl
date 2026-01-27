#!/usr/bin/perl
use warnings;
use strict;

my $fasta_file = shift(@ARGV);
my $packet_prefix = shift(@ARGV);
my $max_number_per_split = shift(@ARGV); #1000

my $fasta = read_fasta($fasta_file);
my %seqs = fasta_parser($fasta);
my $i = 1;
my $final_splits_number = 0;
foreach my $x (sort keys %seqs){
	my $splits_number = 0;
	until((($splits_number * $max_number_per_split)/$i) >= 1){
		$splits_number++;
	}
	if($splits_number > $final_splits_number){
		$final_splits_number = $splits_number;
	}
	open my $out, ">>", "$packet_prefix\_$splits_number\.faa";
        print $out ">$x\n$seqs{$x}\n";
        close $out;
	if(($i % 1000) == 0){
                print STDERR "Printed $i sequences\n";
        }
	$i++;
}
	
print STDERR "Splits --> $final_splits_number\n";

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
