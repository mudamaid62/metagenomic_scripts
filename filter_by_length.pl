#!/usr/bin/perl
use warnings;
use strict;

my $fasta_file = shift(@ARGV);
my $min_len = shift(@ARGV);
my $max_len = shift(@ARGV);

print STDERR "Reading fasta\n";
my $fasta = read_fasta($fasta_file);
my %seqs = fasta_parser($fasta);

print STDERR "Writing filtered fasta\n";
my $kept = 0;
my $discarded = 0;
my $total = scalar(keys %seqs);
foreach my $x(keys %seqs){
	$seqs{$x} =~ s/\*//g;
	my $len = length($seqs{$x});
	if($len >= $min_len and $len <= $max_len){
		print ">$x\n$seqs{$x}\n";
		$kept++;
	}else{
		$discarded++;
	}
}
my $divider = "-" x 30;
my $t_percent = ($total/$total) *100;
my $k_percent = ($kept/$total) * 100;
my $d_percent = ($discarded/$total) *100;
print STDERR "$divider\nSUMMARY\nTotal SEQs-->$total ($t_percent\%)\nKept SEQs-->$kept ($k_percent\%)\nDiscarded SEQs-->$discarded ($d_percent\%)\n$divider\n";

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
			$header =~ s/[ >\:]+/_/g;
                        my $seq = join "",@f_array;
                        $out{$header} = $seq;
                }
        }
        return %out;
}
