#!/usr/bin/perl
use warnings;
use strict;

my $dir = shift(@ARGV);
my $min_len = shift(@ARGV);

my $list = `ls $dir`;
chomp($list);
my @files = split(/\n/,$list);
print "MAG\tNumber_of_contigs\tTotal_length(bp)\tLargest_contig_length(bp)\tShortest_contig_length(bp)\tMedian_contig_length(bp)\tLargest_contig_proportion(\%)\tGC_content(\%)\tN50\tL50\n";

foreach my $f(@files){
	my $fasta = read_fasta("$dir/$f");
	my %seqs = fasta_parser($fasta);
	my $total_len = 0;
	my $gc_len = 0;
	my $contigs = 0;
	my @lens;
	foreach my $x(keys %seqs){
		my $len = length($seqs{$x});
		$seqs{$x} =~ s/[ATU]+//gi;
		my $gc = length($seqs{$x});
		if($len >= $min_len){
			$total_len += $len;
			$gc_len += $gc;
			push @lens, $len;
			$contigs++;
		}
	}
	my $gc_content = ($gc_len/$total_len)*100;
	my @sorted = sort{$b<=>$a}(@lens);
	my $parity = $contigs%2;
	my $median = 0;
	if($parity == 0){
		my $center_1 = $contigs/2;
		my $center_2 = $center_1 - 1;
		$median = ($sorted[$center_1] + $sorted[$center_2])/2;
	}else{
		my $center = int($contigs/2);
		$median = $sorted[$center];
	} 
	my $half = $total_len/2;
	my $current = 0;
	my $i = 0;
	until($current >= $half){
		$current += $sorted[$i];
		$i++;
	}
	my $index = $i - 1;
	my $n50 = $sorted[$index];
	my $largest = $sorted[0];
	my $shortest = pop(@sorted);
	my $prop = ($largest/$total_len)*100;
	print "$f\t$contigs\t$total_len\t$largest\t$shortest\t$median\t$prop\t$gc_content\t$n50\t$i\n";
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
