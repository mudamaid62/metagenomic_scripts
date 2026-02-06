#!/usr/bin/perl
use warnings;
use strict;

my $plass_file = shift(@ARGV);
my $assembly_name = shift(@ARGV);
my $min_len = shift(@ARGV);
my $max_len = shift(@ARGV);

print STDERR "Reading $plass_file\n";
my $fasta = read_fasta($plass_file);
my %seqs = fasta_parser($fasta);
open(COMPLETE, ">$assembly_name\_complete_orfs.faa");
open(SEMI, ">$assembly_name\_semi-partial_orfs.faa");
open(PARTIAL, ">$assembly_name\_partial_orfs.faa");

print STDERR "Filtering $plass_file\n";
my $complete = 0;
my $semi = 0;
my $partial = 0;
my $discarded = 0;
my $total = scalar(keys %seqs);
foreach my $x(keys %seqs){
	my $og_seq = $seqs{$x};
	$seqs{$x} =~ s/[\*]+//g; 
	my $len = length($seqs{$x});
	if($len >= $min_len and $len <= $max_len){
		my $start = substr($og_seq,0,1);
		my $end = substr($og_seq,-1);
		if($start eq "\*" and $end eq "\*"){
			print COMPLETE ">$x\_PLASS_$assembly_name\n$seqs{$x}\n";
			$complete++;
		}elsif($start eq "\*" or $end eq "\*"){
			print SEMI ">$x\_PLASS_$assembly_name\n$seqs{$x}\n";
			$semi++;
		}else{
			print PARTIAL ">$x\_PLASS_$assembly_name\n$seqs{$x}\n";
			$partial++;
		}
	}else{
		$discarded++;
	}
}
my $divider = "-" x 30;
my $t_percent = ($total/$total) *100;
my $c_percent = ($complete/$total) *100;
my $s_percent = ($semi/$total) *100;
my $p_percent = ($partial/$total) *100;
my $d_percent = ($discarded/$total) *100;
print "$divider\nSUMMARY\nTotal ORFs-->$total ($t_percent\%)\nComplete ORFs-->$complete ($c_percent\%)\nSemi-Partial ORFs-->$semi ($s_percent\%)\nPartial ORFs-->$partial ($p_percent\%)\nDiscarded ORFs-->$discarded ($d_percent\%)\n$divider\n";

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
