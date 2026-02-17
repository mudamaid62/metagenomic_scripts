#!/usr/bin/perl
use warnings;
use strict;

my $msa_file = shift(@ARGV);
my $info_file = shift(@ARGV);

my $fasta = read_fasta($msa_file);
my %seqs = fasta_parser($fasta);
my @informative_sites;
my $site_number = 0;
my $informative_number = 0;

open(INFO,"$info_file") or die "$!\n";
while(my $x = <INFO>){
	chomp($x);
	if($x =~ m/^\#/ or $x =~ m/Site/){
		next;
	}	
	my($site,$type) = split(/\t/,$x);
	if($type eq "I"){
		push @informative_sites, $site;
		$informative_number++;
	}
	$site_number++;
}
close INFO;
print STDERR "MSA sites --> $site_number\nInformative sites --> $informative_number\n";

foreach my $f(sort keys %seqs){
	my @full = split(//,$seqs{$f});
	my $informative = "";
	foreach my $y(@informative_sites){
		my $i = $y -  1;
		$informative .= $full[$i];
	}
	print ">$f\n$informative\n";
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
			my @p_array = split(/\t/,$pre);
			my @q_array = split(/\:/,$p_array[0]);
                        $x = ">$q_array[0]";
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
                        $header =~ s/[ \=>\:\[\/\\]+/_/g;
                        $header =~ s/\]//g;
			$header =~ s/.pdb//g;
                        my $seq = join "",@f_array;
                        $out{$header} = $seq;
                }
        }
        return %out;
}

