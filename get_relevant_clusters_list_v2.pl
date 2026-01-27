#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(any);

my $tsv = shift(@ARGV);
my $interest = shift(@ARGV);

print STDERR "Reading tsv\n";
my %mmseqs = read_mm($tsv);
print STDERR "Reading interest file\n";
my @prots = read_p($interest);
my @reps;
my %database;

print "Cluster_representative\tProtein\tDatabase\n";
foreach my $p(@prots){
	$database{$p} = "Reference DB";
	print STDERR "Searching cluster containing $p\n";
	my @matches = grep {$_ eq $p} keys %mmseqs;
	foreach my $y(@matches){
		unless(any {$_ eq $mmseqs{$y}} @reps){
           		push @reps,$mmseqs{$y};
			my $current = scalar(@reps);
			print STDERR "Added $mmseqs{$y} to rep list. Current --> $current reps\n";
          	}
	}
}
foreach my $r(@reps){
	print STDERR "Building $r cluster\n";
	my @r_matches = grep {$mmseqs{$_} eq $r} keys %mmseqs;
	foreach my $z(@r_matches){
		$database{$z} = "Problem DB" if !exists($database{$z});
		print "$mmseqs{$z}\t$z\t$database{$z}\n";      
        }
}

sub read_p {
	my $file = shift;
	my @text;
	open(FILE,"$file") or die "$!\n";
	while(my $x = <FILE>){
		chomp($x);
		$x =~ s/\"//g;
		push @text, $x;
	}
	return @text;
}

sub read_mm {
	my $file = shift;
        my %pairs;
        open(FILE,"$file") or die "$!\n";
        while(my $x = <FILE>){
                chomp($x);
                $x =~ s/\"//g;
		my ($rep,$member) = split(/\t/,$x);
		$pairs{$member} = "$rep";
        }
        return %pairs;
}


