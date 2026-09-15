#!/usr/bin/perl
use warnings;
use strict;

my $mac_table = shift(@ARGV) or die "Usage: $0 <macsyfinder_output.tsv>\n";

# Hierarchy for hit status (higher number = higher priority)
my %status_rank = (
	'mandatory' => 4,
	'accessory' => 3,
	'neutral'   => 2,
	'forbidden' => 1
);

my %best_hit;
my %best_evalue;
my %best_status_rank;

open(my $fh, '<', $mac_table) or die "Could not open file '$mac_table': $!\n";

while (my $line = <$fh>){
	chomp($line);
	next if $line eq '' || $line =~ /^#/ || $line =~ /^replicon/;
	my @fields = split(/\t/, $line);
    	my $site = $fields[0];
	my $hit = $fields[1];
	my $hit_status = $fields[8];
	my $evalue = $fields[10];
	my $name = "$site~$hit";
    	# Fallback to rank 0 if an unknown status string is encountered
	my $current_rank = $status_rank{$hit_status} // 0;
	if(!exists $best_hit{$name}){
		# 1. New record initialization
		save_hit($name, $line, $evalue, $current_rank);
	}elsif($current_rank > $best_status_rank{$name}){
        	# 2. Higher status priority wins strictly regardless of e-value
		save_hit($name, $line, $evalue, $current_rank);
	}elsif($current_rank == $best_status_rank{$name} && $evalue < $best_evalue{$name}){
		# 3. Equal status priority: lower (better) e-value wins
		save_hit($name, $line, $evalue, $current_rank);
	}
}
close($fh);

# Output hits sorted by key
foreach my $key (sort keys %best_hit){
	print "$best_hit{$key}\n";
}

sub save_hit {
	my ($key, $line, $evalue, $rank) = @_;
	$best_hit{$key} = $line;
	$best_evalue{$key} = $evalue;
	$best_status_rank{$key} = $rank;
}
