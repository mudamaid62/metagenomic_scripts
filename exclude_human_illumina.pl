#!/usr/bin/perl
use warnings;
use strict;

my $sam = shift @ARGV;
my ($prefix,$end) = split(/\./,$sam);
system "samtools view -@ 15 -b -e 'flag & 4 || flag & 8' -F 2304 $sam | samtools collate -u -@ 15 -O - | samtools fastq -1 $prefix\_filtered_1.fastq.gz -2 $prefix\_filtered_2.fastq.gz -0 /dev/null -s /dev/null -";
system "rm $sam";
