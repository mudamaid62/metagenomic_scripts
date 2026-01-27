#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

#BECOMES FASTER THAN get_best_hit.pl when the diamond file is about 100000 lines or more

my $file = shift(@ARGV);
my $min_pident = shift(@ARGV);
my $min_qcov = shift(@ARGV);
my $min_scov = shift(@ARGV);
my $min_len_frac = shift(@ARGV);
my $max_len_frac = shift(@ARGV);
my $forks = shift(@ARGV); #64 is the approximate optimal number
my $temp_dir = shift(@ARGV);

my $sep = "-" x 30;
my %hits;
my %queries;
my $query_number = 0;
open(HITS,"$file");
while(my $x = <HITS>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $pair = "$array[0]~$array[2]";
	my $pident = $array[10];
	my $qcov = $array[11];
	my $scov = $array[12];
	my $len_frac = $array[1]/$array[3];
	if($pident >= $min_pident and $qcov >= $min_qcov and $scov >= $min_scov and $len_frac >= $min_len_frac and $len_frac <= $max_len_frac){
		push @array, $len_frac;
		my $register = join "\t",@array;
		$hits{$pair} = $register;
		if(!exists($queries{$array[0]})){
			$queries{$array[0]} = "yes";
			$query_number++;
		}
	}
}
close HITS;
my %best_hits;
my $i = 0;
mkdir $temp_dir;
my $pm = Parallel::ForkManager->new($forks, $temp_dir);
$pm -> run_on_finish(
	sub{
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
		if(defined($data_structure_reference)){
			$i++;
			my $string = ${$data_structure_reference};
			my ($query,$best_result,$log) = split(/~/,$string);
			$best_hits{$query} = "$best_result";
			my $percent = int(($i/$query_number)*100);
        		print STDERR "$sep\nQuery $i / $query_number ($percent\%)\n$log\n";
		}else{
			print STDERR "No message received from child process $pid!\n";
		}
	}
);
GENE_LOOP:

foreach my $q(keys %queries){
	my $pid = $pm->start and next GENE_LOOP;
	my @log_text;
	my %best;
        my %best_score;
        my %best_pident;
        my %best_scov;
        my %best_qcov;
	foreach my $k(sort keys %hits){
		my @pair = split(/~/,$k);
		if($pair[0] eq $q){
			my @register = split(/\t/,$hits{$k});
			my $score = $register[9];
			my $pident = $register[10];
			my $qcov = $register[11];
			my $scov = $register[12];
			my $len_frac = $register[13];
			my $data = join "\t",@register;
			if(!exists($best{$q})){
				$best{$q} = $data;
				$best_score{$q} = $score;
				$best_pident{$q} = $pident;
				$best_scov{$q} = $scov;
				$best_qcov{$q} = $qcov;
				my $log = "No score for $q --> $best_score{$q} for target $pair[1]";
				push @log_text, $log;
			}elsif($score > $best_score{$q}){
				$best{$q} = $data;
				my $log = "Score $score for target $pair[1] better than last $best_score{$q}";
				$best_score{$q} = $score;
				$best_pident{$q} = $pident;
				$best_scov{$q} = $scov;
                                $best_qcov{$q} = $qcov;
				push @log_text, $log;
			}elsif($score == $best_score{$q} and $pident > $best_pident{$q}){
				$best{$q} = $data;
                                my $log = "Score $score for target $pair[1] is the same but ID% $pident is higher than previous $best_pident{$q}";
                                $best_score{$q} = $score;
                                $best_pident{$q} = $pident;
				$best_scov{$q} = $scov;
                                $best_qcov{$q} = $qcov;
				push @log_text, $log;
			}elsif($score == $best_score{$q} and $scov > $best_scov{$q}){
                                $best{$q} = $data;
                                my $log = "Score $score for target $pair[1] is the same but SCOV% $scov is higher than previous $best_scov{$q}";
                                $best_score{$q} = $score;
                                $best_pident{$q} = $pident;
				$best_scov{$q} = $scov;
                                $best_qcov{$q} = $qcov;
				push @log_text, $log;
			}elsif($score == $best_score{$q} and $qcov > $best_qcov{$q}){
                                $best{$q} = $data;
                                my $log = "Score $score for target $pair[1] is the same but QCOV% $qcov is higher than previous $best_qcov{$q}";
                                $best_score{$q} = $score;
                                $best_pident{$q} = $pident;
                                $best_scov{$q} = $scov;
                                $best_qcov{$q} = $qcov;
				push @log_text, $log;
			}
		}
	}
	my $log_line = join "\n", @log_text;
        my $return = "$q~$best{$q}~$log_line";
        $pm->finish(0, \$return);
}
$pm->wait_all_children;
foreach my $w(sort keys %best_hits){
	print "$file\t$best_hits{$w}\n";
}		
rmdir "$temp_dir";
