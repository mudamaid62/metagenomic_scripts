#!/usr/bin/perl
use warnings;
use strict;

my $cons_table = shift(@ARGV);
my $hmm_result = shift(@ARGV);
my $potential_hmm_list = shift(@ARGV);
my $missed_hmm_list = shift(@ARGV);

my %structural_cluster;
my %redundant_protein;
my %ref_class;
my %consensus_class;
my %type;
#in consensus class
my @positives = ('A','B1','B2','B3','C','D');
my @negatives = ('PBP','MBL');

open(CONSENSUS,"$cons_table");
while(my $x = <CONSENSUS>){
	chomp($x);
	my @array = split(/\t/,$x);
	if($array[4] ne "Candidate"){
                next;
        }
	my $protein = $array[0];
	my $sc = $array[1];
	$structural_cluster{$protein} = $sc;
	my $r_class = $array[4];
	$ref_class{$protein} = $r_class;
	my $c_class = $array[5];
	my $confusion_type = "NA";
	$consensus_class{$protein} = $c_class;
	foreach my $p(@positives){
		if($c_class eq $p){
			$confusion_type = "P";
			$type{$protein} = $confusion_type;
			
		}
	}
	foreach my $n(@negatives){
		if($c_class eq $n){
			$confusion_type = "N";
			$type{$protein} = $confusion_type;
		}
	}
	my @redundant;
	if($array[15] ne "-"){
		@redundant = split(/\;/,$array[15]);
		foreach my $z(@redundant){
			$redundant_protein{$z} = "yes";
			$structural_cluster{$z} = $sc;
			$ref_class{$z} = $r_class;
			$consensus_class{$z} = $c_class;
			$type{$z} = $confusion_type;
		}
	}
}	
close CONSENSUS;
my $P_number = 0;
my $N_number = 0;
foreach my $z(keys %type){
	if($type{$z} eq "P" and !exists($redundant_protein{$z})){
		$P_number++;
	}elsif($type{$z} eq "N" and !exists($redundant_protein{$z})){
		$N_number++;
	}
}
print STDERR "Positives $P_number, Negatives $N_number\n";
my %hmm_tscore;
my %hmm_target;
open(HMM_POTENTIAL, ">$potential_hmm_list");
open(HMM_MISSED, ">$missed_hmm_list");
open(HMM,"$hmm_result");
while(my $x = <HMM>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $protein = $array[2];
	my $target = $array[0];
	my $evalue = $array[4];
	my $tvalue = (log($evalue)) * -1;
	$hmm_tscore{$protein} = $tvalue;
	$hmm_target{$protein} = $target; 
}
close HMM;
my @hmm_values = sort{$a<=>$b}(values %hmm_tscore);
my $hmm_value_vector = join "\t",@hmm_values;
print "Method\tThreshold\te-value\tTP\tFP\tTN\tFN\tTPR\tFPR\tFNR\tTNR\tPPV\tFOR\tLR+\tLR-\tACC\tFDR\tNPV\tMK\tDOR\tBA\tF1\tFM\tMCC\tJaccard\tPotential_BLAs\n";
for my $n(1..1000){
	my $threshold = get_subpercentile($n,$hmm_value_vector);
	my $evalue = exp($threshold * -1);
	my %hmm_type;
	my $hmm_TP = 0;
	my $hmm_FP = 0;
	my $hmm_TN = 0;
	my $hmm_FN = 0;
	my $hmm_potential_blas = 0;
	my %hmm_seen;
	for my $protein(sort keys %hmm_tscore){
		if($hmm_tscore{$protein} >= $threshold){
			if(exists($consensus_class{$protein}) and !exists($redundant_protein{$protein})){
				$hmm_seen{$protein} = "seen";
				if($hmm_target{$protein} eq $consensus_class{$protein} and $type{$protein} eq "P"){
					$hmm_type{$protein} = "TP";
					$hmm_TP++;
				}elsif($hmm_target{$protein} ne $consensus_class{$protein} and $type{$protein} eq "P"){
					$hmm_type{$protein} = "FN";
					$hmm_FN++;
				}else{
					$hmm_type{$protein} = "FP";
					#print "$protein\t$target\t$consensus_class{$protein}\t$type{$protein}\n";
					$hmm_FP++;
				}
			}elsif(!exists($redundant_protein{$protein})){
				print HMM_POTENTIAL "$protein\t$threshold\tPotential\n";
				$hmm_potential_blas++;
			}
		}
	}
	foreach my $p(sort keys %type){
		if(!exists($hmm_type{$p})){
			if($type{$p} eq "P" and !exists($redundant_protein{$p})){
				$hmm_type{$p} = "FN";
				$hmm_FN++;
				print HMM_MISSED "$p\t$threshold\tMissed\n";
			}elsif(!exists($redundant_protein{$p})){
				$hmm_type{$p} = "TN";
				$hmm_TN++;
			}
		}
	}
        my $PPV = "-";
        my $FOR = "-";
        my $LR_plus = "-";
        my $LR_minus = "-";
        my $FDR = "-";
        my $NPV = "-";
	my $MK = "-";
        my $DOR = "-";
        my $F1 = "-";
	my $FM = "-";
        my $MCC = "-";
        my $CSI = "-";
	my $TPR = $hmm_TP / $P_number;
	my $FPR = $hmm_FP / $N_number;
	my $FNR = $hmm_FN / $P_number;
	my $TNR = $hmm_TN / $N_number;
	if(($hmm_TP + $hmm_FP) > 0){
		$PPV = $hmm_TP / ($hmm_TP + $hmm_FP);
	}
	if(($hmm_TN + $hmm_FN) > 0){
		$FOR = $hmm_FN / ($hmm_TN + $hmm_FN);
	}
	if($FPR > 0){
		$LR_plus = $TPR / $FPR;
	}
	if($TNR > 0){
		$LR_minus = $FNR / $TNR;
	}
	my $ACC	= ($hmm_TP + $hmm_TN) / ($P_number + $N_number);
	if(($hmm_TP + $hmm_FP) > 0){
		$FDR = $hmm_FP / ($hmm_TP + $hmm_FP);
	}
	if(($hmm_TN + $hmm_FN) > 0){
		$NPV = $hmm_TN / ($hmm_TN + $hmm_FN);
	}
	if($PPV ne "-" and $NPV ne "-"){
		$MK = $PPV + $NPV - 1;
	}
	if($LR_minus ne "-" and $LR_minus > 0){
		$DOR = $LR_plus / $LR_minus;
	}
	my $BA = ($TPR + $TNR) / 2;
	if((2 * $hmm_TP + $hmm_FP + $hmm_FN) > 0){
		$F1 = (2 * $hmm_TP) / (2 * $hmm_TP + $hmm_FP + $hmm_FN);
	}
	if($PPV ne "-" and $TPR ne "-"){
		$FM = ($PPV * $TPR)**(1/2);
	}
	if($TPR ne "-" and $TNR ne "-" and $PPV ne "-" and $NPV ne "-" and $FNR ne "-" and $FPR ne "-" and $FOR ne "-" and $FDR ne "-"){
		$MCC = (($TPR * $TNR * $PPV * $NPV)**(1/2)) - (($FNR * $FPR * $FOR * $FDR)**(1/2));
	}
	if(($hmm_TP + $hmm_FN + $hmm_FP) > 0){
		$CSI = $hmm_TP / ($hmm_TP + $hmm_FN + $hmm_FP);
	}
	print "hmmscan\t$threshold\t$evalue\t$hmm_TP\t$hmm_FP\t$hmm_TN\t$hmm_FN\t$TPR\t$FPR\t$FNR\t$TNR\t$PPV\t$FOR\t$LR_plus\t$LR_minus\t$ACC\t$FDR\t$NPV\t$MK\t$DOR\t$BA\t$F1\t$FM\t$MCC\t$CSI\t$hmm_potential_blas\n";
}
close HMM_POTENTIAL;
close HMM_MISSED;

sub get_subpercentile{
        my ($number,$values) = @_;
        my @array = split(/\t/,$values);
        my $len = scalar(@array);
        my $rank = int(($number/1000)*$len);
        my $index = $rank - 1;
        my $value = $array[$index];
        return $value;
}
