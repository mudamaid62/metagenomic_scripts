#!/usr/bin/perl
use warnings;
use strict;

my $cons_table = shift(@ARGV);
my $mmseqs_result = shift(@ARGV);
my $potential_mmseqs_list = shift(@ARGV);
my $missed_mmseqs_list = shift(@ARGV);

my %structural_cluster;
my %redundant_protein;
my %ref_class;
my %consensus_class;
my %type;
#in consensus class
my @positives = ('A','B1','B2','B3','C','D');
my @negatives = ('PBP','MBL');
my %references;

open(CONSENSUS,"$cons_table");
while(my $x = <CONSENSUS>){
	chomp($x);
	my @array = split(/\t/,$x);
	if($array[4] ne "Candidate"){
                my $class = $array[5];
		$references{$array[0]} = $class;
		my @redundant;
		if($array[15] ne "-"){
                	@redundant = split(/\;/,$array[15]);
                	foreach my $z(@redundant){
                        	$references{$z} = $class;
			}
                }
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
my %mmseqs_tscore;
my %mmseqs_target;
open(MMSEQS_POTENTIAL, ">$potential_mmseqs_list");
open(MMSEQS_MISSED, ">$missed_mmseqs_list");
open(MMSEQS,"$mmseqs_result");
while(my $x = <MMSEQS>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $protein = $array[0];
	my $target = $array[2];
	my $evalue = $array[8];
	my $target_type = $references{$target};
	my $tvalue = (log($evalue)) * -1;
	$mmseqs_tscore{$protein} = $tvalue;
	$mmseqs_target{$protein} = $target_type; 
	#print STDERR "$mmseqs_target{$protein}\t$mmseqs_tscore{$protein}\n";
}
close MMSEQS;
my @mmseqs_values = sort{$a<=>$b}(values %mmseqs_tscore);
my $mmseqs_value_vector = join "\t",@mmseqs_values;
print "Method\tThreshold\te-value\tTP\tFP\tTN\tFN\tTPR\tFPR\tFNR\tTNR\tPPV\tFOR\tLR+\tLR-\tACC\tFDR\tNPV\tMK\tDOR\tBA\tF1\tFM\tMCC\tJaccard\tPotential_BLAs\n";
for my $n(1..1000){
	my $threshold = get_subpercentile($n,$mmseqs_value_vector);
	my $evalue = exp($threshold * -1);
	my %mmseqs_type;
	my $mmseqs_TP = 0;
	my $mmseqs_FP = 0;
	my $mmseqs_TN = 0;
	my $mmseqs_FN = 0;
	my $mmseqs_potential_blas = 0;
	my %mmseqs_seen;
	for my $protein(sort keys %mmseqs_tscore){
		if($mmseqs_tscore{$protein} >= $threshold){
			if(exists($consensus_class{$protein}) and !exists($redundant_protein{$protein})){
				$mmseqs_seen{$protein} = "seen";
				if($mmseqs_target{$protein} eq $consensus_class{$protein} and $type{$protein} eq "P"){
					$mmseqs_type{$protein} = "TP";
					$mmseqs_TP++;
				}elsif($mmseqs_target{$protein} ne $consensus_class{$protein} and $type{$protein} eq "P"){
					$mmseqs_type{$protein} = "FN";
					$mmseqs_FN++;
				}else{
					$mmseqs_type{$protein} = "FP";
					#print "$protein\t$target\t$consensus_class{$protein}\t$type{$protein}\n";
					$mmseqs_FP++;
				}
			}elsif(!exists($redundant_protein{$protein})){
				print MMSEQS_POTENTIAL "$protein\t$threshold\tPotential\n";
				$mmseqs_potential_blas++;
			}
		}
	}
	foreach my $p(sort keys %type){
		if(!exists($mmseqs_type{$p})){
			if($type{$p} eq "P" and !exists($redundant_protein{$p})){
				$mmseqs_type{$p} = "FN";
				$mmseqs_FN++;
				print MMSEQS_MISSED "$p\t$threshold\tMissed\n";
			}elsif(!exists($redundant_protein{$p})){
				$mmseqs_type{$p} = "TN";
				$mmseqs_TN++;
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
	my $TPR = $mmseqs_TP / $P_number;
	my $FPR = $mmseqs_FP / $N_number;
	my $FNR = $mmseqs_FN / $P_number;
	my $TNR = $mmseqs_TN / $N_number;
	if(($mmseqs_TP + $mmseqs_FP) > 0){
		$PPV = $mmseqs_TP / ($mmseqs_TP + $mmseqs_FP);
	}
	if(($mmseqs_TN + $mmseqs_FN) > 0){
		$FOR = $mmseqs_FN / ($mmseqs_TN + $mmseqs_FN);
	}
	if($FPR > 0){
		$LR_plus = $TPR / $FPR;
	}
	if($TNR > 0){
		$LR_minus = $FNR / $TNR;
	}
	my $ACC	= ($mmseqs_TP + $mmseqs_TN) / ($P_number + $N_number);
	if(($mmseqs_TP + $mmseqs_FP) > 0){
		$FDR = $mmseqs_FP / ($mmseqs_TP + $mmseqs_FP);
	}
	if(($mmseqs_TN + $mmseqs_FN) > 0){
		$NPV = $mmseqs_TN / ($mmseqs_TN + $mmseqs_FN);
	}
	if($PPV ne "-" and $NPV ne "-"){
		$MK = $PPV + $NPV - 1;
	}
	if($LR_minus ne "-" and $LR_minus > 0){
		$DOR = $LR_plus / $LR_minus;
	}
	my $BA = ($TPR + $TNR) / 2;
	if((2 * $mmseqs_TP + $mmseqs_FP + $mmseqs_FN) > 0){
		$F1 = (2 * $mmseqs_TP) / (2 * $mmseqs_TP + $mmseqs_FP + $mmseqs_FN);
	}
	if($PPV ne "-" and $TPR ne "-"){
		$FM = ($PPV * $TPR)**(1/2);
	}
	if($TPR ne "-" and $TNR ne "-" and $PPV ne "-" and $NPV ne "-" and $FNR ne "-" and $FPR ne "-" and $FOR ne "-" and $FDR ne "-"){
		$MCC = (($TPR * $TNR * $PPV * $NPV)**(1/2)) - (($FNR * $FPR * $FOR * $FDR)**(1/2));
	}
	if(($mmseqs_TP + $mmseqs_FN + $mmseqs_FP) > 0){
		$CSI = $mmseqs_TP / ($mmseqs_TP + $mmseqs_FN + $mmseqs_FP);
	}
	print "MMSeqs2\t$threshold\t$evalue\t$mmseqs_TP\t$mmseqs_FP\t$mmseqs_TN\t$mmseqs_FN\t$TPR\t$FPR\t$FNR\t$TNR\t$PPV\t$FOR\t$LR_plus\t$LR_minus\t$ACC\t$FDR\t$NPV\t$MK\t$DOR\t$BA\t$F1\t$FM\t$MCC\t$CSI\t$mmseqs_potential_blas\n";
}
close MMSEQS_POTENTIAL;
close MMSEQS_MISSED;

sub get_subpercentile{
        my ($number,$values) = @_;
        my @array = split(/\t/,$values);
        my $len = scalar(@array);
        my $rank = int(($number/1000)*$len);
        my $index = $rank - 1;
        my $value = $array[$index];
        return $value;
}
