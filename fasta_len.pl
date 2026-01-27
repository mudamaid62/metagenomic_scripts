#!/usr/bin/perl
use warnings;
use strict;

my $file = shift(@ARGV);
my $bin_width = shift(@ARGV);

my @lens;
my $total_len = 0;
print "Sequence\tLength\n";
my $current_max_len = 0;
my $current_min_len = 0;
my $current_max = "";
my $current_min = "";
my $current_header  = "";
my $current_seq = "";

open(FILE,"$file");
while(my $x = <FILE>){
	chomp($x);
	if($x =~ m/^>/){
		if($current_seq ne ""){
                	$current_seq =~ s/\*//g;
                	my $len = length($current_seq);
                	push @lens, $len;
			$total_len += $len;
                	print "$current_header\t$len\n";
                	if($current_max eq "" and $current_min eq ""){
                        	$current_max = $current_header;
                        	$current_min = $current_header;
                        	$current_max_len = $len;
                        	$current_min_len = $len;
                	}elsif($len > $current_max_len){
                        	$current_max = $current_header;
                        	$current_max_len = $len;
                	}elsif($len < $current_min_len){
                        	$current_min = $current_header;
                        	$current_min_len = $len;
                	}
			$current_seq = "";
		}
		$x =~ s/>//g;
		$current_header = $x;
	}else{
		$current_seq .= $x;
	}
}
if($current_seq ne ""){
	$current_seq =~ s/\*//g;
	my $len = length($current_seq);
	push @lens, $len;
	$total_len += $len;
	print "$current_header\t$len\n";
	if($current_max eq "" and $current_min eq ""){
		$current_max = $current_header;
		$current_min = $current_header;
		$current_max_len = $len;
		$current_min_len = $len;
	}elsif($len > $current_max_len){
		$current_max = $current_header;
		$current_max_len = $len;
	}elsif($len < $current_min_len){
		$current_min = $current_header;
		$current_min_len = $len;
	}
	$current_seq = "";
}
my @sorted = sort{$a<=>$b} @lens;
my $n = scalar(@sorted);
if($n <= 1){
	print "Total length -> $total_len\n";
	exit;
}
my $min = $sorted[0];
my $max = $sorted[$n - 1];
my $new_n = ($max - $min)/$bin_width;
my $new_n_rest = ($max - $min)%$bin_width;
unless($new_n == 0){
	$new_n = int(($max - $min)/$bin_width) + 1;
}
my $interval = $min + $bin_width;
my @lefts;
my @rights;
my @freqs;
until((scalar(@rights)) == $new_n){
	my $left = $interval - $bin_width;
	my $right = $interval;
	push @lefts, $left;
	push @rights, $right;
	$interval = $interval + $bin_width;
}
my $m = (scalar(@rights)) - 1;
for my $k(0..$m){
	$freqs[$k] = 0;
}
foreach my $y(@sorted){
	for my $i(1..$m){
		if($y > $lefts[$i] and $y <= $rights[$i]){
			$freqs[$i]++;
		}
	}
	if($y >= $lefts[0] and $y <= $rights[0]){
		$freqs[0]++;
	}
}
my @sorted_freqs = sort{$b<=>$a}(@freqs);
my $max_freq = shift(@sorted_freqs);
my $half_freq = int($max_freq/2);
print STDERR "Interval\tfi\tfr\%\tFi\tFr\%\n";
my $total = scalar(@sorted);
my $Fi = 0;
my $Fr = 0;
#my $tick = "\*";
my $tick = "\033[38;5;214m\200\033[38;5;15m";
my $space = " ";
my $max_len = length($max_freq);
my $h_width = (int($max_freq))/10;
my $pre_space = " " x $max_len;
my $h1 = $pre_space."|";
my $h2 = $pre_space."|";
my $h3 = $pre_space."|";
my $h4 = $pre_space."|";
my $h5_len = length($half_freq);
my $h5_rest = " " x ($max_len - $h5_len);
my $h5 = $h5_rest.$half_freq."|";
my $h6 = $pre_space."|";
my $h7 = $pre_space."|";
my $h8 = $pre_space."|";
my $h9 = $pre_space."|";
my $h10 = $max_freq."|";
my $down_rest = $new_n - 2;
my $down2 = "";
if($down_rest > 0){
	$down2 = "-" x $down_rest;
}
my $down3 = "0".$down2."4";
my $down5 = " " x ($max_len - 1);
my $down4 = $down5."Q|".$down3;
for my $j(0..$m){
	$Fi += $freqs[$j];
	my $fr = ($freqs[$j] / $total)*100;
	my $fi = $freqs[$j];
	$Fr += $fr;
	if($j != 0){ 
		print STDERR "\($lefts[$j]\,$rights[$j]\]\t$freqs[$j]\t$fr\t$Fi\t$Fr\n";
	}else{
		print STDERR "\[$lefts[$j]\,$rights[$j]\]\t$freqs[$j]\t$fr\t$Fi\t$Fr\n";
	}
	if($fi >= $h_width){
		$h1 .= $tick;
	}else{
		$h1 .= $space;
	}
	if($fi >= $h_width*2){
                $h2 .= $tick;
        }else{
                $h2 .= $space;
        }	
	if($fi >= $h_width*3){
                $h3 .= $tick;
        }else{
                $h3 .= $space;
        }
	if($fi >= $h_width*4){
                $h4 .= $tick;
        }else{
                $h4 .= $space;
        }
	if($fi >= $h_width*5){
                $h5 .= $tick;
        }else{
                $h5 .= $space;
        }
	if($fi >= $h_width*6){
                $h6 .= $tick;
        }else{
                $h6 .= $space;
        }
        if($fi >= $h_width*7){
                $h7 .= $tick;
        }else{
                $h7 .= $space;
        }
        if($fi >= $h_width*8){
                $h8 .= $tick;
        }else{
                $h8 .= $space;
        }
        if($fi >= $h_width*9){
                $h9 .= $tick;
        }else{
                $h9 .= $space;
        }
        if($fi >= $h_width*10){
                $h10 .= $tick;
        }else{
                $h10 .= $space;
        }
}
print STDERR "$h10\n$h9\n$h8\n$h7\n$h6\n$h5\n$h4\n$h3\n$h2\n$h1\n$down4\n";
print STDERR "min --> $min,max --> $max, Interval_number --> $new_n, Interval_width --> $bin_width\n";
print STDERR "Current_max --> $current_max, Current_min --> $current_min\n";
print STDERR "Total length -> $total_len\n";
