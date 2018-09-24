#!/usr/bin/perl -w
#
# Author: Ruan Jue
#
use strict;
use Getopt::Std;

our ($opt_h, $opt_s);

getopts("hs");

my $total = 0;
my @nums = ();

my $len = 0;
while(<>){
	if(/^>(\S+)/){
		if($opt_s){
			print "\t$len\n" if($len);
			print $1;
		}
		push(@nums, $len),$total+=$len if($len);
		$len = 0;
	} else {
		$len += length($_) - 1;
	}
}
print "\t$len\n" if($opt_s and $len);
push(@nums, $len),$total+=$len if($len);

my $n_seq = @nums;
my $avg = sprintf("%0.2f", $total / $n_seq);

print "Total: $total\n";
print "Count: $n_seq\n";
print "Average: $avg\n";

my @nxxs = ();
for(my $i=0;$i<=10;$i++){
	push(@nxxs, int($total*$i*0.1));
}
push(@nxxs, $total + 1);

my $i = 0;
my $j = 0;

@nums = sort {$b <=> $a} @nums;

my $median = $nums[int($n_seq / 2)];

print "Median: $median\n";

$len = 0;

for(;$i<@nums;$i++){
	$len += $nums[$i];
	while($nxxs[$j] <= $len){
		print "N".$j."0: $nums[$i]\t". ($i + 1) . "\n";
		$j ++;
	}
}

1;
