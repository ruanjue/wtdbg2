#!/usr/bin/perl -w
#
# Author: Ruan Jue
#
use strict;

my $total = 0;
my @nums = ();

my $len = 0;
while(<>){
	next unless(/(\d+)/);
	push(@nums, $1),$total+=$1 if($1);
}

print "Total: $total\n";

my @nxxs = ();
for(my $i=0;$i<=10;$i++){
	push(@nxxs, int($total*$i*0.1));
}
push(@nxxs, $total + 1);

my $i = 0;
my $j = 0;

@nums = sort {$b <=> $a} @nums;

$len = 0;

for(;$i<@nums;$i++){
	$len += $nums[$i];
	while($nxxs[$j] <= $len){
		print "N".$j."0: $nums[$i]\t" . ($i + 1) . "\n";
		$j ++;
	}
}

1;
