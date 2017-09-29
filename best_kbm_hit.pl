#!/usr/bin/perl -w
#
# Author: Jue Ruan
#
use strict;

my $rs = [""];

while(<>){
	my @ts = split;
	if($ts[0] eq $rs->[0]){
		$rs = \@ts if($ts[10] > $rs->[10]);
	} else {
		print join("\t", @{$rs}), "\n" if($rs->[10]);
		$rs = \@ts;
	}
}
print join("\t", @{$rs}), "\n" if($rs->[10]);

1;
