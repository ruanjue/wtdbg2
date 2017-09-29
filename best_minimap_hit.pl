#!/usr/bin/perl -w
#
#Author: Jue Ruan
#
use strict;

my $rs = [""];

while(<>){
	my @ts = split;
	if($ts[0] ne $rs->[0]){
		if(length $rs->[0]){
			print join("\t", @{$rs}), "\n";
		}
		$rs = \@ts;
	} else {
		if($ts[9] > $rs->[9]){
			$rs = \@ts;
		}
	}
}
if(length $rs->[0]){
	print join("\t", @{$rs}), "\n";
}

1;
