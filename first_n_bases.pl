#!/usr/bin/perl -w
#
use strict;

my $n = shift or die("Usage: $0 <num_of_bases> [fasta_file]\n");

my $m = 0;

while(<>){
	if(/^>/){
		last if($m >= $n);
	} else {
		$m += length($_) - 1;
	}
	print;
}

1;
