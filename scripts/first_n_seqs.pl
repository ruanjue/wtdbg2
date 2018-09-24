#!/usr/bin/perl -w
#
use strict;

my $n = shift or die("Usage: $0 <num_of_seqs> [fasta_file]\n");

my $m = 0;

while(<>){
	if(/^>/){
		last if($m == $n);
		$m ++;
	}
	print;
}

1;
