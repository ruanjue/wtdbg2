#!/usr/bin/perl -w
#
#
use strict;

my $cnt = shift or die("Usage: $0 <parts> <index> <fasta_file>\n");
my $idx = shift or die("Usage: $0 <parts> <index> <fasta_file>\n");

$idx --;

my $ns = 0;
my $print = 0;

while(<>){
	if(/^>/){
		if($ns % $cnt == $idx){
			$print = 1;
		} else {
			$print = 0;
		}
		$ns ++;
	}
	print if($print);
}

1;
