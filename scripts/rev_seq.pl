#!/usr/bin/perl -w
#
#AUthor: Ruan Jue
#
use strict;

my $seq = shift;

if($seq){
	rev_seq($seq);
} else {
	$seq = '';
	while(<>){
		if(/^>/){
			&rev_seq($seq);
			print; next;
			$seq = '';
		}
		chomp;
		$seq .= $_;
	}
	&rev_seq($seq);
}

1;

sub rev_seq {
	my $s = shift;
	$s =~tr/ACGTacgt/TGCAtgca/;
	$s = reverse $s;
	while($s=~/(.{1,100})/g){
		print "$1\n";
	}
}
