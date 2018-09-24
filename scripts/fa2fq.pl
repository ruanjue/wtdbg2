#!/usr/bin/perl -w
#
#
#
use strict;

my $Q = '!';

my $tag = '';
my $seq = '';
while(<>){
	if(/^>(.+?)$/){
		if(length $seq){
			print "\@$tag\n";
			print "$seq\n";
			print "+\n";
			print $Q x length($seq);
			print "\n";
		}
		$tag = $1;
		$seq = '';
	} else {
		chomp;
		$seq .= $_;
	}
}

if(length $seq){
	print "\@$tag\n";
	print "$seq\n";
	print "+\n";
	print $Q x length($seq);
	print "\n";
}

1;
