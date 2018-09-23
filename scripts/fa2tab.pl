#!/usr/bin/perl -w
#
use strict;

my $name;
my $seq = '';

while(<>){
	if(/>(\S+)/){
		print "$name\t$seq\n" if(defined $name);
		$name = $1;
		$seq = '';
	} else {
		chomp;
		$seq .= $_;
	}
}
print "$name\t$seq\n" if(defined $name);

1;
