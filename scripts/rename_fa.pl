#!/usr/bin/perl -w
#
# Author: Ruan Jue
#
use strict;
use Getopt::Std;

my $prefix = $ENV{'PARAM_RENAME_FA_PREFIX'} || 'S';
my $suffix = $ENV{'PARAM_RENAME_FA_SUFFIX'} || '';

our ($opt_p, $opt_s, $opt_h, $opt_f, $opt_b, $opt_I);

getopts("hp:s:f:b:I");
die("Usage: $0 [-p name_prefix] [-s name_suffix] [-f trans_file] [-b begin_idx] [-I:include orignial_name] <fasta_file>\n") if($opt_h);
$prefix = $opt_p if(defined $opt_p);
$suffix = $opt_s if(defined $opt_s);
my %hash;
if(defined $opt_f){
	open(IN, "<", $opt_f) or die;
	%hash = ();
	while(<IN>){
		my @ts = split;
		$hash{$ts[0]} = $ts[1];
		#print STDERR "$ts[0]\t$ts[1]\n";
	}
	close IN;
}

my $idx = $opt_b || 0;

while(<>){
	if(/^>(\S+)/){
		my $desc = substr($_, length($1) + 1);
		my $ori = $opt_I? " $1" : "";
		if(%hash){
			if(exists $hash{$1}){
				my $tag = $hash{$1};
				print ">$tag$ori$desc", substr($_, length($1) + 1);
			} else {
				print;
			}
		} else {
			printf(">$prefix%010d$suffix$ori$desc", $idx);
		}
		$idx ++;
	} else {
		print;
	}
}

1;
