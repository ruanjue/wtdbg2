#!/usr/bin/perl -w
#
#Author: Ruan Jue
#
use strict;
use DB_File;

my $fasta_file = shift or die("Usage: $0 <fasta_file>\n");

die("$fasta_file.dbm already exists!!!") if(-e "$fasta_file.dbm");

open(IN, "<", $fasta_file) or die;
my %seqs;

tie %seqs, 'DB_File', "$fasta_file.dbm", O_RDWR | O_CREAT, 0644, $DB_HASH or die "Cannot open $fasta_file.dbm: $!";

my $name = '';
#my $seq = '';

my $off = 0;
my $len = 0;

while(<IN>){
	if(/^>(\S+)/){
		#$seqs{$name} = $seq if($name);
		$seqs{$name} = $off;
		$name = $1;
		$off += $len;
		$len = 0;
		#$seq = '';
	#} else {
		#chomp;
		#$seq .= $_;
	}
	$len += length($_);
}
$seqs{$name} = $off if($name);

untie %seqs;

close IN;

1;
