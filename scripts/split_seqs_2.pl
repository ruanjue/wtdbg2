#!/usr/bin/perl -w
#
#
use strict;

my $cnt = shift or die("Usage: $0 <n_files> <fasta_file> [gzip:0]\n");
my $inf = shift or die("Usage: $0 <n_files> <fasta_file> [gzip:0]\n");
my $zip = shift || 0;
$zip = 0 unless($zip eq '1' or $zip eq 'y');
my $ouf = $inf;
if($ouf eq '-'){
	$ouf = 'unnamed';
} else {
	unshift @ARGV, $inf;
}

my @fhs = ();

for(my $i=1;$i<=$cnt;$i++){
	my $fh;
	if($zip){
		open($fh, "| gzip -c >$ouf.shuffle$i.gz") or die;
	} else {
		open($fh, ">", "$ouf.shuffle$i") or die;
	}
	push(@fhs, $fh);
}

my $n = -1;
while(<>){
	$n ++ if(/^>/);
	my $fh = $fhs[$n % $cnt];
	print $fh $_;
}

foreach my $fh (@fhs) {
	close $fh;
}

1;
