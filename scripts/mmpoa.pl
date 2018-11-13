#!/usr/bin/env perl
#
# Author: Jue Ruan <ruanjue@gmail.com>
#
use Getopt::Std;
use strict;

our ($opt_h, $opt_t, $opt_s, $opt_m, $opt_x, $opt_p);

getopts("ht:s:m:x:p:");

&usage if($opt_h);
&usage(1) if(@ARGV < 2);
my $ref = shift || &usage(1);
my $rst = $ref;
$rst=~s/\.fa$//;
$rst=~s/\.fasta$//;

my $ncpu = 4;
$ncpu = $opt_t if(defined $opt_t);

my $MM = $opt_m || "minimap2";
my $MX = $opt_x || "pb";
my $ST = $opt_s || "samtools";
my $WP = $opt_p || "wtpoa-cns";

my $cmd = '';

$cmd = "$MM -t $ncpu -x map-$MX -a $ref @ARGV | $ST view -Sb - > $rst.bam";
&run($cmd);

$cmd = "$ST sort $rst.bam $rst.srt";
&run($cmd);

$cmd = "$ST view $rst.srt.bam | $WP -t $ncpu -d $ref -i - -fo $rst.mmpoa.fa";
&run($cmd);

1;

sub usage {
	my $ret = shift || 0;
	print qq{Usage: $0 [-t n_cpu:4] [-s samtools] [-m minimap2] [-p wtpoa-cns] [-x pb|ont] <ref.fa> <reads1.fa> [reads2.fa ...]\n};
	exit $ret;
}

sub run {
	my $cmd = shift;
	system("date");
	print "# $cmd\n";
	if(system($cmd)){
		die("$cmd");
	}
	system("date");
}
