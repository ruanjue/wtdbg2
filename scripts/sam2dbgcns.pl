#!/usr/bin/perl -w
#
# Author: Jue Ruan <ruanjue@gmail.com>
#
use strict;
use Getopt::Std;

our ($opt_h, $opt_w, $opt_s, $opt_d, $opt_q, $opt_t);

getopts("hw:s:dq:t");

&usage if($opt_h);
$opt_w = 2000 unless($opt_w);
$opt_s = 1000 unless($opt_s);
$opt_d = 0 unless($opt_d);
$opt_q = 0 unless($opt_q);
$opt_t = 0 unless($opt_t);

&usage if($opt_w <= $opt_s);

my $ctgf = shift or &usage;

my %refs = ();

my $ctag = "";
my $cdes = "";
my $cseq = "";

open(IN, "<", $ctgf) or die;
while(<IN>){
	chomp;
	if(/^>(\S+)/){
		$refs{$ctag} = [$cseq, $cdes] if(length $cseq);
		$ctag = $1;
		$cdes = substr($_, 1 + length($1));
		$cseq = '';
	} else {
		$cseq .= $_;
	}
}
$refs{$ctag} = [$cseq, $cdes] if(length $cseq);
close IN;

$ctag = "";
my $coff = 0;
my $wlen = $opt_w;
my $slen = $opt_s;
my @rseqs = ();
my $ref = undef;
my $rsize = 0;

my $nid = 0;

my $line = 0;

while(<>){
	$line ++;
	next if(/^@/);
	my @ts= split;
	next if($ts[2] eq '*');
	next if($opt_q and $ts[4] < $opt_q);
	my $rtag = $ts[0];
	my $rdir = (($ts[1] >> 4) & 0x01)? "-" : "+";
	my $gtag = $ts[2];
	my $gpos = $ts[3];
	my $rseq = $ts[9];
	if($opt_t){
		my $l = 0;
		my $r = 0;
		if($ts[5]=~/^((\d+)S)/){
			$l = $2;
		}
		if($ts[5]=~/((\d+)S)$/){
			$r = $2;
		}
		$rseq = substr($rseq, $l, length($rseq) - $r - $l) if($l or $r);
	}
	if($gtag ne $ctag){
		while(length $ctag and $coff < $rsize){
			&print_frag;
		}
		$ctag = $gtag;
		#print STDERR "$ctag begin at SAM line $line\n";
		$ref = $refs{$ctag};
		$rsize = length $ref->[0];
		$coff = 0;
		@rseqs = ();
		&print_header;
	} else {
		while($coff + $wlen < $gpos){
			&print_frag;
		}
	}
	my $dup = 0;
	if($opt_d){
		for(my $i=@rseqs-1;$i>=0;$i--){
			last if($rseqs[$i][0] != $gpos);
			if($rseqs[$i][3] eq $rseq){
				$dup = 1;
				last;
			}
		}
	}
	push(@rseqs, [$gpos, $rtag, $rdir, $rseq]) unless($dup);
}

while(length $ctag and $coff < $rsize){
	&print_frag;
}

1;

sub print_header {
	return unless(length $ctag and exists $refs{$ctag});
	print ">$ctag$refs{$ctag}->[1]\n";
}

sub print_frag {
	my $cend = $coff + $wlen;
	my $cnxt = $coff + $slen;
	$cend = $rsize if($cend > $rsize);
	$cnxt = $rsize if($cnxt > $rsize);
	print "E\t$coff\tN$nid\t+\t";
	$nid ++;
	print "N$nid\t+\n";
	$nid ++;
	print "S\t$ctag\_F_$coff\_", ($cend - $coff), "\t+\t$coff\t", ($cend - $coff),"\t", substr($ref->[0], $coff, $cend - $coff), "\n";
	my @nseqs = ();
	foreach my $r (@rseqs){
		next if($r->[0] + length($r->[3]) < $coff);
		if($r->[0] < $cend){
			print "S\t$r->[1]\t$r->[2]\t0\t", length($r->[3]), "\t$r->[3]\n";
			if($r->[0] >= $cnxt){
				push(@nseqs, $r);
			}
		}
	}
	@rseqs = @nseqs;
	$coff = $cnxt;
}

sub usage {
	die(qq{Usage: $0 [options] <ctg_fa_file> [srt_sam_file]
	Option:
	-w <int> window size, [2000]
	-s <int> silding size, [1000]
	-d       remove duplication
	-q <int> filter by mapQ, [0]
	-t <int> trim tailing clip
	});
}
