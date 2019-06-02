#!/usr/bin/perl -w
#
# Author: Jue Ruan
#
use Getopt::Std;
use strict;

our ($opt_h, $opt_l, $opt_f);

getopts("hl:f:");
&usage if($opt_h);

my $map_len = $opt_l || 100;
my $map_cov = $opt_f || 0.75;

my %refs = ();

my %cigar_hash = (M=>[1,1], I=>[1,0], D=>[0,1], N=>[0,1], S=>[1,0], H=>[1,0], P=>[0,0], '='=>[1,1], X=>[1,1]);

my @hits = ();

while(<>){
	if(/^@/){
		print;
		if(/^\@SQ\sSN:(\S+)\sLN:(\d+)/){
			$refs{$1} = $2;
		}
		next;
	}
	my $str = $_;
	my @ts = split;
	my $tag = $ts[0];
	my $flg = $ts[1];
	my $ref = $ts[2];
	next if($ref eq '*');
	my $rlen = $refs{$ref};
	if(not defined $rlen){
		print;
		die("Cannot find '$ref' in SAM Header");
	}
	my $tb = $ts[3];
	my $te = $tb;
	my $qb = 0;
	my $qe = 0;
	my $len = 0;
	my $tmp = 0;
	my $cnt = 0;
	my $cgr = $ts[5];
	while($cgr=~/(\d+)(\D)/g){
		my $mov = $cigar_hash{$2};
		if($mov->[1]){
			$te += $1;
			$qe += $tmp;
			$tmp = 0;
			if($cnt == 0){
				$qb = $qe;
				$cnt = 1;
			}
			if($mov->[0]){
				$qe += $1;
				$len += $1;
			}
		} elsif($mov->[0]){
			$tmp += $1;
			$len += $1;
		}
	}
	if($flg & 0x10){
		$tmp = $len - $qb;
		$qb = $len - $qe;
		$qe = $tmp;
	}
	next unless($qe - $qb >= $map_len and $te - $tb >= $map_len);
	next unless($qe - $qb >= $map_cov * $len or $te - $tb >= $map_cov * $rlen);
	my $hit = [$tag, $len, $qb, $qe, $ref, $rlen, $tb, $te, $str];
	&select_best_hit($hit);
}
&select_best_hit(undef);

1;

sub usage {
	print STDERR qq{Usage: $0 [-l 100:map_len] [-f 0.75:map_cov] <sam_file_in_query_order_with_header>\n};
	exit 1;
}

sub select_best_hit {
	my $hit = shift;
	if(@hits == 0){
		push(@hits, $hit) if($hit);
	} elsif(not defined $hit or $hit->[0] ne $hits[-1][0]){
		@hits = sort {$b->[3] - $b->[2] <=> $a->[3] - $a->[2]} @hits;
		for(my $i=0;$i<@hits;$i++){
			my $pass = 1;
			for(my $j=0;$j<$i;$j++){
				my $x = $hits[$i][2] < $hits[$j][2]? $hits[$j][2] : $hits[$i][2];
				my $y = $hits[$i][3] > $hits[$j][3]? $hits[$j][3] : $hits[$i][3];
				if($y - $x >= (1 - $map_cov) * ($hits[$i][3] - $hits[$i][2])){
					$pass = 0;
					last;
				}
			}
			print $hits[$i][8] if($pass);
		}
		@hits = ();
		push(@hits, $hit) if($hit);
	} else {
		push(@hits, $hit);
	}
}

