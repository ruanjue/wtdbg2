#!/usr/bin/perl -w
#
#Author: Ruan Jue
#
use strict;
use DB_File;

my $dot_file = shift or die("Usage: $0 <dot_file>\n");

die("$dot_file.dbm already exists!!!") if(-e "$dot_file.dbm");

open(IN, "<", $dot_file) or die;
my %hash;

tie %hash, 'DB_File', "$dot_file.dbm", O_RDWR | O_CREAT, 0644, $DB_HASH or die "Cannot open $dot_file.dbm: $!";

my %nodes = ();
my %link = ();

while(<IN>){
	s/^\s+//;
	s/\s+$//;
	if(/^rankdir\s*=\s*(\S+)/){
		$hash{"rankdir"} = $1;
		next;
	}
	my $desc = '';
	while(1){
		if(/\s*(\[[^]]+\]);?$/){
			$_ = substr($_, 0, length($_) - length($1));
			$desc .= $1;
			s/\s+$//;
		} else {
			last;
		}
	}
	my @ts = split;
	if(@ts == 1 and length $desc){
		$nodes{$ts[0]} = $desc;
	} elsif(@ts >= 3 and ($ts[1] eq '->' or $ts[1] eq '-')){
		my ($lnk1, $lnk2) = ("", "");
		if($ts[0]=~/^(\S+?):(\S+)$/){ $ts[0] = $1; $lnk1 = $2; }
		if($ts[2]=~/^(\S+?):(\S+)$/){ $ts[2] = $1; $lnk2 = $2; }
		push(@{$link{$ts[0]}{$ts[2]}}, [$desc, $lnk1, $lnk2]);
		$link{$ts[2]}{$ts[0]} = [] unless(defined $link{$ts[2]}{$ts[0]});
	}
}

close IN;

foreach my $n1 (keys %nodes){
	$hash{$n1} = $nodes{$n1};
}

foreach my $n1 (keys %link){
	my $hx = $link{$n1};
	my $str = (defined $hash{$n1})? $hash{$n1} : "";
	foreach my $n2 (keys %{$hx}){
		my $ls = $hx->{$n2};
		$str .= "\nN\t$n2";
		foreach my $lk (@{$ls}){
			$str .= "\n" . join("\t", "L", @{$lk});
		}
	}
	$hash{$n1} = $str;
}

untie %hash;

1;
