#!/usr/bin/perl -w
#
#Author: Ruan Jue
#
use strict;
use Getopt::Std;
use DB_File;

our ($opt_l, $opt_h);

getopts("l:h");

my $level = $opt_l || 0;
&usage if($opt_h);

my $dbf = shift or &usage;
if($dbf!~/\.dbm$/){
	$dbf .= ".dbm" if(-e "$dbf.dbm");
}

my @nodes = @ARGV;

if(@nodes == 0){
	while(<>){
		chomp;
		push(@nodes, $_);
	}
}

my %hash;

tie %hash, 'DB_File', $dbf, O_RDONLY or die "Cannot open $dbf: $!";

my %levels = map {$_=>0} @nodes;

print "digraph {\n";

print "rank=$hash{rank}\n" if(exists $hash{"rank"});
print "node $hash{node}\n" if(exists $hash{"node"});
print "edge $hash{edge}\n" if(exists $hash{"edge"});

while(@nodes){
	my $nd = shift @nodes;
	my $str = $hash{$nd};
	my @rs = split /\n/, $str;
	print "$nd $rs[0]";
	if($levels{$nd} == 0){
		print " [style=filled fillcolor=yellow]\n"
	} else {
		print "\n";
	}
	my $n2 = '';
	for(my $i=1;$i<@rs;$i++){
		if($rs[$i]=~/^N\s(\S+)/){
			$n2 = $1;
			if($levels{$nd} < $level and not exists $levels{$n2}){
				$levels{$n2} = $levels{$nd} + 1;
				push(@nodes, $n2);
			}
		} elsif($rs[$i]=~/^L\s(.+)$/){
			print " $nd -> $n2 $1\n"
		}
	}
}

print "}\n";

untie %hash;

1;

sub usage {
	die("Usage: $0 [-l TRACE_LEVEL:0] <dot_dbm_file> <node1> ...\n");
}
