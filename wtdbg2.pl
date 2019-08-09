#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (t=>4, m=>"2g");
getopts('t:x:o:a:g:X:M:OPp:k:', \%opts);
die (qq/Usage: wtdbg2.pl [options] <reads.fa>
Options:
  -o STR     output prefix [first input]
  -t INT     number of threads [$opts{t}]
  -x STR     preset: rs, ont, sq, ccs []
  -g NUM     estimated genome size []
  -X FLOAT   coverage cutoff [wtdbg2 default]
  -k INT     k-mer length [wtdbg2 default]
  -p INT     length of HPC k-mers [wtdbg2 default]
  -a STR     additional wtdbg2 options []
  -M STR     minimap2 preset, according to -x when set []
  -m NUM     samtools sort block size [$opts{m}]
  -P         skip minimap2-based polishing
  -O         output without execution
/) if (@ARGV == 0);

$opts{asm} = gwhich("wtdbg2") || die;
$opts{cns} = gwhich("wtpoa-cns") || die;
unless (defined $opts{P}) {
	$opts{smt} = gwhich("samtools") || die;
	$opts{mm2} = gwhich("minimap2") || die;
}

my $prefix = defined($opts{o})? $opts{o} : $ARGV[0];
my $smt_threads = $opts{t} < 4? 4 : $opts{t};

my $asm_opt = "";
$asm_opt .= " -x $opts{x}" if defined($opts{x});
$asm_opt .= " -g $opts{g}" if defined($opts{g}) && $opts{g} =~ /^\d/;
$asm_opt .= " -X $opts{X}" if defined($opts{X});
$asm_opt .= " -k $opts{k}" if defined($opts{k}) && $opts{k} =~ /^\d+$/;
$asm_opt .= " -p $opts{p}" if defined($opts{p}) && $opts{p} =~ /^\d+$/;
$asm_opt .= " $opts{a}" if defined($opts{a});
$asm_opt .= " -t $opts{t} -fo $prefix";
$asm_opt .= " -i $ARGV[$_]" for (0 .. @ARGV-1);

my %map_opts = (
	rs=>'map-pb', rsII=>'map-pb', sq=>'map-pb', sequel=>'map-pb',
	ont=>'map-ont', nanopore=>'map-ont',
	ccs=>'asm20', corrected=>'asm20'
);

if (not defined $opts{M}) {
	if (defined $opts{x}) {
		$opts{M} = $map_opts{lc $opts{x}};
		$opts{M} = 'map-ont' unless(defined $opts{M});
	} else {
		$opts{M} = 'map-ont';
	}
}

my $mm2_input = join(' ', @ARGV);

my @lines = ();
push(@lines, qq/$opts{asm} $asm_opt/);
push(@lines, qq/$opts{cns} -t $opts{t} -i $prefix.ctg.lay.gz -fo $prefix.raw.fa/);
unless (defined $opts{P}) {
	push(@lines, qq/$opts{mm2} -ax $opts{M} -t$opts{t} -r2k $prefix.raw.fa $mm2_input | $opts{smt} sort -m $opts{m} -\@$smt_threads -o $prefix.bam/);
	push(@lines, qq/$opts{smt} view -F0x900 $prefix.bam | $opts{cns} -t $opts{t} -d $prefix.raw.fa -i - -fo $prefix.cns.fa/);
}

if (defined $opts{O}) {
	print "$_\n" for (@lines);
} else {
	for (@lines) {
		system($_) == 0 || die "The following command returns non-zero code:\n  $_\n";
	}
}

sub which {
	my $file = shift;
	my $path = (@_)? shift : $ENV{PATH};
	return if (!defined($path));
	foreach my $x (split(":", $path)) {
		$x =~ s/\/$//;
		return "$x/$file" if (-x "$x/$file") && (-f "$x/$file");
	}
	return;
}

sub gwhich {
    my $progname = shift;
    my $addtional_path = shift if (@_);
    my $dirname = &dirname($0);
    my $tmp;

    chomp($dirname);
    if ($progname =~ /^\// && (-x $progname) && (-f $progname)) {
        return $progname;
    } elsif (defined($addtional_path) && ($tmp = &which($progname, $addtional_path))) {
        return $tmp;
    } elsif (defined($dirname) && (-x "$dirname/$progname") && (-f "$dirname/$progname")) {
        return "$dirname/$progname";
    } elsif ((-x "./$progname") && (-f "./$progname")) {
        return "./$progname";
    } elsif (($tmp = &which($progname))) {
        return $tmp;
    } else {
        return;
    }
}

sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
