#!/usr/bin/perl -w
#
#Author: Ruan Jue
#
use strict;
use DB_File;

my $dbf = shift or die("Usage: $0 <dbm_file> [read_name1 ...]\n");
my $faf = undef;
if($dbf!~/\.dbm$/){
	$dbf .= ".dbm" if(-e "$dbf.dbm");
}

my @names = @ARGV;

my $list_all = 0;

if(@names == 1 and $names[0] eq '#LIST'){
	$list_all = 1;
	@names = ();
} elsif(@names == 0){
	while(<>){
		chomp;
		push(@names, $_);
	}
}

my @tags = ();
foreach my $tag (@names){
	if($tag=~/^(.+?)(\[([+-])(:(-?\d+),(-?\d+))?\])$/){
		push(@tags, [$1, $3 eq '+'? 1:2, (defined $5)? $5:1, (defined $6)? $6:-1, 1]);
	} elsif($tag=~/^(.+?)_([FR])(_(\d+)(_(\d+))?)?$/){
		push(@tags, [$1, $2 eq 'F'? 1:2, (defined $4)? $4:1, (defined $6)? $4 + $6 - 1 : 0, 1]);
	} else {
		push(@tags, [$tag, 1, 1, -1, 0]);
	}
}

foreach my $tag (@tags){
	$tag->[2] = 1 if($tag->[2] < 1);
}

my %seqs;

tie %seqs, 'DB_File', $dbf, O_RDONLY or die "Cannot open $dbf: $!";

if($dbf=~/^(.+?)\.dbm$/){
	$faf = $1;
}

my $fa_file = undef;

foreach my $tag (@tags){
	if(exists $seqs{$tag->[0]}){
		#my $seq = $seqs{$tag->[0]};
		my $seq = &read_fasta($tag->[0], $tag->[4]);
		if($tag->[4]){
			$tag->[3] = length($seq) if($tag->[3] <= $tag->[2]);
			if($tag->[4]){
				print ">", join("_", $tag->[0], $tag->[1] == 1? "F":"R", $tag->[2], $tag->[3] + 1 - $tag->[2]), "\n";
			} else {
				print ">$tag->[0]\n";
			}
			if($tag->[2] < $tag->[3]){
				my $ss = substr($seq, $tag->[2] - 1, $tag->[3] - $tag->[2] + 1);
				if($tag->[1] == 2){
					$ss =~tr/ACGTacgt/TGCAtgca/;
					$ss = reverse $ss;
				}
				while($ss=~/(.{1,100})/g){ print "$1\n"; }
			}
		} else {
			print $seq;
		}
	} else {
		warn("'$tag->[0]' was not found\n");
	}
}

if($list_all){
	&list_all_seqs;
}

untie %seqs;

if($fa_file){
	close $fa_file;
}

1;

sub read_fasta {
	my $tag = shift;
	my $tidy = shift || 0;
	my $obj = $seqs{$tag};
	if($obj!~/^[0-9]/){
		return $obj;
	}
	my $off = $obj;
	if(not defined $fa_file){
		if(not defined $faf){
			die("Cannot find fasta file");
		} else {
			open($fa_file, "<", $faf) or die $!;
		}
	}
	seek($fa_file, $off, 0);
	my $nam = '';
	my $seq = '';
	while(<$fa_file>){
		if(/^>(\S+)/){
			last if(length $nam);
			$nam = $1;
			if(!$tidy){
				$seq .= $_;
			}
		} else {
			if($tidy){
				chomp;
			}
			$seq .= $_;
		}
	}
	if($nam ne $tag){
		die("Broken dbm index, \"$nam\" ne \"$tag\", offset = $off");
	}
	return $seq;
}

sub list_all_seqs {
	while(my ($tag, $seq) = each %seqs){
		print "$tag\t$seq\n";
	}
}

