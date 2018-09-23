#!/usr/bin/perl -w
#
# Author: Jue Ruan <ruanjue@gmail.com>
#
use strict;
use warnings;
use POSIX;
use POSIX ":sys_wait_h";
use Time::HiRes;

my $USER = $ENV{USER};
if(not defined $USER){
	$USER = `whoami`;
	chomp $USER;
}
my $RUNITALL = $ENV{RUNIT_ALL} || 0;

my $cmd = join(" ", @ARGV);
$cmd=~s/^\s+//;
$cmd=~s/\s+$//;
&usage unless($cmd);

my $pagesize = POSIX::sysconf(POSIX::_SC_PAGESIZE);
my $sleep_inv = 1000000; # in mircoseconds

my $maxram = 0; # kb
my $maxcpu = 0;
my $retval = 0;
my $utime = 0;
my $stime = 0;
my $rtime = 0;
my $maxrss = 0;
my $maxvsz = 0;

my %exclusive = ();
if($RUNITALL){
	for my $proc (&get_all_process()){
		$exclusive{$proc->[0]} = 1;
	}
}
my %procs = ();

&get_linux_sys_info();

print STDERR " --------------------------------------------------------------------------------\n";
print STDERR " -- runit.pl is a program launcher and minitor written by Jue Ruan <ruanjue\@gmail.com>\n";
print STDERR " -- RAM : $maxram kB\n";
print STDERR " -- CPU : $maxcpu cores\n";
print STDERR " -- SYS : ", `uname -a`;
print STDERR " -- USER: $USER\n";
print STDERR " -- DATE: ", `date`;
print STDERR " -- PWD : ", `pwd`;
print STDERR " -- CMD : $cmd\n";
print STDERR " --------------------------------------------------------------------------------\n";

my $PID = open(READ, "$cmd |") or die("Cannot invoke commands: $!");

my $begtime = Time::HiRes::time;

if(fork() == 0){
	while(<READ>){
		print;
	}
	exit;
}

while(1){
	my $mrss = 0;
	my $mvsz = 0;
	foreach my $proc (&get_all_process()){
		my ($pid, $cm, $ppid) = @{$proc};
		next if($exclusive{$pid});
		my ($fail, $ut, $st, $rss, $vsz) = &get_linux_proc_info($pid);
		next if($fail);
		if(exists $procs{$pid}){
			$procs{$pid}[0] = $ut;
			$procs{$pid}[1] = $st;
		} else {
			print STDERR " -- RUNIT PID($pid): $ppid\t$cm\n";
			$procs{$pid} = [$ut, $st, $cm, $ppid];
		}
		$mrss += $rss;
		$mvsz += $vsz;
	}
	#print STDERR "$mrss\t$mvsz\n";
	$maxrss = $mrss if($mrss > $maxrss);
	$maxvsz = $mvsz if($mvsz > $maxvsz);
	my $res = waitpid($PID, WNOHANG);
	if($res == -1){
		print STDERR " ** Error ", ($? >> 8), "\n";
		last;
	} elsif($res){
		$retval = $? >> 8;
		last;
	} else {
		Time::HiRes::usleep($sleep_inv);
	}
}
my $endtime = Time::HiRes::time;
$rtime = $endtime - $begtime;

foreach my $pid (sort {$a <=> $b} keys %procs){
	next if($exclusive{$pid});
	print STDERR " -- STAT PID($pid): $procs{$pid}[3]\t$procs{$pid}[0]\t$procs{$pid}[1]\t$procs{$pid}[2]\n";
	$utime += $procs{$pid}[0];
	$stime += $procs{$pid}[1];
}

print STDERR " --------------------------------------------------------------------------------\n";
printf STDERR " -- retval    %16d\n", $retval;
printf STDERR " -- real      %16.3f\n", $rtime;
printf STDERR " -- user      %16.3f\n", $utime;
printf STDERR " -- sys       %16.3f\n", $stime;
printf STDERR " -- maxrss    %16.3f kB\n", $maxrss;
printf STDERR " -- maxvsz    %16.3f kB\n", $maxvsz;
print STDERR " --------------------------------------------------------------------------------\n";

1;

sub usage {
	print qq{Launch program and minitor the cputime and ram usage of it and its childs\n};
	print qq{Usage: $0 \$'commands'\n};
	exit 1;
}

sub get_linux_sys_info {
	open(IN, "/proc/meminfo") or die;
	while(<IN>){
		my @ts = split;
		if($ts[0]=~/^MemTotal/){
			$maxram = $ts[1];
			last;
		}
	}
	close IN;
	open(IN, "/proc/cpuinfo") or die;
	$maxcpu = 0;
	while(<IN>){
		if(/^processor/){
			$maxcpu ++;
		}
	}
	close IN;
}

sub get_all_process {
	my %ps = ();
	my @pids = ();
	if(my $psid = open(IN, "ps -o ppid,pid,cmd --no-headers --user $USER 2>/dev/null |")){
		while(<IN>){
			chomp;
			my @ts = split /\s+/, $_, 3;
			next if($ts[1] == $psid);
			if($RUNITALL){
				push(@pids, [$ts[1], $ts[2], $ts[0]]);
			} else {
				$ps{$ts[1]}[0] = $ts[2];
				$ps{$ts[1]}[1] = $ts[0];
				$ps{$ts[1]}[2] = [] if(not defined $ps{$ts[1]}[2]);
				push(@{$ps{$ts[0]}[2]}, $ts[1]);
			}
		}
		close IN;
	}
	if($RUNITALL){
	} else {
		my @stack = ($PID);
		while(@stack){
			my $pid = pop @stack;
			next unless(exists $ps{$pid});
			my $p = $ps{$pid};
			push(@pids, [$pid, $p->[0], $p->[1]]);
			push(@stack, @{$p->[2]});
		}
	}
	return @pids;
}

sub get_linux_proc_info {
	my $pid = shift;
	my ($fail, $ut, $st, $rss, $vsz);
	if(open(IN, "/proc/$pid/stat")){
		if($_ = <IN>){
			my @ts = split /\s/, $_;
			#print STDERR join("|", @ts), "\n";
			#for(my $i=0;$i<@ts;$i++){
			#	print STDERR ($i+1), "\t", $ts[$i], "\n";
			#}
			$ut = $ts[13] / 100;
			$st = $ts[14] / 100;
			$vsz = $ts[22] / 1024;
			$rss = $ts[23] * $pagesize / 1024;
			$fail = 0;
		} else {
			$fail = 1;
		}
		close IN;
	} else {
		$fail = 1;
	}
	if($fail){
		print STDERR " ** FAIL STAT($pid)\n";
	}
	return ($fail, $ut, $st, $rss, $vsz);
}

