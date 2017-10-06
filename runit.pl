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
my $RUNITALL = $ENV{RUNIT_ALL} || 0;

my $cmd = join(" ", @ARGV);
$cmd=~s/^\s+//;
$cmd=~s/\s+$//;
&usage unless($cmd);

my $pagesize = POSIX::sysconf(POSIX::_SC_PAGESIZE);
my $sleep_inv = 1000; # mircoseconds

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
	for my $proc (&get_user_all_process()){
		$exclusive{$proc->[0]} = 1;
	}
}
my %procs = ();

&get_linux_sys_info();

print STDERR " -- Total memory $maxram kB\n";
print STDERR " -- $maxcpu cores\n";
print STDERR " -- CMD: $cmd\n";
print STDERR " --\n";

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
	foreach my $proc ($RUNITALL? &get_user_all_process() : &get_child_process($PID)){
		my ($pid, $cm) = @{$proc};
		next if($exclusive{$pid});
		my ($fail, $ut, $st, $rss, $vsz) = &get_linux_proc_info($pid);
		next if($fail);
		$procs{$pid} = [$ut, $st, $cm];
		$mrss += $rss;
		$mvsz += $vsz;
	}
	#print STDERR "$mrss\t$mvsz\n";
	$maxrss = $mrss if($mrss > $maxrss);
	$maxvsz = $mvsz if($mvsz > $maxvsz);
	my $res = waitpid($PID, WNOHANG);
	if($res == -1){
		print STDERR "** Error ", ($? >> 8), "\n";
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
	print STDERR " -- STAT PID($pid): $procs{$pid}[0]\t$procs{$pid}[1]\t$procs{$pid}[2]\n";
	$utime += $procs{$pid}[0];
	$stime += $procs{$pid}[1];
}

printf STDERR " -- retval    %16d\n", $retval;
printf STDERR " -- real      %16.3f\n", $rtime;
printf STDERR " -- user      %16.3f\n", $utime;
printf STDERR " -- sys       %16.3f\n", $stime;
printf STDERR " -- maxrss    %16.3f kB\n", $maxrss;
printf STDERR " -- maxvsz    %16.3f kB\n", $maxvsz;

1;

sub usage {
	print qq{$0 \$'commands'\n};
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

sub get_child_process {
	my @pids = ();
	if(open(IN, "ps -o pid,cmd --no-headers --ppid $PID 2>/dev/null |")){
		while(<IN>){
			my @ts = split;
			push(@pids, \@ts);
		}
		close IN;
	}
	push(@pids, [$PID, $cmd]);
	return @pids;
}

sub get_user_all_process {
	my @pids = ();
	if(open(IN, "ps -o pid,cmd --no-headers --user $USER 2>/dev/null |")){
		while(<IN>){
			my @ts = split;
			push(@pids, \@ts);
		}
		close IN;
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
	return ($fail, $ut, $st, $rss, $vsz);
}

