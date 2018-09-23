#!/usr/bin/perl -w
my $count=0;
while(<>){
	if($count%4==0){
		print ">", substr($_, 1);
	} elsif($count%4==1){
		print $_;
	}
	$count++;
}

1;
