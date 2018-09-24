#!/bin/bash
#
# Author: Jue Ruan
#
REF=$1
CTG=$2
# change this path to your mummer
MUM=/public/software/mummer-323-intel

if [ -z $REF ] || [ -z $CTG ]; then
	echo "Usage: $0 <ref> <ctg>"
	exit
fi

echo "REF:$REF"
echo "CTG:$REF"

$MUM/nucmer --mumreference -l 100 -c 1000 -d 10 --banded -D 5 $REF $CTG

$MUM/delta-filter -i 95 -o 95 out.delta > out.best.delta

$MUM/dnadiff -d out.best.delta

$MUM/mummerplot out.best.delta --fat -f -png

