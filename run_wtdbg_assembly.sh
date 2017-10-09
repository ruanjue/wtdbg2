#!/bin/bash
#
# Author: Jue Ruan <ruanjue@gmail.com>
#


print_usage(){
cat <<EOF
Usage: ${0} [options]
Options:
 -i <string> Reads file, [wt.fa or wt.fa.gz]
 -o <string> Prefix of output files, [dbg]
 -t <int>    Number of threads, [96]
 -S <int>    See wtdbg -h, [4]
 -k <int>    See wtdbg -h, [0]
 -p <int>    See wtdbg -h, [21]
 -T          Print commands and exit
 -f          force overwrite
 -P <string> Path of wtdbg and kbm
EOF
}

OPTIND=1

PVER=1.2.8
WT_S=4
WT_k=0
WT_P=21
NCPU=96
if [ -e wt.fa ]; then
	READS=wt.fa
elif [ -e wt.fa.gz ]; then
	READS=wt.fa.gz
else
	READS=seqs.fa
fi
PREFIX=dbg
OVERWRITE=0
STEP=0
EXEC_CMD=1

while getopts "h?fTi:o:t:S:k:p:P:" opt; do
	case "$opt" in
	h|\?)
		print_usage
		exit 0
		;;
	i)
		READS=$OPTARG
		;;
	o)
		PREFIX=$OPTARG
		;;
	t)
		NCPU=$OPTARG
		;;
	f)
		OVERWRITE=1
		;;
	S)
		WT_S=$OPTARG
		;;
	k)
		WT_k=$OPTARG
		;;
	p)
		WT_P=$OPTARG
		;;
	T)
		EXEC_CMD=0
		;;
	P)
		export PATH="$OPTARG:$PATH"
		;;
	*)
		echo "Unknown paramemter"
		print_usage
		exit 0
		;;
	esac
done

echo "#!/bin/bash"
echo "#WTDBG, a de novo assembler for long noisy reads implemented fuzzy bruijn graph"
echo "#It is ultra-efficient in both CPU time and peak memory, and able to assemble huge genomes"
echo "#Author: Jue Ruan <ruanjue@gmail.com>"
echo "#NOTICE: If you are familiar with wtdbg, please add more options to get better results"
echo "#        First round of error correction is good enough for short reads mapping\n"

echo "### checking"
echo -n "#"
which wtdbg-$PVER || exit
echo -n "#"
which wtdbg-cns || exit
echo -n "#"
which kbm-$PVER || exit
echo -n "#"
which best_kbm_hit.pl || exit
echo -n "#"
which awk || exit
echo -n "#"
which tee || exit
echo -n "#"
which date || exit
echo -n "#"
which map2dbgcns || exit
echo -n "#"
ls -l $READS || exit

if [ -e $PREFIX.ctg.lay ]; then
	STEP=1
	if [ -e $PREFIX.ctg.lay.fa ]; then
		STEP=2
		if [ -e $PREFIX.map ]; then
			STEP=3
			if [ -e $PREFIX.map.lay ]; then
				STEP=4
				if [ -e $PREFIX.map.fa ]; then
					STEP=5
				fi
			fi
		fi
	fi
fi

if [[ $STEP > 0 ]] && [[ $OVERWRITE > 0 ]] ; then
	echo "#[WARNNING] Force overwrite!!!"
	STEP=0
fi

if [[ $EXEC_CMD > 0 ]]; then
	echo "Pipeline will start in 5 seconds, 'kill -9 $$' to cancel"
	sleep 5
fi

echo "### assembling"
CMD="wtdbg-$PVER -t $NCPU -i $READS --tidy-reads 5000 -fo $PREFIX -k $WT_k -p $WT_P -S $WT_S --rescue-low-cov-edges"
if [[ $STEP < 1 ]]; then
echo $CMD
if [[ $EXEC_CMD > 0 ]] ; then
	date
	eval $CMD
fi
else
	echo "# SKIP!"
fi

echo "### first round of correction"
CMD="wtdbg-cns -t $NCPU -i $PREFIX.ctg.lay -fo $PREFIX.ctg.lay.fa -c 0"
if [[ $STEP < 2 ]]; then
echo $CMD
if [ $EXEC_CMD -gt 0 ] ; then
	date
	eval $CMD
fi
else
	echo "# SKIP!"
fi

echo "### mapping"
CMD="kbm-$PVER -t $NCPU -d $PREFIX.ctg.lay.fa -i $READS -k $WT_K -p $WT_P -S $WT_S -O 0 | best_kbm_hit.pl | awk '{print \$6\"\t\"\$9\"\t\"\$10\"\t\"\$1\"\t\"\$2\"\t\"\$4\"\t\"\$5}' >$PREFIX.map"
if [[ $STEP < 3 ]]; then
echo $CMD
if [ $EXEC_CMD -gt 0 ] ; then
	date
	eval $CMD
fi
else
	echo "# SKIP!"
fi

echo "### generating new layout"
CMD="map2dbgcns $PREFIX.ctg.lay.fa $READS $PREFIX.map >$PREFIX.map.lay"
if [[ $STEP < 4 ]]; then
echo $CMD
if [ $EXEC_CMD -gt 0 ] ; then
	date
	eval $CMD
fi
else
	echo "# SKIP!"
fi

echo "### second round of correction"
CMD="wtdbg-cns -t $NCPU -i $PREFIX.map.lay -fo $PREFIX.map.fa -k 13 -c 3"
if [[ $STEP < 5 ]]; then
echo $CMD
if [ $EXEC_CMD -gt 0 ] ; then
	date
	eval $CMD
fi
else
	echo "# SKIP!"
fi

echo "### Finished"
if [ $EXEC_CMD -gt 0 ] ; then
	date
fi
