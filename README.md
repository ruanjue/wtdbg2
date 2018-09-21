# NEWS
* |2018-09-21| Rename wtdbg-1.2.8 to wtdbg2
* |2018-09-19| GFA supports
```sh
wtdbg-dot2gfa.pl dbg.3.dot >dbg.3.gfa
```
* |2018-09-01| New consensus module: POACNS <br>
```sh
wtpoa-cns -t 0 -i dbg.ctg.lay -fo dbg.ctg.poacns.fa
```
1, wtpoa-cns implements POA to generate MSA from reads fragments within an edge. Used SSE instructions, but seems still can be further improved in speed <br>
2, performs realignment on MSA. <br>
3, recalibrates homopolymers, and make consensus sequence for an edge. Most reamining small errors come from homopolymer, more efforts are needed <br>
4, joins edges' sequences into a contig. <br>

Welcome to test POACNS and feedback.

# WTDBG

A fuzzy Bruijn graph (FBG) approach to long noisy reads assembly

# Introduction
A challenge in assembling long noisy reads from third generation sequencing (TGS) is reducing its requirement of computing resource, especially for large genomes.
To address this issue, I developed a novel sequence alignment algorithm and a new assembly graph for efficiently assembling large genomes using TGS data.

* Alignment <br>
KBM: Kmer-BIN-Mapping.<br>
KBM groups k-mers from each non-overlapped sliding 256 bp fragments in long reads into bins.
Bins of which most k-mers are high frequency, are filtered as highly repetitive ones.
Then, KBM searches synteny of matched bin pairs in sequences in a dynamic programming way.
A matched bin pair in two sequences is defined as two bins different by original but share a set of k-mers.
The result of alignments in KBM have the same features of traditional sequence alignment, excepting the unit of KBM alignments is 256 bp bin instead of single base.

* Assembly <br>
FBG: Fuzzy Bruijn Graph. <br>
FBG is composed of vertices in length of 1024 bp from reads, and edges connecting vertices in their order on read paths.
Comparing with DBG, the size of vertices in FBG are much bigger, thus won¡¯t be sensitive to small repeat.
To tolerate high sequencing errors, FBG's vertices are found using gapped sequence alignments from KBM or other aligners, comparing with searching identical k-mers in DBG.

* Why choose wtdbg in genome assembly <br>
There are many assemblers for long noisy reads assembly, e.g. FALCON, CANU, miniasm, and SMARTdenovo (progenitor of wtdbg). If you have a genome of 10G bp or bigger in size,
wtdbg is your first or even the only option. For small but complicated genomes (< 3 G), wtdbg was often reported to yield better assembly by my friends.
Besides, KBM is easy to use when you are setting up a web-server for long reads mapping (see Example 2).

* Limitation <br>
Max read length is 0x0003FFFFU (256 Kb), longer reads will be split. <br>
Max number of reads is 0x03FFFFFFU (64 M). If your data volume exceeds, please filter relative shorter reads. <br>
In KBM, max read length is 0xFFFFFFFFU (4 Gb), max number of reads is 0x0FFFFFFFU (256 M). <br>
Max number of threads is 4096. <br>
Cannot parallelly run in multiple nodes. But you can implement it simplely using `kbm` and `wtdbg --load-alignments` <br>
Developed and tested in Linux-GCC only. <br>
Only accepts fasta/fastq format for input, '.gz' suffixed files will be piped by `gzip -dc`.

# Installation
```sh
git clone https://github.com/ruanjue/wtdbg2.git
cd wtdbg2
make
```

# Long reads mapping

Supposes you have `hg19.fa` as reference sequences, and `reads.fq.gz` as query sequences.
## Example 1
```sh
kbm -t 64 -d hg19.fa -i reads.fa.gz -o reads.kbmap

```
### output format
* COL1  `qry_name`
* COL2  `qry_strand`
* COL3  `qry_length`
* COL4  `qry_beg`
* COL5  `qry_end`
* COL6  `ref_name`
* COL7  `ref_strand` (always equals `+`)
* COL8  `ref_length`
* COL9  `ref_beg`
* COL10 `ref_end`
* COL11 `match_len` (length of matched k-mers)
* COL12 `align_len` (length of aligned)
* COL13 `#kcnt`     (number of matched k-mers)
* COL14 `#gap`      (number of gapped BINs)
* COL15 `cigar`     (256 x SAM's cigar)

## Example 2
Suitable for online tools, or frequently used references
### Build KBM-INDEX
```sh
kbm -t 64 -d hg19.fa -i /dev/null -W hg19.kbmidx
```
### Set up KBM server
```sh
kbm -R hg19.kbmidx start
```
Now, hg19.kbmidx is cached in memory for further call.

### Mapping with KBM-INDEX
```sh
kbm -R hg19.kbmidx -t 64 -i reads.fa -o reads.kbmap
```
Please note that, hg19.kbmidx can be multilple called by any processes in the same computer. <br>

### Shutdown KBM server
```sh
kbm -R hg19.kbmidx stop
```

# Long reads assembling

## Quick Start
```sh
echo "manual"
run_wtdbg_assembly.sh -h
echo "generating shell script"
run_wtdbg_assembly.sh -t 0 -i reads.fa.gz -o dbg -T >run.sh
```
`run.sh` is ready for invoked

## Play with wtdbg
```sh
wtdbg2 -h
wtdbg2 --help
```
### options
```sh
 -t <int>    Number of threads, 0: all cores, [0]
 -i <string> Long reads sequences file, + *
 -I <string> Error-free sequences file, +
 -o <string> Prefix of output files, *
 -f          Force overwrite
 -k <int>    Kmer fsize, 0 <= k <= 25, [0]
 -p <int>    Kmer psize, 0 <= p <= 25, [21]
             k + p <= 25, seed is <k-mer>+<p-homopolymer-compressed>
 -K <float>  Filter high frequency kmers, maybe repetitive, [1000]
             if K >= 1, take the integer value as cutoff
             else, mask the top fraction part high frequency kmers
 -E <int>    Min kmer frequency, [2]
 -F          Filter low frequency kmers by a 4G-bytes array (max_occ=3 2-bits). Here, -E must greater than 1
 -S <int>    Subsampling kmers, 1/(<-S>) kmers are indexed, [4]
             -S is very useful in saving memeory and speeding up
             please note that subsampling kmers will have less matched length
 -X <int>    Max number of bin(256bp) in one gap, [4]
 -Y <int>    Max number of bin(256bp) in one deviation, [4]
 -x <int>    penalty for BIN gap, [-7]
 -y <int>    penalty for BIN deviation, [-21]
 -l <float>  Min length of alignment, [2048]
 -m <float>  Min matched, [200]
 -s <float>  Max length variation of two aligned fragments, [0.2]
 -q          Quiet
 -v          Verbose, +
 --help      Show more options
```

### For higher error rate long sequences
Decrease `-p`. Try `-p 19` or `-p 17` <br>
Decrease `-S`. Try `-S 2` or `-S 1` <br>
Both will increase computing time.

### For very high coverage
Increase `--edge-min`. Try `--edge-min 4`, or higher.

### For low coverage
Decrease `--edge-min`. Try `--edge-min 2 --rescue-low-cov-edges`.

### Filter reads
`--tidy-reads 5000`. Will filtered shorter sequences. If names in format of `\/\d+_\d+$`, will selected the longest subread.

### output
Suppose the prefix is `dbg`
* dbg.1.dot <br>
DOT file for initialized graph
* dbg.1.nodes <br>
nodes and their positions in reads
* dbg.1.reads <br>
reads and their nodes
* dbg.2.dot <br>
DOT file after transitive reduction
* dbg.3.dot <br>
DOT file after merging bubble and remove tips
* dbg.alignments <br>
KBMAP file, all vs all alignments
* dbg.binkmer <br>
Distribution of number of k-mers in a BIN
* dbg.closed\_bins <br>
Filtered BINs
* dbg.clps
Reads clip information. <br>
COL1 read\_name <br>
COL2 read\_length <br>
COL3 keep\_offset <br>
COL4 keep\_length
* dbg.ctg.dot <br>
DOT file for contigs
* dbg.ctg.lay <br>
Contigs layout file. Will be read by `wtdbg-cns`. This file is the main result of `wtdbg` <br>
**Format**:<br>
```
>ctg(\d+) nodes=(\d+) len=(\d+)
E <OFFSET> <NODE1> <STRAND1> <NODE2> <STRAND2>
S <READ_NAME> <STRAND> <REG_OFFSET> <REG_LENGTH> <REG_SEQ>
S <READ_NAME> <STRAND> <REG_OFFSET> <REG_LENGTH> <REG_SEQ>
S <READ_NAME> <STRAND> <REG_OFFSET> <REG_LENGTH> <REG_SEQ>
...
E ...
...

```
One contig contains many edges (starting with 'E'), each edge contains many regions inside reads. <br>
Please note that one read often contains many REGs.

* dbg.events <br>
Log file of graph simplification
* dbg.frg.dot <br>
DOT file for unitigs
* dbg.frg.nodes <br>
unitigs and their nodes
* dbg.kmerdep <br>
Distribution of k-mer depth
* STDERR stream <br>
wtdbg print runtime information on progrom's STDERR stream. `--quiet` to disiable it

## Consensus
```sh
wtdbg-cns -t 64 -i dbg.ctg.lay -o dbg.ctg.lay.fa
```
The output file `dbg.ctg.lay.fa` is ready for further polished by `PILON` or `QUIVER`.

```sh
wtpoa-cns -t 64 -i dbg.ctg.lay -o dbg.ctg.lay.fa
```
wtpoa-cns is slower than wtdbg-cns, but offer more accurate consensus sequences.
I will update it in following development.

# Performance
## Human (3G) CHM1 PacBio P5C3 dataset, 65.5 core.hours
* Data Source
http://datasets.pacb.com/2014/Human54x/fasta.html
* Command
```sh
wtdbg2 -t 96 -i pb.fa -fo dbg --tidy-reads 5000 --edge-min 2 --rescue-low-cov-edges
```
* Contigs
`TOT 2978536704, CNT 8752, AVG 340327, MAX 11662848, N50 1925120, L50 453, N90 400128, L90 1727, Min 5120`
* Runtime
`real 6131.803 sec, user 201836.200 sec, sys 33956.790 sec, maxrss 117281672.0 kB, maxvsize 202422172.0 kB`

## Human (3G) CHM1 PacBio P6C4 dataset, 211.3 core.hours
* Data Source
http://www.ebi.ac.uk/ena/data/view/PRJNA246220
* Command
```sh
wtdbg2 -t 96 -i wt.fa -fo dbg --tidy-reads 5000 --edge-min 4 --rescue-low-cov-edges
```
* Contigs
`TOT 2964872448, CNT 1909, AVG 1553103, MAX 105310208, N50 23586816, L50 34, N90 3326976, L90 158, Min 5120`
* Runtime
`real 16806.534 sec, user 681278.770 sec, sys 79371.630 sec, maxrss 264956752.0 kB, maxvsize 443356532.0 kB`

## Axolotl (32G) PacBio dataset, 32 X, 3053 core.hours
* Command
```sh
wtdbg2 -t 96 -i ../rawdata/pacbio.fa.gz -p 21 -S 2 --aln-noskip --rescue-low-cov-edges --tidy-reads 5000 -fo axolotl
```
* Contigs
`TOT 27375160576, CNT 115355, AVG 237313, MAX 7812608, N50 606976, L50 12527, N90 144896, L90 47295, Min 5120`
* Runtime
`real 190237.591 sec, user 10994200.800 sec, sys 488715.030 sec, maxrss 1671005352.0 kB, maxvsize 2365400208.0 kB`

## Human (3G) NA12878 ONT dataset, 197.5 core.hours
* Data Source
https://github.com/nanopore-wgs-consortium/NA12878
* Command
```sh
wtdbg2 -t 64 -i NA12878-ONT.fa.gz -fo dbg -S 2 --edge-min 2 --rescue-low-cov-edges
```
* Contigs
`TOT 2827644928, CNT 19473, AVG 145209, MAX 31366400, N50 4540672, L50 162, N90 172800, L90 1111, Min 5120`

* Runtime
`real 14992.925 sec, user 649202.270 sec, sys 61638.300 sec, maxrss 256840096.0 kB, maxvsize 356668088.0 kB`

# Citation
To be published. <br>
URL **https://github.com/ruanjue/wtdbg2/** <br>

# Contact
Jue Ruan <ruanjue@gmail.com> <br>
Jue Ruan <ruanjue@caas.cn>
