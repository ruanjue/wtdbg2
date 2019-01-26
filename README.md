## <a name="start"></a>Getting Started
```sh
git clone https://github.com/ruanjue/wtdbg2
cd wtdbg2 && make
# assemble long reads
./wtdbg2 -x rs -g 4.6m -i reads.fa.gz -t16 -fo prefix
# derive consensus
./wtpoa-cns -t16 -i prefix.ctg.lay.gz -fo prefix.ctg.fa

# polish consensus, not necessary if you want to polish the assemblies using other tools
minimap2 -t 16 -x map-pb -a prefix.ctg.fa reads.fa.gz | samtools view -Sb - >prefix.ctg.map.bam
samtools sort prefix.ctg.map.bam prefix.ctg.map.srt
samtools view prefix.ctg.map.srt.bam | ./wtpoa-cns -t 16 -d prefix.ctg.fa -i - -fo prefix.ctg.2nd.fa

# polish contigs using short reads
bwa mem -t 16 prefix.ctg.fa sr.1.fa sr.2.fa | samtools view -Sb - >sr.bam
samtools sort sr.bam sr.srt
samtools view sr.srt.bam | ./wtpoa-cns -t 16 -x sam-sr -d prefix.ctg.fa -i - -fo prefix.ctg.3rd.fa
```

## <a name="intro"></a>Introduction

Wtdbg2 is a *de novo* sequence assembler for long noisy reads produced by
PacBio or Oxford Nanopore Technologies (ONT). It assembles raw reads without
error correction and then builds the consensus from intermediate assembly
output. Wtdbg2 is able to assemble the human and even the 32Gb
[Axolotl][Axolotl] genome at a speed tens of times faster than [CANU][canu] and
[FALCON][falcon] while producing contigs of comparable base accuracy.

During assembly, wtdbg2 chops reads into 1024bp segments, merges similar
segments into a vertex and connects vertices based on the segment adjacency on
reads. The resulting graph is called fuzzy Bruijn graph (FBG). It is akin to De
Bruijn graph but permits mismatches/gaps and keeps read paths when collapsing
k-mers. The use of FBG distinguishes wtdbg2 from the majority of long-read
assemblers.

## <a name="install"></a>Installation

Wtdbg2 only works on 64-bit Linux. To compile, please type `make` in the source
code directory. You can then copy `wtdbg2` and `wtpoa-cns` to your `PATH`.

Wtdbg2 also comes with an approxmimate read mapper `kbm`, a faster but less
accurate consesus tool `wtdbg-cns` and many auxiliary scripts in the `scripts`
directory.

## <a name="use"></a>Usage

Wtdbg2 has two key components: an assembler **wtdbg2** and a consenser
**wtpoa-cns**. Executable **wtdbg2** assembles raw reads and generates the
contig layout and edge sequences in a file "*prefix*.ctg.lay.gz". Executable
**wtpoa-cns** takes this file as input and produces the final consensus in
FASTA. A typical workflow looks like this:
```sh
./wtdbg2 -x rs -g 4.6m -t 16 -i reads.fa.gz -fo prefix
./wtpoa-cns -t 16 -i prefix.ctg.lay.gz -fo prefix.ctg.fa
```
where `-g` is the estimated genome size and `-x` specifies the sequencing
technology, which could take value "rs" for PacBio RSII, "sq" for PacBio
Sequel, "ccs" for PacBio CCS reads and "ont" for Oxford Nanopore. This option
sets multiple parameters and should be **applied before other parameters**.
When you are unable to get a good assembly, you may need to tune other
parameters as follows.

Wtdbg2 combines normal k-mers and homopolymer-compressed (HPC) k-mers to find
read overlaps. Option `-k` specifies the length of normal k-mers, while `-p`
specifies the length of HPC k-mers. By default, wtdbg2 samples a fourth of all
k-mers by their hashcodes. For data of relatively low coverage, you may
increase this sampling rate by reducing `-S`. This will greatly increase the
peak memory as a cost, though. Option `-e`, which defaults to 3, specifies the
minimum read coverage of an edge in the assembly graph. You may adjust this
option according to the overall sequencing depth, too. Option `-A` also helps
relatively low coverage data at the cost of performance. For PacBio data,
`-L5000` often leads to better assemblies emperically, so is recommended.
Please run `wtdbg2 --help` for a complete list of available options or consult
[README-ori.md](README-ori.md) for more help.

The following table shows various command lines and their resource usage for
the assembly step:

|Dataset                 |GSize |Cov     |Asm options        |CPU asm |CPU cns |Real tot|     RAM|
|:-----------------------|-----:|-------:|:------------------|-------:|-------:|-------:|-------:|
|[E. coli][pbcr]         |4.6Mb |PB x20  |-x rs -g4.6m -t16  |     53s|   8m54s|     42s|    1.0G|
|[C. elegans][ce]        |100Mb |PB x80  |-x rs -g100m -t32  |   1h07m|   5h06m|  13m42s|   11.6G|
|[D. melanogaster A4][dm2]| 144m|PB x120 |-x rs -g144m -t32  |   2h06m|   5h11m|  26m17s|   19.4G|
|[D. melanogaster ISO1][dm1]|144m|ONT x32|-xont -g144m -t32  |   5h12m|   4h30m|  25m59s|   17.3G|
|[A. thaliana][at]       |125Mb |PB x75  |-x sq -g125m -t32  |  11h26m|   4h57m|  49m35s|   25.7G|
|[Human NA12878][na12878]|3Gb   |ONT x36 |-x ont -g3g -t31   | 793h11m|  97h46m|  31h03m|  221.8G|
|[Human NA19240][na19240]|3Gb   |ONT x35 |-x ont -g3g -t31   | 935h31m|  89h17m|  35h20m|  215.0G|
|[Human HG00733][hg00733]|3Gb   |PB x93  |-x sq -g3g -t47    |2114h26m| 152h24m|  52h22m|  338.1G|
|[Human NA24385][na24385]|3Gb   |CCS x28 |-x ccs -g3g -t31   | 231h25m|  58h48m|  10h14m|  112.9G|
|[Human CHM1][chm1]      |3Gb   |PB x60  |-x rs -g3g -t96    | 105h33m| 139h24m|   5h17m|  225.1G|
|[Axolotl][axosra]       |32Gb  |PB x32  |-x rs -g32g -t96   |2806h40m|1456h13m| 110h16m| 1788.1G|

The timing was obtained on three local servers with different hardware
configurations. There are also run-to-run fluctuations. Exact timing on your
machines may differ. The assembled contigs can be found at the following FTP:
```txt
ftp://ftp.dfci.harvard.edu/pub/hli/wtdbg/
```

## Limitations

* For Nanopore data, wtdbg2 may produce an assembly smaller than the true
  genome.

* When inputing multiple files of both fasta and fastq format, please put fastq first, then fasta.
  Otherwise, program cannot find '>' in fastq, and append all fastq in one read.

## Citing wtdbg2

If you use wtdbg2, please cite:

> Ruan, J. and Li, H. (2019) Fast and accurate long-read assembly with wtdbg2. *bioRxiv*. doi:10.1101/530972


## Getting Help

Please use the [GitHub's Issues page][issue] if you have questions. You may
also directly contact Jue Ruan at ruanjue@gmail.com.

[miniasm]: https://github.com/lh3/miniasm
[canu]: https://github.com/marbl/canu
[falcon]: https://github.com/PacificBiosciences/FALCON
[Axolotl]: https://www.nature.com/articles/nature25458
[chm1]: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP044331
[na12878]: https://github.com/nanopore-wgs-consortium/NA12878/blob/master/rel5.md
[na19240]: https://www.ebi.ac.uk/ena/data/view/PRJEB26791
[pbcr]: http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData.tar.gz
[ce]: https://github.com/PacificBiosciences/DevNet/wiki/C.-elegans-data-set
[axosra]: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA378970
[issue]: https://github.com/ruanjue/wtdbg2/issues
[at]: https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/
[dm1]: https://www.ebi.ac.uk/ena/data/view/SRR6702603
[dm2]: https://www.ebi.ac.uk/ena/data/view/SRR5439404
[hg00733]: https://www.ebi.ac.uk/ena/data/view/SRR7615963
[na24385]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/
