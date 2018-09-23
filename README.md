## <a name="start"></a>Getting Started
```sh
git clone https://github.com/ruanjue/wtdbg2
cd wtdbg2 && make
# assemble PacBio reads
./wtdbg2 -t 16 -i pacbio.fa.gz -fo prefix -L 5000
# assemble Nanopore reads
./wtdbg2 -t 16 -i ont.fa.gz -fo prefix
# derive consensus
./wtpoa-cns -t 16 -i prefix.ctg.lay > prefix.ctg.lay.fa
```

## <a name="intro"></a>Introduction

Wtdbg2 is a *de novo* sequence assembler for long noisy reads produced by
PacBio or Oxford Nanopore Technologies (ONT). It assembles raw reads without
error correction and then builds the consensus from intermediate assembly
output. Wtdbg2 is able to assemble the human and even the 32Gb
[Axolotl][Axolotl] genome at a speed tens of times faster than [CANU][canu] and
[FALCON][falcon] while producing contigs of comparable base accuracy. Using 96
CPU cores, it can assemble a mammlian genome in several wall-clock hours.

## <a name="install"></a>Installation

Wtdbg2 only works on 64-bit Linux. It depends on [zlib][zlib] for direct
access to gzip'd files. To compile, please type `make` in the source code
directory. You can copy `wtdbg2` and `wtpoa-cns` to your `PATH`.

Wtdbg2 also comes with an approxmimate read mapper `kbm`, a faster but less
accurate consesus tool `wtdbg-cns` and many auxiliary scripts in the `scripts`
directory.

## <a name="use"></a>Usage

Wtdbg2 has two key components: an assembler **wtdg2** and a consenser
**wtpoa-cns**. Executable **wtdbg2** assembles raw reads and generates the
layout of contigs and sequences on the edges/nodes in a file
"*prefix*.ctg.lay".  Executable **wtpoa-cns** takes this file as input and
produces the final consensus in FASTA. A typical workflow looks like this:
```sh
./wtdbg2 -t 16 -i reads.fa.gz -fo prefix
./wtpoa-cns -t 16 -i prefix.ctg.lay > prefix.ctg.lay.fa
```
where `-t` specifies the number of CPU cores. When the default doesn't work
well, you may need to apply more options.

Wtdbg2 combines normal k-mers and homopolymer-compressed (HPC) k-mers to find
read overlaps. Option `-k` specifies the length of normal k-mers, while `-p`
specifies the length of HPC k-mers. By default, wtdbg2 chooses one every four
consecutive k-mers. For data of relatively low coverage, you may increase this
sampling rate by reducing `-S`. This will greatly increase the peak memory as a
cost. Reducing `-e` to 2 also helps lower coverage. For PacBio data, option
`-L5000` is usually recommended. Please run `wtdbg2 --help` for a complete list
of available options.

The following table shows various command lines we used and their resource
usage for the assembly step (not including the consensus step):

|Dataset                 |Genome|Coverage  |Wtdbg2 options|CPU hours|Peak RAM|
|:-----------------------|-----:|---------:|:-------------|--------:|-------:|
|[E. coli][pbcr]         |4.6Mb |PacBio x20|              |         |        |
|[C. elegans][ce]        |100Mb |PacBio x80|              |         |        |
|[Human CHM1][chm1]      |3Gb   |PacBio x60|-L5000 -e4    |378.5    |  252.7G|
|[Human NA12878][na12878]|3Gb   |ONT x30   |-S2 -e2       |197.4    |  244.9G|
|[Axolotl][axosra]       |32Gb  |PacBio x32|-L5000 -AS2   |3189.7   | 1593.6G|

## Getting Help

Please use the [GitHub's Issues page][issue] if you have questions. You may
also directly contact Jue Ruan at ruanjue@gmail.com.

[miniasm]: https://github.com/lh3/miniasm
[canu]: https://github.com/marbl/canu
[falcon]: https://github.com/PacificBiosciences/FALCON
[Axolotl]: https://www.nature.com/articles/nature25458
[zlib]: http://zlib.net
[chm1]: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA246220
[na12878]: https://github.com/nanopore-wgs-consortium/NA12878
[pbcr]: http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData.tar.gz
[ce]: https://github.com/PacificBiosciences/DevNet/wiki/C.-elegans-data-set
[axosra]: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA378970
[issue]: https://github.com/ruanjue/wtdbg2/issues
