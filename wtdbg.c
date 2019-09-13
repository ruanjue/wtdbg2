/*
 *
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "wtdbg.h"
#include "wtdbg-graph.h"
#include <getopt.h>
#include <regex.h>

#ifndef VERSION
#define VERSION 0.0
#endif
#ifndef RELEASE
#define RELEASE 19830203
#endif

static struct option prog_opts[] = {
	{"cpu",                              1, 0, 't'},
	{"input",                            1, 0, 'i'},
	{"err-free-seq",                     1, 0, 'I'},
	{"force",                            0, 0, 'f'},
	{"prefix",                           1, 0, 'o'},
	{"preset",                           1, 0, 'x'},
	{"kmer-fsize",                       1, 0, 'k'},
	{"kmer-psize",                       1, 0, 'p'},
	{"kmer-depth-max",                   1, 0, 'K'},
	{"kmer-depth-min",                   1, 0, 'E'},
	{"genome-size",                      1, 0, 'g'},
	{"rdcov-cutoff",                     1, 0, 'X'},
	{"rdname-filter",                    1, 0, 3007},
	{"rdname-includeonly",               1, 0, 3008},
	{"rdcov-filter",                     1, 0, 2009},
	//{"kmer-depth-min-filter",          0, 0, 'F'},
	{"kmer-subsampling",                 1, 0, 'S'},
	{"kbm-parts",                        1, 0, 1035},
	{"dp-max-gap",                       1, 0, 2005},
	{"dp-max-var",                       1, 0, 2006},
	{"dp-penalty-gap",                   1, 0, 2007},
	{"dp-penalty-var",                   1, 0, 2008},
	{"aln-min-length",                   1, 0, 'l'},
	{"aln-min-match",                    1, 0, 'm'},
	{"aln-min-similarity",               1, 0, 's'},
	{"aln-max-var",                      1, 0, 2004},
	{"realign",                          0, 0, 'R'},
	{"realn-kmer-psize",                 1, 0, 3001},
	{"realn-kmer-subsampling",           1, 0, 3002},
	{"realn-min-length",                 1, 0, 3003},
	{"realn-min-match",                  1, 0, 3004},
	{"realn-min-similarity",             1, 0, 3005},
	{"realn-max-var",                    1, 0, 3006},
	//{"corr-mode",                        1, 0, 2010},
	//{"corr-min",                         1, 0, 2012},
	//{"corr-max",                         1, 0, 2013},
	//{"corr-cov",                         1, 0, 2014},
	//{"corr-block-size",                  1, 0, 2015},
	//{"corr-block-step",                  1, 0, 2016},
	{"keep-multiple-alignment-parts",    1, 0, 2011},
	{"verbose",                          0, 0, 'v'},
	{"quiet",                            0, 0, 'q'},
	{"version",                          0, 0, 'V'},
	{"help",                             0, 0, 1000}, // detailed document
	{"tidy-reads",                       1, 0, 'L'},
	{"tidy-name",                        0, 0, 1001},
	{"err-free-nodes",                   0, 0, 1002},
	{"limit-input",                      1, 0, 1003},
	{"node-len",                         1, 0, 1004},
	{"node-ovl",                         1, 0, 1005},
	{"node-drop",                        1, 0, 1006},
	{"edge-min",                         1, 0, 'e'},
	{"edge-max-span",                    1, 0, 3009},
	{"node-min",                         1, 0, 1007},
	{"node-max",                         1, 0, 1008},
	{"ttr-cutoff-depth",                 1, 0, 1009},
	{"ttr-cutoff-ratio",                 1, 0, 1010},
	{"dump-seqs",                        1, 0, 1036},
	{"dump-kbm",                         1, 0, 1011},
	{"load-seqs",                        1, 0, 2002},
	{"load-kbm",                         1, 0, 1012},
	{"load-alignments",                  1, 0, 1013},
	{"load-nodes",                       1, 0, 2000},
	{"load-clips",                       1, 0, 2001},
	{"aln-strand",                       1, 0, 1014},
	{"bubble-step",                      1, 0, 1015},
	{"tip-step",                         1, 0, 1016},
	{"ctg-min-length",                   1, 0, 1017},
	{"ctg-min-nodes",                    1, 0, 1018},
	{"minimal-output",                   0, 0, 1019},
	{"bin-complexity-cutoff",            1, 0, 1020},
	{"aln-dovetail",                     1, 0, 1021},
	{"no-local-graph-analysis",          0, 0, 1022},
	{"no-read-length-sort",              0, 0, 1023},
	{"keep-isolated-nodes",              0, 0, 1024},
	{"no-read-clip",                     0, 0, 1025},
	{"no-chainning-clip",                0, 0, 1026},
	{"aln-bestn",                        1, 0, 1027},
	{"aln-maxhit",                       1, 0, 1028},
	{"aln-kmer-sampling",                1, 0, 1029},
	{"aln-noskip",                       0, 0, 'A'},
	{"node-matched-bins",                1, 0, 1031},
	{"rescue-low-cov-edges",             0, 0, 1032},
	{"drop-low-cov-edges",               0, 0, 1033},
	{"mem-stingy",                       0, 0, 1034},
	{0, 0, 0, 0}
};

int usage(int level){
	printf(
	"WTDBG: De novo assembler for long noisy sequences\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: %s (%s)\n"
	"Usage: wtdbg2 [options] -i <reads.fa> -o <prefix> [reads.fa ...]\n"
	"Options:\n"
	" -i <string> Long reads sequences file (REQUIRED; can be multiple), []\n"
	" -o <string> Prefix of output files (REQUIRED), []\n"
	" -t <int>    Number of threads, 0 for all cores, [4]\n"
	" -f          Force to overwrite output files\n"
	" -x <string> Presets, comma delimited, []\n"
	"            preset1/rsII/rs: -p 21 -S 4 -s 0.05 -L 5000\n"
	"                    preset2: -p 0 -k 15 -AS 2 -s 0.05 -L 5000\n"
	"                    preset3: -p 19 -AS 2 -s 0.05 -L 5000\n"
	"                  sequel/sq\n"
	"               nanopore/ont:\n"
	"            (genome size < 1G: preset2) -p 0 -k 15 -AS 2 -s 0.05 -L 5000\n"
	"            (genome size >= 1G: preset3) -p 19 -AS 2 -s 0.05 -L 5000\n"
	"      preset4/corrected/ccs: -p 21 -k 0 -AS 4 -K 0.05 -s 0.5\n"
	" -g <number> Approximate genome size (k/m/g suffix allowed) [0]\n"
	" -X <float>  Choose the best <float> depth from input reads(effective with -g) [50.0]\n"
	" -L <int>    Choose the longest subread and drop reads shorter than <int> (5000 recommended for PacBio) [0]\n"
	"             Negative integer indicate tidying read names too, e.g. -5000.\n"
	" -k <int>    Kmer fsize, 0 <= k <= 23, [0]\n"
	" -p <int>    Kmer psize, 0 <= p <= 23, [21]\n"
	"             k + p <= 25, seed is <k-mer>+<p-homopolymer-compressed>\n"
	" -K <float>  Filter high frequency kmers, maybe repetitive, [1000.05]\n"
	"             >= 1000 and indexing >= (1 - 0.05) * total_kmers_count\n"
	" -S <float>  Subsampling kmers, 1/(<-S>) kmers are indexed, [4.00]\n"
	"             -S is very useful in saving memeory and speeding up\n"
	"             please note that subsampling kmers will have less matched length\n"
	" -l <float>  Min length of alignment, [2048]\n"
	" -m <float>  Min matched length by kmer matching, [200]\n"
	" -R          Enable realignment mode\n"
	" -A          Keep contained reads during alignment\n"
	" -s <float>  Min similarity, calculated by kmer matched length / aligned length, [0.05]\n"
	" -e <int>    Min read depth of a valid edge, [3]\n"
	" -q          Quiet\n"
	" -v          Verbose (can be multiple)\n"
	" -V          Print version information and then exit\n"
	" --help      Show more options\n"
	, TOSTR(VERSION), TOSTR(RELEASE)
	);
	if(level > 0){
		printf(
	" ** more options **\n"
	" --cpu <int>\n"
	"   See -t 0, default: all cores\n"
	" --input <string> +\n"
	"   See -i\n"
	//" --err-free-seq <string> +\n"
	//"   See -I. Error-free sequences will be firstly token for nodes, if --err-free-nodes is specified, only select nodes from those sequences\n"
	" --force\n"
	"   See -f\n"
	" --prefix <string>\n"
	"   See -o\n"
	" --preset <string>\n"
	"   See -x\n"
	" --kmer-fsize <int>\n"
	"   See -k 0\n"
	" --kmer-psize <int>\n"
	"   See -p 21\n"
	" --kmer-depth-max <float>\n"
	"   See -K 1000.05\n"
	" -E, --kmer-depth-min <int>\n"
	"   Min kmer frequency, [2]\n"
	//" --kmer-depth-min-filter\n"
	//"   See -F\n"
	//"   `wtdbg` uses a 4 Gbytes array to counting the occurence (0-3) of kmers in the way of counting-bloom-filter. It will reduce memory space largely\n"
	//"    Orphaned kmers won't appear in building kbm-index\n"
	" --kmer-subsampling <float>\n"
	"   See -S 4.0\n"
	" --kbm-parts <int>\n"
	"   Split total reads into multiple parts, index one part by one to save memory, [1]\n"
	" --aln-kmer-sampling <int>\n"
	"   Select no more than n seeds in a query bin, default: 256\n"
	" --dp-max-gap <int>\n"
	"   Max number of bin(256bp) in one gap, [4]\n"
	" --dp-max-var <int>\n"
	"   Max number of bin(256bp) in one deviation, [4]\n"
	" --dp-penalty-gap <int>\n"
	"   Penalty for BIN gap, [-7]\n"
	" --dp-penalty-var <int>\n"
	"   Penalty for BIN deviation, [-21]\n"
	" --aln-min-length <int>\n"
	"   See -l 2048\n"
	" --aln-min-match <int>\n"
	"   See -m 200. Here the num of matches counting basepair of the matched kmer's regions\n"
	" --aln-min-similarity <float>\n"
	"   See -s 0.05\n"
	" --aln-max-var <float>\n"
	"   Max length variation of two aligned fragments, default: 0.25\n"
	" --aln-dovetail <int>\n"
	"   Retain dovetail overlaps only, the max overhang size is <--aln-dovetail>, the value should be times of 256, -1 to disable filtering, default: 256\n"
	" --aln-strand <int>\n"
	"   1: forward, 2: reverse, 3: both. Please don't change the deault vaule 3, unless you exactly know what you are doing\n"
	" --aln-maxhit <int>\n"
	"   Max n hits for each read in build graph, default: 1000\n"
	" --aln-bestn <int>\n"
	"   Use best n hits for each read in build graph, 0: keep all, default: 500\n"
	"   <prefix>.alignments always store all alignments\n"
	" -R, --realign\n"
	"   Enable re-alignment, see --realn-kmer-psize=15, --realn-kmer-subsampling=1, --realn-min-length=2048, --realn-min-match=200, --realn-min-similarity=0.1, --realn-max-var=0.25\n"
	" --realn-kmer-psize <int>\n"
	"   Set kmer-psize in realignment, (kmer-ksize always eq 0), default:15\n"
	" --realn-kmer-subsampling <int>\n"
	"   Set kmer-subsampling in realignment, default:1\n"
	" --realn-min-length <int>\n"
	"   Set aln-min-length in realignment, default: 2048\n"
	" --realn-min-match <int>\n"
	"   Set aln-min-match in realignment, default: 200\n"
	" --realn-min-similarity <float>\n"
	"   Set aln-min-similarity in realignment, default: 0.1\n"
	" --realn-max-var <float>\n"
	"   Set aln-max-var in realignment, default: 0.25\n"
	" -A, --aln-noskip\n"
	"   Even a read was contained in previous alignment, still align it against other reads\n"
	//" --corr-mode <float>\n"
	//"   Default: 0.0. If set > 0 and set --g <genome_size>, will turn on correct-align mode.\n"
	//"   Wtdbg will select <genome_size> * <corr-mode> bases from reads of middle length, and align them aginst all reads.\n"
	//"   Then, wtdbg will correct them using POACNS, and query corrected sequences against all reads again\n"
	//"   In correct-align mode, --aln-bestn = unlimited, --no-read-clip, --no-chaining-clip. Will support those features in future\n"
	//" --corr-min <int>\n"
	//" --corr-max <int>\n"
	//"   For each read to be corrected, uses at least <corr-min> alignments, and at most <corr-max> alignments\n"
	//"   Default: --corr_min = 5, --corr-max = 10\n"
	//" --corr-cov <float>\n"
	//"   Default: 0.75. When aligning reads to be corrected, the alignments should cover at least <corr-cov> of read length\n"
	//" --corr-block-size <int>\n"
	//"   Default: 2048. MUST be times of 256bp. Used in POACNS\n"
	//" --corr-block-step <int>\n"
	//"   Default: 1536. MUST be times of 256bp. Used in POACNS\n"
	" --keep-multiple-alignment-parts\n"
	"   By default, wtdbg will keep only the best alignment between two reads after chainning. This option will disable it, and keep multiple\n"
	" --verbose +\n"
	"   See -v. -vvvv will display the most detailed information\n"
	" --quiet\n"
	"   See -q\n"
	" --limit-input <int>\n"
	"   Limit the input sequences to at most <int> M bp. Usually for test\n"
	" -L <int>, --tidy-reads <int>\n"
	"   Default: 0. Pick longest subreads if possible. Filter reads less than <--tidy-reads>. Please add --tidy-name or set --tidy-reads to nagetive value\n" 
	"   if want to rename reads. Set to 0 bp to disable tidy. Suggested value is 5000 for pacbio RSII reads\n"
	" --tidy-name\n"
	"   Rename reads into 'S%%010d' format. The first read is named as S0000000001\n"
	" --rdname-filter <string>\n"
	"   A file contains lines of reads name to be discarded in loading. If you want to filter reads by yourself, please also set -X 0\n"
	" --rdname-includeonly <string>\n"
	"   Reverse manner with --rdname-filter\n"
	//" --keep-name\n"
	//"   Keep orignal read names even with --tidy-reads, '-L 5000 --keep-name' equals '-L -5000'\n"
	" -g <number>, --genome-size <number>\n"
	"   Provide genome size, e.g. 100.4m, 2.3g. In this version, it is used with -X/--rdcov-cutoff in selecting reads just after readed all.\n"
	" -X <float>, --rdcov-cutoff <float>\n"
	"   Default: 50.0. Retaining 50.0 folds of genome coverage, combined with -g and --rdcov-filter.\n"
	" --rdcov-filter [0|1]\n"
	"   Default 0. Strategy 0: retaining longest reads. Strategy 1: retaining medain length reads. \n"
	" --err-free-nodes\n"
	"   Select nodes from error-free-sequences only. E.g. you have contigs assembled from NGS-WGS reads, and long noisy reads.\n"
	"   You can type '--err-free-seq your_ctg.fa --input your_long_reads.fa --err-free-nodes' to perform assembly somehow act as long-reads scaffolding\n"
	" --node-len <int>\n"
	"   The default value is 1024, which is times of KBM_BIN_SIZE(always equals 256 bp). It specifies the length of intervals (or call nodes after selecting).\n"
	"   kbm indexs sequences into BINs of 256 bp in size, so that many parameter should be times of 256 bp. There are: --node-len, --node-ovl, --aln-min-length, --aln-dovetail ."
	"   Other parameters are counted in BINs, --dp-max-gap, --dp-max-var .\n"
	" --node-matched-bins <int>\n"
	"   Min matched bins in a node, default:1\n"
	" --node-ovl <int>\n"
	"   Default: 256. Max overlap size between two adjacent intervals in any read. It is used in selecting best nodes representing reads in graph\n"
	" --node-drop <float>\n"
	"   Default: 0.25. Will discard an node when has more this ratio intervals are conflicted with previous generated node\n"
	" -e <int>, --edge-min=<int>\n"
	"   Default: 3. The minimal depth of a valid edge is set to 3. In another word, Valid edges must be supported by at least 3 reads\n"
	"   When the sequence depth is low, have a try with --edge-min 2. Or very high, try --edge-min 4\n"
	" --edge-max-span <int>\n"
	"   Default: 1024 BINs. Program will build edges of length no large than 1024\n"
	" --drop-low-cov-edges\n"
	"   Don't attempt to rescue low coverage edges\n"
	" --node-min <int>\n"
	"   Min depth of an interval to be selected as valid node. Defaultly, this value is automaticly the same with --edge-min.\n"
	" --node-max <int>\n"
	"   Nodes with too high depth will be regarded as repetitive, and be masked. Default: 200, more than 200 reads contain this node\n"
	" --ttr-cutoff-depth <int>, 0\n"
	" --ttr-cutoff-ratio <float>, 0.5\n"
	"   Tiny Tandom Repeat. A node located inside ttr will bring noisy in graph, should be masked. The pattern of such nodes is:\n"
	"   depth >= <--ttr-cutoff-depth>, and none of their edges have depth greater than depth * <--ttr-cutoff-ratio 0.5>\n"
	"   set --ttr-cutoff-depth 0 to disable ttr masking\n"
	" --dump-kbm <string>\n"
	"   Dump kbm index into file for loaded by `kbm` or `wtdbg`\n"
	" --dump-seqs <string>\n"
	"   Dump kbm index (only sequences, no k-mer index) into file for loaded by `kbm` or `wtdbg`\n"
	"   Please note: normally load it with --load-kbm, not with --load-seqs\n"
	" --load-kbm <string>\n"
	"   Instead of reading sequences and building kbm index, which is time-consumed, loading kbm-index from already dumped file.\n"
	"   Please note that, once kbm-index is mmaped by kbm -R <kbm-index> start, will just get the shared memory in minute time.\n"
	"   See `kbm` -R <your_seqs.kbmidx> [start | stop]\n"
	" --load-seqs <string>\n"
	"   Similar with --load-kbm, but only use the sequences in kbmidx, and rebuild index in process's RAM.\n"
	" --load-alignments <string> +\n"
	"   `wtdbg` output reads' alignments into <--prefix>.alignments, program can load them to fastly build assembly graph. Or you can offer\n"
	"   other source of alignments to `wtdbg`. When --load-alignment, will only reading long sequences but skip building kbm index\n"
	"   You can type --load-alignments <file> more than once to load alignments from many files\n"
	" --load-clips <string>\n"
	"   Combined with --load-nodes. Load reads clips. You can find it in `wtdbg`'s <--prefix>.clps\n"
	" --load-nodes <sting>\n"
	"   Load dumped nodes from previous execution for fast construct the assembly graph, should be combined with --load-clips. You can find it in `wtdbg`'s <--prefix>.1.nodes\n"
	" --bubble-step <int>\n"
	"   Max step to search a bubble, meaning the max step from the starting node to the ending node. Default: 40\n"
	" --tip-step <int>\n"
	"   Max step to search a tip, 10\n"
	" --ctg-min-length <int>\n"
	"   Min length of contigs to be output, 5000\n"
	" --ctg-min-nodes <int>\n"
	"   Min num of nodes in a contig to be ouput, 3\n"
	" --minimal-output\n"
	"   Will generate as less output files (<--prefix>.*) as it can\n"
	" --bin-complexity-cutoff <int>\n"
	"   Used in filtering BINs. If a BIN has less indexed valid kmers than <--bin-complexity-cutoff 2>, masks it.\n"
	" --no-local-graph-analysis\n"
	"   Before building edges, for each node, local-graph-analysis reads all related reads and according nodes, and builds a local graph to judge whether to mask it\n"
	"   The analysis aims to find repetitive nodes\n"
	" --no-read-length-sort\n"
	"   Defaultly, `wtdbg` sorts input sequences by length DSC. The order of reads affects the generating of nodes in selecting important intervals\n"
	" --keep-isolated-nodes\n"
	"   In graph clean, `wtdbg` normally masks isolated (orphaned) nodes\n"
	" --no-read-clip\n"
	"   Defaultly, `wtdbg` clips a input sequence by analyzing its overlaps to remove high error endings, rolling-circle repeats (see PacBio CCS), and chimera.\n"
	"   When building edges, clipped region won't contribute. However, `wtdbg` will use them in the final linking of unitigs\n"
	" --no-chainning-clip\n"
	"   Defaultly, performs alignments chainning in read clipping\n"
	"   ** If '--aln-bestn 0 --no-read-clip', alignments will be parsed directly, and less RAM spent on recording alignments\n"
	"\n"
		);
	}
	return (level < 0)? 1 : 0;
}

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

int main(int argc, char **argv){
	Graph *g;
	KBMPar *par, *rpar;
	KBM *kbm;
	FileReader *fr;
	BioSequence *seqs[2], *seq;
	chash *rdtaghash[2];
	cplist *pbs, *ngs, *pws;
	FILE *evtlog;
	char *prefix, *rdtag_filter[2], *dump_seqs, *load_seqs, *dump_kbm, *load_kbm, *load_nodes, *load_clips;
	char regtag[14];
	int len, tag_size, asyn_read, preset;
	u8i tot_bp, cnt, bub, tip, rep, yarn, max_bp, max_idx_bp, nfix, opt_flags;
	uint32_t i, j, k;
	int c, opt_idx, ncpu, only_fix, realign, node_cov, max_node_cov, exp_node_cov, min_bins, edge_cov, edge_span, store_low_cov_edge, reglen, regovl, bub_step, tip_step, rep_step;
	int frgtip_len, ttr_n_cov;
	int quiet, tidy_reads, filter_rd_strategy, tidy_rdtag, less_out, tip_like, cut_tip, rep_filter, out_alns, cnn_filter, log_rep, rep_detach, del_iso, rdclip, chainning, uniq_hit, bestn, rescue_low_edges;
	int min_ctg_len, min_ctg_nds, max_trace_end, max_overhang, overwrite, node_order, fast_mode, corr_min, corr_max, corr_bsize, corr_bstep, mem_stingy, num_index;
	double genome_size, genome_depx;
	float node_drop, node_mrg, ttr_e_cov, fval, cut_low_edges, corr_mode, corr_cov;
	pbs = init_cplist(4);
	ngs = init_cplist(4);
	pws = init_cplist(4);
	asyn_read = 1;
	ncpu = 4;
	mem_stingy = 0;
	tidy_reads = 0;
	tidy_rdtag = -1;
	preset = 0;
	genome_size = 0;
	genome_depx = 50.0;
	num_index = 1;
	filter_rd_strategy = 0;
	fast_mode = 0;
	// no longer supports  corr-mode 
	corr_mode = 0;
	corr_min = 5;
	corr_max = 10;
	corr_cov = 0.75;
	corr_bsize = 2048;
	corr_bstep = 2048 - 512;
	// -------
	max_bp = 0;
	max_idx_bp = 0LLU * 1000 * 1000 * 1000; // unlimited
	rdtag_filter[0] = NULL;
	rdtag_filter[1] = NULL;
	rdtaghash[0] = NULL;
	rdtaghash[1] = NULL;
	reglen = 1024;
	regovl = 256;
	node_drop = 0.25;
	node_mrg = 0.9;
	only_fix = 0;
	node_cov = 0; // will equal edge_cov, if no --node-cov
	max_node_cov = 200;
	exp_node_cov = 40;
	min_bins = 1;
	edge_cov = 0; // will be set to 3, if no genome_size available and no -e
	edge_span = 1024;
	rdclip = 1;
	chainning = 1;
	uniq_hit = 1;
	bestn = 500;
	ttr_n_cov = 0;
	ttr_e_cov = 0.5;
	dump_seqs = NULL;
	load_seqs = NULL;
	dump_kbm = NULL;
	load_kbm = NULL;
	load_clips = NULL;
	load_nodes = NULL;
	store_low_cov_edge = 1;
	cut_low_edges = 0.0;
	rescue_low_edges = 1;
	bub_step = 40;
	tip_step = 10;
	rep_step = 0;
	max_trace_end = 5;
	frgtip_len = 50000;
	prefix = NULL;
	overwrite = 0;
	less_out = 0;
	quiet = 0;
	rep_filter = 1;
	tip_like = 0;
	cut_tip = 1;
	cnn_filter = 1;
	log_rep = 1;
	rep_detach = 0;
	del_iso = 1;
	max_overhang = 256;
	min_ctg_len = 5000;
	min_ctg_nds = 3;
	node_order = 0;
	out_alns = 1;
	par = init_kbmpar();
	par->ksize = 0;
	par->psize = 21;
	par->kmer_mod = KBM_N_HASH * 4;
	par->kmin = 2;
	par->max_bgap = 4;
	par->max_bvar = 4;
	par->self_aln = 1; // won't perform B->A when existing A->B
	par->rd_len_order = 1;
	par->min_aln = 2048 / KBM_BIN_SIZE;
	par->min_mat = 200;
	par->min_sim = 0.05;
	par->aln_var = 0.25;
	realign = 0;
	rpar = init_kbmpar();
	rpar->ksize = 0;
	rpar->psize = 15;
	rpar->kmer_mod = KBM_N_HASH;
	rpar->kmin = 1;
	rpar->max_bgap = 4;
	rpar->max_bvar = 4;
	rpar->self_aln = 0; // won't perform B->A when existing A->B
	rpar->rd_len_order = 0;
	rpar->min_aln = 2048 / KBM_BIN_SIZE;
	rpar->min_mat = 200;
	rpar->min_sim = 0.1;
	rpar->aln_var = 0.25;
	opt_flags = 0;
	while((c = getopt_long(argc, argv, "ht:i:fo:x:E:k:p:K:S:l:m:s:RvqVe:L:Ag:X:", prog_opts, &opt_idx)) != -1){
		switch(c){
			case 't': ncpu = atoi(optarg); break;
			case 'i': push_cplist(pbs, optarg); break;
			//case 'I': push_cplist(ngs, optarg); par->rd_len_order = 0; break;
			case 'f': overwrite = 1; break;
			case 'o': prefix = optarg; break;
			case 'x':
					{
						char *ptr, *beg;
						beg = optarg;
						do {
							ptr = index(beg, ',');
							if(ptr) *ptr = 0;
							if(KBM_LOG){
								fprintf(KBM_LOGF, " -- Preset: '%s' --", beg); fflush(KBM_LOGF);
							}
							if(strcasecmp(beg, "preset1") == 0 || strcasecmp(beg, "rs") == 0 || strcasecmp(beg, "rsII") == 0){
								preset = 1;
							} else if(strcasecmp(beg, "preset2") == 0){
								preset = 2;
							} else if(strcasecmp(beg, "preset3") == 0){
								preset = 3;
							} else if(strcasecmp(beg, "sq") == 0 || strcasecmp(beg, "sequel") == 0){
								preset = -1;
							} else if(strcasecmp(beg, "ont") == 0 || strcasecmp(beg, "nanopore") == 0){
								preset = -1;
							} else if(strcasecmp(beg, "preset4") == 0 || strcasecmp(beg, "ccs") == 0 || strcasecmp(beg, "corrected") == 0){
								preset = 4;
							} else {
								fprintf(stderr, " ** ERROR: cannot recognize '%s' in '-x %s'\n", beg, optarg);
								exit(1);
							}
							if(KBM_LOG){
								fprintf(KBM_LOGF, "\n"); fflush(KBM_LOGF);
							}
							if(ptr){
								*ptr = ',';
								beg = ptr + 1;
							} else {
								break;
							}
						} while(1);
					}
					break;
			case 'k': par->ksize = atoi(optarg); opt_flags |= (1 << 1); break;
			case 'p': par->psize = atoi(optarg); opt_flags |= (1 << 0); break;
			case 'K': fval = atof(optarg); par->kmax = fval; par->ktop = fval - par->kmax; opt_flags |= (1 << 6); break;
			case 'E': par->kmin = atoi(optarg); break;
			case 'S': par->kmer_mod = UInt(atof(optarg) * KBM_N_HASH); opt_flags |= (1 << 2);break;
			case 'g': genome_size = mm_parse_num(optarg); break;
			case 'X': genome_depx = atof(optarg); break;
			case 3007: rdtag_filter[0] = optarg; break;
			case 3008: rdtag_filter[1] = optarg; break;
			case 2009: filter_rd_strategy = atoi(optarg); break;
			case 2005: par->max_bgap = atoi(optarg); break;
			case 2006: par->max_bvar = atoi(optarg); break;
			case 2007: par->pgap = atoi(optarg); break;
			case 2008: par->pvar = atoi(optarg); break;
			case 'l': par->min_aln = atoi(optarg) / KBM_BIN_SIZE; break;
			case 'm': par->min_mat = atoi(optarg); break;
			case 2004: par->aln_var = atof(optarg); break;
			case 's': par->min_sim = atof(optarg); opt_flags |= (1 << 3); break;
			//case 2010: corr_mode = atof(optarg); break;
			//case 2012: corr_min = atoi(optarg); break;
			//case 2013: corr_max = atoi(optarg); break;
			//case 2014: corr_cov = atof(optarg); break;
			//case 2015: corr_bsize = atoi(optarg); break;
			//case 2016: corr_bstep = atoi(optarg); break;
			case 2011: uniq_hit = 0; break;
			case 'v': KBM_LOG ++; break;
			case 'q': quiet = 1; break;
			case 'h': return usage(0);
			case 1000: return usage(1);
			case 'L':  tidy_reads = atoi(optarg); opt_flags |= (1 << 4); break;
			case 1001: tidy_rdtag = 1; break;
			case 1002: only_fix = 1; break;
			case 1003: max_bp = atol(optarg); break;
			case 1035: num_index = atoi(optarg); break;
			case 1004: reglen = atoi(optarg); break;
			case 1005: regovl = atoi(optarg); break;
			case 1006: node_drop = atof(optarg); break;
			case 'e':  edge_cov = atoi(optarg); break;
			case 3009: edge_span = atoi(optarg); break;
			case 1007: node_cov = atoi(optarg); break;
			case 1008: max_node_cov = atoi(optarg); break;
			case 1009: ttr_n_cov = atoi(optarg); break;
			case 1010: ttr_e_cov = atof(optarg); break;
			case 1036: dump_seqs = optarg; break;
			case 2002: load_seqs = optarg; break;
			case 1011: dump_kbm = optarg; break;
			case 1012: load_kbm = optarg; break;
			case 2000: load_nodes = optarg; break;
			case 2001: load_clips = optarg; break;
			case 1013: push_cplist(pws, optarg); break;
			case 1014: par->strand_mask = atoi(optarg); break;
			case 1015: bub_step = atoi(optarg); break;
			case 1016: tip_step = atoi(optarg); break;
			case 1017: min_ctg_len = atoi(optarg); break;
			case 1018: min_ctg_nds = atoi(optarg); break;
			case 1019: less_out = 1; break;
			case 1020: par->min_bin_degree = atoi(optarg); break;
			case 1021: max_overhang = atoi(optarg); break;
			case 1022: cnn_filter = 0; break;
			case 1023: par->rd_len_order = 0; break;
			case 1024: del_iso = 0; break;
			case 1025: rdclip = 0; break;
			case 1026: chainning = 0; break;
			case 1027: bestn = atoi(optarg); break;
			case 1028: par->max_hit = atoi(optarg); break;
			case 1029: par->ksampling = atoi(optarg); break;
			case 'R': realign = 1; break;
			case 3001: rpar->psize = atoi(optarg); break;
			case 3002: rpar->kmer_mod = UInt(atof(optarg) * KBM_N_HASH); break;
			case 3003: rpar->min_aln = atoi(optarg) / KBM_BIN_SIZE; break;
			case 3004: rpar->min_mat = atoi(optarg); break;
			case 3005: rpar->min_sim = atof(optarg); break;
			case 3006: rpar->aln_var = atof(optarg); break;
			case 'A':  par->skip_contained = 0; opt_flags |= (1 << 5); break;
			case 1031: min_bins = atoi(optarg); break;
			case 1032: rescue_low_edges = 1; break;
			case 1033: rescue_low_edges = 0; break;
			case 'V': fprintf(stdout, "wtdbg2 %s\n", TOSTR(VERSION)); return 0;
			case 1034: mem_stingy = 1; break;
			default: return usage(-1);
		}
	}
	if(optind == 1) return usage(-1);
	if(optind < argc){
		fprintf(stderr, "WARNING: unused command-line arguments. For multiple input files, please apply multiple -i.\n");
		fprintf(stderr, "WARNING: try to recognize and add to input files list\n");
		for(c=optind;c<argc;c++){
			if(file_exists(argv[c])){
				fprintf(stderr, " * \"%s\" exists, added.\n", argv[c]);
				push_cplist(pbs, argv[c]);
			}
		}
	}
	if(prefix == NULL) {
		fprintf(stderr, "ERROR: please specify the output prefix with -o\n");
		return 1;
	}
	if(load_seqs == NULL && load_kbm == NULL && pbs->size + ngs->size == 0) {
		fprintf(stderr, "ERROR: please specify the input with -i/--load-seqs/--load-kbm\n");
		return 1;
	}
	if((reglen % KBM_BIN_SIZE)){
		reglen = ((reglen + KBM_BIN_SIZE - 1) / KBM_BIN_SIZE) * KBM_BIN_SIZE;
		fprintf(stderr, " ** Adjust -j to %d\n", reglen);
	}
	if(!overwrite && file_exists(prefix)){
		fprintf(stderr, "File exists! '%s'\n\n", prefix);
		return usage(-1);
	}
	if(max_idx_bp == 0) max_idx_bp = 0xFFFFFFFFFFFFFFFFLLU;
	if(preset == -1){
		if(genome_size && genome_size < 1000000000LLU){
			preset = 2;
		} else {
			preset = 3;
		}
	}
	switch(preset){
		case 1:
			if(!(opt_flags & (1 << 1))) par->ksize = 0;
			if(!(opt_flags & (1 << 0))) par->psize = 21;
			if(!(opt_flags & (1 << 2))) par->kmer_mod = 4 * KBM_N_HASH;
			if(!(opt_flags & (1 << 3))) par->min_sim = 0.05;
			if(!(opt_flags & (1 << 5))) par->skip_contained = 1;
			if(!(opt_flags & (1 << 4))) tidy_reads = 5000;
			break;
		case 2:
			if(!(opt_flags & (1 << 1))) par->ksize = 15;
			if(!(opt_flags & (1 << 0))) par->psize = 0;
			if(!(opt_flags & (1 << 2))) par->kmer_mod = 2 * KBM_N_HASH;
			if(!(opt_flags & (1 << 3))) par->min_sim = 0.05;
			if(!(opt_flags & (1 << 5))) par->skip_contained = 0;
			if(!(opt_flags & (1 << 4))) tidy_reads = 5000;
			break;
		case 3:
			if(!(opt_flags & (1 << 1))) par->ksize = 0;
			if(!(opt_flags & (1 << 0))) par->psize = 19;
			if(!(opt_flags & (1 << 2))) par->kmer_mod = 2 * KBM_N_HASH;
			if(!(opt_flags & (1 << 3))) par->min_sim = 0.05;
			if(!(opt_flags & (1 << 5))) par->skip_contained = 0;
			if(!(opt_flags & (1 << 4))) tidy_reads = 5000;
			break;
		case 4:
			if(!(opt_flags & (1 << 1))) par->ksize = 0;
			if(!(opt_flags & (1 << 0))) par->psize = 21;
			if(!(opt_flags & (1 << 2))) par->kmer_mod = 4 * KBM_N_HASH;
			if(!(opt_flags & (1 << 3))) par->min_sim = 0.5;
			if(!(opt_flags & (1 << 5))) par->skip_contained = 0;
			if(!(opt_flags & (1 << 6))){ par->kmax = 0; par->ktop = 0.05; }
			//if(!(opt_flags & (1 << 4))) tidy_reads = 5000;
			break;
	}
	if(par->ksize + par->psize > KBM_MAX_KSIZE){
		fprintf(stderr, " -- Invalid kmer size %d+%d=%d > %d in %s -- %s:%d --\n", par->ksize, par->psize, par->ksize + par->psize, KBM_MAX_KSIZE,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(quiet){
		int devnull;
		devnull = open("/dev/null", O_WRONLY);
		dup2(devnull, STDERR_FILENO);
	}
	if(tidy_rdtag == -1){
		if(tidy_reads >= 0){
			tidy_rdtag = 0;
		} else {
			tidy_rdtag = 1;
		}
	}
	if(tidy_reads < 0) tidy_reads = - tidy_reads;
	max_bp *= 1000000;
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	if(ncpu <= 0 && _sig_proc_deamon) ncpu = _sig_proc_deamon->ncpu;
	if(ncpu <= 0){
		fprintf(stderr, " -- Invalid cpu number '%d' in %s -- %s:%d --\n", ncpu, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(load_kbm){
		fprintf(KBM_LOGF, "[%s] loading kbm index from %s\n", date(), load_kbm);
		if((kbm = mem_find_obj_file(&kbm_obj_desc, load_kbm, NULL, NULL, NULL, NULL, 0)) == NULL){
			fprintf(KBM_LOGF, " -- cannot find mmap object %s --\n", load_kbm);
			fprintf(KBM_LOGF, " -- try read from file --\n");
			kbm = mem_read_obj_file(&kbm_obj_desc, load_kbm, NULL, NULL, NULL, NULL);
		}
		nfix = 0;
		tot_bp = 0;
		for(i=0;i<kbm->reads->size;i++) tot_bp += kbm->reads->buffer[i].bincnt * KBM_BIN_SIZE;
		fprintf(KBM_LOGF, "[%s] Done. %u sequences, %llu bp, parameter('-S %d')\n", date(), (u4i)kbm->reads->size, tot_bp, kbm->par->kmer_mod / KBM_N_HASH);
		{
			// check KBMPar
			if((opt_flags >> 0) & 0x01){
				if(kbm->par->psize != par->psize){
					fprintf(KBM_LOGF, " ** -p is different, %d != %d\n", kbm->par->psize, par->psize); exit(1);
				}
			} else {
				par->psize = kbm->par->psize;
			}
			if((opt_flags >> 1) & 0x01){
				if(kbm->par->ksize != par->ksize){
					fprintf(KBM_LOGF, " ** -k is different, %d != %d\n", kbm->par->ksize, par->ksize); exit(1);
				}
			} else {
				par->ksize = kbm->par->ksize;
			}
			if((opt_flags >> 2) & 0x01){
				if(kbm->par->kmer_mod != par->kmer_mod){
					fprintf(KBM_LOGF, " ** -S is different, %d != %d\n", kbm->par->kmer_mod / KBM_N_HASH, par->kmer_mod / KBM_N_HASH); exit(1);
				}
			} else {
				par->kmer_mod = kbm->par->kmer_mod;
			}
			if((opt_flags >> 3) & 0x01){
				if(kbm->par->rd_len_order != par->rd_len_order){
					fprintf(KBM_LOGF, " ** par->rd_len_order is different, %d != %d\n", kbm->par->rd_len_order, par->rd_len_order); exit(1);
				}
			} else {
				par->rd_len_order = kbm->par->rd_len_order;
			}
		}
	} else if(load_seqs){
		fprintf(KBM_LOGF, "[%s] loading kbm index from %s\n", date(), load_seqs);
		if((kbm = mem_find_obj_file(&kbm_obj_desc, load_seqs, NULL, NULL, NULL, NULL, 0)) == NULL){
			fprintf(KBM_LOGF, " -- cannot find mmap object %s --\n", load_seqs);
			fprintf(KBM_LOGF, " -- try read from file --\n");
			kbm = mem_read_obj_file(&kbm_obj_desc, load_seqs, NULL, NULL, NULL, NULL);
		}
		kbm = clone_seqs_kbm(kbm, par);
		nfix = 0;
		tot_bp = 0;
		for(i=0;i<kbm->reads->size;i++) tot_bp += kbm->reads->buffer[i].bincnt * KBM_BIN_SIZE;
		fprintf(KBM_LOGF, "[%s] Done. %u sequences, %llu bp\n", date(), (u4i)kbm->reads->size, tot_bp);
	} else {
		kbm = init_kbm(par);
		for(i=0;i<2;i++){
			if(rdtag_filter[i]){
				rdtaghash[i] = init_chash(1023);
				fr = open_filereader(rdtag_filter[i], 0);
				while(readline_filereader(fr)){
					char *str = strdup(fr->line->string);
					put_chash(rdtaghash[i], str);
				}
				close_filereader(fr);
				tidy_reads = 0;
			} else {
				rdtaghash[i] = NULL;
			}
		}
		fprintf(KBM_LOGF, "[%s] loading reads\n", date());
		tot_bp = 0;
		nfix = 0;
		seqs[0] = init_biosequence();
		seqs[1] = init_biosequence();
		regex_t reg;
		regmatch_t mats[3];
		int z;
		z = regcomp(&reg, "^(.+?)/[0-9]+_[0-9]+$", REG_EXTENDED);
		if(z){
			regerror(z, &reg, regtag, 13);
			fprintf(stderr, " -- REGCOMP: %s --\n", regtag); fflush(stderr);
			return 1;
		}
		for(j=0;j<2;j++){
			if(j == 0){
				if(ngs->size == 0){
					continue;
				} else {
					fr = open_all_filereader(ngs->size, ngs->buffer, asyn_read);
				}
			} else {
				if(pbs->size == 0){
					continue;
				} else {
					fr = open_all_filereader(pbs->size, pbs->buffer, asyn_read);
				}
			}
			k = 0;
			reset_biosequence(seqs[0]);
			reset_biosequence(seqs[1]);
			while(1){
				int has = readseq_filereader(fr, seqs[k]);
				if(tidy_reads){
					if(has){
						if((z = regexec(&reg, seqs[k]->tag->string, 3, mats, 0)) == 0){
							trunc_string(seqs[k]->tag, mats[1].rm_eo);
						} else if(z != REG_NOMATCH){
							regerror(z, &reg, regtag, 13);
							fprintf(stderr, " -- REGEXEC: %s --\n", regtag); fflush(stderr);
						}
						//fprintf(stderr, "1: %s len=%d\n", seqs[k]->tag->string, seqs[k]->seq->size); fflush(stderr);
						//fprintf(stderr, "2: %s len=%d\n", seqs[!k]->tag->string, seqs[!k]->seq->size); fflush(stderr);
						if(seqs[k]->tag->size == seqs[!k]->tag->size && strcmp(seqs[k]->tag->string, seqs[!k]->tag->string) == 0){
							if(seqs[k]->seq->size > seqs[!k]->seq->size){
								k = !k;
							}
							continue;
						} else {
							seq = seqs[!k];
							k = !k;
						}
					} else {
						seq = seqs[!k];
					}
					if(seq->seq->size < tidy_reads){
						if(has) continue;
						else break;
					}
					if(tidy_rdtag){
						sprintf(regtag, "S%010llu", (u8i)kbm->reads->size);
						clear_string(seq->tag);
						append_string(seq->tag, regtag, 11);
					}
				} else {
					if(has == 0) break;
					seq = seqs[k];
					if(rdtaghash[0]){
						if(exists_chash(rdtaghash[0], seq->tag->string)){
							continue;
						}
					}
					if(rdtaghash[1]){
						if(!exists_chash(rdtaghash[1], seq->tag->string)){
							continue;
						}
					}
				}
				tag_size = seq->tag->size;
				for(i=0;(int)i<seq->seq->size;i+=WT_MAX_RDLEN){
					len = num_min(seq->seq->size - i, WT_MAX_RDLEN);
					if(i){
						append_string(seq->tag, "_V", 2);
						add_int_string(seq->tag, i / WT_MAX_RDLEN);
					}
					if(!KBM_LOG && (kbm->reads->size % 10000) == 0){ fprintf(KBM_LOGF, "\r%u", (u4i)kbm->reads->size); fflush(KBM_LOGF); }
					//fprintf(stderr, " -- %s len=%d in %s -- %s:%d --\n", seq->tag->string, seq->seq->size, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					if(kbm->reads->size >= WT_MAX_RD){
						fprintf(stderr, " -- Read Number Out of Range: %u --\n", (u4i)kbm->reads->size); fflush(stderr);
						break;
					}
					push_kbm(kbm, seq->tag->string, seq->tag->size, seq->seq->string + i, len);
					if(i){ seq->tag->size = tag_size; seq->tag->string[tag_size] = '\0'; }
					if(j == 0) nfix ++;
				}
				tot_bp += seq->seq->size;
				if(max_bp && tot_bp >= max_bp){ break; }
				if(has == 0) break;
				if(kbm->reads->size >= WT_MAX_RD){
					fprintf(stderr, " -- Read Number Out of Range: %u --\n", (u4i)kbm->reads->size); fflush(stderr);
					break;
				}
			}
			close_filereader(fr);
			for(i=0;i<2;i++){
				if(rdtaghash[i]){
					char **str;
					reset_iter_chash(rdtaghash[i]);
					while((str = ref_iter_chash(rdtaghash[i]))){
						free(*str);
					}
					free_chash(rdtaghash[i]);
				}
			}
		}
		regfree(&reg);
		free_biosequence(seqs[0]);
		free_biosequence(seqs[1]);
		if(!KBM_LOG){ fprintf(KBM_LOGF, "\r%u reads", (unsigned)kbm->reads->size); fflush(KBM_LOGF); }
		{
			if(par->rd_len_order && genome_size > 0 && genome_depx > 0){
				cnt = genome_size * genome_depx;
				if(cnt < tot_bp){
					fprintf(KBM_LOGF, "\n[%s] filtering from %u reads (>=%u bp), %llu bp. Try selecting %llu bp", date(), (unsigned)kbm->reads->size, tidy_reads, tot_bp, cnt); fflush(KBM_LOGF);
					tot_bp = filter_reads_kbm(kbm, cnt, filter_rd_strategy);
				}
			}
			ready_kbm(kbm);
			fprintf(KBM_LOGF, "\n[%s] Done, %u reads (>=%u bp), %llu bp, %u bins\n", date(), (unsigned)kbm->reads->size, tidy_reads, tot_bp, (u4i)kbm->bins->size); fflush(KBM_LOGF);
			if(dump_seqs){
				FILE *dump;
				fprintf(KBM_LOGF, "[%s] dump kbm-index (only seqs) to %s ...", date(), dump_seqs); fflush(KBM_LOGF);
				dump = open_file_for_write(dump_seqs, NULL, 1);
				mem_dump_obj_file(kbm, 1, &kbm_obj_desc, 1, 0, dump);
				fclose(dump);
				fprintf(KBM_LOGF, " Done\n"); fflush(KBM_LOGF);
			}
		}
	}
	print_proc_stat_info(0);
	if(edge_cov <= 0){
		if(genome_size > 0){
			float dep;
			dep = tot_bp / genome_size;
			if(dep <= 40){
				edge_cov = 2;
			} else if(dep >= 80){
				edge_cov = 4;
			} else {
				edge_cov = 3;
			}
		} else {
			edge_cov = 3;
		}
		fprintf(KBM_LOGF, "[%s] Set --edge-cov to %d\n", date(), edge_cov); fflush(KBM_LOGF);
	}
	//if(genome_size <= 0 && corr_mode > 0){
		//fprintf(KBM_LOGF, "[%s] MUST set -g <?> with --corr-mode %f\n", date(), corr_mode); fflush(KBM_LOGF);
		//return 1;
	//}
	if(node_cov == 0) node_cov = edge_cov;
	fprintf(KBM_LOGF, "KEY PARAMETERS: -k %d -p %d -K %f %s-S %f -s %f -g %llu -X %f -e %d -L %d\n",
		par->ksize, par->psize, par->kmax + par->ktop, par->skip_contained? "" : "-A ", ((double)par->kmer_mod) / KBM_N_HASH, par->min_sim, (u8i)genome_size, genome_depx, edge_cov, tidy_reads);
	g = init_graph(kbm);
	{
		g->rpar = realign? rpar : NULL;
		g->genome_size = genome_size;
		g->num_index = num_index;
		//g->corr_mode = (corr_mode > 0 && genome_size > 0)? 1 : 0;
		//g->corr_gcov = corr_mode;
		//g->corr_min = corr_min;
		//g->corr_max = corr_max;
		//g->corr_cov = corr_cov;
		//g->corr_bsize = corr_bsize;
		//g->corr_bstep = corr_bstep;
		g->node_order = node_order;
		g->mem_stingy = mem_stingy;
		g->reglen = reglen / KBM_BIN_SIZE;
		g->regovl = regovl / KBM_BIN_SIZE;
		g->max_overhang = max_overhang / KBM_BIN_SIZE;
		g->node_max_conflict = node_drop;
		g->node_merge_cutoff = node_mrg;
		g->min_node_cov = node_cov;
		g->max_node_cov_sg = node_cov;
		g->max_node_cov = max_node_cov;
		g->exp_node_cov = exp_node_cov;
		g->min_node_mats = min_bins;
		g->min_edge_cov = edge_cov;
		g->max_edge_span = edge_span;
		g->max_sg_end = max_trace_end;
		g->store_low_cov_edge = store_low_cov_edge;
		g->bub_step = bub_step;
		g->tip_step = tip_step;
		g->rep_step = rep_step;
		g->min_ctg_len = min_ctg_len;
		g->min_ctg_nds = min_ctg_nds;
		g->n_fix = nfix;
		g->only_fix = only_fix;
		g->rep_filter = rep_filter;
		g->rep_detach = rep_detach;
		g->cut_tip = cut_tip;
		g->chainning_hits = chainning;
		g->uniq_hit = uniq_hit;
		g->bestn = bestn;
		g->minimal_output = less_out;
	}
	g->par = par;
	if(log_rep && !less_out){
		evtlog = open_file_for_write(prefix, ".events", 1);
	} else evtlog = NULL;
	if(load_nodes && load_clips){
		fprintf(KBM_LOGF, "[%s] loading nodes from %s ... ", date(), load_nodes); fflush(KBM_LOGF);
		FileReader *clp = open_filereader(load_clips, asyn_read);
		FileReader *nds = open_filereader(load_nodes, asyn_read);
		load_nodes_graph(g, clp, nds);
		close_filereader(clp);
		close_filereader(nds);
		fprintf(KBM_LOGF, " %llu nodes\n", (u8i)g->nodes->size);
		print_proc_stat_info(0);
	} else if(pws->size){
		fprintf(KBM_LOGF, "[%s] loading alignments from ", date());
		for(i=0;i<pws->size;i++){
			if(i){
				fprintf(KBM_LOGF, ",\"%s\"", pws->buffer[i]);
			} else {
				fprintf(KBM_LOGF, "\"%s\"", pws->buffer[i]);
			}
		}
		fprintf(KBM_LOGF, "\n");
		fr = open_all_filereader(pws->size, pws->buffer, asyn_read);
		build_nodes_graph(g, max_idx_bp, ncpu, fr, rdclip, prefix, NULL);
		close_filereader(fr);
		fprintf(KBM_LOGF, "[%s] Done, %llu nodes\n", date(), (unsigned long long)g->nodes->size);
	} else {
		fprintf(KBM_LOGF, "[%s] generating nodes, %d threads\n", date(), ncpu);
		build_nodes_graph(g, max_idx_bp, ncpu, NULL, rdclip, prefix, dump_kbm);
		fprintf(KBM_LOGF, "[%s] Done, %llu nodes\n", date(), (unsigned long long)g->nodes->size);
	}
	if(load_nodes == NULL || strlen(load_nodes) != strlen(prefix) + strlen(".1.nodes") || strncmp(load_nodes, prefix, strlen(prefix)) || strcmp(load_nodes + strlen(prefix), ".1.nodes")){
		generic_print_graph(g, print_nodes_graph, prefix, ".1.nodes");
	}
	if(1){
		estimate_genome_size(g, tot_bp, KBM_LOGF);
		cnt = mask_nodes_by_cov_graph(g, evtlog);
		fprintf(KBM_LOGF, "[%s] masked %llu high coverage nodes (>%d or <%d)\n", date(), (unsigned long long)cnt, max_node_cov, node_cov);
	}
	if(cnn_filter){
		cnt = mask_nodes_by_connectivity_graph(g, ncpu, evtlog);
		fprintf(KBM_LOGF, "[%s] masked %llu repeat-like nodes by local subgraph analysis\n", date(), (unsigned long long)cnt);
	}
	if(tip_like){
		cnt = mask_possible_tip_nodes_graph(g);
		fprintf(KBM_LOGF, "[%s] masked %llu tip-like nodes\n", date(), (unsigned long long)cnt);
	}
	fprintf(KBM_LOGF, "[%s] generating edges\n", date());
	build_edges_graph(g, ncpu, evtlog);
	fprintf(KBM_LOGF, "[%s] Done, %llu edges\n", date(), (unsigned long long)g->edges->size);
	if(ttr_n_cov){
		//print_node_edges_cov_graph(g, evtlog);
		cnt = mask_nodes_by_edge_cov_graph(g, ttr_n_cov, ttr_e_cov, evtlog);
		fprintf(KBM_LOGF, "[%s] deleted %llu nodes, might be tandom repeats\n", date(), (unsigned long long)cnt);
	}
	if(!less_out) generic_print_graph(g, print_reads_graph, prefix, ".1.reads");
	if(!less_out) generic_print_graph(g, print_dot_full_graph,   prefix, ".1.dot.gz");
	fprintf(KBM_LOGF, "[%s] graph clean\n", date()); fflush(KBM_LOGF);
	if(0){
		cnt = mask_read_weak_regs_graph(g, ncpu);
		fprintf(KBM_LOGF, "[%s] masked %llu regions(%d bp) as unreliable, total regs %llu\n", date(), (unsigned long long)cnt, reglen, (u8i)g->regs->size);
	}
	if(cut_low_edges){
		cnt = cut_relative_low_cov_edges_graph(g, cut_low_edges);
		fprintf(KBM_LOGF, "[%s] cut %llu low cov edges\n", date(), (unsigned long long)cnt);
	}
	if(rescue_low_edges){
		//cnt = rescue_low_cov_tip_edges_graph(g);
		//cnt = rescue_low_cov_edges_graph(g);
		cnt = rescue_mercy_edges_graph(g);
		fprintf(KBM_LOGF, "[%s] rescued %llu low cov edges\n", date(), (unsigned long long)cnt);
	}
	cnt = cut_binary_edges_graph(g);
	fprintf(KBM_LOGF, "[%s] deleted %llu binary edges\n", date(), (unsigned long long)cnt);
	if(!g->rep_detach && del_iso){
		cnt = del_isolated_nodes_graph(g, evtlog);
		fprintf(KBM_LOGF, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
	}
	//cnt = reduce_transitive_edges_graph(g);
	cnt = myers_transitive_reduction_graph(g, 1.2f);
	set_init_ends_graph(g);
	fprintf(KBM_LOGF, "[%s] cut %llu transitive edges\n", date(), (unsigned long long)cnt);
	if(del_iso){
		cnt = del_isolated_nodes_graph(g, evtlog);
		if(cnt){
			fprintf(KBM_LOGF, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
		}
	}
	if(!less_out) generic_print_graph(g, print_dot_graph,   prefix, ".2.dot.gz");
	{
		bub = tip = rep = yarn = 0;
		u8i high = 0;
		int safe = 1;
		do {
			c = 0;
			do {
				cnt = trim_tips_graph(g, tip_step, bub > 0);
				tip += cnt;
				if(cnt) c = 1;
			} while(cnt);
			do {
				cnt = pop_bubbles_graph(g, bub_step, safe);
				bub += cnt;
				if(cnt) c = 1;
			} while(cnt);
			do {
				cnt = trim_blunt_tips_graph(g);
				tip += cnt;
				if(cnt) c = 1;
			} while(cnt);
			do {
				cnt = pop_bubbles_graph(g, bub_step, safe);
				bub += cnt;
				if(cnt) c = 1;
			} while(cnt);
			if(c) continue;
			if(safe == 1){
				safe = 0;
				c = 1;
				continue;
			}
			do {
				cnt = resolve_yarns_graph(g, bub_step * 5);
				yarn += cnt;
				if(cnt) c = 1;
			} while(cnt);
		} while(c);
		do {
			c = 0;
			cnt = rescue_high_cov_edges_graph(g, 2, 20);
			if(cnt){
				high += cnt;
				c = 1;
			}
		} while(c);
		if(bub + tip + high){ fprintf(KBM_LOGF, "[%s] %llu bubbles; %llu tips; %llu yarns; rescued %llu high edges\n", date(), bub, tip, yarn, high); fflush(KBM_LOGF); }
		//if(bub + tip + yarn){ fprintf(KBM_LOGF, "[%s] %llu bubbles; %llu tips; %llu yarns\n", date(), bub, tip, yarn); fflush(KBM_LOGF); }
	}
	if(del_iso){
		cnt = del_isolated_nodes_graph(g, evtlog);
		fprintf(KBM_LOGF, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
	}
	if(!less_out) generic_print_graph(g, print_dot_graph,   prefix, ".3.dot.gz");
	rep = mask_all_branching_nodes_graph(g);
	fprintf(KBM_LOGF, "[%s] cut %llu branching nodes\n", date(), rep);
	if(del_iso){
		cnt = del_isolated_nodes_graph(g, evtlog);
		fprintf(KBM_LOGF, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
	}
	fprintf(KBM_LOGF, "[%s] building unitigs\n", date());
	gen_unitigs_graph(g);
	//fprintf(KBM_LOGF, "[%s] trimming and extending unitigs by local assembly, %d threads\n", date(), ncpu);
	unitigs2frgs_graph(g, ncpu);
	if(!less_out) generic_print_graph(g, print_frgs_nodes_graph, prefix, ".frg.nodes");
	fprintf(KBM_LOGF, "[%s] generating links\n", date());
	cnt = gen_lnks_graph(g, ncpu, evtlog);
	fprintf(KBM_LOGF, "[%s] generated %llu links\n", date(), cnt);
	if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.dot.gz");
	if(1){
		cnt = rescue_weak_tip_lnks_graph(g);
		fprintf(KBM_LOGF, "[%s] rescue %llu weak links\n", date(), (unsigned long long)cnt);
	}
	cnt = cut_binary_lnks_graph(g, evtlog);
	fprintf(KBM_LOGF, "[%s] deleted %llu binary links\n", date(), (unsigned long long)cnt);
	//cnt = reduce_transitive_lnks_graph(g);
	cnt = myers_transitive_reduction_frg_graph(g, 10000.1f / KBM_BIN_SIZE);
	fprintf(KBM_LOGF, "[%s] cut %llu transitive links\n", date(), (unsigned long long)cnt);
	cnt = remove_boomerangs_frg_graph(g, 30 * 1000 / KBM_BIN_SIZE);
	fprintf(KBM_LOGF, "[%s] remove %llu boomerangs\n", date(), (unsigned long long)cnt);
	cnt = cut_weak_branches_frg_graph(g);
	fprintf(KBM_LOGF, "[%s] remove %llu weak branches\n", date(), (unsigned long long)cnt);
	//cnt = cut_low_cov_lnks_graph(g, 1);
	//fprintf(KBM_LOGF, "[%s] deleted %llu low cov links\n", date(), (unsigned long long)cnt);
	//if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.2.dot");
	cnt = trim_frgtips_graph(g, frgtip_len);
	fprintf(KBM_LOGF, "[%s] cut %llu tips\n", date(), (unsigned long long)cnt);
	//if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.3.dot");
	bub = 0;
	do {
		cnt = pop_frg_bubbles_graph(g, bub_step);
		bub += cnt;
	} while(cnt);
	fprintf(KBM_LOGF, "[%s] pop %llu bubbles\n", date(), bub);
	{
		cnt = 0;
		while(1){
			u8i c;
			if((c = detach_repetitive_frg_graph(g, 100 * 1000 / KBM_BIN_SIZE)) == 0){
				break;
			}
			cnt += c;
		}
		fprintf(KBM_LOGF, "[%s] detached %llu repeat-associated paths\n", date(), (unsigned long long)cnt);
	}
	cnt = trim_frgtips_graph(g, frgtip_len);
	fprintf(KBM_LOGF, "[%s] cut %llu tips\n", date(), (unsigned long long)cnt);
	if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".ctg.dot.gz");
	fprintf(KBM_LOGF, "[%s] building contigs\n", date());
	cnt = gen_contigs_graph(g, evtlog);
	fprintf(KBM_LOGF, "[%s] searched %llu contigs\n", date(), (unsigned long long)cnt);
	if(0){
		cnt = gen_complex_contigs_graph(g);
		u8i sum;
		seqletv *qs;
		sum = 0;
		for(i=g->major_nctg;i<g->ctgs->size;i++){
			qs = (seqletv*)get_vplist(g->ctgs, i);
			sum += qs->buffer[qs->size - 1].off + qs->buffer[qs->size - 1].len;
		}
		fprintf(KBM_LOGF, "[%s] added %llu unsolved repetitive contigs, %llu bp\n", date(), (unsigned long long)cnt, sum);
	}
	n50_stat_contigs_graph(g);
	//cnt = generic_print_graph(g, print_isolated_nodes_dot_graph, prefix, ".4.dot");
	//fprintf(KBM_LOGF, "[%s] %llu nodes not in contigs\n", date(), (unsigned long long)cnt);
	//cnt = count_isolated_reads_graph(g);
	//fprintf(KBM_LOGF, "[%s] %llu reads not in contigs\n", date(), (unsigned long long)cnt);
	cnt = print_ctgs_graph(g, 0, 0, g->major_nctg, prefix, ".ctg.lay.gz", ncpu, evtlog);
	if(0){
		fprintf(KBM_LOGF, "[%s] outputing reptigs\n", date());
		cnt = print_ctgs_graph(g, cnt, g->major_nctg, g->ctgs->size, prefix, ".rtg.lay", ncpu, evtlog);
	}
	if(evtlog) fclose(evtlog);
	free_cplist(pbs);
	free_cplist(ngs);
	free_cplist(pws);
	if(load_kbm == NULL) free_kbm(kbm);
	free_kbmpar(par);
	free_kbmpar(rpar);
	free_graph(g);
	fprintf(KBM_LOGF, "[%s] Program Done\n", date());
	END_STAT_PROC_INFO(stderr);
	return 0;
}
