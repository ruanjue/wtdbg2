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

#include "wtpoa.h"

#ifndef VERSION
#define VERSION 0.0
#endif
#ifndef RELEASE
#define RELEASE 19830203
#endif

int usage(){
	printf(
	"WTPOA-CNS: Consensuser for wtdbg using PO-MSA\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: %s (%s)\n"
	"Usage: wtpoa-cns [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [4]\n"
	" -d <string> Reference sequences for SAM input, will invoke sorted-SAM input mode\n"
	" -u          Only process reference regions present in/between SAM alignments\n"
	" -r          Force to use reference mode\n"
	" -p <string> Similar with -d, but translate SAM into wtdbg layout file\n"
	" -i <string> Input file(s) *.ctg.lay from wtdbg, +, [STDIN]\n"
	"             Or sorted SAM files when having -d\n"
	" -o <string> Output files, [STDOUT]\n"
	" -f          Force overwrite\n"
	" -j <int>    Expected max length of node, or say the overlap length of two adjacent units in layout file, [1500] bp\n"
	" -b <int>    Bonus for tri-bases match, [0]\n"
	" -M <int>    Match score, [2]\n"
	" -X <int>    Mismatch score, [-5]\n"
	" -I <int>    Insertion score, [-2]\n"
	" -D <int>    Deletion score, [-4]\n"
	" -H <float>  Homopolymer merge score used in dp-call-cns mode, [-3]\n"
	" -B <int>    Bandwidth, [96]\n"
	" -W <int>    Window size in the middle of the first read for fast align remaining reads, [200]\n"
	"             If $W is negative, will disable fast align, but use the abs($W) as Band align score cutoff\n"
	" -w <int>    Min size of aligned size in window, [$W * 0.5]\n"
	"             In sorted-SAM input mode, -w is the sliding window size [2000]\n"
	" -A          Abort TriPOA when any read cannot be fast aligned, then try POA\n"
	" -S <int>    Shuffle mode, 0: don't shuffle reads, 1: by shared kmers, 2: subsampling. [1]\n"
	" -R <int>    Realignment bandwidth, 0: disable, [16]\n"
	" -c <int>    Consensus mode: 0, run-length; 1, dp-call-cns, [0]\n"
	" -C <int>    Min count of bases to call a consensus base, [3]\n"
	" -F <float>  Min frequency of non-gap bases to call a consensus base, [0.5]\n"
	" -N <int>    Max number of reads in PO-MSA [20]\n"
	"             Keep in mind that I am not going to generate high accurate consensus sequences here\n"
	" -x <string> Presets, []\n"
	"             sam-sr: polishs contigs from short reads mapping, accepts sorted SAM files\n"
	"                     shorted for '-w 200 -j 150 -R 0 -b 1 -c 1 -N 50 -rS 2'\n"
	" -v          Verbose\n"
	" -V          Print version information and then exit\n"
	"\n", TOSTR(VERSION), TOSTR(RELEASE));
	return 1;
}

int main(int argc, char **argv){
	POGPar par;
	CTGCNS *cc;
	SeqBank *refs;
	FileReader *fr, *db;
	cplist *infs, *dbfs;
	FILE *out;
	char *outf;
	u4i i;
	int reglen, shuffle, winlen, winmin, fail_skip;
	int seqmax, wsize, print_lay, sam_present;
	int c, ncpu, overwrite;
	par = DEFAULT_POG_PAR;
	ncpu = 4;
	shuffle = 1;
	seqmax = 20;
	winlen = 200;
	winmin = 0;
	fail_skip = 1;
	reglen = 1500;
	wsize = 2000; // used for SAM input
	infs = init_cplist(4);
	dbfs = init_cplist(4);
	outf = NULL;
	overwrite = 0;
	print_lay = 0;
	sam_present = 0;
	while((c = getopt(argc, argv, "hvVt:d:rp:ui:o:fj:S:B:W:w:Ab:M:X:I:D:E:H:R:c:C:F:N:x:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'p': print_lay = 1;
			case 'd': push_cplist(dbfs, optarg); break;
			case 'u': sam_present = 1; break;
			case 'r': par.refmode = 1; break;
			case 'i': push_cplist(infs, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'j': reglen = atoi(optarg); break;
			case 'S': shuffle = atoi(optarg); break;
			case 'B': par.W = atoi(optarg); break;
			case 'W': winlen = atoi(optarg); break;
			case 'w': wsize = winmin = atoi(optarg); break;
			case 'A': fail_skip = 0; break;
			case 'b': par.tribase = atoi(optarg); break;
			case 'M': par.M = atoi(optarg); break;
			case 'X': par.X = atoi(optarg); break;
			case 'I': par.I = atoi(optarg); break;
			case 'D': par.D = atoi(optarg); break;
			case 'E': par.E = atoi(optarg); break;
			case 'H': par.H = atof(optarg); break;
			case 'R': par.rW = atoi(optarg); break;
			case 'c': par.cnsmode = atoi(optarg); break;
			case 'C': par.msa_min_cnt = atoi(optarg); break;
			case 'F': par.msa_min_freq = atof(optarg); break;
			case 'N': seqmax = atoi(optarg); break;
			case 'x':
				if(strcasecmp(optarg, "sam-sr") == 0){
					wsize = winmin = 200;
					reglen = 150;
					par.rW = 0;
					par.tribase = 1;
					par.cnsmode = 1;
					seqmax = 50;
					par.refmode = 1;
					shuffle = 2;
				} else {
					fprintf(stderr, "Unknown preset[%s]\n", optarg);
					return 1;
				}
				break;
			case 'v': cns_debug ++; break;
			case 'V': fprintf(stdout, "wtpoa-cns %s\n", TOSTR(VERSION)); return 0;
			default: return usage();
		}
	}
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	//SET_PROC_LIMIT(10 * 1024 * 1024 * 1024ULL, 0); // TODO: remove it after debug
	if(winlen < 0){
		par.W_score = - winlen;
		winlen = 0;
	}
	if(winmin <= 0){
		winmin = winlen * 0.5;
	}
	if(seqmax <= 0 || seqmax >= POG_RDCNT_MAX){
		seqmax = POG_RDCNT_MAX;
	}
	if(ncpu <= 0 && _sig_proc_deamon) ncpu = _sig_proc_deamon->ncpu;
	if(ncpu <= 0){
		fprintf(stderr, " -- Invalid cpu number '%d' in %s -- %s:%d --\n", ncpu, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(outf && !overwrite && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	if(infs->size) fr = open_all_filereader(infs->size, infs->buffer, 1);
	else fr = open_filereader(NULL, 1);
	if(outf){
		out = open_file_for_write(outf, NULL, 1);
	} else out = stdout;
	if(shuffle == 2){
		srand48(13);
	}
	if(dbfs->size == 0){
		WTLAYBlock *wb;
		wb = init_wtlayblock(fr);
		cc = init_ctgcns(wb, iter_wtlayblock, info_wtlayblock, ncpu, shuffle, seqmax, winlen, winmin, fail_skip, reglen, &par);
		cc->print_progress = 100;
		if(print_lay){
			print_lays_ctgcns(cc, out);
		} else {
			while(iter_ctgcns(cc)){
				fprintf(out, ">%s len=%d\n", cc->tag->string, (u4i)cc->cns->size);
				for(i=0;i<cc->cns->size;i+=100){
					if(i + 100 <= cc->cns->size){
						println_seq_basebank(cc->cns, i, 100, out);
					} else {
						println_seq_basebank(cc->cns, i, cc->cns->size - i, out);
					}
				}
				fflush(out);
			}
		}
		free_ctgcns(cc);
		free_wtlayblock(wb);
	} else {
		SAMBlock *sb;
		BioSequence *seq;
		if(reglen > wsize || 2 * reglen < wsize){
			fprintf(stderr, " -- SAM Input mode: -w wsize(%d), -j reglen(%d), MUST has reglen <= wsize and 2 * reglen >= wsize in %s -- %s:%d --\n", wsize, reglen, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			return 1;
		}
		refs = init_seqbank();
		db = open_all_filereader(dbfs->size, dbfs->buffer, 1);
		seq = init_biosequence();
		while(readseq_filereader(db, seq)){
			push_seqbank(refs, seq->tag->string, seq->tag->size, seq->seq->string, seq->seq->size);
		}
		free_biosequence(seq);
		close_filereader(db);
		sb = init_samblock(refs, fr, wsize, reglen, sam_present);
		cc = init_ctgcns(sb, iter_samblock, info_samblock, ncpu, shuffle, seqmax, 0, 0, fail_skip, UInt((wsize - reglen) * 1.2 + 100), &par);
		cc->print_progress = 100;
		if(print_lay){
			print_lays_ctgcns(cc, out);
		} else {
			while(iter_ctgcns(cc)){
				fprintf(out, ">%s len=%d\n", cc->tag->string, (u4i)cc->cns->size);
				for(i=0;i<cc->cns->size;i+=100){
					if(i + 100 <= cc->cns->size){
						println_seq_basebank(cc->cns, i, 100, out);
					} else {
						println_seq_basebank(cc->cns, i, cc->cns->size - i, out);
					}
				}
				fflush(out);
			}
		}
		free_ctgcns(cc);
		free_samblock(sb);
		free_seqbank(refs);
	}
	close_filereader(fr);
	if(outf) fclose(out);
	free_cplist(infs);
	free_cplist(dbfs);
	END_STAT_PROC_INFO(stderr);
	return 0;
}
