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

int usage(){
	printf(
	"WTPOA-CNS: Consensuser for wtdbg using PO-MSA\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.2\n"
	"Usage: wtpoa-cns [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [4]\n"
	" -d <string> Reference sequences for SAM input, will invoke sorted-SAM input mode and auto set '-j 100 -W 0 -w 1000'\n"
	" -u          Only process reference regions present in/between SAM alignments\n"
	" -r          Force to use reference mode\n"
	" -p <string> Similar with -d, but translate SAM into wtdbg layout file\n"
	" -i <string> Input file(s) *.ctg.lay from wtdbg, +, [STDIN]\n"
	"             Or sorted SAM files when having -d\n"
	" -o <string> Output files, [STDOUT]\n"
	" -f          Force overwrite\n"
	" -j <int>    Expected max length of node, or say the overlap length of two adjacent units in layout file, [1500] bp\n"
	" -M <int>    Match score, [2]\n"
	" -X <int>    Mismatch score, [-5]\n"
	" -I <int>    Insertion score, [-2]\n"
	" -D <int>    Deletion score, [-4]\n"
	" -B <int>    Bandwidth, [96]\n"
	" -W <int>    Window size in the middle of the first read for fast align remaining reads, [200]\n"
	"             If $W is negative, will disable fast align, but use the abs($W) as Band align score cutoff\n"
	" -w <int>    Min size of aligned size in window, [$W * 0.5]\n"
	"             In sorted-SAM input mode, -w is the sliding window size [2000]\n"
	" -A          Abort TriPOA when any read cannot be fast aligned, then try POA\n"
	//" -S <int>    SSE speedup in graph alignment, [1]\n"
	" -R <int>    Realignment bandwidth, 0: disable, [16]\n"
	" -C <int>    Min count of bases to call a consensus base, [3]\n"
	" -F <float>  Min frequency of non-gap bases to call a consensus base, [0.5]\n"
	" -N <int>    Max number of reads in PO-MSA [20]\n"
	"             Keep in mind that I am not going to generate high accurate consensus sequences here\n"
	" -v          Verbose\n"
	" -V          Print version information and then exit\n"
	"\n");
	return 1;
}

int main(int argc, char **argv){
	CTGCNS *cc;
	SeqBank *refs;
	FileReader *fr, *db;
	cplist *infs, *dbfs;
	FILE *out;
	char *outf;
	u4i i;
	int reglen, use_sse, refmode, bandwidth, rW, winlen, winmin, fail_skip, M, X, I, D, E, mincnt, seqmax, wsize, print_lay, sam_present;
	float minfreq;
	int c, ncpu, overwrite;
	ncpu = 4;
	use_sse = 2;
	refmode = 0;
	seqmax = 20;
	bandwidth = 96;
	winlen = 200;
	winmin = 0;
	fail_skip = 1;
	reglen = 1500;
	wsize = 2000; // used for SAM input
	M = 2;
	X = -5;
	I = -2;
	D = -4;
	E = -1;
	rW = 16;
	mincnt = 3;
	minfreq = 0.5;
	infs = init_cplist(4);
	dbfs = init_cplist(4);
	outf = NULL;
	overwrite = 0;
	print_lay = 0;
	sam_present = 0;
	while((c = getopt(argc, argv, "hvVt:d:rp:ui:o:fj:S:B:W:w:AM:X:I:D:E:R:C:F:N:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'p': print_lay = 1;
			case 'd': push_cplist(dbfs, optarg); reglen = 100; winlen = 0; winmin = 1000; break;
			case 'u': sam_present = 1; break;
			case 'r': refmode = 1; break;
			case 'i': push_cplist(infs, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'j': reglen = atoi(optarg); break;
			case 'S': use_sse = atoi(optarg); break;
			case 'B': bandwidth = atoi(optarg); break;
			case 'W': winlen = atoi(optarg); break;
			case 'w': wsize = winmin = atoi(optarg); break;
			case 'A': fail_skip = 0; break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'I': I = atoi(optarg); break;
			case 'D': D = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'R': rW = atoi(optarg); break;
			case 'C': mincnt = atoi(optarg); break;
			case 'F': minfreq = atof(optarg); break;
			case 'N': seqmax = atoi(optarg); break;
			case 'v': cns_debug ++; break;
			case 'V': fprintf(stdout, "wtpoa-cns 1.2\n"); return 0;
			default: return usage();
		}
	}
	BEG_STAT_PROC_INFO(stderr, argc, argv);
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
	if(dbfs->size == 0){
		WTLAYBlock *wb;
		wb = init_wtlayblock(fr);
		cc = init_ctgcns(wb, iter_wtlayblock, info_wtlayblock, ncpu, refmode, seqmax, winlen, winmin, fail_skip, bandwidth, M, X, I, D, -1, rW, mincnt, minfreq, reglen);
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
		refs = init_seqbank();
		db = open_all_filereader(dbfs->size, dbfs->buffer, 1);
		seq = init_biosequence();
		while(readseq_filereader(db, seq)){
			push_seqbank(refs, seq->tag->string, seq->tag->size, seq->seq->string, seq->seq->size);
		}
		free_biosequence(seq);
		close_filereader(db);
		sb = init_samblock(refs, fr, wsize, reglen, sam_present);
		cc = init_ctgcns(sb, iter_samblock, info_samblock, ncpu, 1, seqmax, 0, 0, fail_skip, bandwidth, M, X, I, D, -1, rW, mincnt, minfreq, UInt((wsize - reglen) * 1.2 + 100));
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
