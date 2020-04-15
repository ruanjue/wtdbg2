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
	" -u <int>    XORed flags to handle SAM input. [0]\n"
	"             0x1: Only process reference regions present in/between SAM alignments\n"
	"             0x2: Don't fileter secondary/supplementary SAM records with flag (0x100 | 0x800)\n"
	" -r          Force to use reference mode\n"
	" -p <string> Similar with -d, but translate SAM into wtdbg layout file\n"
	" -i <string> Input file(s) *.ctg.lay from wtdbg, +, [STDIN]\n"
	"             Or sorted SAM files when having -d/-p\n"
	" -o <string> Output files, [STDOUT]\n"
	" -e <string> Output the coordinates of orignal edges in consensus sequences, [NULL]\n"
	" -f          Force overwrite\n"
	" -j <int>    Expected max length of node, or say the overlap length of two adjacent units in layout file, [1500] bp\n"
	"             overlap will be default to 1500(or 150 for sam-sr) when having -d/-p, block size will be 2.5 * overlap\n"
	" -b <int>    Bonus for tri-bases match, [0]\n"
	" -M <int>    Match score, [2]\n"
	" -X <int>    Mismatch score, [-5]\n"
	" -I <int>    Insertion score, [-2]\n"
	" -D <int>    Deletion score, [-4]\n"
	" -H <float>  Homopolymer merge score used in dp-call-cns mode, [-3]\n"
	" -B <expr>   Bandwidth in POA, [Wmin[,Wmax[,mat_rate]]], mat_rate = matched_bases/total_bases [64,1024,0.92]\n"
	"             Program will double bandwidth from Wmin to Wmax when mat_rate is lower than setting\n"
	" -W <int>    Window size in the middle of the first read for fast align remaining reads, [200]\n"
	"             If $W is negative, will disable fast align, but use the abs($W) as Band align score cutoff\n"
	" -w <int>    Min size of aligned size in window, [$W * 0.5]\n"
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
	"                     shorted for '-j 50 -W 0 -R 0 -b 1 -c 1 -N 50 -rS 2'\n"
	" -q          Quiet\n"
	" -v          Verbose\n"
	" -V          Print version information and then exit\n"
	"\n", TOSTR(VERSION), TOSTR(RELEASE));
	return 1;
}

int main(int argc, char **argv){
	POGPar par;
	CTGCNS *cc;
	ctg_cns_t *ctg;
	SeqBank *refs;
	FileReader *fr, *db;
	cplist *infs, *dbfs;
	FILE *out, *map;
	char *outf, *mapf;
	u4i i;
	int reglen, shuffle, winlen, winmin, fail_skip;
	int seqmax, wsize, print_lay, flags;
	int c, ncpu, overwrite, quiet;
	par = DEFAULT_POG_PAR;
	ncpu = 4;
	quiet = 0;
	shuffle = 1;
	seqmax = 20;
	winlen = 200;
	winmin = 0;
	fail_skip = 1;
	reglen = 0;
	infs = init_cplist(4);
	dbfs = init_cplist(4);
	outf = NULL;
	mapf = NULL;
	map = NULL;
	overwrite = 0;
	print_lay = 0;
	flags = 0;
	while((c = getopt(argc, argv, "hqvVt:d:rp:u:i:o:e:fj:S:B:W:w:Ab:M:X:I:D:E:H:R:c:C:F:N:x:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'p': print_lay = 1; // fall through
			case 'd': push_cplist(dbfs, optarg); break;
			case 'u': flags = atoi(optarg); break;
			case 'r': par.refmode = 1; break;
			case 'i': push_cplist(infs, optarg); break;
			case 'o': outf = optarg; break;
			case 'e': mapf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'j': reglen = atoi(optarg); break;
			case 'S': shuffle = atoi(optarg); break;
			case 'B':
				{
					char *ptr;
					par.W = strtol(optarg, &ptr, 10);
					if(ptr && ptr[0] == ','){
						par.Wmax = strtol(ptr + 1, &ptr, 10);
						if(ptr && ptr[0] == ','){
							par.W_mat_rate = atof(ptr + 1);
						}
					}
				}
				break;
			case 'W': winlen = atoi(optarg); break;
			case 'w': winmin = atoi(optarg); break;
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
					winlen = 0;
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
			case 'q': quiet = 1; break;
			case 'v': cns_debug ++; break;
			case 'V': fprintf(stdout, "wtpoa-cns %s\n", TOSTR(VERSION)); return 0;
			default: return usage();
		}
	}
	if(quiet){
		int devnull;
		devnull = open("/dev/null", O_WRONLY);
		dup2(devnull, STDERR_FILENO);
	}
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	//SET_PROC_LIMIT(10 * 1024 * 1024 * 1024ULL, 0); // TODO: remove it after debug
	if(reglen == 0){
		reglen = 1500;
	}
	wsize = 2.5 * reglen;
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
	if(mapf){
		map = open_file_for_write(mapf, NULL, 1);
		fprintf(map, "#ctg ctg_off edge edge_full_len edge_off edge_len\n");
	}
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
			while((ctg = iter_ctgcns(cc))){
				fprintf(out, ">%s len=%d\n", ctg->tag->string, (u4i)ctg->cns->size);
				for(i=0;i<ctg->cns->size;i+=100){
					if(i + 100 <= ctg->cns->size){
						println_seq_basebank(ctg->cns, i, 100, out);
					} else {
						println_seq_basebank(ctg->cns, i, ctg->cns->size - i, out);
					}
				}
				fflush(out);
				if(mapf){
					edge_cns_t *edge;
					int coff = 0;
					for(i=0;i<ctg->rs->size;i++){
						edge = ref_edgecnsv(ctg->rs, i);
						fprintf(map, "%s\t%d\tE%u\t%d\t%d\t%d\n", ctg->tag->string, coff, i, edge->slen, edge->beg, edge->end);
						if(edge->end > edge->beg){
							coff += edge->end - edge->beg;
						} else if(coff + edge->end > edge->beg){
							coff += edge->end - edge->beg;
						} else {
							coff = 0;
						}
					}
				}
				if(cc->cycs->size){ // keep only one for reuse
					free_ctg(ctg);
				} else {
					repay_ctgcns(cc, ctg);
				}
			}
		}
		free_ctgcns(cc);
		free_wtlayblock(wb);
	} else {
		SAMBlock *sb;
		BioSequence *seq;
		if(2 * reglen >= wsize){
			fprintf(stderr, " -- SAM Input mode: -w wsize(%d), -j reglen(%d), MUST has 2 * reglen >= wsize in %s -- %s:%d --\n", wsize, reglen, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
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
		sb = init_samblock(refs, fr, wsize, reglen * 2 / 3, flags);
		cc = init_ctgcns(sb, iter_samblock, info_samblock, ncpu, shuffle, seqmax, winlen, winmin, fail_skip, reglen, &par);
		cc->print_progress = 100;
		if(print_lay){
			print_lays_ctgcns(cc, out);
		} else {
			while((ctg = iter_ctgcns(cc))){
				fprintf(out, ">%s len=%d\n", ctg->tag->string, (u4i)ctg->cns->size);
				for(i=0;i<ctg->cns->size;i+=100){
					if(i + 100 <= ctg->cns->size){
						println_seq_basebank(ctg->cns, i, 100, out);
					} else {
						println_seq_basebank(ctg->cns, i, ctg->cns->size - i, out);
					}
				}
				fflush(out);
				if(mapf){
					edge_cns_t *edge;
					int coff = 0;
					for(i=0;i<ctg->rs->size;i++){
						edge = ref_edgecnsv(ctg->rs, i);
						fprintf(map, "%s\t%d\tE%u\t%d\t%d\t%d\n", ctg->tag->string, coff, i, edge->slen, edge->beg, edge->end);
						if(edge->end > edge->beg){
							coff += edge->end - edge->beg;
						} else if(coff + edge->end > edge->beg){
							coff += edge->end - edge->beg;
						} else {
							coff = 0;
						}
					}
				}
				if(cc->cycs->size){ // keep only one for reuse
					free_ctg(ctg);
				} else {
					repay_ctgcns(cc, ctg);
				}
			}
		}
		free_ctgcns(cc);
		free_samblock(sb);
		free_seqbank(refs);
	}
	close_filereader(fr);
	if(outf) fclose(out);
	if(mapf) fclose(map);
	free_cplist(infs);
	free_cplist(dbfs);
	END_STAT_PROC_INFO(stderr);
	return 0;
}
