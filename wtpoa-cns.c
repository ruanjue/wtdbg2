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

#include "tripoa.h"
#include "kswx.h"
#include "filereader.h"

typedef struct {
	u8i idx;
	u4i node1, node2, soff, slen;
	int beg, end;
} edge_cns_t;
define_list(edgecnsv, edge_cns_t);

typedef struct {
	edgecnsv *heap;
	edgecnsv *rs;
	String *tag;
	BaseBank *seq;
	u8i widx;
	u4i cidx, eidx;
	int M, X, O, E, reglen;
	FILE *out;
} CTGCNS;

static inline CTGCNS* init_ctgcns(int M, int X, int O, int E, int reglen, FILE *out){
	CTGCNS *cc;
	cc = malloc(sizeof(CTGCNS));
	cc->heap = init_edgecnsv(64);
	cc->rs   = init_edgecnsv(64);
	cc->tag = init_string(64);
	cc->seq = init_basebank();
	cc->widx = 0;
	cc->cidx = 0;
	cc->eidx = 0;
	cc->M = M;
	cc->X = X;
	cc->O = O;
	cc->E = E;
	cc->reglen = reglen;
	cc->out  = out;
	return cc;
}

static inline void free_ctgcns(CTGCNS *cc){
	free_edgecnsv(cc->heap);
	free_edgecnsv(cc->rs);
	free_string(cc->tag);
	free_basebank(cc->seq);
	free(cc);
}

int revise_joint_point(u32list *cigars, int *qe, int *te){
	u4i i, op, ln, max;
	int qq, q, tt, t;
	q = t = 0;
	qq = tt = 0;
	max = 0;
	for(i=1;i<=cigars->size;i++){
		op = cigars->buffer[cigars->size - i] & 0xF;
		ln = cigars->buffer[cigars->size - i] >> 4;
		if(op == 0){
			if(ln > max){
				qq = q; tt = t;
				max = ln;
			}
			q += ln;
			t += ln;
		} else if(op == 1){
			q += ln;
		} else {
			t += ln;
		}
	}
	*qe -= qq;
	*te -= tt;
	return 1;
}

thread_beg_def(mcns);
CTGCNS *cc;
TriPOG *g;
edge_cns_t edge;
thread_end_def(mcns);

thread_beg_func(mcns);
CTGCNS *cc;
TriPOG *g;
edge_cns_t *edge;
u1v *seq1, *seq2;
kswx_t XX, *xs[2];
u8list *mem_cache;
u32list *cigars[2];
u4i eidx;
int qb, qe, tb, te, b, e;
cc = mcns->cc;
g = mcns->g;
seq1 = init_u1v(1024);
seq2 = init_u1v(1024);
mem_cache = init_u8list(1024);
cigars[0] = init_u32list(64);
cigars[1] = NULL;
xs[0] = &XX;
xs[1] = NULL;
thread_beg_loop(mcns);
if(g->seqs->nseq){
	end_tripog(g);
}
eidx = MAX_U4;
thread_beg_syn(mcns);
if(cc->eidx + 1 < cc->rs->size){
	eidx = cc->eidx;
	edge = ref_edgecnsv(cc->rs, eidx);
	clear_and_inc_u1v(seq1, edge->slen);
	bitseq_basebank(cc->seq, edge->soff, edge->slen, seq1->buffer);
	edge = ref_edgecnsv(cc->rs, eidx + 1);
	clear_and_inc_u1v(seq2, edge->slen);
	bitseq_basebank(cc->seq, edge->soff, edge->slen, seq2->buffer);
	cc->eidx ++;
}
thread_end_syn(mcns);
if(eidx != MAX_U4){
	qb = 0; qe = seq1->size;
	tb = 0; te = seq2->size;
	if(qe > cc->reglen) qb = qe - cc->reglen;
	if(te > cc->reglen) te = cc->reglen;
	kswx_overlap_align_core(xs, cigars, qe - qb, seq1->buffer + qb, te - tb, seq2->buffer + tb, 1, cc->M, cc->X, cc->O, cc->O, cc->E, mem_cache);
	XX.qb += qb;
	XX.qe += qb;
	XX.tb += tb;
	XX.te += tb;
	//if(cns_debug){
		//fprintf(stderr, "#%llu_%s\t%d\t%d\t%d", edge->idx - 1, ctg->string, (int)cseqs->size, XX.qb, XX.qe);
		//fprintf(stderr, "\t%llu_%s\t%d\t%d\t%d", edge->idx ,   ctg->string, (int)edge->len,   XX.tb, XX.te);
		//fprintf(stderr, "\t%d\t%d\t%d\t%d\t%d\n", XX.aln, XX.mat, XX.mis, XX.ins, XX.del);
	//}
	b = XX.qe;
	e = XX.te;
	revise_joint_point(cigars[0], &b, &e);
	thread_beg_syn(mcns);
	if(cns_debug){
		fprintf(stderr, "JOINT\t%llu\tqe = %d -> %d\tte = %d -> %d\n", ref_edgecnsv(cc->rs, eidx)->idx, XX.qe, b, XX.te, e);
		fflush(stderr);
	}
	ref_edgecnsv(cc->rs, eidx    )->end = b;
	ref_edgecnsv(cc->rs, eidx + 1)->beg = e;
	thread_end_syn(mcns);
}
thread_end_loop(mcns);
free_u1v(seq1);
free_u1v(seq2);
free_u8list(mem_cache);
free_u32list(cigars[0]);
thread_end_func(mcns);

int run_cns(FileReader *fr, u4i ncpu, int use_sse, u4i seqmax, int winlen, int winmin, int fail_skip, int W, int M, int X, int I, int D, int E, int rW, int mincnt, float minfreq, int reglen, FILE *out){
	CTGCNS *cc;
	BaseBank *seq;
	edge_cns_t *edge;
	u8i ridx;
	u4i nrun, eidx, meths[2];
	u4i i;
	int c, j, sl, state, next, flag;
	char *ss;
	thread_preprocess(mcns);
	cc = init_ctgcns(M, X, (I + D) / 2, E, reglen, out);
	seq = init_basebank();
	thread_beg_init(mcns, ncpu);
	mcns->cc = cc;
	mcns->g = init_tripog(winlen, winmin, fail_skip, M, X, I, D, W, use_sse, rW, mincnt, minfreq);
	ZEROS(&(mcns->edge));
	thread_end_init(mcns);
	meths[0] = meths[1] = 0;
	clear_string(cc->tag);
	append_string(cc->tag, "anonymous", 9);
	thread_wait_one(mcns);
	cc->widx = 1;
	ridx = 1;
	state = 1;
	next = 0;
	c = 0;
	nrun = flag = 0;
	while(state){
		if(state == 1){
			c = readtable_filereader(fr);
			if(c == -1 || fr->line->string[0] == '>'){
				if(c != -1) cc->cidx ++;
				if(mcns->edge.idx) thread_wake(mcns);
				nrun = ncpu + 1;
				state = 2;
			} else if(fr->line->string[0] == 'E'){
				if(mcns->edge.idx) thread_wake(mcns);
				nrun = 1;
				state = 3;
			} else if(fr->line->string[0] == 'S' || fr->line->string[0] == 's'){
				if(mcns->g->seqs->nseq >= seqmax) continue;
				ss = get_col_str(fr, 5);
				sl = get_col_len(fr, 5);
				//if(UInt(sl) > POG_RDLEN_MAX){
					//sl = POG_RDLEN_MAX;
				//}
				push_tripog(mcns->g, ss, sl);
			}
		} else if(state == 2){
			if(nrun){ // wait for all edge consensus
				thread_wait_next(mcns);
				nrun --;
				state = 4;
				next = 2;
			} else { // end of one contig consensus
				if(cc->rs->size){
					while(1){
						thread_beg_syn(mcns);
						eidx = cc->eidx;
						thread_end_syn(mcns);
						if(eidx + 1 >= cc->rs->size) break;
						thread_wait_one(mcns);
						thread_wake(mcns);
					}
					thread_wait_all(mcns);
					clear_basebank(seq);
					for(i=0;i<cc->rs->size;i++){
						edge = ref_edgecnsv(cc->rs, i);
						if(edge->end > edge->beg){
							fast_fwdbits2basebank(seq, cc->seq->bits, edge->soff + edge->beg, edge->end - edge->beg);
						} else if(Int(seq->size) + edge->end > edge->beg){
							seq->size = seq->size + edge->end - edge->beg;
							normalize_basebank(seq);
						} else {
							clear_basebank(seq);
						}
					}
					fprintf(out, ">%s len=%d\n", cc->tag->string, (u4i)seq->size);
					for(i=0;i<seq->size;i+=100){
						if(i + 100 <= seq->size){
							println_seq_basebank(seq, i, 100, out);
						} else {
							println_seq_basebank(seq, i, seq->size - i, out);
						}
					}
					fflush(out);
				}
				if(c == -1){
					state = 0;
				} else {
					state = 1;
					clear_string(cc->tag);
					ss = get_col_str(fr, 0) + 1;
					sl = get_col_len(fr, 0) - 1;
					for(j=0;j<sl;j++){
						if(ss[j] == ' ') break;
					}
					append_string(cc->tag, ss, j);
					clear_basebank(cc->seq);
					clear_edgecnsv(cc->heap);
					clear_edgecnsv(cc->rs);
					cc->eidx = 0;
				}
			}
		} else if(state == 3){
			if(nrun){
				thread_wait_one(mcns);
				next = 3;
				state = 4;
				nrun --;
			} else {
				beg_tripog(mcns->g);
				mcns->edge.idx = ridx ++;
				mcns->edge.node1 = atoll(get_col_str(fr, 2) + 1);
				mcns->edge.node2 = atoll(get_col_str(fr, 4) + 1);
				mcns->edge.soff = 0;
				mcns->edge.slen = 0;
				mcns->edge.beg = 0;
				mcns->edge.end = 0;
				state = 1;
			}
		} else if(state == 4){ // collect edge consensus
			if(mcns->edge.idx){
				meths[mcns->g->is_tripog] ++;
				if(cns_debug){
					thread_beg_syn(mcns);
					fprintf(stderr, "%llu_%s_N%u_N%u\t%d\t%d\t", mcns->edge.idx, cc->tag->string, mcns->edge.node1, mcns->edge.node2, mcns->g->seqs->nseq, (u4i)mcns->g->cns->size);
					println_seq_basebank(mcns->g->cns, 0, mcns->g->cns->size, stderr);
					thread_end_syn(mcns);
				}
				mcns->edge.soff = cc->seq->size;
				mcns->edge.slen = mcns->g->cns->size;
				mcns->edge.beg = 0;
				mcns->edge.end = mcns->g->cns->size;
				fast_fwdbits2basebank(cc->seq, mcns->g->cns->bits, 0, mcns->g->cns->size);
				array_heap_push(cc->heap->buffer, cc->heap->size, cc->heap->cap, edge_cns_t, mcns->edge, num_cmp(a.idx, b.idx));
				ZEROS(&mcns->edge);
			}
			thread_beg_syn(mcns);
			while(cc->heap->size){
				edge = head_edgecnsv(cc->heap);
				if(edge->idx == cc->widx){
					if(!cns_debug && (cc->widx % 100) == 0){
						fprintf(stderr, "\r%u contigs %llu edges", cc->cidx, cc->widx); fflush(stderr);
					}
					push_edgecnsv(cc->rs, *edge);
					array_heap_remove(cc->heap->buffer, cc->heap->size, cc->heap->cap, edge_cns_t, 0, num_cmp(a.idx, b.idx));
					cc->widx ++;
				} else {
					break;
				}
			}
			thread_end_syn(mcns);
			state = next;
		}
	}
	if(!cns_debug){
		fprintf(stderr, "\r%u contigs %llu edges\n", cc->cidx, cc->widx); fflush(stderr);
	}
	fprintf(stderr, " -- TRIPOA: not=%u yes=%u --\n", meths[0], meths[1]); fflush(stderr);
	thread_beg_close(mcns);
	free_tripog(mcns->g);
	thread_end_close(mcns);
	free_basebank(seq);
	free_ctgcns(cc);
	return 0;
}

int usage(){
	printf(
	"WTPOA-CNS: Consensuser for wtdbg using PO-MSA\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.1\n"
	"Usage: wtpoa-cns [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [4]\n"
	" -i <string> Input file(s) *.ctg.lay from wtdbg, +, [STDIN]\n"
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
	FileReader *fr;
	cplist *infs;
	FILE *out;
	char *outf;
	int reglen, use_sse, bandwidth, rW, winlen, winmin, fail_skip, M, X, I, D, E, mincnt, seqmax;
	float minfreq;
	int c, ncpu, overwrite;
	ncpu = 4;
	use_sse = 2;
	seqmax = 20;
	bandwidth = 96;
	winlen = 200;
	winmin = 0;
	fail_skip = 1;
	reglen = 1500;
	M = 2;
	X = -5;
	I = -2;
	D = -4;
	E = -1;
	rW = 16;
	mincnt = 3;
	minfreq = 0.5;
	infs = init_cplist(4);
	outf = NULL;
	overwrite = 0;
	while((c = getopt(argc, argv, "hvVt:i:o:fj:S:B:W:w:AM:X:I:D:E:R:C:F:N:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'i': push_cplist(infs, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'j': reglen = atoi(optarg); break;
			case 'S': use_sse = atoi(optarg); break;
			case 'B': bandwidth = atoi(optarg); break;
			case 'W': winlen = atoi(optarg); break;
			case 'w': winmin = atoi(optarg); break;
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
			case 'V': fprintf(stdout, "wtpoa-cns 1.1\n"); return 0;
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
	run_cns(fr, ncpu, use_sse, seqmax, winlen, winmin, fail_skip, bandwidth, M, X, I, D, -1, rW, mincnt, minfreq, reglen, out);
	close_filereader(fr);
	if(outf) fclose(out);
	free_cplist(infs);
	END_STAT_PROC_INFO(stderr);
	return 0;
}
