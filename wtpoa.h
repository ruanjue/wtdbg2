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

#ifndef __WTDBG_POA_CNS_RJ_H
#define __WTDBG_POA_CNS_RJ_H

#include "tripoa.h"
#include "kswx.h"
#include "filereader.h"

typedef struct {
	char *reftag;
	u4i reflen, refoff;
	u4i node1, node2;
} lay_blk_t;

typedef struct {
	u4i chridx, bidx; // block idx
	u4i rdidx, rdoff:31, rddir:1;
	char *rdtag;
	BaseBank *seq;
	u4i rbeg, rend;
} lay_seq_t;

static inline void _init_lay_seq_t(lay_seq_t *sq){
	ZEROS(sq);
	sq->seq = init_basebank();
}

static inline void _free_lay_seq_t(lay_seq_t *sq){
	free_basebank(sq->seq);
}

define_recycle_list(layseqr, lay_seq_t, u4i, _init_lay_seq_t(a), _free_lay_seq_t(a));

typedef lay_seq_t* (*iter_cns_block)(void *obj);
typedef void       (*info_cns_block)(void *obj, lay_seq_t *sq, lay_blk_t *bk);

typedef struct {
	u4i cidx, eidx;
	u4i node1, node2, soff, slen;
	int beg, end;
} edge_cns_t;
define_list(edgecnsv, edge_cns_t);

typedef struct {
	u4i cidx;
	edgecnsv *rs;
	String *tag;
	BaseBank *seq, *cns;
	u4i cnt;
} ctg_cns_t;
define_list(ctgcnsv, ctg_cns_t*);

static inline ctg_cns_t* init_ctg(u4i cidx){
	ctg_cns_t *ctg;
	ctg = malloc(sizeof(ctg_cns_t));
	ctg->cidx = cidx;
	ctg->rs  = init_edgecnsv(2);
	ctg->tag = init_string(32);
	ctg->seq = init_basebank(2048);
	ctg->cns = init_basebank(2048);
	ctg->cnt = 0;
	return ctg;
}

static inline void clear_ctg(ctg_cns_t *ctg){
	ctg->cidx = MAX_U4;
	ctg->cnt  = 0;
	clear_edgecnsv(ctg->rs);
	clear_string(ctg->tag);
	clear_basebank(ctg->seq);
	clear_basebank(ctg->cns);
}

static inline void free_ctg(ctg_cns_t *ctg){
	free_edgecnsv(ctg->rs);
	free_string(ctg->tag);
	free_basebank(ctg->seq);
	free_basebank(ctg->cns);
	free(ctg);
}

typedef struct {
	u4i cidx, eidx;
	int type; // 0, none; 1, edge cns; 2, edges aligning; 3, merging into cns
} cns_evt_t;
define_list(cnsevtv, cns_evt_t);

struct CTGCNS;

thread_beg_def(mcns);
struct CTGCNS *cc;
TriPOG *g;
cns_evt_t task;
edge_cns_t edges[2];
u1v *seq1, *seq2;
thread_end_def(mcns);

typedef struct CTGCNS {
	void *obj;
	iter_cns_block itercns;
	info_cns_block infocns;
	lay_seq_t *sq;
	lay_blk_t BK;
	u4i ncpu;
	ctgcnsv *ctgs, *cycs;
	u4i nctg;
	cnsevtv *evts;
	u8i ridx, bases;
	u4i chridx, bidx;
	int state;
	u4i seqmax, inpmax;
	int winlen, winmin, fail_skip;
	int M, X, I, D, E, reglen;
	u8i tri_rets[2];
	u4i print_progress;
	thread_def_shared_vars(mcns);
} CTGCNS;

static inline int revise_joint_point(u32list *cigars, int *qe, int *te){
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

thread_beg_func(mcns);
struct CTGCNS *cc;
TriPOG *g;
ctg_cns_t *ctg;
edge_cns_t *edge;
u1v *seq1, *seq2;
kswx_t XX, *xs[2];
u8list *mem_cache;
u32list *cigars[2];
u4i i;
int qb, qe, tb, te, b, e, ol;
cc = mcns->cc;
g = mcns->g;
seq1 = mcns->seq1;
seq2 = mcns->seq2;
mem_cache = init_u8list(1024);
cigars[0] = init_u32list(64);
cigars[1] = NULL;
xs[0] = &XX;
xs[1] = NULL;
thread_beg_loop(mcns);
if(mcns->task.type == 1){
	if(g->seqs->nseq){
		end_tripog(g);
	}
} else if(mcns->task.type == 2){
	qb = 0; qe = seq1->size;
	tb = 0; te = seq2->size;
	if(qe > cc->reglen) qb = qe - cc->reglen;
	if(te > cc->reglen) te = cc->reglen;
	ol = num_min(qe -qb, te - tb);
	kswx_overlap_align_core(xs, cigars, qe - qb, seq1->buffer + qb, te - tb, seq2->buffer + tb, 1, cc->M, cc->X, (cc->I + cc->D) / 2, (cc->I + cc->D) / 2, cc->E, mem_cache);
	if(XX.aln < Int(0.5 * cc->reglen) || XX.mat < Int(XX.aln * 0.9)){
		// full length alignment
		int maxl;
		maxl = num_min(seq1->size, seq2->size);
		maxl = num_min(maxl, cc->reglen * 4);
		qb = 0; qe = seq1->size;
		tb = 0; te = seq2->size;
		if(qe > maxl) qb = qe - maxl;
		if(te > maxl) te = maxl;
		ol = num_min(qe -qb, te - tb);
		kswx_overlap_align_core(xs, cigars, qe - qb, seq1->buffer + qb, te - tb, seq2->buffer + tb, 1, cc->M, cc->X, (cc->I + cc->D) / 2, (cc->I + cc->D) / 2, cc->E, mem_cache);
	}
	XX.qb += qb;
	XX.qe += qb;
	XX.tb += tb;
	XX.te += tb;
	b = XX.qe;
	e = XX.te;
	revise_joint_point(cigars[0], &b, &e);
	mcns->edges[0].end = b;
	mcns->edges[1].beg = e;
	if(cns_debug){
		thread_beg_syn(mcns);
		fprintf(stderr, "JOINT\tC%u\tE%u\tqe = %d -> %d\tte = %d -> %d\t[%5d,%5d]\t[%5d,%5d,%d,%d,%d]\n", mcns->edges[0].cidx, mcns->edges[0].eidx, XX.qe, b, XX.te, e, (int)seq1->size, (int)seq2->size, XX.aln, XX.mat, XX.mis, XX.ins, XX.del);
		fflush(stderr);
		thread_end_syn(mcns);
	}
} else if(mcns->task.type == 3){
	ctg = get_ctgcnsv(cc->ctgs, mcns->task.cidx);
	clear_basebank(ctg->cns);
	for(i=0;i<ctg->rs->size;i++){
		edge = ref_edgecnsv(ctg->rs, i);
		if(edge->end > edge->beg){
			fast_fwdbits2basebank(ctg->cns, ctg->seq->bits, edge->soff + edge->beg, edge->end - edge->beg);
		} else if(Int(ctg->cns->size) + edge->end > edge->beg){
			ctg->cns->size = ctg->cns->size + edge->end - edge->beg;
			normalize_basebank(ctg->cns);
		} else {
			clear_basebank(ctg->cns);
		}
	}
}
thread_end_loop(mcns);
free_u8list(mem_cache);
free_u32list(cigars[0]);
thread_end_func(mcns);

static inline CTGCNS* init_ctgcns(void *obj, iter_cns_block itercns, info_cns_block infocns, u4i ncpu, int shuffle_rds, u4i seqmax, int winlen, int winmin, int fail_skip, int reglen, POGPar *par){
	CTGCNS *cc;
	thread_prepare(mcns);
	cc = malloc(sizeof(CTGCNS));
	cc->obj = obj;
	cc->itercns = itercns;
	cc->infocns = infocns;
	cc->sq = NULL;
	ZEROS(&cc->BK);
	cc->ncpu = ncpu;
	cc->ctgs = init_ctgcnsv(8);
	cc->cycs = init_ctgcnsv(8);
	cc->evts = init_cnsevtv(64);
	cc->nctg = 0;
	cc->ridx = 0;
	cc->bases = 0;
	cc->state = 1;
	cc->seqmax = seqmax;
	cc->inpmax = seqmax * 5;
	cc->winlen = winlen;
	cc->winmin = winmin;
	cc->fail_skip = fail_skip;
	cc->M = par->M;
	cc->X = par->X;
	cc->I = par->I;
	cc->D = par->D;
	cc->E = par->E;
	cc->reglen = reglen;
	thread_beg_init(mcns, ncpu);
	mcns->cc = cc;
	mcns->g = init_tripog(seqmax, shuffle_rds, winlen, winmin, fail_skip, par);
	ZEROS(&mcns->task);
	mcns->seq1 = init_u1v(1024);
	mcns->seq2 = init_u1v(1024);
	ZEROS(&(mcns->edges[0]));
	ZEROS(&(mcns->edges[1]));
	mcns->edges[0].eidx = MAX_U4;
	mcns->edges[1].eidx = MAX_U4;
	thread_end_init(mcns);
	cc->tri_rets[0] = 0;
	cc->tri_rets[1] = 0;
	cc->print_progress = 0;
	cc->chridx = MAX_U4;
	cc->bidx = MAX_U4;
	thread_export(mcns, cc);
	return cc;
}

static inline void reset_ctgcns(CTGCNS *cc, void *obj, iter_cns_block itercns, info_cns_block infocns){
	ctg_cns_t *ctg;
	u4i i;
	cc->obj = obj;
	cc->itercns = itercns;
	cc->infocns = infocns;
	cc->sq = NULL;
	ZEROS(&cc->BK);
	cc->ridx = 0;
	cc->state = 1;
	cc->chridx = MAX_U4;
	cc->bidx = MAX_U4;
	for(i=0;i<cc->ctgs->size;i++){
		ctg = get_ctgcnsv(cc->ctgs, i);
		if(ctg == NULL) continue;
		free_ctg(ctg);
	}
	clear_ctgcnsv(cc->ctgs);
	cc->nctg = 0;
	clear_cnsevtv(cc->evts);
	cc->bases = 0;
}

static inline void free_ctgcns(CTGCNS *cc){
	ctg_cns_t *ctg;
	u4i i;
	thread_prepare(mcns);
	thread_import(mcns, cc);
	thread_beg_close(mcns);
	free_tripog(mcns->g);
	free_u1v(mcns->seq1);
	free_u1v(mcns->seq2);
	thread_end_close(mcns);
	for(i=0;i<cc->ctgs->size;i++){
		ctg = get_ctgcnsv(cc->ctgs, i);
		if(ctg == NULL) continue;
		free_ctg(ctg);
	}
	free_ctgcnsv(cc->ctgs);
	for(i=0;i<cc->cycs->size;i++){
		ctg = get_ctgcnsv(cc->cycs, i);
		if(ctg == NULL) continue;
		free_ctg(ctg);
	}
	free_ctgcnsv(cc->cycs);
	free_cnsevtv(cc->evts);
	free(cc);
}

static inline void print_lays_ctgcns(CTGCNS *cc, FILE *out){
	lay_blk_t *bk;
	bk = &cc->BK;
	while((cc->sq = cc->itercns(cc->obj))){
		if(cc->sq->chridx != cc->chridx){
			cc->chridx = cc->sq->chridx;
			cc->bidx = MAX_U4;
			cc->infocns(cc->obj, cc->sq, bk);
			fflush(out);
			fprintf(out, ">%s len=%u\n", bk->reftag, bk->reflen);
		}
		if(cc->sq->bidx != cc->bidx){
			cc->bidx = cc->sq->bidx;
			cc->infocns(cc->obj, cc->sq, bk);
			fprintf(out, "E\t%u\tN%u\t+\tN%u\t+\n", bk->refoff, bk->node1, bk->node2);
		}
		{
			if(cc->sq->rdtag){
				fprintf(out, "S\t%s\t%c\t%u\t%u\t", cc->sq->rdtag, "+-"[cc->sq->rddir], cc->sq->rdoff, (u4i)cc->sq->seq->size);
			} else {
				fprintf(out, "S\tR%llu\t%c\t%u\t%u\t", (u8i)cc->sq->rdidx, "+-"[cc->sq->rddir], cc->sq->rdoff, (u4i)cc->sq->seq->size);
			}
			print_seq_basebank(cc->sq->seq, 0, cc->sq->seq->size, out);
			fprintf(out, "\t%u\t%u\n", cc->sq->rbeg, cc->sq->rend);
		}
	}
}

static inline void repay_ctgcns(CTGCNS *cc, ctg_cns_t *ctg){
	clear_ctg(ctg);
	push_ctgcnsv(cc->cycs, ctg);
}

static inline ctg_cns_t* iter_ctgcns(CTGCNS *cc){
	ctg_cns_t *ctg, *ret;
	edge_cns_t *edge;
	cns_evt_t EVT;
	lay_blk_t *bk;
	u8i eidx;
	thread_prepare(mcns);
	thread_import(mcns, cc);
	bk = &cc->BK;
	ret = NULL;
	while(1){
		cc->sq = cc->itercns(cc->obj);
		if(cc->sq){
			if(cc->sq->chridx != cc->chridx){
				if(cc->cycs->size){
					pop_ctgcnsv(cc->cycs, &ctg);
					ctg->cidx = cc->ctgs->size;
				} else {
					ctg = init_ctg(cc->ctgs->size);
				}
				push_ctgcnsv(cc->ctgs, ctg);
				cc->chridx = cc->sq->chridx;
				cc->bidx = MAX_U4;
				cc->infocns(cc->obj, cc->sq, bk);
				if(bk->reftag){
					append_string(ctg->tag, bk->reftag, strlen(bk->reftag));
				} else {
					append_string(ctg->tag, "anonymous", 10);
				}
			}
			if(cc->sq->bidx != cc->bidx){
				cc->ridx ++;
				if(!cns_debug && cc->print_progress && (cc->ridx % cc->print_progress) == 0){
					fprintf(stderr, "\r%u contigs %llu edges %llu bases", (u4i)cc->ctgs->size, cc->ridx, cc->bases); fflush(stderr);
				}
			} else {
				if(mcns->g->seqs->nseq < cc->inpmax){
					fwdbitpush_tripog(mcns->g, cc->sq->seq->bits, 0, cc->sq->seq->size, cc->sq->rbeg, cc->sq->rend);
				}
				continue;
			}
		}
		if(mcns->task.type != 0){
			thread_wake(mcns);
		}
		while(1){
			thread_wait_one(mcns);
			if(mcns->task.type == 2){
				ctg = get_ctgcnsv(cc->ctgs, mcns->task.cidx);
				ctg->rs->buffer[mcns->edges[0].eidx].end = mcns->edges[0].end;
				ctg->rs->buffer[mcns->edges[1].eidx].beg = mcns->edges[1].beg;
				mcns->edges[0].eidx = MAX_U4;
				mcns->edges[1].eidx = MAX_U4;
				mcns->task.type = 0;
				ctg->cnt ++;
				if(ctg->cnt + 1 ==  ctg->rs->size && (ctg->cidx + 1 < cc->ctgs->size || cc->sq == NULL)){
					EVT.cidx = ctg->cidx;
					EVT.eidx = MAX_U4;
					EVT.type = 3;
					array_heap_push(cc->evts->buffer, cc->evts->size, cc->evts->cap, cns_evt_t, EVT, num_cmpxx(b.type, a.type, a.cidx, b.cidx, a.eidx, b.eidx));
				}
			} else if(mcns->task.type == 1){
				ctg = get_ctgcnsv(cc->ctgs, mcns->task.cidx);
				cc->tri_rets[mcns->g->is_tripog] ++;
				if(cns_debug){
					thread_beg_syn(mcns);
					fprintf(stderr, "%u_%s_N%u_N%u\t%d\t%d\t", mcns->task.eidx, ctg->tag->string, mcns->edges[0].node1, mcns->edges[0].node2, mcns->g->seqs->nseq, (u4i)mcns->g->cns->size);
					println_seq_basebank(mcns->g->cns, 0, mcns->g->cns->size, stderr);
					thread_end_syn(mcns);
				}
				mcns->edges[0].soff = ctg->seq->size;
				mcns->edges[0].slen = mcns->g->cns->size;
				mcns->edges[0].beg = 0;
				mcns->edges[0].end = mcns->g->cns->size;
				fast_fwdbits2basebank(ctg->seq, mcns->g->cns->bits, 0, mcns->g->cns->size);
				eidx = mcns->task.eidx;
				ctg->rs->buffer[eidx] = mcns->edges[0];
				mcns->edges[0].eidx = MAX_U4;
				mcns->task.type = 0;
				if(eidx && ctg->rs->buffer[eidx - 1].eidx != MAX_U4){
					EVT.cidx = ctg->cidx;
					EVT.eidx = eidx - 1;
					EVT.type = 2;
					array_heap_push(cc->evts->buffer, cc->evts->size, cc->evts->cap, cns_evt_t, EVT, num_cmpxx(b.type, a.type, a.cidx, b.cidx, a.eidx, b.eidx));
				}
				if(eidx + 1 < ctg->rs->size && ctg->rs->buffer[eidx + 1].eidx != MAX_U4){
					EVT.cidx = ctg->cidx;
					EVT.eidx = eidx;
					EVT.type = 2;
					array_heap_push(cc->evts->buffer, cc->evts->size, cc->evts->cap, cns_evt_t, EVT, num_cmpxx(b.type, a.type, a.cidx, b.cidx, a.eidx, b.eidx));
				}
				if(ctg->rs->size == 1 && (ctg->cidx + 1 < cc->ctgs->size || cc->sq == NULL)){
					EVT.cidx = ctg->cidx;
					EVT.eidx = MAX_U4;
					EVT.type = 3;
					array_heap_push(cc->evts->buffer, cc->evts->size, cc->evts->cap, cns_evt_t, EVT, num_cmpxx(b.type, a.type, a.cidx, b.cidx, a.eidx, b.eidx));
				}
			} else if(mcns->task.type == 3){
				ctg = get_ctgcnsv(cc->ctgs, mcns->task.cidx);
				cc->bases += ctg->cns->size;
				cc->nctg ++;
				ret = ctg;
				set_ctgcnsv(cc->ctgs, mcns->task.cidx, NULL);
				mcns->task.type = 0;
				break;
			}
			if(cc->evts->size){
				EVT = cc->evts->buffer[0];
				array_heap_remove(cc->evts->buffer, cc->evts->size, cc->evts->cap, cns_evt_t, 0, num_cmpxx(b.type, a.type, a.cidx, b.cidx, a.eidx, b.eidx));
				mcns->task = EVT;
				if(EVT.type == 2){
					ctg = get_ctgcnsv(cc->ctgs, EVT.cidx);
					mcns->edges[0] = ctg->rs->buffer[EVT.eidx];
					mcns->edges[0].eidx = EVT.eidx;
					clear_and_inc_u1v(mcns->seq1, mcns->edges[0].slen);
					bitseq_basebank(ctg->seq, mcns->edges[0].soff, mcns->edges[0].slen, mcns->seq1->buffer);
					mcns->edges[1] = ctg->rs->buffer[EVT.eidx + 1];
					mcns->edges[1].eidx = EVT.eidx + 1;
					clear_and_inc_u1v(mcns->seq2, mcns->edges[1].slen);
					bitseq_basebank(ctg->seq, mcns->edges[1].soff, mcns->edges[1].slen, mcns->seq2->buffer);
				}
				thread_wake(mcns);
			} else {
				break;
			}
		}
		if(cc->sq){
			ctg = get_ctgcnsv(cc->ctgs, cc->ctgs->size - 1);
			cc->bidx = cc->sq->bidx;
			cc->infocns(cc->obj, cc->sq, bk);
			mcns->task.cidx = ctg->cidx;
			mcns->task.eidx = ctg->rs->size;
			mcns->task.type = 1;
			mcns->edges[1].eidx = MAX_U4;
			edge = next_ref_edgecnsv(ctg->rs);
			edge->cidx  = ctg->cidx;
			edge->eidx  = MAX_U4;
			edge->node1 = bk->node1;
			edge->node2 = bk->node2;
			edge->soff = 0;
			edge->slen = 0;
			edge->beg = 0;
			edge->end = 0;
			mcns->edges[0] = *edge;
			mcns->edges[0].eidx = offset_edgecnsv(ctg->rs, edge);
			beg_tripog(mcns->g);
			fwdbitpush_tripog(mcns->g, cc->sq->seq->bits, 0, cc->sq->seq->size, cc->sq->rbeg, cc->sq->rend);
		}
		if(ret){
			break;
		} else if(cc->sq == NULL){
			if(cc->nctg == cc->ctgs->size){ // all finished
				if(!cns_debug && cc->print_progress){
					fprintf(stderr, "\r%u contigs %llu edges %llu bases\n", (u4i)cc->ctgs->size, cc->ridx, cc->bases); fflush(stderr);
				}
				break;
			}
		}
	}
	thread_export(mcns, cc);
	return ret;
}

/*
static inline int iter_ctgcns2(CTGCNS *cc){
	edge_cns_t *edge;
	lay_blk_t *bk;
	u8i i, eidx;
	int next, waitall;
	thread_prepare(mcns);
	thread_import(mcns, cc);
	next = 0;
	waitall = 0;
	bk = &cc->BK;
	while(cc->state){
		if(cc->state == 1){
			if(cc->sq == NULL){
				cc->sq = cc->itercns(cc->obj);
			}
			if(cc->sq == NULL || cc->sq->chridx != cc->chridx){
				if(cc->chridx != MAX_U4){
					if(mcns->edges[0].idx != MAX_U8){
						thread_wake(mcns);
						cc->erun ++;
					}
					cc->state = 2;
					next = 4;
					waitall = 1;
				} else if(cc->sq){
					cc->state = 5;
				} else {
					thread_export(mcns, cc);
					return 0;
				}
			} else if(cc->sq->bidx != cc->bidx){
				cc->ridx ++;
				if(mcns->edges[0].idx != MAX_U8){
					thread_wake(mcns);
					cc->erun ++;
				}
				cc->state = 2;
				next = 3;
				waitall = 0;
				if(!cns_debug && cc->print_progress && (cc->ridx % cc->print_progress) == 0){
					fprintf(stderr, "\r%u contigs %llu edges %llu bases", (u4i)cc->ctgs->size, cc->ridx, cc->bases); fflush(stderr);
				}
			} else {
				if(0){
					if(cc->sq->rdtag){
						fprintf(stdout, "S\t%s\t%c\t%u\t%u\t", cc->sq->rdtag, "+-"[cc->sq->rddir], cc->sq->rdoff, (u4i)cc->sq->seq->size);
					} else {
						fprintf(stdout, "S\tR%llu\t%c\t%u\t%u\t", (u8i)cc->sq->rdidx, "+-"[cc->sq->rddir], cc->sq->rdoff, (u4i)cc->sq->seq->size);
					}
					print_seq_basebank(cc->sq->seq, 0, cc->sq->seq->size, stdout);
					fprintf(stdout, "\t%u\t%u\n", cc->sq->rbeg, cc->sq->rend);
				}
				if(mcns->g->seqs->nseq < cc->inpmax){
					fwdbitpush_tripog(mcns->g, cc->sq->seq->bits, 0, cc->sq->seq->size, cc->sq->rbeg, cc->sq->rend);
				}
				cc->sq = NULL;
			}
		} else if(cc->state == 2){
			while(1){
				if(waitall){
					if(cc->tasks->size == 0){
						thread_wait_done(mcns);
					} else {
						thread_wait_one(mcns);
					}
				} else {
					thread_wait_one(mcns);
				}
				if(mcns->edges[1].idx != MAX_U8){
					cc->erun --;
					cc->rs->buffer[mcns->edges[0].idx].end = mcns->edges[0].end;
					cc->rs->buffer[mcns->edges[1].idx].beg = mcns->edges[1].beg;
					mcns->edges[0].idx = MAX_U8;
					mcns->edges[1].idx = MAX_U8;
				} else if(mcns->edges[0].idx != MAX_U8){
					cc->erun --;
					cc->tri_rets[mcns->g->is_tripog] ++;
					if(cns_debug){
						thread_beg_syn(mcns);
						fprintf(stderr, "%llu_%s_N%u_N%u\t%d\t%d\t", mcns->edges[0].idx, cc->tag->string, mcns->edges[0].node1, mcns->edges[0].node2, mcns->g->seqs->nseq, (u4i)mcns->g->cns->size);
						println_seq_basebank(mcns->g->cns, 0, mcns->g->cns->size, stderr);
						thread_end_syn(mcns);
					}
					mcns->edges[0].soff = cc->seq->size;
					mcns->edges[0].slen = mcns->g->cns->size;
					mcns->edges[0].beg = 0;
					mcns->edges[0].end = mcns->g->cns->size;
					fast_fwdbits2basebank(cc->seq, mcns->g->cns->bits, 0, mcns->g->cns->size);
					eidx = mcns->edges[0].idx;
					cc->rs->buffer[eidx] = mcns->edges[0];
					mcns->edges[0].idx = MAX_U8;
					if(eidx && cc->rs->buffer[eidx - 1].idx != MAX_U8){
						push_u8v(cc->tasks, eidx - 1);
					}
					if(eidx + 1 < cc->rs->size && cc->rs->buffer[eidx + 1].idx != MAX_U8){
						push_u8v(cc->tasks, eidx);
					}
				}
				if(pop_u8v(cc->tasks, &eidx)){
					mcns->edges[0] = cc->rs->buffer[eidx];
					mcns->edges[0].idx = eidx;
					clear_and_inc_u1v(mcns->seq1, mcns->edges[0].slen);
					bitseq_basebank(cc->seq, mcns->edges[0].soff, mcns->edges[0].slen, mcns->seq1->buffer);
					mcns->edges[1] = cc->rs->buffer[eidx + 1];
					mcns->edges[1].idx = eidx + 1;
					clear_and_inc_u1v(mcns->seq2, mcns->edges[1].slen);
					bitseq_basebank(cc->seq, mcns->edges[1].soff, mcns->edges[1].slen, mcns->seq2->buffer);
					thread_wake(mcns);
					cc->erun ++;
				} else {
					if(waitall){
						if(cc->erun == 0){
							break;
						}
					} else {
						break;
					}
				}
			}
			cc->state = next;
		} else if(cc->state == 3){
			cc->bidx = cc->sq->bidx;
			cc->infocns(cc->obj, cc->sq, bk);
			mcns->edges[1].idx = MAX_U8;
			edge = next_ref_edgecnsv(cc->rs);
			edge->idx = MAX_U8;
			edge->node1 = bk->node1;
			edge->node2 = bk->node2;
			edge->soff = 0;
			edge->slen = 0;
			edge->beg = 0;
			edge->end = 0;
			mcns->edges[0] = *edge;
			mcns->edges[0].idx = offset_edgecnsv(cc->rs, edge);
			beg_tripog(mcns->g);
			cc->state = 1;
			if(0){
				fprintf(stdout, "E\t%u\tN%u\t+\tN%u\t+\n", bk->refoff, bk->node1, bk->node2);
			}
		} else if(cc->state == 4){
			clear_basebank(cc->cns);
			for(i=0;i<cc->rs->size;i++){
				edge = ref_edgecnsv(cc->rs, i);
				if(edge->end > edge->beg){
					fast_fwdbits2basebank(cc->cns, cc->seq->bits, edge->soff + edge->beg, edge->end - edge->beg);
				} else if(Int(cc->cns->size) + edge->end > edge->beg){
					cc->cns->size = cc->cns->size + edge->end - edge->beg;
					normalize_basebank(cc->cns);
				} else {
					clear_basebank(cc->cns);
				}
			}
			cc->bases += cc->cns->size;
			if(cc->sq == NULL){
				if(!cns_debug && cc->print_progress){
					fprintf(stderr, "\r%u contigs %llu edges %llu bases\n", (u4i)cc->ctgs->size, cc->ridx, cc->bases); fflush(stderr);
				}
				cc->state = 0;
			} else {
				cc->state = 5;
			}
			thread_export(mcns, cc);
			return 1;
		} else if(cc->state == 5){
			cc->chridx = cc->sq->chridx;
			//cc->cidx ++;
			cc->bidx = MAX_U4;
			cc->infocns(cc->obj, cc->sq, bk);
			clear_string(cc->tag);
			if(bk->reftag){
				append_string(cc->tag, bk->reftag, strlen(bk->reftag));
			} else {
				append_string(cc->tag, "anonymous", 10);
			}
			clear_basebank(cc->seq);
			clear_u8v(cc->tasks);
			clear_edgecnsv(cc->rs);
			cc->state = 1;
			if(0){
				fflush(stdout);
				fprintf(stdout, ">%s len=%u\n", bk->reftag, bk->reflen);
			}
		}
	}
	thread_export(mcns, cc);
	return 0;
}
*/

typedef struct {
	u4i     chridx;
	u4v     *chrs;
	SeqBank *refs;
	BitVec  *vsts;
	FileReader *fr;
	layseqr *seqs;
	u4v     *heap;
	u8i     rdidx;
	u4i     cidx, bidx, sidx, lidx;
	u2i     bsize, bstep; // block size, and slide steps
	int     flags; // if flags & 0x1, only polish reference presented in SAM lines. 0x2: Don't filter secondary/supplementary alignments
	u4i     bidx2;
} SAMBlock;

static inline SAMBlock* init_samblock(SeqBank *refs, FileReader *fr, u2i bsize, u2i bovlp, int flags){
	SAMBlock *sb;
	lay_seq_t *stop;
	u2i bstep;
	bstep = bsize - bovlp;
	assert(bstep <= bsize && 2 * bstep >= bsize);
	sb = malloc(sizeof(SAMBlock));
	sb->chridx = 0;
	sb->chrs = init_u4v(32);
	push_u4v(sb->chrs, MAX_U4);
	sb->refs = refs;
	sb->vsts = init_bitvec(sb->refs->nseq);
	sb->fr   = fr;
	sb->seqs = init_layseqr(1024);
	sb->heap = init_u4v(1024);
	sb->rdidx = 0;
	sb->cidx = 1;
	sb->bidx = 0;
	sb->bidx2 = MAX_U4;
	sb->bsize = bsize;
	sb->bstep = bstep;
	sb->sidx = MAX_U4; // last output lay_seq
	sb->lidx = MAX_U4; // last input lay_seq
	sb->flags = flags;
	stop = pop_layseqr(sb->seqs);
	stop->chridx = refs->nseq + 1;
	stop->bidx   = 0;
	return sb;
}

static inline void free_samblock(SAMBlock *sb){
	free_u4v(sb->chrs);
	free_bitvec(sb->vsts);
	free_layseqr(sb->seqs);
	free_u4v(sb->heap);
	free(sb);
}

static inline lay_seq_t* _push_padding_ref_samblock(SAMBlock *sb, lay_seq_t *sq){
	lay_seq_t *st;
	u8i off;
	u4i sqidx, len, i;
	if(sb->cidx > sb->refs->nseq) return sq;
	sqidx = offset_layseqr(sb->seqs, sq);
	if((sb->flags & 0x1) && sb->bidx2 == MAX_U4){
		sb->bidx = sq->bidx;
	}
	sb->bidx2 = sb->bidx;
	while(sb->cidx < sq->chridx || (sb->cidx == sq->chridx && sb->bidx <= sq->bidx)){
		st = pop_layseqr(sb->seqs);
		st->chridx = sb->cidx;
		st->bidx  = sb->bidx;
		st->rdidx = 0;
		st->rddir = 0;
		st->rdoff = st->bidx * sb->bstep;
		st->rbeg  = st->rdoff;
		if(sb->chrs->size <= sb->cidx){
			for(i=0;i<sb->refs->nseq;i++){
				if(get_bitvec(sb->vsts, i) == 0) break;
			}
			if(i == sb->refs->nseq) break;
			push_u4v(sb->chrs, i);
			one_bitvec(sb->vsts, i);
		}
		off = sb->refs->rdoffs->buffer[sb->chrs->buffer[sb->cidx]];
		len = sb->refs->rdlens->buffer[sb->chrs->buffer[sb->cidx]];
		st->rend = num_min(st->rbeg + sb->bsize, len);
		clear_basebank(st->seq);
		fast_fwdbits2basebank(st->seq, sb->refs->rdseqs->bits, off + st->rdoff, st->rend - st->rbeg);
		array_heap_push(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, offset_layseqr(sb->seqs, st),
			num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
				ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
		if(st->rend >= len){
			sb->cidx ++;
			sb->bidx = 0;
			sb->bidx2 = sb->bidx;
			if(sb->cidx >= sb->refs->nseq) break;
		} else {
			sb->bidx ++;
		}
		sq = ref_layseqr(sb->seqs, sqidx);
	}
	sq = ref_layseqr(sb->seqs, sqidx);
	return sq;
}

static inline lay_seq_t* iter_samblock(void *obj){
	SAMBlock *sb;
	lay_seq_t *sc, *sl;
	u4i chr, chridx, scidx, off, minlen, rddir, rdoff, nxt, val, len, op;
	char *ptr, *str;
	int c, samflag;
	sb = (SAMBlock*)obj;
	if(sb->sidx != MAX_U4){
		recyc_layseqr(sb->seqs, sb->sidx);
		sb->sidx = MAX_U4;
	}
	do {
		if(sb->heap->size == 0){
			sb->lidx = MAX_U4;
		}
		while(sb->lidx == MAX_U4){
			c = readtable_filereader(sb->fr);
			if(c == -1){
				if(sb->flags & 0x1){
				} else {
					sc = ref_layseqr(sb->seqs, 0);
					sc = _push_padding_ref_samblock(sb, sc);
				}
				break;
			}
			if(sb->fr->line->string[0] == '@') continue;
			if(c < 11) continue;
			samflag = atoi(get_col_str(sb->fr, 1));
			if(samflag & 0x004) continue;
			//if(get_col_str(sb->fr, 9)[0] == '*') continue;
			if(!(sb->flags & 0x2)){
				if(samflag & (0x100 | 0x800)) continue; // filter secondary/supplementary alignments
			}
			// chr
			if((chr = getval_cuhash(sb->refs->rdhash, get_col_str(sb->fr, 2))) == MAX_U4){
				continue;
			}
			if(chr != sb->chrs->buffer[sb->chridx]){
				chridx = ++ sb->chridx;
				push_u4v(sb->chrs, chr);
				one_bitvec(sb->vsts, chr);
			} else {
				chridx = sb->chridx;
			}
			off = atol(get_col_str(sb->fr, 3)) - 1;
			ptr = get_col_str(sb->fr, 5);
			if((*ptr) == '*') continue;
			sb->rdidx ++;
			rdoff = 0;
			minlen = num_min(get_col_len(sb->fr, 9), sb->bsize) / 2;
			rddir = (atol(get_col_str(sb->fr, 1)) & 0x10) >> 4;
			str = get_col_str(sb->fr, 9);
			{
				encap_layseqr(sb->seqs, 2);
				sc = pop_layseqr(sb->seqs);
				sc->chridx = chridx;
				sc->bidx = (off / sb->bstep);
				sc->rdidx = sb->rdidx;
				sc->rddir = rddir;
				sc->rdoff = rdoff;
				sc->rbeg = off;
				sc->rend = 0;
				clear_basebank(sc->seq);
			}
			sl = NULL;
			// parse SAM cigar
			{
				nxt = (off / sb->bstep) * sb->bstep;
				if(nxt && off - nxt < UInt(sb->bsize - sb->bstep)){
					if(sc->bidx){
						scidx = offset_layseqr(sb->seqs, sc);
						sl = pop_layseqr(sb->seqs);
						sc = ref_layseqr(sb->seqs, scidx);
						sl->chridx = chridx;
						sl->bidx = sc->bidx - 1;
						sl->rdidx = sb->rdidx;
						sl->rddir = rddir;
						sl->rdoff = rdoff;
						sl->rbeg = off;
						sl->rend = 0;
						clear_basebank(sl->seq);
					}
					nxt += sb->bsize - sb->bstep;
				} else {
					nxt += sb->bstep;
				}
			}
			val = 0;
			while(*ptr){
				op = MAX_U4;
				switch(*ptr){
					case '0' ... '9': val = val * 10 + (*ptr) - '0'; break;
					case 'X':
					case '=': 
					case 'M': op = 0b111; break;
					case 'I': op = 0b110; break;
					case 'D': op = 0b101; break;
					case 'N': op = 0b001; break;
					case 'S': op = 0b010; break;
					case 'H':
					case 'P': op = 0b000; break;
					default:
						fprintf(stderr, " -- %s\n", get_col_str(sb->fr, 5));
						fprintf(stderr, " -- Bad cigar %d '%c' in %s -- %s:%d --\n", (int)(ptr - get_col_str(sb->fr, 5)), *ptr, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
				}
				ptr ++;
				if(op == MAX_U4) continue;
				if(val == 0) val = 1;
				while(val){
					if(op & 0b001){
						len = num_min(val, nxt - off);
						off += len;
					} else {
						len = val;
					}
					val -= len;
					if(op & 0b010){
						rdoff += len;
						if(op & 0b100){
							if(sc){
								if(sc->seq->size == 0){
									sc->rbeg = (op & 0b001)? off - len : off;
								}
								seq2basebank(sc->seq, str, len);
								sc->rend = off;
							}
							if(sl){
								if(sl->seq->size == 0){
									sl->rbeg = (op & 0b001)? off - len : off;
								}
								seq2basebank(sl->seq, str, len);
								sl->rend = off;
							}
						}
						str += len;
					}
					if(off == nxt){
						if(sl){
							if(sl->rend >= sl->rbeg + minlen){
								scidx = offset_layseqr(sb->seqs, sc);
								sl = _push_padding_ref_samblock(sb, sl);
								sc = ref_layseqr(sb->seqs, scidx);
								if(sb->lidx == MAX_U4){
									sb->lidx = offset_layseqr(sb->seqs, sl);
								}
								array_heap_push(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, offset_layseqr(sb->seqs, sl),
									num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
										ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
							} else {
								push_layseqr(sb->seqs, sl);
							}
							sl = NULL;
							nxt += 2 * sb->bstep - sb->bsize;
						} else {
							scidx = offset_layseqr(sb->seqs, sc);
							encap_layseqr(sb->seqs, 1);
							sl = ref_layseqr(sb->seqs, scidx);
							sc = pop_layseqr(sb->seqs);
							sc->chridx = chridx;
							sc->bidx = sl->bidx + 1;
							sc->rdidx = sb->rdidx;
							sc->rddir = rddir;
							sc->rdoff = rdoff;
							sc->rbeg = off;
							sc->rend = 0;
							clear_basebank(sc->seq);
							nxt += sb->bsize - sb->bstep;
						}
					}
				}
			}
			if(sl){
				if(sl->rend >= sl->rbeg + minlen){
					u4i scidx;
					scidx = sc? offset_layseqr(sb->seqs, sc) : MAX_U4;
					sl = _push_padding_ref_samblock(sb, sl);
					sc = scidx == MAX_U4? NULL : ref_layseqr(sb->seqs, scidx);
					if(sb->lidx == MAX_U4){
						sb->lidx = offset_layseqr(sb->seqs, sl);
					}
					array_heap_push(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, offset_layseqr(sb->seqs, sl),
						num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
							ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
				} else {
					push_layseqr(sb->seqs, sl);
				}
			}
			if(sc){
				if(sc->rend >= sc->rbeg + minlen){
					scidx = sl? offset_layseqr(sb->seqs, sl) : MAX_U4;
					sc = _push_padding_ref_samblock(sb, sc);
					sl = scidx == MAX_U4? NULL : ref_layseqr(sb->seqs, scidx);
					if(sb->lidx == MAX_U4){
						sb->lidx = offset_layseqr(sb->seqs, sc);
					}
					array_heap_push(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, offset_layseqr(sb->seqs, sc),
						num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
							ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
				} else {
					push_layseqr(sb->seqs, sc);
				}
			}
		}
		if(sb->heap->size == 0) break;
		if(sb->lidx != MAX_U4){
			sl = ref_layseqr(sb->seqs, sb->lidx);
			sb->sidx = sb->heap->buffer[0];
			sc = ref_layseqr(sb->seqs, sb->sidx);
			if(sc->chridx == sl->chridx && sc->bidx + 1 >= sl->bidx){
				sb->lidx = MAX_U4;
				continue;
			}
			array_heap_remove(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, 0,
				num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
					ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
			sc->rbeg -= sc->bidx * sb->bstep;
			sc->rend -= sc->bidx * sb->bstep;
		} else {
			sb->sidx = sb->heap->buffer[0];
			array_heap_remove(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, 0,
				num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
					ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
			sc = ref_layseqr(sb->seqs, sb->sidx);
			sc->rbeg -= sc->bidx * sb->bstep;
			sc->rend -= sc->bidx * sb->bstep;
		}
		return sc;
	} while(1);
	sb->sidx = MAX_U4;
	return NULL;
}

static inline void info_samblock(void *obj, lay_seq_t *sq, lay_blk_t *bk){
	SAMBlock *sb;
	sb = (SAMBlock*)obj;
	if(sq == NULL) return;
	bk->node1 = sq->bidx;
	bk->node2 = sq->bidx + 1;
	bk->reftag = sb->refs->rdtags->buffer[sb->chrs->buffer[sq->chridx]];
	bk->reflen = sb->refs->rdlens->buffer[sb->chrs->buffer[sq->chridx]];
	bk->refoff = sq->bidx * sb->bstep;
}

typedef struct {
	String *tag;
	FileReader *fr;
	u4i chridx, bidx;
	u8i rdidx;
	lay_seq_t *key;
	lay_blk_t *blk;
} WTLAYBlock;

static inline WTLAYBlock* init_wtlayblock(FileReader *fr){
	WTLAYBlock *wb;
	wb = malloc(sizeof(WTLAYBlock));
	wb->tag = init_string(32);
	append_string(wb->tag, "annonymous", 10);
	wb->fr = fr;
	wb->chridx = 0;
	wb->bidx = 0;
	wb->rdidx = 0;
	wb->key = calloc(1, sizeof(lay_seq_t));
	wb->key->seq = init_basebank();
	wb->blk = calloc(1, sizeof(lay_blk_t));
	return wb;
}

static inline void free_wtlayblock(WTLAYBlock *wb){
	free_string(wb->tag);
	free_basebank(wb->key->seq);
	free(wb->key);
	free(wb->blk);
	free(wb);
}

static inline lay_seq_t* iter_wtlayblock(void *obj){
	WTLAYBlock *wb;
	char *ss;
	int c, sl;
	wb = (WTLAYBlock*)obj;
	clear_basebank(wb->key->seq);
	while((c = readtable_filereader(wb->fr)) != -1){
		switch(wb->fr->line->string[0]){
			case '>':
				wb->chridx ++;
				wb->bidx = 0;
				clear_string(wb->tag);
				ss = get_col_str(wb->fr, 0) + 1;
				for(sl=0;ss[sl] && ss[sl] != ' ' && ss[sl] != '\t';sl++){
				}
				append_string(wb->tag, ss, sl);
				wb->blk->node1 = 0;
				wb->blk->node2 = 0;
				ss = strstr(ss + sl, "len=");
				if(ss){
					wb->blk->reflen = atol(ss);
				} else {
					wb->blk->reflen = 0;
				}
				wb->blk->reftag = wb->tag->string;
				break;
			case 'E':
				wb->bidx ++;
				wb->blk->refoff = atoll(get_col_str(wb->fr, 1));
				wb->blk->node1 = atoll(get_col_str(wb->fr, 2) + 1);
				wb->blk->node2 = atoll(get_col_str(wb->fr, 4) + 1);
				wb->blk->reftag = wb->tag->string;
				break;
			case 'S':
			case 's':
				ss = get_col_str(wb->fr, 5);
				sl = get_col_len(wb->fr, 5);
				wb->key->chridx = wb->chridx;
				wb->key->bidx   = wb->bidx;
				wb->key->rdidx  = wb->rdidx ++;
				wb->key->rddir  = (get_col_str(wb->fr, 2)[0] == '-');
				wb->key->rdoff  = atoi(get_col_str(wb->fr, 3));
				if(c >= 8){
					wb->key->rbeg = atoi(get_col_str(wb->fr, 6));
					wb->key->rend = atoi(get_col_str(wb->fr, 7));
				} else {
					wb->key->rbeg = 0;
					wb->key->rend = 0;
				}
				seq2basebank(wb->key->seq, ss, sl);
				return wb->key;
		}
	}
	return NULL;
}

static inline void info_wtlayblock(void *obj, lay_seq_t *sq, lay_blk_t *bk){
	if(sq == NULL) return;
	memcpy(bk, ((WTLAYBlock*)obj)->blk, sizeof(lay_blk_t));
}

#endif
