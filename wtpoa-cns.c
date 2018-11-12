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
		fprintf(stderr, "JOINT\t%llu\tqe = %d -> %d\tte = %d -> %d\t[%5d,%5d]\t[%5d,%5d,%d,%d,%d]\n", ref_edgecnsv(cc->rs, eidx)->idx, XX.qe, b, XX.te, e, (int)seq1->size, (int)seq2->size, XX.aln, XX.mat, XX.mis, XX.ins, XX.del);
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

typedef struct {
	char *reftag;
	u4i reflen, refoff;
	u4i node1, node2;
} lay_blk_t;

typedef struct {
	u4i chridx, bidx; // block idx
	u4i rdidx, rdoff;
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

typedef struct {
	u4i     chridx;
	u4v     *chrs;
	SeqBank *refs;
	FileReader *fr;
	layseqr *seqs;
	u4v     *heap;
	u8i     rdidx;
	u4i     cidx, bidx, sidx, lidx;
	u2i     bsize, bstep; // block size, and slide steps
	int     sam_present; // if 1, only polish reference presented in SAM lines
	u4i     bidx2;
} SAMBlock;

static inline SAMBlock* init_samblock(SeqBank *refs, FileReader *fr, u2i bsize, u2i bstep, int sam_present){
	SAMBlock *sb;
	lay_seq_t *stop;
	assert(bstep <= bsize && 2 * bstep >= bsize);
	sb = malloc(sizeof(SAMBlock));
	sb->chridx = 0;
	sb->chrs = init_u4v(32);
	push_u4v(sb->chrs, MAX_U4);
	sb->refs = refs;
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
	sb->sam_present = sam_present;
	stop = pop_layseqr(sb->seqs);
	stop->chridx = refs->nseq + 1;
	stop->bidx   = 0;
	return sb;
}

static inline void free_samblock(SAMBlock *sb){
	free_u4v(sb->chrs);
	free_layseqr(sb->seqs);
	free_u4v(sb->heap);
	free(sb);
}

static inline lay_seq_t* _push_padding_ref_samblock(SAMBlock *sb, lay_seq_t *sq){
	lay_seq_t *st;
	u8i off;
	u4i sqidx, len, i, j;
	if(sb->cidx > sb->refs->nseq) return sq;
	sqidx = offset_layseqr(sb->seqs, sq);
	if(sb->sam_present && sb->bidx2 == MAX_U4){
		sb->bidx = sq->bidx;
	}
	sb->bidx2 = sb->bidx;
	while(sb->cidx < sq->chridx || (sb->cidx == sq->chridx && sb->bidx <= sq->bidx)){
		st = pop_layseqr(sb->seqs);
		st->chridx = sb->cidx;
		st->bidx = sb->bidx;
		st->rdidx = 0;
		st->rdoff = st->bidx * sb->bstep;
		st->rbeg = st->rdoff;
		if(sb->chrs->size <= sb->cidx){
			for(i=0;i<sb->refs->nseq;i++){
				for(j=1;j<sb->chrs->size;j++){
					if(sb->chrs->buffer[j] == i) break;
				}
				if(j == sb->chrs->size) break;
			}
			if(i == sb->refs->nseq) break;
			push_u4v(sb->chrs, i);
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
	u4i chr, chridx, off, rdoff, nxt, val, len, op;
	char *ptr, *str;
	int c;
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
				if(sb->sam_present){
				} else {
					sc = ref_layseqr(sb->seqs, 0);
					sc = _push_padding_ref_samblock(sb, sc);
				}
				break;
			}
			if(sb->fr->line->string[0] == '@') continue;
			if(c < 11) continue;
			if(get_col_str(sb->fr, 9)[0] == '*') continue;
			// chr
			if((chr = getval_cuhash(sb->refs->rdhash, get_col_str(sb->fr, 2))) == MAX_U4){
				continue;
			}
			if(chr != sb->chrs->buffer[sb->chridx]){
				chridx = ++ sb->chridx;
				push_u4v(sb->chrs, chr);
			} else {
				chridx = sb->chridx;
			}
			off = atol(get_col_str(sb->fr, 3));
			ptr = get_col_str(sb->fr, 5);
			if((*ptr) == '*') continue;
			sb->rdidx ++;
			rdoff = 0;
			str = get_col_str(sb->fr, 9);
			{
				encap_layseqr(sb->seqs, 2);
				sc = pop_layseqr(sb->seqs);
				sc->chridx = chridx;
				sc->bidx = (off / sb->bstep);
				sc->rdidx = sb->rdidx;
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
						sl = pop_layseqr(sb->seqs);
						sl->chridx = chridx;
						sl->bidx = sc->bidx - 1;
						sl->rdidx = sb->rdidx;
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
							if(sl->rend > sl->rbeg){
								sl = _push_padding_ref_samblock(sb, sl);
								if(sb->lidx == MAX_U4){
									sb->lidx = offset_layseqr(sb->seqs, sl);
								}
								array_heap_push(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, offset_layseqr(sb->seqs, sl),
									num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
										ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
							}
							sl = NULL;
							nxt += 2 * sb->bstep - sb->bsize;
						} else {
							u4i scidx;
							scidx = offset_layseqr(sb->seqs, sc);
							encap_layseqr(sb->seqs, 1);
							sl = ref_layseqr(sb->seqs, scidx);
							sc = pop_layseqr(sb->seqs);
							sc->chridx = chridx;
							sc->bidx = sl->bidx + 1;
							sc->rdidx = sb->rdidx;
							sc->rdoff = rdoff;
							sc->rbeg = off;
							sc->rend = 0;
							clear_basebank(sc->seq);
							nxt += sb->bsize - sb->bstep;
						}
					}
				}
			}
			if(sl && sl->rend > sl->rbeg){
				sl = _push_padding_ref_samblock(sb, sl);
				if(sb->lidx == MAX_U4){
					sb->lidx = offset_layseqr(sb->seqs, sl);
				}
				array_heap_push(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, offset_layseqr(sb->seqs, sl),
					num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
						ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
			}
			if(sc && sc->rend > sc->rbeg){
				sc = _push_padding_ref_samblock(sb, sc);
				if(sb->lidx == MAX_U4){
					sb->lidx = offset_layseqr(sb->seqs, sc);
				}
				array_heap_push(sb->heap->buffer, sb->heap->size, sb->heap->cap, u4i, offset_layseqr(sb->seqs, sc),
					num_cmpxx(ref_layseqr(sb->seqs, a)->chridx, ref_layseqr(sb->seqs, b)->chridx, ref_layseqr(sb->seqs, a)->bidx,
						ref_layseqr(sb->seqs, b)->bidx, ref_layseqr(sb->seqs, a)->rdidx, ref_layseqr(sb->seqs, b)->rdidx));
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
				wb->key->rdoff  = 0;
				wb->key->rbeg   = 0;
				wb->key->rend   = 0;
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

typedef lay_seq_t* (*iter_cns_block)(void *obj);
typedef void       (*info_cns_block)(void *obj, lay_seq_t *sq, lay_blk_t *bk);

int run_cns(void *obj, iter_cns_block itercns, info_cns_block infocns, u4i ncpu, int refmode, u4i seqmax, int winlen, int winmin, int fail_skip, int W, int M, int X, int I, int D, int E, int rW, int mincnt, float minfreq, int reglen, int print_lay, FILE *out){
	CTGCNS *cc;
	BaseBank *seq;
	edge_cns_t *edge;
	lay_seq_t *sq;
	lay_blk_t *bk;
	u8i ridx;
	u4i nrun, eidx, meths[2];
	u4i i, chridx, bidx;
	int state, next, flag;
	thread_preprocess(mcns);
	cc = init_ctgcns(M, X, (I + D) / 2, E, reglen, out);
	seq = init_basebank();
	thread_beg_init(mcns, ncpu);
	mcns->cc = cc;
	mcns->g = init_tripog(seqmax, refmode, winlen, winmin, fail_skip, M, X, I, D, W, 2, rW, mincnt, minfreq);
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
	nrun = flag = 0;
	sq = NULL;
	chridx = MAX_U4;
	bidx = MAX_U4;
	bk = alloca(sizeof(lay_blk_t));
	if(print_lay){
		while((sq = itercns(obj))){
			if(sq->chridx != chridx){
				chridx = sq->chridx;
				bidx = MAX_U4;
				infocns(obj, sq, bk);
				fflush(out);
				fprintf(out, ">%s len=%u\n", bk->reftag, bk->reflen);
			}
			if(sq->bidx != bidx){
				bidx = sq->bidx;
				infocns(obj, sq, bk);
				fprintf(out, "E\t%u\tN%u\t+\tN%u\t+\n", bk->refoff, bk->node1, bk->node2);
			}
			{
				fprintf(out, "S\tR%llu\t+\t%u\t%u\t", (u8i)sq->rdidx, sq->rdoff, (u4i)sq->seq->size);
				print_seq_basebank(sq->seq, 0, sq->seq->size, out);
				fprintf(out, "\t%u\t%u\n", sq->rbeg, sq->rend);
			}
		}
	} else {
		while(state){
			if(state == 1){
				if(sq == NULL){
					sq = itercns(obj);
				}
				if(sq == NULL || sq->chridx != chridx){
					if(sq){
						chridx = sq->chridx;
						cc->cidx ++;
						bidx = MAX_U4;
						infocns(obj, sq, bk);
					}
					if(mcns->edge.idx) thread_wake(mcns);
					nrun = ncpu + 1;
					state = 2;
				} else if(sq->bidx != bidx){
					if(mcns->edge.idx) thread_wake(mcns);
					bidx = sq->bidx;
					infocns(obj, sq, bk);
					nrun = 1;
					state = 3;
				} else {
					fwdbitpush_tripog(mcns->g, sq->seq->bits, 0, sq->seq->size, sq->rbeg, sq->rend);
					sq = NULL;
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
					if(sq == NULL){
						state = 0;
					} else {
						state = 1;
						clear_string(cc->tag);
						if(bk->reftag){
							append_string(cc->tag, bk->reftag, strlen(bk->reftag));
						} else {
							append_string(cc->tag, "anonymous", 10);
						}
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
					mcns->edge.node1 = bk->node1;
					mcns->edge.node2 = bk->node2;
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
	}
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
	"Version: 1.2\n"
	"Usage: wtpoa-cns [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [4]\n"
	" -d <string> Reference sequences for SAM input, will invoke sorted-SAM input mode and auto set '-j 100 -W 0 -w 1000'\n"
	" -u          Only process reference regions present in/between SAM alignments\n"
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
	SeqBank *refs;
	FileReader *fr, *db;
	cplist *infs, *dbfs;
	FILE *out;
	char *outf;
	int reglen, use_sse, bandwidth, rW, winlen, winmin, fail_skip, M, X, I, D, E, mincnt, seqmax, wsize, print_lay, sam_present;
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
	while((c = getopt(argc, argv, "hvVt:d:p:ui:o:fj:S:B:W:w:AM:X:I:D:E:R:C:F:N:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'p': print_lay = 1;
			case 'd': push_cplist(dbfs, optarg); reglen = 100; winlen = 0; winmin = 1000; break;
			case 'u': sam_present = 1; break;
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
		run_cns(wb, iter_wtlayblock, info_wtlayblock, ncpu, 0, seqmax, winlen, winmin, fail_skip, bandwidth, M, X, I, D, -1, rW, mincnt, minfreq, reglen, print_lay, out);
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
		sb = init_samblock(refs, fr, winmin, winmin - reglen, sam_present);
		run_cns(sb, iter_samblock, info_samblock, ncpu, 1, seqmax, 0, 0, fail_skip, bandwidth, M, X, I, D, -1, rW, mincnt, minfreq, reglen, print_lay, out);
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
