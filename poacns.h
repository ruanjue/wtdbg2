/*
 * 
 * Copyright (c) 2018, Jue Ruan <ruanjue@gmail.com>
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

#ifndef PO_MSA_CNS_RJ_H
#define PO_MSA_CNS_RJ_H

#include "bit2vec.h"
#include "dna.h"
#include "chararray.h"
#include "list.h"
#include "hashset.h"
#include <emmintrin.h>
#include <tmmintrin.h>

#if __BYTE_ORDER == 1234
//#pragma message(" ** " __FILE__ " has been tested in LITTLE_ENDIAN **\n")
#else
#pragma message(" ** " __FILE__ " hasn't been tested in BIG_ENDIAN **\n")
#endif
static int cns_debug = 0;

#define POG_RDLEN_MAX	0x7FF8
#define POG_RDCNT_MAX	0x3FFF
#define POG_BASEWIDTH	3
#define POG_BASEMAX		0x3F
#define POG_HEAD_NODE	0
#define POG_TAIL_NODE	1
#define POG_DP_BT_M	0
#define POG_DP_BT_I	1
#define POG_DP_BT_D	2
#define POG_ALNMODE_OVERLAP	0
#define POG_ALNMODE_GLOBAL	1

#define POG_VST_MAX	MAX_U2

#define POG_SCORE_MIN	(-(MAX_B2 >> 1))

typedef struct {
	u2i rid, pos, base;
	u2i rbeg, rend, rmax;
	u2i nin, vst;
	u4i edge, erev;
	u4i aligned;
	union {
		u4i aux;
		u4i bpos; // backbone pos
	};
	u4i coff;
	union {
		u8i flag;
		struct { u4i roff, voff; };
	};
} pog_node_t;
define_list(pognodev, pog_node_t);

typedef struct {
	u4i node;
	u4i cov:31, is_aux:1;
	u4i next;
} pog_edge_t;
define_list(pogedgev, pog_edge_t);

typedef struct {
	int cnsmode; // 0: gen_cns_pog, 1: dp_call_cns_pog
	int refmode;
	int alnmode;
	int tribase;
	int sse;
	int near_dialog;
	int W_score;
	int W, rW, cW; // 64, 16, 8
	int Wmax; // 1024
	float W_mat_rate; // 0.92
	int M, X, I, D, E;
	float H; // homopolymer merge
	int T; // bonus for end
	u4i msa_min_cnt;
	float msa_min_freq;
} POGPar;

static const POGPar DEFAULT_POG_PAR = (POGPar){0, 0, POG_ALNMODE_OVERLAP, 0, 1, 0, 0, 64, 16, 8, 1024, 0.92, 2, -5, -2, -4, -1, -3, 20, 3, 0.5};

typedef struct {
	u4i coff:29, bt:3;
	float max;
} pog_cns_t;

typedef struct {
	SeqBank  *seqs;
	u2v      *sbegs, *sends; // suggested [beg, end) on ref(1st seq) in reference-based mode
	pognodev *nodes;
	pogedgev *edges;
	POGPar *par;
	b2v *qprof;
	b2v *rows;
	u4v *rowr;
	b2v *btds;
	u2v *btvs;
	u4v *btxs;
	u4v *stack;
	u2i backbone;
	u4i msa_len;
	u1v *msa;
	u2v *bcnts[7];
	u2v *hcovs;
	u1v *cbts;
	BaseBank *cns;
	u8i ncall;
} POG;

static inline POG* init_pog(POGPar par){
	POG *g;
	g = malloc(sizeof(POG));
	g->seqs = init_seqbank();
	g->sbegs = init_u2v(32);
	g->sends = init_u2v(32);
	g->nodes = init_pognodev(16 * 1024);
	g->edges = init_pogedgev(16 * 1024);
	g->par   = malloc(sizeof(POGPar));
	memcpy(g->par, &par, sizeof(POGPar));
	g->qprof = init_b2v(4 * 1024);
	g->rows  = init_b2v(16 * 1024);
	g->rowr  = init_u4v(64);
	g->btds  = init_b2v(16 * 1024);
	g->btvs  = init_u2v(8 * 1024 * 1024);
	g->btxs  = init_u4v(1024);
	g->stack = init_u4v(32);
	g->backbone = 0;
	g->msa_len = 0;
	g->msa = init_u1v(16 * 1024);
	g->bcnts[0] = init_u2v(4 * 1024);
	g->bcnts[1] = init_u2v(4 * 1024);
	g->bcnts[2] = init_u2v(4 * 1024);
	g->bcnts[3] = init_u2v(4 * 1024);
	g->bcnts[4] = init_u2v(4 * 1024);
	g->bcnts[5] = init_u2v(4 * 1024);
	g->bcnts[6] = init_u2v(4 * 1024);
	g->hcovs = init_u2v(4 * 1024);
	g->cbts  = init_u1v(5 * 2 * 1024);
	g->cns = init_basebank();
	g->ncall = 0;
	return g;
}

static inline void renew_pog(POG *g){
	free_seqbank(g->seqs); g->seqs = init_seqbank();
	renew_u2v(g->sbegs, 32);
	renew_u2v(g->sends, 32);
	renew_pognodev(g->nodes, 16 * 1024);
	renew_pogedgev(g->edges, 16 * 1024);
	renew_b2v(g->qprof, 4 * 1024);
	renew_b2v(g->rows, 16 * 1024);
	renew_u4v(g->rowr, 64);
	renew_b2v(g->btds, 16 * 1024);
	renew_u2v(g->btvs, 8 * 1024 * 1024);
	renew_u4v(g->btxs, 16 * 1024);
	renew_u4v(g->stack, 32);
	renew_u1v(g->msa, 16 * 1024);
	renew_u2v(g->bcnts[0], 4 * 1024);
	renew_u2v(g->bcnts[1], 4 * 1024);
	renew_u2v(g->bcnts[2], 4 * 1024);
	renew_u2v(g->bcnts[3], 4 * 1024);
	renew_u2v(g->bcnts[4], 4 * 1024);
	renew_u2v(g->bcnts[5], 4 * 1024);
	renew_u2v(g->bcnts[6], 4 * 1024);
	renew_u2v(g->hcovs, 4 * 1024);
	renew_u1v(g->cbts, 5 * 2 * 1024);
	free_basebank(g->cns); g->cns = init_basebank();
}

static inline void free_pog(POG *g){
	free_seqbank(g->seqs);
	free_u2v(g->sbegs);
	free_u2v(g->sends);
	free_pognodev(g->nodes);
	free_pogedgev(g->edges);
	free_b2v(g->qprof);
	free_b2v(g->rows);
	free_u4v(g->rowr);
	free_b2v(g->btds);
	free_u2v(g->btvs);
	free_u4v(g->btxs);
	free_u4v(g->stack);
	free_u1v(g->msa);
	free_u2v(g->bcnts[0]);
	free_u2v(g->bcnts[1]);
	free_u2v(g->bcnts[2]);
	free_u2v(g->bcnts[3]);
	free_u2v(g->bcnts[4]);
	free_u2v(g->bcnts[5]);
	free_u2v(g->bcnts[6]);
	free_u2v(g->hcovs);
	free_u1v(g->cbts);
	free_basebank(g->cns);
	free(g->par);
	free(g);
}

static inline void push_pog_core(POG *g, char *seq, u4i len, u2i refbeg, u2i refend){
	if(g->seqs->nseq < POG_RDCNT_MAX){
		len = num_min(len, POG_RDLEN_MAX);
		push_seqbank(g->seqs, NULL, 0, seq, len);
		push_u2v(g->sbegs, refbeg);
		push_u2v(g->sends, refend);
	}
}

static inline void push_pog(POG *g, char *seq, u1i len){ push_pog_core(g, seq, len, 0, 0); }

static inline void fwdbitpush_pog_core(POG *g, u8i *bits, u8i off, u4i len, u2i refbeg, u2i refend){
	if(g->seqs->nseq < POG_RDCNT_MAX){
		len = num_min(len, POG_RDLEN_MAX);
		fwdbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
		push_u2v(g->sbegs, refbeg);
		push_u2v(g->sends, refend);
	}
}

static inline void fwdbitpush_pog(POG *g, u8i *bits, u8i off, u4i len){ fwdbitpush_pog_core(g, bits, off, len, 0, 0); }

static inline void revbitpush_pog_core(POG *g, u8i *bits, u8i off, u4i len, u2i refbeg, u2i refend){
	if(g->seqs->nseq < POG_RDCNT_MAX){
		len = num_min(len, POG_RDLEN_MAX);
		revbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
		push_u2v(g->sbegs, refbeg);
		push_u2v(g->sends, refend);
	}
}

static inline void revbitpush_pog(POG *g, u8i *bits, u8i off, u4i len){ revbitpush_pog_core(g, bits, off, len, 0, 0); }

static inline void print_dot_pog(POG *g, FILE *out){
	pog_node_t *n;
	pog_edge_t *e;
	u4i nidx, eidx;
	fprintf(out, "digraph {\n");
	fprintf(out, "rankdir=LR\n");
	fprintf(out, "N0 [label=\"BEG\"]\n");
	fprintf(out, "N1 [label=\"END\"]\n");
	for(nidx=POG_TAIL_NODE+1;nidx<g->nodes->size;nidx++){
		n = ref_pognodev(g->nodes, nidx);
		fprintf(out, "N%u [label=R%u_%u_%c]\n", nidx, n->rid, n->pos, bit_base_table[(n->base) & 0x03]);
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		n = ref_pognodev(g->nodes, nidx);
		if(n->aligned != nidx){
			fprintf(out, "N%u -> N%u [color=magenta style=dashed]\n", nidx, n->aligned);
		}
		eidx = n->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
#if DEBUG
			if(e->next == eidx){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			eidx = e->next;
			if(e->is_aux) continue;
			fprintf(out, "N%u -> N%u [label=%u]\n", nidx, e->node, e->cov);
		}
	}
	fprintf(out, "}\n");
}

static inline void fprint_dot_pog(POG *g, char *prefix, char *suffix){
	FILE *out;
	out = open_file_for_write(prefix, suffix, 1);
	print_dot_pog(g, out);
	fclose(out);
}

static inline void print_vstdot_pog(POG *g, char *fname){
	FILE *out;
	pog_node_t *n;
	pog_edge_t *e;
	u4i nidx, eidx;
	out = open_file_for_write(fname, NULL, 1);
	fprintf(out, "digraph {\n");
	for(nidx=0;nidx<g->nodes->size;nidx++){
		n = ref_pognodev(g->nodes, nidx);
		if(nidx && n->vst == 0) continue;
		fprintf(out, "N%u [label=\"N%u:%u:%u:%d\"]\n", nidx, nidx, n->aux, n->nin, n->vst);
		eidx = n->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			if(e->is_aux) continue;
			fprintf(out, "N%u -> N%u [label=%u]\n", nidx, e->node, e->cov);
		}
	}
	fprintf(out, "}\n");
	fclose(out);
}

static inline void print_seqs_pog(POG *g, char *prefix, char *suffix){
	FILE *out;
	u4i i;
	out = open_file_for_write(prefix, suffix, 1);
	for(i=0;i<g->seqs->nseq;i++){
		fprintf(out, ">S%u len=%u\n", i, g->seqs->rdlens->buffer[i]);
		println_fwdseq_basebank(g->seqs->rdseqs, g->seqs->rdoffs->buffer[i], g->seqs->rdlens->buffer[i], out);
	}
	fclose(out);
}

static inline void print_msa_pog(POG *g, FILE *out){
	char *str;
	u1i *cns;
	u4i i, j, b, e, c, n;
	str = malloc(g->msa_len + 1);
	fprintf(out, "[POS] ");
	for(i=j=0;i<g->msa_len;i++){
		if((i % 10) == 0){
			fprintf(out, "|%04u", i + 1);
			j += 5;
		} else if(i >= j){
			putc(' ', out);
			j ++;
		}
	}
	fputc('\n', out);
	for(i=0;i<g->seqs->nseq+1;i++){
		if(i == g->seqs->nseq){
			fprintf(out, "[CNS] ");
		} else {
			fprintf(out, "[%03u] ", i);
		}
		b = i * g->msa_len;
		e = b + g->msa_len;
		n = 0;
		for(j=b;j<e;j++){
			c = g->msa->buffer[j];
			if(c < 4) n ++;
			str[j-b] = "ACGT-acgt*"[c];
			//fputc("ACGT-"[c], out);
		}
		str[e-b] = 0;
		fputs(str, out);
		if(i == g->seqs->nseq){
			fprintf(out, "\t%u\t%u\n", (u4i)g->cns->size, n);
		} else {
			fprintf(out, "\t%u\t%u\n", g->seqs->rdlens->buffer[i], n);
		}
	}
	fprintf(out, "[POS] ");
	cns = ref_u1v(g->msa, g->msa_len * g->seqs->nseq);
	for(i=j=b=0;i<g->msa_len;i++){
		if(cns[i] < 4){
			j ++;
			if((j % 10) == 1){
				while(b < i){
					fputc(' ', out);
					b ++;
				}
				fprintf(out, "|%04u", j);
				b += 5;
			}
		}
	}
	fprintf(out, "\n");
	if(1){
		u4i divs[5];
		u2i *hcovs;
		hcovs = g->hcovs->buffer;
		divs[0] = 10000;
		divs[1] = 1000;
		divs[2] = 100;
		divs[3] = 10;
		divs[4] = 1;
		for(i=0;i<4;i++){
			for(j=0;j<g->msa_len;j++){
				if(hcovs[j] < divs[i + 1]){
					str[j] = ' ';
				} else {
					str[j] = '0' + ((hcovs[j] % divs[i]) / divs[i + 1]);
				}
			}
			str[j] = 0;
			fprintf(out, "[%c  ] %s\n", "HCOV"[i], str);
		}
	}
	free(str);
}

static inline pog_node_t* add_node_pog(POG *g, u2i rid, u2i pos, u2i base){
	pog_node_t *u;
	u = next_ref_pognodev(g->nodes);
	ZEROS(u);
	u->rid  = rid;
	u->pos  = pos;
	u->base = base;
	u->aligned = offset_pognodev(g->nodes, u);
	return u;
}

static inline pog_edge_t* add_edge_pog(POG *g, pog_node_t *u, pog_node_t *v, int cov, int is_aux){
	pog_edge_t *e, *f, *p;
	u4i eidx;
#if DEBUG
	if(u == v){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	if(u->edge == 0){
		e = next_ref_pogedgev(g->edges);
		u->edge = offset_pogedgev(g->edges, e);
		e->node = offset_pognodev(g->nodes, v);
		e->cov = cov;
		e->is_aux = is_aux;
		e->next = 0;
		v->nin ++;
		return e;
	} else {
		e = f = p = NULL;
		encap_pogedgev(g->edges, 1);
		eidx = u->edge;
		while(eidx){
			f = ref_pogedgev(g->edges, eidx);
			eidx = f->next;
			if(f->node == offset_pognodev(g->nodes, v)){
				e = f;
				break;
			}
			p = f;
		}
		if(e == NULL){
			e = next_ref_pogedgev(g->edges);
			e->node = offset_pognodev(g->nodes, v);
			e->cov = cov;
			e->is_aux = is_aux;
			e->next = 0;
			v->nin ++;
		} else {
			e->is_aux &= is_aux;
			e->cov += cov;
			if(cov == 0){
				return e;
			} else if(p == NULL){
				return e;
			}
			// detach e from u->edge
			p->next = e->next;
			e->next = 0;
		}
		f = ref_pogedgev(g->edges, u->edge);
		if(e->cov > f->cov){
			e->next = u->edge;
			u->edge = offset_pogedgev(g->edges, e);
		} else {
			while(f->next){
				if(ref_pogedgev(g->edges, f->next)->cov < e->cov){
					break;
				}
				f = ref_pogedgev(g->edges, f->next);
			}
			e->next = f->next;
			f->next = offset_pogedgev(g->edges, e);
		}
		return e;
	}
}

static inline pog_node_t* merge_node_pog(POG *g, pog_node_t *x, pog_node_t *u){
	pog_node_t *v, *w;
	pog_edge_t *e;
	u4i xidx, eidx;
	xidx = offset_pognodev(g->nodes, x);
	do {
		v = ref_pognodev(g->nodes, xidx);
		if(v->nin && v->base == u->base){
			eidx = u->edge;
			while(eidx){
				e = ref_pogedgev(g->edges, eidx);
				eidx = e->next;
				w = ref_pognodev(g->nodes, e->node);
				add_edge_pog(g, v, w, e->cov, e->is_aux);
				w->nin --; // will delete e
			}
			u->edge = 0;
			break;
		}
		v = NULL;
	} while(xidx != offset_pognodev(g->nodes, x));
	u->aligned = x->aligned;
	x->aligned = offset_pognodev(g->nodes, u);
	return v? v : u;
}

static inline void beg_pog_core(POG *g, BaseBank *cns, u8i off, u2i len, int clear_all){
	pog_node_t *head, *tail, *u, *v;
	pog_edge_t *e;
	u4i i, bb, bk;
	g->ncall ++;
	if((g->ncall % 16) == 0){
		renew_pog(g);
	}
	clear_and_encap_pognodev(g->nodes, 2 + len);
	clear_pogedgev(g->edges);
	ZEROS(next_ref_pogedgev(g->edges));
	head = next_ref_pognodev(g->nodes);
	ZEROS(head);
	if(g->par->tribase){
		head->base = POG_BASEMAX + 1;
	} else {
		head->base = 4;
	}
	head->aligned = POG_HEAD_NODE;
	tail = next_ref_pognodev(g->nodes);
	ZEROS(tail);
	if(g->par->tribase){
		tail->base = POG_BASEMAX + 1;
	} else {
		tail->base = 4;
	}
	tail->aligned = POG_TAIL_NODE;
	u = head;
	bk = 0;
	// add backbone nodes
	g->backbone = len;
	for(i=0;i<len;i++){
		bb = get_basebank(cns, off + i);
		if(g->par->tribase){
			bk = ((bk << 2) | bb) & POG_BASEMAX;
		} else {
			bk = bb;
		}
		v = add_node_pog(g, g->seqs->nseq, i, bk);
		v->bpos = i;
		e = add_edge_pog(g, u, v, 0, 1);
		u = v;
	}
	e = add_edge_pog(g, u, tail, 0, 1);
	if(clear_all){
		clear_seqbank(g->seqs);
		clear_u2v(g->sbegs);
		clear_u2v(g->sends);
		clear_basebank(g->cns);
	}
}

static inline void beg_pog(POG *g){
	beg_pog_core(g, NULL, 0, 0, 1);
}

static inline void prepare_rd_align_pog(POG *g, u2i rid){
	pog_node_t *u, *v;
	pog_edge_t *e;
	b2i *row, *btds;
	u8i rmax;
	u4i seqoff, seqlen, seqlex, slen, seqinc;
	u4i nidx, eidx, i, j, bb;
	seqoff = g->seqs->rdoffs->buffer[rid];
	seqlen = g->seqs->rdlens->buffer[rid];
	seqlex = roundup_times(seqlen, 8); // 128 = 8 * sizeof(b2i)
	slen = seqlex >> 3;
	seqinc = seqlex + 8; // seqlex, max_idx, beg, end, and paddings
	// init graph
	encap_pognodev(g->nodes, seqlen + 2);
	encap_pogedgev(g->edges, seqlen + 2);
	// estimate cap of g->rows
	for(i=0;i<g->nodes->size;i++){
		v = ref_pognodev(g->nodes, i);
		v->roff = 0;
		v->voff = 0;
		v->coff = 0;
		v->vst  = 0;
		v->aux  = 0;
	}
	clear_u4v(g->stack);
	push_u4v(g->stack, POG_HEAD_NODE);
	rmax = 0;
	bb = 1;
	while(pop_u4v(g->stack, &nidx)){
		u = ref_pognodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
#if DEBUG
			if(eidx == e->next){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			eidx = e->next;
			v = ref_pognodev(g->nodes, e->node);
			if(v->vst == 0){
				bb ++;
				if(bb > rmax){
					rmax = bb;
				}
			}
			v->vst ++;
			if(v->vst == v->nin){
				push_u4v(g->stack, e->node);
			}
		}
		bb --;
	}
#if DEBUG
	if(bb){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	rmax += 2;
	head_sl_b2v(g->rows, g->rows->n_head);
	clear_and_encap_b2v(g->rows, rmax * seqinc + 8);
	head_sr_b2v(g->rows, (16 - (((u8i)g->rows->buffer) & 0xF)) >> 1);
	head_sl_b2v(g->btds, g->btds->n_head);
	clear_and_encap_b2v(g->btds, rmax * seqinc + 8);
	head_sr_b2v(g->btds, (16 - (((u8i)g->btds->buffer) & 0xF)) >> 1);
	inc_b2v(g->rows, seqinc);
	inc_b2v(g->btds, seqinc);
	clear_u4v(g->rowr);
	bb = 0;
	for(i=0;i<g->nodes->size;i++){
		v = ref_pognodev(g->nodes, i);
		v->vst = 0;
		v->erev = bb;
		bb += v->nin;
		if(i < 2) bb ++; // in case of add mnode
	}
	clear_u2v(g->btvs);
	inc_u2v(g->btvs, 8);
	clear_and_inc_u4v(g->btxs, bb + 1);
	//g->btxs->buffer[0] = 0;
	memset(g->btxs->buffer, 0, (bb + 1) * sizeof(u4i));
	u = ref_pognodev(g->nodes, POG_HEAD_NODE);
	if(g->rowr->size){
		u->coff = u->roff = g->rowr->buffer[--g->rowr->size];
	} else {
		u->roff = g->rows->size;
		inc_b2v(g->rows, seqinc);
		u->coff = g->btds->size;
		inc_b2v(g->btds, seqinc);
	}
	row = ref_b2v(g->rows, u->roff);
	btds = ref_b2v(g->btds, u->coff);
	if(g->par->alnmode == POG_ALNMODE_OVERLAP){ // overlap alignment
		memset(row,  0, seqinc * sizeof(b2i));
		//row[0] = g->par->T;
	} else { // global
		__m128i MIN;
		MIN = _mm_set1_epi16(POG_SCORE_MIN);
		for(j=0;j<slen;j++){
			_mm_store_si128(((__m128i*)row) + j, MIN);
		}
		//row[0] = g->par->T;
	}
	memset(btds, 0, seqinc * sizeof(b2i));
	if(g->par->near_dialog && g->par->W_score <= 0){
		row[seqlex] = g->par->W;
		row[seqlex + 1] = 0;
		row[seqlex + 2] = num_min(Int(seqlex), 2 * g->par->W + 1);
	} else {
		row[seqlex] = 0;
		row[seqlex + 1] = 0;
		row[seqlex + 2] = seqlex;
	}
}

// SSE
// BANDED, Auto fit the previous-row's max score in center
// OVERLAP
static inline void sse_band_row_rdaln_pog(POG *g, u4i nidx1, u4i nidx2, u4i seqlen, u2i vst, u4i coff1, u4i coff2, u4i roff1, u4i roff2, b2i *qp, int W, int center){
	__m128i I, D, H, E, F, S, MAX, MIN, CMP;
	__m128i BT, BT1, BT2, BT1_MASK, BT2_MASK, BTX_MASK;
	b2i *row1, *row2, *btds;
	u4i i, slen, seqlex, beg, end;
	int lsth, lstf, msk, mi;
	UNUSED(coff1);
	slen = (seqlen + 7) / 8;
	seqlex = slen * 8;
	I = _mm_set1_epi16(g->par->I);
	D = _mm_set1_epi16(g->par->D);
	MIN = _mm_set1_epi16(POG_SCORE_MIN);
	BT1_MASK = _mm_set1_epi16(0b10);
	BT2_MASK = _mm_set1_epi16(0b01);
	BTX_MASK = _mm_set1_epi16(vst << 2);
	row1 = ref_b2v(g->rows, roff1);
	row2 = ref_b2v(g->rows, roff2);
	btds = ref_b2v(g->btds, coff2);
	if(W){
		if(center < 0){
			if(row1[row1[seqlex]] >= g->par->W_score){
				if(g->par->near_dialog){
					if(row1[seqlex + 1] > 0 || row1[seqlex + 2] < Int(seqlex)){
						if(row1[seqlex] == (row1[seqlex + 1] + row1[seqlex + 2]) / 2){
							center = (row1[seqlex + 1] + row1[seqlex + 2]) / 2 + 1;
						} else if(row1[seqlex] > (row1[seqlex + 1] + row1[seqlex + 2]) / 2){
							center = (row1[seqlex + 1] + row1[seqlex + 2]) / 2 + 2;
						} else {
							center = (row1[seqlex + 1] + row1[seqlex + 2]) / 2;
						}
					} else { // first set center
						center = row1[seqlex];
					}
				} else {
					center = row1[seqlex];
				}
				beg = num_max(center - W, row1[seqlex + 1]);
				end = num_min(center + 1 + W, row1[seqlex + 2] + 1);
				if(end > seqlex) end = seqlex;
			} else {
				beg = row1[seqlex + 1];
				end = num_min(row1[seqlex + 2] + 1, Int(seqlex));
			}
		} else {
			beg = num_max(center - W, row1[seqlex + 1]);
			end = num_min(center + 1 + W, row1[seqlex + 2] + 1);
			if(end > seqlex) end = seqlex;
		}
	} else {
		beg = 0;
		end = seqlex;
	}
	if(beg & 0x7U) beg = beg & (~0x7U);
	if(end & 0x7U) end = (end + 8) & (~0x7U);
	beg = beg / 8;
	end = (end + 7) / 8;
	if(row1[seqlex + 1] >= row1[seqlex + 2]){ // happening in realign mode
		for(i=beg*8;i<end*8;i+=8){
			_mm_store_si128(((__m128i*)(row1 + i)), MIN);
		}
	} else {
		if(Int(beg * 8) < row1[seqlex + 1]){
			for(i=beg*8;Int(i)<row1[seqlex + 1];i+=8){
				_mm_store_si128(((__m128i*)(row1 + i)), MIN);
			}
		}
		if(Int(end * 8) > row1[seqlex + 2]){
			for(i=row1[seqlex+2];i<end*8;i+=8){
				_mm_store_si128(((__m128i*)(row1 + i)), MIN);
			}
		}
	}
	if(nidx2 == POG_TAIL_NODE){
		while(end * 8 < seqlex){
			_mm_store_si128(((__m128i*)(row1)) + end, MIN);
			end ++;
		}
	}
	MAX = _mm_set1_epi16(POG_SCORE_MIN);
	if(beg){
		lstf = lsth = POG_SCORE_MIN;
	} else {
		if(row1[seqlex + 1] >= row1[seqlex + 2]){
			lstf = POG_SCORE_MIN;
			lsth = 0; // auto global alignment
		} else {
			lstf = (g->par->alnmode == POG_ALNMODE_OVERLAP)? 0 : POG_SCORE_MIN;
			if(nidx1){
				lsth = (g->par->alnmode == POG_ALNMODE_OVERLAP)? 0 : POG_SCORE_MIN;;
			} else {
				lsth = g->par->T;
			}
		}
	}
	for(i=beg;i<end;i++){
		E = _mm_load_si128(((__m128i*)(row1)) + i);
		H = _mm_slli_si128(E, 2);
		H = _mm_insert_epi16(H, lsth, 0);
		H = _mm_adds_epi16(H, _mm_load_si128(((__m128i*)qp) + i));
		lsth = _mm_extract_epi16(E, 7);
		E = _mm_adds_epi16(E, D);
		S = _mm_max_epi16(H, E);
		BT1 = _mm_cmpgt_epi16(S, H);
		BT1 = _mm_and_si128(BT1, BT1_MASK);
		F = _mm_slli_si128(S, 2);
		F = _mm_insert_epi16(F, lstf, 0);
		F = _mm_adds_epi16(F, I);
		F = _mm_max_epi16(S, F);
		while(1){
			H = _mm_slli_si128(F, 2);
			H = _mm_insert_epi16(H, lstf, 0);
			H = _mm_adds_epi16(H, I);
			CMP = _mm_cmpgt_epi16(H, F);
			msk = _mm_movemask_epi8(CMP);
			if(msk == 0){
				break;
			} else {
				F = _mm_max_epi16(H, F);
			}
		}
		MAX = _mm_max_epi16(MAX, F);
		lstf = _mm_extract_epi16(F, 7);
		BT2 = _mm_cmpgt_epi16(F, S);
		BT2 = _mm_and_si128(BT2, BT2_MASK);
		BT  = _mm_xor_si128(BT1, BT2); // (0, (1, 3), 2)
		BT  = _mm_xor_si128(BT, BTX_MASK); // (0, (1, 3), 2) | (vst << 2)
		_mm_store_si128(((__m128i*)row2) + i, F);
		_mm_stream_si128(((__m128i*)btds) + i, BT);
	}

	if(1){
		b2i ary[8];
		u4i j, k;
		_mm_store_si128((__m128i*)ary, MAX);
		k = 0;
		for(j=1;j<8;j++){
			if(ary[j] > ary[k]){
				k = j;
			}
		}
		mi = beg * 8 + k;
		for(i=beg+1;i<end;i++){
			if(row2[i * 8 + k] > row2[mi]){
				mi = i * 8 + k;
			}
		}
	} else {
		mi = beg * 8;
		for(i=mi+1;i<end*8;i++){
			if(row2[i] > row2[mi]){
				mi = i;
			}
		}
	}
	row2[seqlex] = mi;
	row2[seqlex + 1] = beg * 8;
	row2[seqlex + 2] = end * 8;
}

static inline void merge_row_rdaln_pog(POG *g, u4i seqlen, u4i coff1, u4i coff2, u4i roff1, u4i roff2){
	b2i *row1, *row2, *btd1, *btd2;
	u4i i, seqlex, beg[3], end[3], sz, b;
	row1 = ref_b2v(g->rows, roff1);
	row2 = ref_b2v(g->rows, roff2);
	btd1 = ref_b2v(g->btds, coff1);
	btd2 = ref_b2v(g->btds, coff2);
	seqlex = roundup_times(seqlen, 8);
	beg[0] = (row1[seqlex + 1]);
	beg[1] = (row2[seqlex + 1]);
	end[0] = (row1[seqlex + 2]);
	end[1] = (row2[seqlex + 2]);
	beg[2] = num_max(beg[0], beg[1]);
	end[2] = num_min(end[0], end[1]);
	if(1){
		__m128i R1, R2, D1, D2, MK;
		for(i=beg[2];i&0x7;i++){
			if(row2[i] < row1[i]){
				row2[i] = row1[i];
				btd2[i] = btd1[i];
			}
		}
		for(;i+8<=end[2];i+=8){
			R1 = _mm_load_si128((__m128i*)(row1 + i));
			R2 = _mm_load_si128((__m128i*)(row2 + i));
			D1 = _mm_load_si128((__m128i*)(btd1 + i));
			D2 = _mm_load_si128((__m128i*)(btd2 + i));
			MK = _mm_cmpgt_epi16(R1, R2);
			R2 = _mm_max_epi16(R1, R2);
			_mm_store_si128((__m128i*)(row2 + i), R2);
			D2 = _mm_or_si128(_mm_and_si128(MK, D1), _mm_andnot_si128(MK, D2));
			_mm_store_si128((__m128i*)(btd2 + i), D2);
		}
		for(;i<end[2];i++){
			if(row2[i] < row1[i]){
				row2[i] = row1[i];
				btd2[i] = btd1[i];
			}
		}
	} else {
		for(i=beg[2];i<end[2];i++){
			if(row2[i] < row1[i]){
				row2[i] = row1[i];
				btd2[i] = btd1[i];
			}
		}
	}
	if(0){
		if(beg[0] < beg[2]){
			sz = num_min(beg[2], end[0]) - beg[0];
			memcpy(row2 + beg[0], row1 + beg[0], sz * sizeof(b2i));
			memcpy(btd2 + beg[0], btd1 + beg[0], sz * sizeof(b2i));
		}
		for(i=end[2];i<beg[2];i++){ // in case of two independent regions
			row2[i] = POG_SCORE_MIN;
		}
		if(end[0] > end[2]){
			sz = num_min(end[0], seqlex) - end[2];
			memcpy(row2 + end[2], row1 + end[2], sz * sizeof(b2i));
			memcpy(btd2 + end[2], btd1 + end[2], sz * sizeof(b2i));
		}
	} else {
		if(beg[0] < beg[1]){
			b = beg[0];
			if(end[0] >= beg[1]){
				sz = beg[1] - beg[0];
				memcpy(row2 + b, row1 + b, sz * sizeof(b2i));
				memcpy(btd2 + b, btd1 + b, sz * sizeof(b2i));
			} else {
				__m128i F;
				F = _mm_set1_epi16(POG_SCORE_MIN);
				sz = end[0] - beg[0];
				memcpy(row2 + b, row1 + b, sz * sizeof(b2i));
				memcpy(btd2 + b, btd1 + b, sz * sizeof(b2i));
				for(i=end[0];i<beg[1];i+=8){
					_mm_store_si128(((__m128i*)(row2 + i)), F);
				}
			}
		}
		if(end[0] > end[1]){
			if(beg[0] <= end[1]){
				b = end[1];
				sz = end[0] - b;
				memcpy(row2 + b, row1 + b, sz * sizeof(b2i));
				memcpy(btd2 + b, btd1 + b, sz * sizeof(b2i));
			} else {
				__m128i F;
				F = _mm_set1_epi16(POG_SCORE_MIN);
				b = beg[0];
				sz = end[0] - b;
				memcpy(row2 + b, row1 + b, sz * sizeof(b2i));
				memcpy(btd2 + b, btd1 + b, sz * sizeof(b2i));
				for(i=end[1];i<beg[0];i+=8){
					_mm_store_si128(((__m128i*)(row2 + i)), F);
				}
			}
		}
	}
	row2[seqlex + 1] = num_min(beg[0], beg[1]);
	row2[seqlex + 2] = num_max(end[0], end[1]);
	if(row1[(row1[seqlex])] > row2[(row2[seqlex])]){
		row2[seqlex] = row1[seqlex];
	}
}

static inline void set_rd_query_prof(POG *g, u4i rid){
	b2i *qp;
	u4i seqoff, seqlen, seqlex, slen, seqinc;
	u4i i, j, bb, bk;
	seqoff = g->seqs->rdoffs->buffer[rid];
	seqlen = g->seqs->rdlens->buffer[rid];
	seqlex = roundup_times(seqlen, 8);
	slen = seqlex >> 3;
	seqinc = seqlex + 8;
	head_sl_b2v(g->qprof, g->qprof->n_head);
	if(g->par->tribase){
		clear_and_encap_b2v(g->qprof, seqlex * (POG_BASEMAX + 2) + 8);
		head_sr_b2v(g->qprof, (16 - (((u8i)g->qprof->buffer) & 0xF)) >> 1);
		for(i=0;i<=POG_BASEMAX;i++){
			qp = g->qprof->buffer + i * seqlex;
			bk = 0;
			for(j=0;j<seqlen;j++){
				bb = get_basebank(g->seqs->rdseqs, seqoff + j);
				bk = ((bk << 2) | bb) & POG_BASEMAX;
				if(bb == (i & 0x03)){
					if(bk == i){
						qp[j] = g->par->M + g->par->tribase;
					} else {
						qp[j] = g->par->M;
					}
				} else {
					qp[j] = g->par->X;
				}
			}
			for(;j<seqlex;j++){
				qp[j] = POG_SCORE_MIN;
			}
		}
	} else {
		clear_and_encap_b2v(g->qprof, seqlex * (4 + 1) + 8);
		head_sr_b2v(g->qprof, (16 - (((u8i)g->qprof->buffer) & 0xF)) >> 1);
		for(i=0;i<4;i++){
			qp = g->qprof->buffer + i * seqlex;
			for(j=0;j<seqlen;j++){
				bb = get_basebank(g->seqs->rdseqs, seqoff + j);
				if(bb == (i & 0x03)){
					qp[j] = g->par->M;
				} else {
					qp[j] = g->par->X;
				}
			}
			for(;j<seqlex;j++){
				qp[j] = POG_SCORE_MIN;
			}
		}
	}
	{
		__m128i X;
		qp = g->qprof->buffer + i * seqlex;
		X = _mm_set1_epi16(g->par->X);
		for(j=0;j<slen;j++){
			_mm_stream_si128(((__m128i*)qp) + j, X);
		}
	}
}

static inline u2i get_rdbase_pog(POG *g, u4i rid, u4i pos){
	u4i seqoff;
	seqoff = g->seqs->rdoffs->buffer[rid];
	if(g->par->tribase == 0){
		return get_basebank(g->seqs->rdseqs, seqoff + pos);
	}
	if(pos >= POG_BASEWIDTH){
		return sub32seqbits(g->seqs->rdseqs->bits, seqoff + pos + 1 - POG_BASEWIDTH) >> (64 - (POG_BASEWIDTH << 1));
	} else {
		return (sub32seqbits(g->seqs->rdseqs->bits, seqoff)) >> (64 - ((pos + 1) << 1));
	}
}

static inline int _cal_matches_alignment_pog(POG *g, u4i rid, int *xb, int xe, int *badtail){
	pog_node_t *u, *v;
	u4i nidx, seqoff, btx, bt, vst;
	int x, seqlen, mat, score, x0, x1, flag;
	seqlen = g->seqs->rdlens->buffer[rid];
	seqoff = g->seqs->rdoffs->buffer[rid];
	x = xe;
	score = 0;
	x0 = x1 = x;
	flag = 0;
	nidx = POG_TAIL_NODE;
	v = ref_pognodev(g->nodes, nidx);
	btx = 1;
	mat = 0;
	while(x >= 0){
		u  = ref_pognodev(g->nodes, nidx);
		bt = g->btvs->buffer[u->voff + x - u->rbeg];
		vst = bt >> 2;
		bt  = bt & 0x03;
		if(bt == POG_DP_BT_M){
			if(flag) score += ((u->base & 0x03) == get_basebank(g->seqs->rdseqs, seqoff + x))? g->par->M : g->par->X;
			else { flag = 1; x1 = x; }
			mat ++;
			x --;
		} else if(bt & POG_DP_BT_I){
			if(flag) score += g->par->I;
			x --;
		} else {
			if(flag) score += g->par->D;
			else { flag = 1; x1 = x; }
		}
		if(score == 10 * g->par->M){
			x0 = x + 1;
		}
		if(bt & POG_DP_BT_I){
		} else {
			nidx = g->btxs->buffer[g->nodes->buffer[nidx].erev + vst];
		}
		if(x < 0){
			break;
		}
		u = ref_pognodev(g->nodes, nidx);
		if(x < u->rbeg || x >= u->rend){
			break;
		}
	}
	*xb = x;
	*badtail = x1 > x0? x1 - x0 : 0;
	return mat;
}

static inline int _alignment2graph_pog(POG *g, u4i rid, int xb, int xe){
	pog_node_t *u, *v;
	pog_edge_t *e;
	u4i nidx, i, seqoff, btx, bt, vst;
	int x, seqlen;
	seqlen = g->seqs->rdlens->buffer[rid];
	seqoff = g->seqs->rdoffs->buffer[rid];
	x = xe;
	nidx = POG_TAIL_NODE;
	v = ref_pognodev(g->nodes, nidx);
	int my_print = (cns_debug > 3);
	if(x + 1 < (int)seqlen){
		for(i=seqlen-1;(int)i>x;i--){
			if(my_print){
				fprintf(stderr, "BT[%u] y=N%u_R%u_%u_%c x=%d:%c z=%d\n", rid, nidx, v->rid, v->pos, bit_base_table[(v->base) & 0x03], i, bit_base_table[get_basebank(g->seqs->rdseqs, seqoff + i)], -1);
				fflush(stderr);
			}
			u = add_node_pog(g, rid, i, get_rdbase_pog(g, rid, i));
			add_edge_pog(g, u, v, 1, 0);
			v = u;
		}
		v->coff = ref_pognodev(g->nodes, POG_TAIL_NODE)->coff;
	}
	btx = 1;
	while(1){
		if(x >= 0){
			u  = ref_pognodev(g->nodes, nidx);
			bt = g->btvs->buffer[u->voff + x - u->rbeg];
		} else {
			bt = 0;
		}
		vst = bt >> 2;
		bt  = bt & 0x03;
		if(my_print){
			pog_node_t *w;
			w = ref_pognodev(g->nodes, nidx);
			fprintf(stderr, "BT[%u] y=N%u_R%u_%u_%c x=%d:%c z=%d\n", rid, nidx, w->rid, w->pos, bit_base_table[(w->base) & 0x03], x, x>=0? bit_base_table[get_basebank(g->seqs->rdseqs, seqoff + x)] : '*', bt);
			fflush(stderr);
		}
		if(bt == POG_DP_BT_M){
			if(btx){
				u = add_node_pog(g, rid, x, get_rdbase_pog(g, rid, x));
				e = add_edge_pog(g, u, v, 1, 0);
				x --;
				if(nidx > POG_TAIL_NODE){
					u = merge_node_pog(g, ref_pognodev(g->nodes, nidx), u);
				}
				v = u;
			} else {
				xb = x;
				if(nidx == POG_HEAD_NODE || nidx + 1 >= g->nodes->size || g->nodes->buffer[nidx].rid != g->nodes->buffer[nidx + 1].rid){
					nidx = POG_HEAD_NODE;
				} else {
					nidx ++;
				}
				u = ref_pognodev(g->nodes, nidx);
				e = add_edge_pog(g, u, v, 1, 0);
				if(nidx != POG_HEAD_NODE){
					v = u;
					u = ref_pognodev(g->nodes, POG_HEAD_NODE);
					//e = add_edge_pog(g, u, v, 1, 0);
					e = add_edge_pog(g, u, v, 0, 1);
				}
				break;
			}
		} else if(bt & POG_DP_BT_I){
			u = add_node_pog(g, rid, x, get_rdbase_pog(g, rid, x));
			e = add_edge_pog(g, u, v, 1, 0);
			x --;
			v = u;
		} else {
		}
		if(x < 0){ // && nidx
			if(bt == POG_DP_BT_I){
			} else {
				nidx = 0;
			}
			btx = 0;
			continue;
		}
		if(bt & POG_DP_BT_I){
		} else {
			nidx = g->btxs->buffer[g->nodes->buffer[nidx].erev + vst];
		}
#if DEBUG
		if(nidx >= g->nodes->size){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
#endif
		u = ref_pognodev(g->nodes, nidx);
		if(x < u->rbeg || x >= u->rend){
			xb = x;
			while(x >= 0){
				u = add_node_pog(g, rid, x, get_rdbase_pog(g, rid, x));
				e = add_edge_pog(g, u, v, 1, 0);
				x --;
				v = u;
			}
			{
				u = ref_pognodev(g->nodes, POG_HEAD_NODE);
				e = add_edge_pog(g, u, v, 1, 0);
				v = u;
			}
			break;
		}
	}
	return xb;
}

static inline int align_rd_pog_core(POG *g, u2i rid, int W, int *xe){
	pog_node_t *u, *v;
	pog_edge_t *e;
	u4i seqoff, seqlen, seqlex, slen, seqinc;
	b4i score, x;
	b2i *qp, *row;
	__m128i SMASK;
	u4i nidx, eidx, coff, roff, mnode;
	u4i i, j;
	seqoff = g->seqs->rdoffs->buffer[rid];
	seqlen = g->seqs->rdlens->buffer[rid];
	seqlex = roundup_times(seqlen, 8); // 128 = 8 * sizeof(b2i)
	slen = seqlex >> 3;
	seqinc = seqlex + 8; // seqlex, max_idx, beg, end, and paddings
#if __BYTE_ORDER == 1234
	SMASK = _mm_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80);
#else
	// TODO: need to test in BIG_ENDIAN
	SMASK = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80);
#endif
	// set query prof
	// fit buffer to 16 bytes aligned
	set_rd_query_prof(g, rid);
	prepare_rd_align_pog(g, rid);
	clear_u4v(g->stack);
	push_u4v(g->stack, POG_HEAD_NODE);
	score = POG_SCORE_MIN;
	mnode = POG_TAIL_NODE; // point to POG_TAIL
	x = seqlen - 1;
	int my_print = (cns_debug > 3);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_pognodev(g->nodes, nidx);
		u->rmax = g->rows->buffer[u->roff + seqlex];
		u->rbeg = g->rows->buffer[u->roff + seqlex + 1];
		u->rend = g->rows->buffer[u->roff + seqlex + 2];
		if(my_print){
			fprintf(stderr, "NODEALN[%u:R%u:%u] %d-%d:%d:%d %d %d\n", nidx, u->rid, u->pos, u->rbeg, u->rend, u->rmax,
				g->rows->buffer[u->roff + u->rmax], g->btds->buffer[u->coff + u->rmax] >> 2, g->btds->buffer[u->coff + u->rmax] & 0x03);
		}
		if(g->par->alnmode == POG_ALNMODE_OVERLAP){
			if(u->rend >= seqlen &&  g->rows->buffer[u->roff + x] > score){ // overlap alignment
				score = g->rows->buffer[u->roff + x];
				mnode = nidx;
			}
		}
		eidx = u->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_pognodev(g->nodes, e->node);
			qp = g->qprof->buffer + v->base * seqlex;
			if(g->rowr->size){
				coff = roff = g->rowr->buffer[-- g->rowr->size];
			} else {
#if 0
				if(g->rows->size + seqinc + g->rows->n_head > g->rows->cap){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
#endif
				roff = g->rows->size;
				inc_b2v(g->rows, seqinc);
				coff = g->btds->size;
				inc_b2v(g->btds, seqinc);
			}
			g->btxs->buffer[v->erev + v->vst] = nidx; // save backtrace nidx
			sse_band_row_rdaln_pog(g, nidx, e->node, seqlen, v->vst, u->coff, coff, u->roff, roff, qp, W, -1);
			if(v->vst){
				merge_row_rdaln_pog(g, seqlen, coff, v->coff, roff, v->roff);
				push_u4v(g->rowr, roff);
			} else {
				v->coff = coff;
				v->roff = roff;
			}
			v->vst ++;
			if(v->vst == v->nin){
				push_u4v(g->stack, e->node);
#if DEBUG
			} else if(v->vst > v->nin){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
#endif
			}
		}
		push_u4v(g->rowr, u->roff);
		{ // compress-copy btds into btvs
			u->voff = g->btvs->size;
			inc_u2v(g->btvs, u->rend - u->rbeg);
			if(0){
				for(i=u->rbeg;i<u->rend;i++){
					g->btvs->buffer[u->voff + i - u->rbeg] = (u1i)g->btds->buffer[u->coff + i];
				}
			} else {
				memcpy(g->btvs->buffer + u->voff, g->btds->buffer + u->coff + u->rbeg, (u->rend - u->rbeg) * sizeof(u2i));
			}
		}
	}
	v = ref_pognodev(g->nodes, POG_TAIL_NODE);
#if DEBUG
	if(v->roff == 0){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	if(v->voff == 0){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
#endif
	if(W == 0){
		row = ref_b2v(g->rows, v->roff);
		j = 0;
		for(i=1;i<seqlen;i++){
			if(row[i] > row[j]){
				j = i;
			}
		}
		v->rmax = row[seqlex] = j;
		v->rbeg = row[seqlex + 1] = 0;
		v->rend = row[seqlex + 2] = seqlex;
	}
	if(g->par->alnmode == POG_ALNMODE_OVERLAP){
		if(g->rows->buffer[v->roff + v->rmax] + g->par->T > score){
			score = g->rows->buffer[v->roff + (g->rows->buffer[v->roff + seqlex])];
			mnode = POG_TAIL_NODE;
			x = (g->rows->buffer[v->roff + seqlex]);
		}
	} else {
		score = g->rows->buffer[v->roff + seqlen - 1];
		mnode = POG_TAIL_NODE;
		x = seqlen - 1;
	}
	if(x == Int(seqlen - 1) && mnode != POG_TAIL_NODE){
		v = ref_pognodev(g->nodes, POG_TAIL_NODE);
		g->btvs->buffer[v->voff + x - v->rbeg] = (2 | (v->vst << 2));
		g->btxs->buffer[v->erev + v->vst] = mnode;
		v->vst ++;
	}
	*xe = x;
	return score;
}

static inline int align_rd_pog(POG *g, u2i rid){
	int xb, xe, rlen, bad, xlen, score, W, mat;
	rlen = g->seqs->rdlens->buffer[rid];
	xlen = rlen / 8 * 8;
	W = g->par->W? g->par->W : xlen;
	xlen = num_min(xlen, g->par->Wmax);
	mat = 0;
	// try increase W when align score is low
	while(1){
		score = align_rd_pog_core(g, rid, W, &xe);
		mat = _cal_matches_alignment_pog(g, rid, &xb, xe, &bad);
		if(cns_debug > 1){
			fprintf(stderr, "ALIGN[%03d] len=%u ref=%d,%d band=%d aligned=%d,%d tail=%d mat=%d,%0.3f score=%d\n", rid, g->seqs->rdlens->buffer[rid], g->sbegs->buffer[rid], g->sends->buffer[rid], W, xb + 1, xe + 1, bad, mat, 1.0 * mat / rlen, score);
		}
		if(rid == 0){
			break;
		}
		if(mat >= rlen * g->par->W_mat_rate && bad < 50){
			break;
		}
		if(W >= xlen) break;
		W = num_min(W * 2, xlen);
	}
	xb = 0;
	xb = _alignment2graph_pog(g, rid, xb, xe);
	return score;
}

static inline u4i realign_msa_pog_core(POG *g, u4i ridx, int W){
	u2i *bcnts[7], *hcovs, *bts, *bs, *bs1;
	u1i *r, *s;
	int wins[20], winl, winy, winx;
	u4i i, *dps, f, h, e, *row1, *row2, max, bt;
	int j, beg, end, off, x, y;
	hcovs = g->hcovs->buffer;
	for(i=0;i<7;i++){
		clear_and_encap_u2v(g->bcnts[i], g->msa_len);
		bcnts[i] = g->bcnts[i]->buffer;
		memset(bcnts[i], 0, g->msa_len * sizeof(u2i));
	}
	for(i=0;i<g->seqs->nseq;i++){
		if(i == ridx) continue;
		r = ref_u1v(g->msa, g->msa_len * i);
		for(j=0;j<(int)g->msa_len;j++){
			bcnts[r[j]][j] ++;
		}
	}
	for(i=0;i<g->msa_len;i++){
		max = 0;
		for(j=1;j<4;j++){
			if(bcnts[j][i] > bcnts[max][i]){
				max = j;
			}
		}
		bcnts[5][i] = max;
		bcnts[6][i] = bcnts[0][i] + bcnts[1][i] + bcnts[2][i] + bcnts[3][i];
	}
	s = ref_u1v(g->msa, g->msa_len * ridx);
	winl = 20;
	winx = 0;
	winy = 0;
	beg = -1;
	memset(wins, 0, winl * sizeof(int));
	if(0){
		for(i=0;i<g->msa_len;i++){
			max = bcnts[5][i];
			winy -= wins[winx];
			if(s[i] < 4 && ((hcovs[i] >= 4 && bcnts[s[i]][i] == 0))){
				wins[winx] = 1;
			} else {
				wins[winx] = 0;
			}
			winy += wins[winx];
			//fprintf(stderr, " -- hcovs = %d s[i] = %c:%d i = %u winy = %d in %s -- %s:%d --\n", hcovs[i], "ACGT-"[s[i]], bcnts[s[i]][i], i, winy, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			if(Int(i) > winl){
				if(winy >= Int(0.9 * winl)){
					if(beg < 0){
						beg = i - winl;
					}
				} else {
					if(beg >= 0){
						end = i;
						if(cns_debug > 1){
							fprintf(stderr, " -- remove low quality fragment %d - %d: ", beg, end);
							for(j=beg;j<end;j++){
								fputc("ACGT-"[s[j]], stderr);
							}
							fputc('\n', stderr);
							fflush(stderr);
						}
						for(j=beg;j<end;j++){
							s[j] = 4;
						}
					}
					beg = -1;
				}
			}
			winx = (winx + 1) % winl;
		}
	}
	clear_and_encap_b2v(g->btds, (g->msa_len + 1) * W * 2);
	bts = (u2i*)g->btds->buffer;
	clear_and_encap_u4v(g->btxs, (g->msa_len + 1) * W * 2);
	dps = g->btxs->buffer;
	memset(bts, 0, 2 * W * sizeof(u2i));
	memset(dps, 0, 2 * W * sizeof(u4i));
	bts += W * 2;
	dps += W * 2;
	for(i=0;i<g->msa_len;i++){
		bs1 = bts + 2 * W * (Int(i) - 1);
		bs  = bts + 2 * W * i;
		row1 = dps + 2 * W * (Int(i) - 1);
		row2 = dps + 2 * W * i;
		f = 0;
		off = Int(i) - W;
		if(Int(i) <= W){
			beg = W - i;
		} else {
			beg = 0;
		}
		if(i + W - 1 >= g->msa_len){
			end = W + g->msa_len - i;
		} else {
			end = W * 2 - 1;
		}
		if(s[i] == 4){
			for(j=beg;j<end;j++){
				e = row1[j + 1];
				//h = row1[j] + 0 + 1;
				h = row1[j];
				if(h >= e){
					bs[j] = 0;
				} else {
					h = e;
					bs[j] = 1;
				}
				if(h >= f){
					f = h;
				} else {
					h = f;
					bs[j] = 2;
				}
				row2[j] = h;
			}
		} else {
			for(j=beg;j<end;j++){
				e = row1[j + 1];
				//h = row1[j] + (s[i] == bcnts[5][off + j]? bcnts[s[i]][off + j] : 0);
				if(Int(bcnts[s[i]][off + j] + 1) >= Int(bcnts[6][off + j] - bcnts[s[i]][off + j])){
					h = row1[j] + 2 * bcnts[s[i]][off + j] - bcnts[6][off + j] + 1;
				} else {
					h = row1[j] + 1;
				}
				// bonus for putting homo together
				if(i && bs1[j] == 0 && s[i] == s[i - 1]){
					h += 1;
				}
				if(h >= e){
					bs[j] = 0;
				} else {
					h = e;
					bs[j] = 1;
				}
				if(h >= f){
					f = h;
				} else {
					h = f;
					bs[j] = 2;
				}
				row2[j] = h;
			}
		}
		if(beg){
			row2[beg - 1] = 0;
			bs[beg - 1]   = 0;
		}
		row2[end] = 0;
		bs[end] = 0;
		if(cns_debug > 3){
			fprintf(stderr, "ROW[%u] '%c' %03d-%03d:", i, "ACGT-"[s[i]], beg, end);
			for(j=beg;j<end;j++){
				fprintf(stderr, " %5u[%d]", row2[j], bs[j]);
			}
			fprintf(stderr, "\n");
		}
	}
	row2 = dps + 2 * W * (g->msa_len - 1);
	max = row2[W];
	r = ref_u1v(g->msa, g->msa_len * g->seqs->nseq);
	x = y = g->msa_len - 1;
	while(y >= 0 && x >= 0){
		off = x - W;
		bs  = bts + 2 * W * x;
		bt  = bs[y - off];
		switch(bt){
			case 0: r[y] = s[x]; x --; y --; break;
			case 1: x --; break;
			case 2: r[y] = 4; y --; break;
		}
		if(cns_debug > 3 && bt){
			fprintf(stderr, " -- x = %d y = %d bt = %d in %s -- %s:%d --\n", x, y, bt, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
	}
	if(cns_debug > 1){
		fprintf(stderr, " -- max = %d x = %d y = %d in %s -- %s:%d --\n", max, x, y, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	}
	while(y >= 0){
		r[y] = 4;
		y --;
	}
	return max;
}

static inline void realign_msa_pog(POG *g){
	u4i ridx;
	if(g->par->rW <= 0 || g->seqs->nseq < 3) return;
	if(cns_debug > 1){
		fprintf(stderr, "RAW MSA\n");
		print_msa_pog(g, stderr);
	}
	for(ridx=0;ridx<g->seqs->nseq;ridx++){
		if(cns_debug > 1){
			fprintf(stderr, "REALIGN[%u]\n", ridx);
		}
		realign_msa_pog_core(g, ridx, g->par->rW);
		if(cns_debug > 1){
			print_msa_pog(g, stderr);
		}
		// update ridx alignment
		memcpy(g->msa->buffer + g->msa_len * ridx, g->msa->buffer + g->msa_len * g->seqs->nseq, g->msa_len);
	}
}

static const float homo_merging_cmp_norm[20] = {
	0.95, 0.90, 0.80, 0.75, 0.70,
	0.70, 0.65, 0.60, 0.60, 0.55,
	0.50, 0.45, 0.40, 0.30, 0.20,
	0.20, 0.20, 0.20, 0.20, 0.20
};

/*
static inline void dp_gen_cns_pog(POG *g){
	pog_cns_t dps[2 * 5], *dp1, *dp2;
	b2i *row1, *row2, *row3;
	u1i *s, *r;
	u4i mcnt;
	int bcnts[5];
	u2i *hcovs, *rexs;
	u4i rid, beg, end, x, val, fidx, col, dep;
	u4i mx, mi, i, j;
	float DP_MIN, score0, score1;
	DP_MIN = - (MAX_B4 >> 1);
	mcnt = num_min(g->par->msa_min_cnt, g->seqs->nseq);
	clear_and_encap_u1v(g->cbts, g->msa_len * 5);
	// scan read end
	clear_u4v(g->btxs);
	clear_u2v(g->btvs); // whether read cov current column
	hcovs = g->hcovs->buffer;
	memset(hcovs, 0, g->msa_len * sizeof(u2i));
	for(rid=0;rid<g->seqs->nseq;rid++){
		push_u2v(g->btvs, MAX_U2);
		r = ref_u1v(g->msa, g->msa_len * rid);
		beg = 0;
		while(beg < g->msa_len && r[beg] == 4) beg ++;
		push_u4v(g->btxs, (0U << 31) | (rid << 16) | beg);
		end = g->msa_len;
		while(end > beg && r[end - 1] == 4) end --;
		push_u4v(g->btxs, (1U << 31) | (rid << 16) | end);
		for(i=beg;i<end;i++) hcovs[i] ++;
	}
	sort_array(g->btxs->buffer, g->btxs->size, u4i, num_cmpgt(a & 0xFFFF, b & 0xFFFF));
	mx = 0;
	rexs = g->btvs->buffer;
	memset(dps, 0, 2 * 5 * sizeof(pog_cns_t));
	fidx = 0;
	clear_b2v(g->rows);
	clear_u4v(g->rowr);
	dp1 = dps;
	inc_b2v(g->rows, g->par->cW * g->seqs->nseq);
	for(i=0;i<5;i++){
		dp1[i].coff = g->rows->size;
		inc_b2v(g->rows, g->par->cW * g->seqs->nseq);
		row1 = g->rows->buffer + dp1[i].coff;
		memset(row1, 0, g->par->cW * g->seqs->nseq * sizeof(b2i));
	}
	dep = 0;
	for(col=0;col<g->msa_len;col++){
		fidx = !fidx;
		dp1 = dps + fidx * 5;
		dp2 = dps + (!fidx) * 5;
		// collect read coverage
		while(mx < g->btxs->size){
			val = g->btxs->buffer[mx];
			if((val & 0xFFFF) <= col){
				rid = (val << 1) >> 17;
				if(val & (1U << 31)){
#if DEBUG
					if(rexs[rid] == MAX_U2){
						fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
#endif
					rexs[rid] = MAX_U2;
					dep --;
				} else {
#if DEBUG
					if(rexs[rid] < MAX_U2){
						fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
#endif
					rexs[rid] = 0;
					dep ++;
				}
				mx ++;
			} else {
				break;
			}
		}
		for(i=0;i<5;i++){
			dp2[i].max = DP_MIN;
			dp2[i].bt  = 4;
			if(g->rowr->size){
				dp2[i].coff = pop_u4v(g->rowr);
			} else {
				dp2.coff = g->rows->size;
				inc_b2v(g->rows, g->par->cW * g->seqs->nseq);
			}
		}
		row2 = g->rows->buffer; // temprary
		if(dep < mcnt && !(g->par->refmode && rexs[0] != MAX_U2)){
			//TODO: coding from here
			mi = 4;
			for(i=0;i<4;i++){
				if(dp1[i] > dp1[mi]){
					mi = i;
				}
			}
			dp2[4] = dp1[mi];
			bts[4] = mi;
		} else {
			for(i=0;i<=4;i++){
				score0 = dp1[i];
				for(j=0;j<=4;j++){
					score1 = score0;
					for(x=0;x<=4;x++){
						if(j == 4){
							if(x != 4){
								score1 += bcnts[x] * g->par->I;
							}
						} else {
							if(x == 4){
								if(i == j){
									score1 += bcnts[x] * g->par->H;
								} else {
									score1 += bcnts[x] * g->par->D;
								}
							} else if(x == j){
								score1 += bcnts[x] * g->par->M;
							} else {
								score1 += bcnts[x] * g->par->X;
							}
						}
					}
					if(score1 > dp2[j]){
						dp2[j] = score1;
						bts[j] = i;
					}
				}
			}
			// Prior to the first read
			for(rid=0;rid<g->seqs->nseq;rid++){
				if(rexs[i] == 0){
					continue;
				}
				r = ref_u1v(g->msa, g->msa_len * rid);
				dp2[r[col]] += 0.1;
				break;
			}
		}
		for(i=0;i<=4;i++){
			push_u1v(g->cbts, bts[i]);
		}
		if(cns_debug > 2){
			fprintf(stderr, "[DPCNS%03u]\t[%2d, %2d, %2d, %2d, %2d]", col, bcnts[0], bcnts[1], bcnts[2], bcnts[3], bcnts[4]);
			//fprintf(stderr, "\t[%0.1f, %0.1f, %0.1f, %0.1f, %0.1f]", dp1[0], dp1[1], dp1[2], dp1[3], dp1[4]);
			fprintf(stderr, "\t[%0.1f:%d, %0.1f:%d, %0.1f:%d, %0.1f:%d, %0.1f:%d]\n", dp2[0], bts[0], dp2[1], bts[1], dp2[2], bts[2], dp2[3], bts[3], dp2[4], bts[4]);
		}
	}
	// Backtrace consensus sequence
	dp2 = dps + (!fidx) * 5;
	mx = 4;
	for(i=0;i<4;i++){
		if(dp2[i] > dp2[mx]){
			mx = i;
		}
	}
	x = g->msa_len;
	s = ref_u1v(g->msa, g->msa_len * g->seqs->nseq);
	while(x){
		x --;
		s[x] = mx;
		mx = g->cbts->buffer[x * 5 + mx];
	}
	for(i=0;i<g->msa_len;i++){
		if(s[i] < 4){
			bit2basebank(g->cns, s[i]);
		}
	}
}
*/

static inline void dp_call_cns_pog(POG *g){
	u1i *s, *r, bts[5];
	u4i mcnt;
	int bcnts[5];
	u2i *hcovs, *rexs;
	u4i rid, beg, end, x, val, fidx, col, dep;
	u4i mx, mi, i, j;
	float dps[2 * 5], *dp1, *dp2, DP_MIN, score0, score1;
	DP_MIN = - (MAX_B4 >> 1);
	mcnt = num_min(g->par->msa_min_cnt, g->seqs->nseq);
	clear_and_encap_u1v(g->cbts, g->msa_len * 5);
	// scan read end
	clear_u4v(g->btxs);
	clear_u2v(g->btvs); // whether read cov current column
	hcovs = g->hcovs->buffer;
	memset(hcovs, 0, g->msa_len * sizeof(u2i));
	for(rid=0;rid<g->seqs->nseq;rid++){
		push_u2v(g->btvs, 0);
		r = ref_u1v(g->msa, g->msa_len * rid);
		beg = 0;
		while(beg < g->msa_len && r[beg] == 4) beg ++;
		push_u4v(g->btxs, (0U << 31) | (rid << 16) | beg);
		end = g->msa_len;
		while(end > beg && r[end - 1] == 4) end --;
		push_u4v(g->btxs, (1U << 31) | (rid << 16) | end);
		for(i=beg;i<end;i++) hcovs[i] ++;
	}
	sort_array(g->btxs->buffer, g->btxs->size, u4i, num_cmpgt(a & 0xFFFF, b & 0xFFFF));
	memset(dps, 0, 2 * 5 * sizeof(float));
	fidx = 0;
	dep = 0;
	mx = 0;
	rexs = g->btvs->buffer;
	for(col=0;col<g->msa_len;col++){
		fidx = !fidx;
		dp1 = dps + fidx * 5;
		dp2 = dps + (!fidx) * 5;
		// collect read coverage
		while(mx < g->btxs->size){
			val = g->btxs->buffer[mx];
			if((val & 0xFFFF) <= col){
				rid = (val << 1) >> 17;
				if(val & (1U << 31)){
#if DEBUG
					if(rexs[rid] == 0){
						fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
#endif
					rexs[rid] = 0;
					dep --;
				} else {
#if DEBUG
					if(rexs[rid] == 1){
						fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						abort();
					}
#endif
					rexs[rid] = 1;
					dep ++;
				}
				mx ++;
			} else {
				break;
			}
		}
		bts[0] = bts[1] = bts[2] = bts[3] = bts[4] = 4;
		dp2[0] = DP_MIN;
		dp2[1] = DP_MIN;
		dp2[2] = DP_MIN;
		dp2[3] = DP_MIN;
		dp2[4] = DP_MIN;
		memset(bcnts, 0, sizeof(int) * 5);
		for(rid=0;rid<g->seqs->nseq;rid++){
			if(rexs[rid] == 0){
				continue;
			}
			r = ref_u1v(g->msa, g->msa_len * rid);
			bcnts[r[col]] ++;
		}
		if(dep < mcnt && !(g->par->refmode && rexs[0])){
			mi = 4;
			for(i=0;i<4;i++){
				if(dp1[i] > dp1[mi]){
					mi = i;
				}
			}
			dp2[4] = dp1[mi];
			bts[4] = mi;
		} else {
			for(i=0;i<=4;i++){
				score0 = dp1[i];
				for(j=0;j<=4;j++){
					score1 = score0;
					for(x=0;x<=4;x++){
						if(j == 4){
							if(x != 4){
								score1 += bcnts[x] * g->par->I;
							}
						} else {
							if(x == 4){
								if(i == j){
									score1 += bcnts[x] * g->par->H;
								} else {
									score1 += bcnts[x] * g->par->D;
								}
							} else if(x == j){
								score1 += bcnts[x] * g->par->M;
							} else {
								score1 += bcnts[x] * g->par->X;
							}
						}
					}
					if(score1 > dp2[j]){
						dp2[j] = score1;
						bts[j] = i;
					}
				}
			}
			// Prior to the first read
			for(rid=0;rid<g->seqs->nseq;rid++){
				if(rexs[i] == 0){
					continue;
				}
				r = ref_u1v(g->msa, g->msa_len * rid);
				dp2[r[col]] += 0.1;
				break;
			}
		}
		for(i=0;i<=4;i++){
			push_u1v(g->cbts, bts[i]);
		}
		if(cns_debug > 2){
			fprintf(stderr, "[DPCNS%03u]\t[%2d, %2d, %2d, %2d, %2d]", col, bcnts[0], bcnts[1], bcnts[2], bcnts[3], bcnts[4]);
			//fprintf(stderr, "\t[%0.1f, %0.1f, %0.1f, %0.1f, %0.1f]", dp1[0], dp1[1], dp1[2], dp1[3], dp1[4]);
			fprintf(stderr, "\t[%0.1f:%d, %0.1f:%d, %0.1f:%d, %0.1f:%d, %0.1f:%d]\n", dp2[0], bts[0], dp2[1], bts[1], dp2[2], bts[2], dp2[3], bts[3], dp2[4], bts[4]);
		}
	}
	// Backtrace consensus sequence
	dp2 = dps + (!fidx) * 5;
	mx = 4;
	for(i=0;i<4;i++){
		if(dp2[i] > dp2[mx]){
			mx = i;
		}
	}
	x = g->msa_len;
	s = ref_u1v(g->msa, g->msa_len * g->seqs->nseq);
	while(x){
		x --;
		s[x] = mx;
		mx = g->cbts->buffer[x * 5 + mx];
	}
	for(i=0;i<g->msa_len;i++){
		if(s[i] < 4){
			bit2basebank(g->cns, s[i]);
		}
	}
}

static inline void gen_cns_pog(POG *g){
	u1i *s, *r;
	u4i idx, lst, lsv, lsx, rid, max, i, mcnt, min_cnt, max_cnt, runlen, cnl, beg, end, cbeg, cend;
	u4i freqs[5][11], fmax1, fmax2, tot, cut, *hcnts;
	u2i cnts[5], *bls, cl, bx, bl, bm, dl, dm, *hcovs, corr;
	//u2i *vsts;
	float fl, fx1, fx2, norm, flush;
	mcnt = num_min(g->par->msa_min_cnt, g->seqs->nseq);
	min_cnt = num_max(mcnt, UInt(g->par->msa_min_freq * g->seqs->nseq));
	max_cnt = g->seqs->nseq - min_cnt;
	s = ref_u1v(g->msa, g->msa_len * g->seqs->nseq);
	memset(s, 4, g->msa_len);
	clear_basebank(g->cns);
	hcovs = g->hcovs->buffer;
	bls   = hcovs + 8 + g->msa_len; // bls is all zeros, see end_pog
	if(g->par->rW){
		memset(hcovs, 0, g->msa_len * sizeof(u2i));
		for(rid=0;rid<g->seqs->nseq;rid++){
			r = ref_u1v(g->msa, g->msa_len * rid);
			beg = 0;
			while(beg < g->msa_len && r[beg] == 4) beg ++;
			end = g->msa_len;
			while(end > beg && r[end - 1] == 4) end --;
			for(i=beg;i<end;i++) hcovs[i] ++;
		}
	}
	// revise mcnt
	if(g->par->refmode){
	} else if(0){
		tot = 0;
		clear_and_encap_u4v(g->btxs, mcnt + 1);
		hcnts = g->btxs->buffer;
		//hcnts = alloca((mcnt + 1) * sizeof(u4i));
		memset(hcnts, 0, (mcnt + 1) * sizeof(u4i));
		for(i=0;i<g->msa_len;i++){
			tot += (hcovs[i] + 1) * (hcovs[i]) / 2;
			hcnts[num_min(hcovs[i], mcnt)] += (hcovs[i] + 1) * (hcovs[i]) / 2;
		}
		cut = tot * 0.8;
		tot = 0;
		for(i=mcnt;i>1;i--){
			tot += hcnts[i];
			if(tot >= cut){
				break;
			}
		}
		if(i < mcnt){
			if(cns_debug > 1){
				fprintf(stderr, " Revise mcnt %u -> %u\n", mcnt, i); fflush(stderr);
			}
			mcnt = i;
			min_cnt = i;
			max_cnt = g->seqs->nseq - min_cnt;
		}
	}
	// gen cns
	cbeg = 0;
	cend = g->msa_len;
	if(g->par->refmode){
		rid = 0;
		{
			r = ref_u1v(g->msa, g->msa_len * rid);
			cbeg = 0;
			while(cbeg < g->msa_len && r[cbeg] == 4){ hcovs[cbeg] = 0; cbeg ++; }
			end = g->msa_len;
			while(cend > cbeg && r[cend - 1] == 4){ cend --; hcovs[cend] = 0; }
		}
	} else {
		cbeg = 0;
		while(cbeg < g->msa_len && hcovs[cbeg] < min_cnt) cbeg ++;
		cend = g->msa_len;
		while(cend > cbeg && hcovs[cend - 1] < min_cnt) cend --;
	}
	fmax1 = 5;
	fmax2 = 10;
	if(cns_debug > 1){
		for(i=0;i<fmax1;i++){
			memset(freqs[i], 0, (fmax2 + 1) * sizeof(u4i));
		}
	}
	lst = lsv = MAX_U4;
	lsx = 0;
	runlen = 0;
	cl = 0;
	cnl = 0;
	end = 0;
	flush = 0;
	for(idx=cbeg;idx<=cend;idx++){
		if(idx < cend){
			flush = 0;
			memset(cnts, 0, 5 * 2);
			if(g->par->refmode){
				if((max = g->msa->buffer[g->msa_len * 0 + idx]) == 4){
					max = 0;
				} else {
					if(hcovs[idx] < min_cnt){
						flush = 1;
					}
				}
			} else {
				max = 0;
			}
			for(rid=0;rid<g->seqs->nseq;rid++){
				cnts[g->msa->buffer[g->msa_len * rid + idx]] ++;
			}
			for(i=0;i<4;i++){
				if(cnts[i] > cnts[max]){
					max = i;
				}
			}
			if(flush){
				s[idx] = max;
				end = idx;
			} else if(hcovs[idx] >= min_cnt && cnts[max] >= mcnt && cnts[4] - (g->seqs->nseq - hcovs[idx]) < (g->seqs->nseq - cnts[4])){
				s[idx] = max;
				end = idx;
				flush = 1;
			}
		} else {
			if(end + 1 == idx){
				end = idx;
			}
			max = 4;
			flush = 1;
		}
		if(flush){
			if(lst == MAX_U4){ lst = idx; lsv = idx; }
			if(idx == cend || s[lst] != max){
				runlen = 0;
				for(rid=0;rid<g->seqs->nseq;rid++){
					r = ref_u1v(g->msa, g->msa_len * rid);
					bl = 0;
					for(i=lsv;i<end;i++){
						if(r[i] == s[lst]) bl ++;
					}
					bls[rid] = bl;
					runlen += bl;
				}
				sort_array(bls, g->seqs->nseq, u2i, num_cmpgt(b, a));
				bx = 1;
				bl = bls[0];
				bm = 1;
				dl = 0;
				dm = 0;
				for(rid=1;rid<=hcovs[lst];rid++){
					if(rid < hcovs[lst] && bls[rid] == bls[rid - 1]){
						bx ++;
					} else {
						if(cns_debug > 1 && cl <= fmax1 && bls[rid - 1] <= fmax2){
							freqs[cl - 1][bls[rid - 1]] += bx;
						}
						if(bx > bm){
							if(dm < bm){
								dm = bm;
								dl = bl;
							}
							bm = bx;
							bl = bls[rid - 1];
						} else if(bx > dm){
							dm = bx;
							dl = bls[rid - 1];
						}
						bx = 1;
					}
				}
				fl = runlen * 1.0 / hcovs[lst];
				if(cns_debug > 1){
					fprintf(stderr, "HOMO[%3u] %4u\t%c\t%u\t%u:%u\t%u:%u\t%u\t%u\t%0.2f", lst, cnl, bit_base_table[s[lst]], cl, bl, bm, dl, dm, runlen, hcovs[lst], fl);
					if(cl != bl){
						fprintf(stderr, " *\n");
					} else {
						fprintf(stderr, "\n");
					}
				}
				if(cl != bl){
					if(bl > cl + 1 && bl > Int(fl) + 1) bl = Int(fl) + 1; // Don't correct too much
					if(bl + 1 < cl && (float)bl < fl) bl = cl - 1;
					if(bl < cl && cl == dl && dm >= Int(0.8 * bm)){ // prefer increase runlen
						norm = 1.0 / homo_merging_cmp_norm[num_min(cl, 19)];
					} else {
						norm = homo_merging_cmp_norm[num_min(bl, 19)];
					}
					fx1 = num_abs(cl - fl);
					fx2 = num_abs(bl - fl);
					if(norm * fx2 <= fx1){
						corr = bl;
						if(cns_debug > 1){
							fprintf(stderr, "CORR[%03u] %4u\t%c\t%u->%u\n", lst, cnl, bit_base_table[s[lst]], cl, bl);
						}
					} else {
						corr = cl;
					}
				} else {
					corr = cl;
				}
				cnl += cl;
				for(i=0;i<corr;i++){
					bit2basebank(g->cns, s[lst]);
					if(i){
						s[lst + i] = s[lst];
					}
				}
				for(i+=lst;i<end;i++){
					s[i] = 4;
				}
				runlen = 0;
				lst = idx;
				lsv = lsx;
				lsx = idx;
				cl = 1;
			} else {
				cl ++;
				lsx = idx;
			}
		}
	}
	if(cns_debug > 1){
		for(cl=0;cl<fmax1;cl++){
			fprintf(stderr, "RUNLEN[%u]", cl + 1);
			for(bl=0;bl<=fmax2;bl++){
				fprintf(stderr, "%10u,", freqs[cl][bl]);
			}
			fprintf(stderr, "\n");
		}
	}
}

static inline void check_dup_edges_pog(POG *g){
	pog_node_t *u;
	pog_edge_t *e;
	u32hash *hash;
	u4i nidx, eidx, *t;
	int exists;
	hash = init_u32hash(13);
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_pognodev(g->nodes, nidx);
		clear_u32hash(hash);
		eidx = u->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			t = prepare_u32hash(hash, e->node, &exists);
			if(exists){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			} else {
				*t = e->node;
			}
		}
	}
	free_u32hash(hash);
}

static inline void end_pog(POG *g){
	pog_node_t *u, *v, *x;
	pog_edge_t *e;
	u1i *r;
	u4i i, ridx, nidx, eidx, xidx, moff, ready, beg, end, reflen, margin;
	int score;
	clear_basebank(g->cns);
	if(g->seqs->nseq == 0) return;
	score = align_rd_pog(g, 0);
	if(g->par->refmode){
		// add edges to refbeg and refend, to make sure reads can be aligned quickly and correctly within small bandwidth
		u = ref_pognodev(g->nodes, POG_HEAD_NODE);
		reflen = g->seqs->rdlens->buffer[0];
		margin = 5;
		for(i=0;i<g->sbegs->size;i++){
			if(g->sbegs->buffer[i] < margin || g->sbegs->buffer[i] + margin >= reflen) continue;
			nidx = POG_TAIL_NODE + 1 + (reflen - 1 - g->sbegs->buffer[i]);
			v = ref_pognodev(g->nodes, nidx);
			e = add_edge_pog(g, u, v, 0, 1);
		}
		v = ref_pognodev(g->nodes, POG_TAIL_NODE);
		for(i=0;i<g->sends->size;i++){
			if(g->sends->buffer[i] < margin || g->sends->buffer[i] + margin >= reflen) continue;
			nidx = POG_TAIL_NODE + 1 + (reflen - g->sends->buffer[i]);
			u = ref_pognodev(g->nodes, nidx);
			e = add_edge_pog(g, u, v, 0, 1);
		}
	}
	for(ridx=1;ridx<g->seqs->nseq;ridx++){
		score = align_rd_pog(g, ridx);
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_pognodev(g->nodes, nidx);
		u->vst = 0;
		u->aux = 0;
		u->coff = 0;
		u->erev = 0;
	}
	// calcuate msa_len
	clear_u4v(g->stack);
	push_u4v(g->stack, POG_HEAD_NODE);
	nidx = POG_HEAD_NODE;
	while(pop_u4v(g->stack, &nidx)){
		u = ref_pognodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_pognodev(g->nodes, e->node);
			if(u->coff + 1 > v->coff){
				v->coff = u->coff + 1;
			}
			v->vst ++;
			if(v->vst > v->nin){
				print_vstdot_pog(g, "1.dot");
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				check_dup_edges_pog(g);
				abort();
			}
		}
		eidx = u->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_pognodev(g->nodes, e->node);
			if(v->aux) continue; // already pushed
			if(v->vst > v->nin){
				print_vstdot_pog(g, "1.dot");
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				check_dup_edges_pog(g);
				abort();
			}
			if(v->vst == v->nin){
				ready = 1;
				{
					xidx = v->aligned;
					moff = v->coff;
					while(xidx != e->node){
						x = ref_pognodev(g->nodes, xidx);
						if(x->nin > x->vst){
							ready = 0;
							break;
						}
						if(x->coff > moff){
							moff = x->coff;
						}
						xidx = x->aligned;
					}
				}
				if(ready){
					v->coff = moff;
					v->aux  = 1;
					push_u4v(g->stack, e->node);
					xidx = v->aligned;
					while(xidx != e->node){
						x = ref_pognodev(g->nodes, xidx);
						x->coff = moff;
						if(x->edge){
							push_u4v(g->stack, xidx);
							x->aux = 1;
						}
						xidx = x->aligned;
					}
				}
			}
		}
	}
	if(nidx != POG_TAIL_NODE){
		fprint_dot_pog(g, "1.dot", NULL);
		print_seqs_pog(g, "1.seqs.fa", NULL);
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	// generate MSA
	u = ref_pognodev(g->nodes, POG_TAIL_NODE);
	g->msa_len = u->coff;
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_pognodev(g->nodes, nidx);
		u->vst = 0;
	}
	clear_and_encap_u1v(g->msa, g->msa_len * (g->seqs->nseq + 2));
	g->msa->size = g->msa_len * (g->seqs->nseq + 1);
	memset(g->msa->buffer, 4, g->msa->size);
	clear_u4v(g->stack);
	push_u4v(g->stack, POG_HEAD_NODE);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_pognodev(g->nodes, nidx);
		eidx = u->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_pognodev(g->nodes, e->node);
			v->vst ++;
			if(v->vst >= v->nin){
				ready = 1;
				xidx = v->aligned;
				while(xidx != e->node){
					x = ref_pognodev(g->nodes, xidx);
					if(x->nin > x->vst){
						ready = 0;
						break;
					}
					xidx = x->aligned;
				}
				if(ready){
					xidx = e->node;
					do {
						x = ref_pognodev(g->nodes, xidx);
						if(xidx != POG_TAIL_NODE){
							set_u1v(g->msa, x->rid * g->msa_len + x->coff, (x->base) & 0x03);
						}
						if(x->edge){
							push_u4v(g->stack, xidx);
						}
						xidx = x->aligned;
					} while(xidx != e->node);
				}
			}
		}
	}
	clear_and_encap_u2v(g->hcovs, g->msa_len + 8 + g->seqs->nseq);
	memset(g->hcovs->buffer, 0, (g->msa_len + 8 + g->seqs->nseq) * sizeof(u2i));
	for(ridx=0;ridx<g->seqs->nseq;ridx++){
		r = ref_u1v(g->msa, g->msa_len * ridx);
		beg = 0;
		while(beg < g->msa_len && r[beg] == 4) beg ++;
		end = g->msa_len;
		while(end && r[end - 1] == 4) end --;
		for(i=beg;i<end;i++) g->hcovs->buffer[i] ++;
	}
	// realignment
	realign_msa_pog(g);
	//// realignment again
	//realign_msa_pog(g);
	// gen consensus line
	if(g->par->cnsmode == 1){
		dp_call_cns_pog(g);
	} else {
		gen_cns_pog(g);
	}
	if(cns_debug > 1){
		print_msa_pog(g, stderr);
	}
	if(0){
		fprintf(stderr, " -- seqs\t%llu --\n", (u8i)g->seqs->rdseqs->cap); fflush(stderr);
		fprintf(stderr, " -- nodes\t%llu --\n", (u8i)g->nodes->cap); fflush(stderr);
		fprintf(stderr, " -- edges\t%llu --\n", (u8i)g->edges->cap); fflush(stderr);
		fprintf(stderr, " -- rows\t%llu/%llu --\n", (u8i)g->rows->size, (u8i)g->rows->cap); fflush(stderr);
		fprintf(stderr, " -- btds\t%llu/%llu --\n", (u8i)g->rows->size, (u8i)g->btds->cap); fflush(stderr);
		fprintf(stderr, " -- btxs\t%llu/%llu --\n", (u8i)g->btxs->size, (u8i)g->btxs->cap); fflush(stderr);
		fprintf(stderr, " -- cns\t%llu --\n", (u8i)g->cns->cap); fflush(stderr);
	}
}

#endif
