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
#include "string.h"
#include "list.h"
#include "hashset.h"
#include <emmintrin.h>

static int cns_debug = 0;

#define POG_RDLEN_MAX	MAX_B2
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
	u4i rid:15, pos:14, base:3;
	u2i nin, vst;
	u4i edge;
	u4i aligned;
	u4i btx;
	u4i coff;
} pog_node_t;
define_list(pognodev, pog_node_t);

typedef struct {
	u4i node;
	u4i cov;
	u4i next;
} pog_edge_t;
define_list(pogedgev, pog_edge_t);

typedef struct {
	SeqBank  *seqs;
	pognodev *nodes;
	pogedgev *edges;
	int sse, W_score, near_dialog, aln_mode;
	int M, X, I, D, T, W, rW;
	u4i msa_min_cnt;
	float msa_min_freq;
	b2v *qprof;
	b2v *rows;
	b2v *btds;
	u4v *btxs;
	u4v *stack;
	u4i msa_len;
	u1v *msa;
	u2v *bcnts[7];
	u2v *hcovs;
	BaseBank *cns;
} POG;

static inline POG* init_pog(int M, int X, int I, int D, int W, int Wscore_cutoff, int use_sse, int rW, u4i min_cnt, float min_freq){
	POG *g;
	g = malloc(sizeof(POG));
	g->seqs = init_seqbank();
	g->nodes = init_pognodev(1024);
	g->edges = init_pogedgev(1024);
	g->sse = use_sse;
	g->W_score = Wscore_cutoff;
	g->near_dialog = 0;
	g->aln_mode = POG_ALNMODE_OVERLAP; // 0: overlap, 1: global
	g->W = W;
	g->M = M;
	g->X = X;
	g->I = I;
	g->D = D;
	g->T = 10 * M;
	g->rW = rW;
	g->msa_min_cnt = min_cnt;
	g->msa_min_freq = min_freq;
	g->qprof = init_b2v(1024);
	g->rows  = init_b2v(1024);
	g->btds  = init_b2v(1024);
	g->btxs  = init_u4v(1024);
	g->stack = init_u4v(32);
	g->msa_len = 0;
	g->msa = init_u1v(1024);
	g->bcnts[0] = init_u2v(1024);
	g->bcnts[1] = init_u2v(1024);
	g->bcnts[2] = init_u2v(1024);
	g->bcnts[3] = init_u2v(1024);
	g->bcnts[4] = init_u2v(1024);
	g->bcnts[5] = init_u2v(1024);
	g->bcnts[6] = init_u2v(1024);
	g->hcovs = init_u2v(1024);
	g->cns = init_basebank();
	return g;
}

static inline void free_pog(POG *g){
	free_seqbank(g->seqs);
	free_pognodev(g->nodes);
	free_pogedgev(g->edges);
	free_b2v(g->qprof);
	free_b2v(g->rows);
	free_b2v(g->btds);
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
	free_basebank(g->cns);
	free(g);
}

static inline void print_dot_pog(POG *g, FILE *out){
	pog_node_t *n;
	pog_edge_t *e;
	u4i nidx, eidx;
	fprintf(out, "digraph {\n");
	//fprintf(out, "rankdir=LR\n");
	fprintf(out, "N0 [label=\"BEG\"]\n");
	fprintf(out, "N1 [label=\"END\"]\n");
	for(nidx=POG_TAIL_NODE+1;nidx<g->nodes->size;nidx++){
		n = ref_pognodev(g->nodes, nidx);
		fprintf(out, "N%u [label=R%u_%u_%c]\n", nidx, n->rid, n->pos, bit_base_table[n->base]);
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		n = ref_pognodev(g->nodes, nidx);
		if(n->aligned != nidx){
			fprintf(out, "N%u -> N%u [color=magenta style=dashed]\n", nidx, n->aligned);
		}
		eidx = n->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			if(nidx == POG_HEAD_NODE && e->node == POG_TAIL_NODE) continue;
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
		fprintf(out, "N%u [label=\"N%u:%u:%u:%d\"]\n", nidx, nidx, n->btx, n->nin, n->vst);
		eidx = n->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			if(nidx == POG_HEAD_NODE && e->node == POG_TAIL_NODE) continue;
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

static inline void beg_pog(POG *g){
	pog_node_t *head, *tail;
	pog_edge_t *e;
	clear_seqbank(g->seqs);
	clear_pognodev(g->nodes);
	clear_pogedgev(g->edges);
	if(1){
		renew_b2v(g->rows, 1024);
		renew_b2v(g->btds, 1024);
		renew_u4v(g->btxs, 1024);
	} else {
		clear_b2v(g->rows);
		clear_b2v(g->btds);
		clear_u4v(g->btxs);
	}
	clear_basebank(g->cns);
	head = next_ref_pognodev(g->nodes);
	ZEROS(head);
	tail = next_ref_pognodev(g->nodes);
	ZEROS(tail);
	tail->base = 4;
	tail->aligned = POG_TAIL_NODE;
	e = next_ref_pogedgev(g->edges);
	ZEROS(e);
	e = next_ref_pogedgev(g->edges);
	e->node = POG_TAIL_NODE;
	e->cov  = 1;
	e->next = 0;
	tail->nin ++;
	head->edge = offset_pogedgev(g->edges, e);
}

static inline void v1_sse_band_row_rdaln_pog(POG *g, u4i nidx1, u4i nidx2, u4i seqlen, u4i coff1, u4i coff2, b2i *qp){
	__m128i D, H, E, S, BT_MASK, BT, MIN, MAX;
	b2i *row1, *row2, *btds;
	u4i i, slen, seqlex, beg, end, *btxs;
	int f, mi, center;
	slen = (seqlen + 7) / 8;
	seqlex = slen * 8;
	D = _mm_set1_epi16(g->D);
	MIN = _mm_set1_epi16(POG_SCORE_MIN);
	BT_MASK = _mm_set1_epi16(0b10);
	row1 = ref_b2v(g->rows, coff1);
	row2 = ref_b2v(g->rows, coff2);
	btds = ref_b2v(g->btds, coff2);
	btxs = ref_u4v(g->btxs, coff2);
	if(row1[row1[seqlex]] >= g->W_score){
		if(g->near_dialog){
			if(row1[seqlex + 1] > 0 || row1[seqlex + 2] < Int(seqlen)){
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
		beg = num_max(center - g->W, row1[seqlex + 1]);
		end = num_min(center + 1 + g->W, row1[seqlex + 2] + 1);
		if(end > seqlen) end = seqlen;
	} else {
		beg = row1[seqlex + 1];
		end = num_min(row1[seqlex + 2] + 1, Int(seqlen));
	}
	if(nidx2 == POG_TAIL_NODE){
		end = row1[seqlex + 2];
		while(end & 0x07){
			row1[end] = POG_SCORE_MIN;
			end ++;
		}
		while(end < seqlex){
			_mm_store_si128(((__m128i*)(row1 + end)), MIN);
			end += 8;
		}
	}
	beg = beg / 8;
	end = (end + 7) / 8;
	if(beg == 0){
		H = _mm_load_si128((__m128i*)row1);
		H = _mm_slli_si128(H, 2);
	} else {
		if(1){
			_mm_storeu_si128(((__m128i*)row2) + beg - 1, MIN);
		} else {
			for(i=0;i<beg;i++){
				_mm_store_si128(((__m128i*)row2) + i, MIN);
			}
		}
		H = _mm_loadu_si128((__m128i*)(row1 + beg * 8 - 1));
	}
	for(i=beg;i<end;i++){
		H = _mm_adds_epi16(H, _mm_load_si128(((__m128i*)qp) + i));
		E = _mm_load_si128(((__m128i*)row1) + i);
		E = _mm_adds_epi16(E, D);
		S = _mm_max_epi16(H, E);
		BT = _mm_cmpgt_epi16(S, H);
		BT = _mm_and_si128(BT, BT_MASK);
		_mm_store_si128(((__m128i*)row2) + i, S);
		_mm_store_si128(((__m128i*)btds) + i, BT);
		H = _mm_loadu_si128((__m128i*)(row1 + i * 8 + 7));
	}
	{
		__m128i F, I, CMP;
		u4i j, msk;
		b2i *r, *d;
		u4i *x;
		{
			F = _mm_set1_epi32(nidx1);
			for(i=beg;i<end<<1;i++){
				_mm_store_si128(((__m128i*)btxs) + i, F);
			}
		}
		I = _mm_set1_epi16(g->I);
		for(i=beg;i<end;i++){
			H = _mm_load_si128(((__m128i*)row2) + i);
			if(i > beg){
				F = _mm_loadu_si128((__m128i*)(row2 + i * 8 - 1));
			} else {
				F = _mm_slli_si128(H, 2);
				if(i){
					F = _mm_insert_epi16(F, POG_SCORE_MIN, 0);
				}
			}
			F = _mm_adds_epi16(F, I);
			CMP = _mm_cmpgt_epi16(F, H);
			msk = _mm_movemask_epi8(CMP);
			if(msk){
				r = row2 + i * 8;
				d = btds + i * 8;
				x = btxs + i * 8;
				j = __builtin_ctz(msk) >> 1;
				f = (j || i)? (r + j - 1)[0] : POG_SCORE_MIN;
				for(;j<8;j++){
					f = f + g->I;
					if(r[j] < f){
						r[j] = f;
						d[j] = 1;
						x[j] = nidx2;
					} else {
						f = r[j];
					}
				}
			}
		}
	}
	if(1){
		b2i ary[8];
		u4i j, k;
		MAX = _mm_set1_epi16(POG_SCORE_MIN);
		for(i=beg;i<end;i++){
			H = _mm_load_si128(((__m128i*)row2) + i);
			MAX = _mm_max_epi16(MAX, H);
		}
		_mm_storeu_si128((__m128i*)ary, MAX);
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
	if(1){
		if(end < slen){
			_mm_store_si128(((__m128i*)row2) + end, MIN);
		}
	} else {
		for(i=end;i<slen;i++){
			_mm_store_si128(((__m128i*)row2) + i, MIN);
		}
	}
	row2[seqlex] = mi;
	row2[seqlex + 1] = beg * 8;
	row2[seqlex + 2] = num_min(end * 8, seqlen);
}

// SSE
// BANDED, Auto fit the previous-row's max score in center
// OVERLAP
static inline void sse_band_row_rdaln_pog(POG *g, u4i nidx1, u4i nidx2, u4i seqlen, u4i coff1, u4i coff2, b2i *qp){
	__m128i I, D, H, E, F, S, MAX, MIN, CMP;
	__m128i BT, BT1, BT2, BT1_MASK, BT2_MASK;
	b2i *row1, *row2, *btds;
	u4i i, slen, seqlex, beg, end, *btxs;
	int lsth, lstf, msk, mi, center;
	slen = (seqlen + 7) / 8;
	seqlex = slen * 8;
	I = _mm_set1_epi16(g->I);
	D = _mm_set1_epi16(g->D);
	MIN = _mm_set1_epi16(POG_SCORE_MIN);
	BT1_MASK = _mm_set1_epi16(0b10);
	BT2_MASK = _mm_set1_epi16(0b01);
	row1 = ref_b2v(g->rows, coff1);
	row2 = ref_b2v(g->rows, coff2);
	btds = ref_b2v(g->btds, coff2);
	btxs = ref_u4v(g->btxs, coff2);
	if(g->W){
		if(row1[row1[seqlex]] >= g->W_score){
			if(g->near_dialog){
				if(row1[seqlex + 1] > 0 || row1[seqlex + 2] < Int(seqlen)){
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
			beg = num_max(center - g->W, row1[seqlex + 1]);
			end = num_min(center + 1 + g->W, row1[seqlex + 2] + 1);
			if(end > seqlen) end = seqlen;
		} else {
			beg = row1[seqlex + 1];
			end = num_min(row1[seqlex + 2] + 1, Int(seqlen));
		}
	} else {
		beg = 0;
		end = seqlen;
	}
	if(nidx2 == POG_TAIL_NODE){
		end = row1[seqlex + 2];
		while(end & 0x07){
			row1[end] = POG_SCORE_MIN;
			end ++;
		}
		while(end < seqlex){
			_mm_store_si128(((__m128i*)(row1 + end)), MIN);
			end += 8;
		}
	}
	beg = beg / 8;
	end = (end + 7) / 8;
	{
		F = _mm_set1_epi32(nidx1);
		for(i=beg;i<(end<<1);i++){
			_mm_stream_si128(((__m128i*)btxs) + i, F);
		}
	}
	if(beg){
		_mm_store_si128(((__m128i*)row2) + beg - 1, MIN);
	}

	MAX = _mm_set1_epi16(POG_SCORE_MIN);
	if(g->aln_mode == POG_ALNMODE_OVERLAP){
		lstf = lsth = beg? POG_SCORE_MIN : 0;
	} else {
		lstf = lsth = POG_SCORE_MIN;
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
		BT  = _mm_xor_si128(BT1, BT2); // 0, (1, 3), 2
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
	if(end < slen){
		_mm_store_si128(((__m128i*)row2) + end, MIN);
	}
	row2[seqlex] = mi;
	row2[seqlex + 1] = beg * 8;
	row2[seqlex + 2] = num_min(end * 8, seqlen);
}

static inline void merge_row_rdaln_pog(POG *g, u4i seqlen, u4i coff1, u4i coff2){
	b2i *row1, *row2, *btd1, *btd2;
	u4i i, *btx1, *btx2;
	u4i seqlex, beg[3], end[3], sz;
	int delta;
	row1 = ref_b2v(g->rows, coff1);
	row2 = ref_b2v(g->rows, coff2);
	btd1 = ref_b2v(g->btds, coff1);
	btd2 = ref_b2v(g->btds, coff2);
	btx1 = ref_u4v(g->btxs, coff1);
	btx2 = ref_u4v(g->btxs, coff2);
	seqlex = roundup_times(seqlen, 8);
	beg[0] = (row1[seqlex + 1]);
	beg[1] = (row2[seqlex + 1]);
	end[0] = (row1[seqlex + 2]);
	end[1] = (row2[seqlex + 2]);
	beg[2] = num_max(beg[0], beg[1]);
	end[2] = num_min(end[0], end[1]);
	if(1){
		__m128i R1, R2, D1, D2, X1, X2, MK, ML, MH;
		//b2i *ckd;
		//u4i *ckx;
		//ckd = calloc(end[2] - beg[2], 2);
		//ckx = calloc(end[2] - beg[2], 4);
		//for(i=beg[2];i<end[2];i++){
			//ckd[i - beg[2]] = row2[i] >= row1[i]? btd2[i] : btd1[i];
			//ckx[i - beg[2]] = row2[i] >= row1[i]? btx2[i] : btx1[i];
		//}
		for(i=beg[2];i&0x7;i++){
			if(row2[i] < row1[i]){
				row2[i] = row1[i];
				btd2[i] = btd1[i];
				btx2[i] = btx1[i];
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
			//ML = _mm_unpacklo_epi16(MK, _mm_setzero_si128());
			//MH = _mm_unpackhi_epi16(MK, _mm_setzero_si128());
			ML = _mm_unpacklo_epi16(MK, MK);
			MH = _mm_unpackhi_epi16(MK, MK);
			X1 = _mm_load_si128((__m128i*)(btx1 + i));
			X2 = _mm_load_si128((__m128i*)(btx2 + i));
			X2  = _mm_or_si128(_mm_and_si128(ML, X1), _mm_andnot_si128(ML, X2));
			_mm_store_si128((__m128i*)(btx2 + i), X2);
			X1 = _mm_load_si128((__m128i*)(btx1 + i + 4));
			X2 = _mm_load_si128((__m128i*)(btx2 + i + 4));
			X2  = _mm_or_si128(_mm_and_si128(MH, X1), _mm_andnot_si128(MH, X2));
			_mm_store_si128((__m128i*)(btx2 + i + 4), X2);
		}
		for(;i<end[2];i++){
			if(row2[i] < row1[i]){
				row2[i] = row1[i];
				btd2[i] = btd1[i];
				btx2[i] = btx1[i];
			}
		}
		//for(i=beg[2];i<end[2];i++){
			//if(ckd[i-beg[2]] != btd2[i]){
				//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				//abort();
			//}
			//if(ckx[i-beg[2]] != btx2[i]){
				//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				//abort();
			//}
		//}
		//free(ckd);
		//free(ckx);
	} else {
		for(i=beg[2];i<end[2];i++){
			if(row2[i] < row1[i]){
				row2[i] = row1[i];
				btd2[i] = btd1[i];
				btx2[i] = btx1[i];
			}
		}
	}
	if(beg[0] < beg[2]){
		if(beg[0] > 8){
			delta = 8;
		} else {
			delta = beg[0];
		}
		sz = num_min(beg[2], end[0]) - beg[0] + delta;
		memcpy(row2 + beg[0] - delta, row1 + beg[0] - delta, sz * sizeof(b2i));
		memcpy(btd2 + beg[0] - delta, btd1 + beg[0] - delta, sz * sizeof(b2i));
		memcpy(btx2 + beg[0] - delta, btx1 + beg[0] - delta, sz * sizeof(u4i));
	}
	for(i=end[2];i<beg[2];i++){ // in case of two independent regions
		row2[i] = POG_SCORE_MIN;
	}
	if(end[0] > end[2]){
		sz = num_min(end[0] + 8, seqlex) - end[2];
		memcpy(row2 + end[2], row1 + end[2], sz * sizeof(b2i));
		memcpy(btd2 + end[2], btd1 + end[2], sz * sizeof(b2i));
		memcpy(btx2 + end[2], btx1 + end[2], sz * sizeof(u4i));
	}
	row2[seqlex + 1] = num_min(beg[0], beg[1]);
	row2[seqlex + 2] = num_max(end[0], end[1]);
	if(row1[(row1[seqlex])] > row2[(row2[seqlex])]){
		row2[seqlex] = row1[seqlex];
	}
}

static inline void add_edge_core_pog(POG *g, pog_node_t *v, pog_edge_t *e){
	pog_edge_t *f;
	if(0){
		e->next = v->edge;
		v->edge = offset_pogedgev(g->edges, e);
		return;
	}
	if(v->edge == 0){
		v->edge = offset_pogedgev(g->edges, e);
		e->next = 0;
	} else {
		f = ref_pogedgev(g->edges, v->edge);
		if(e->cov > f->cov){
			v->edge = offset_pogedgev(g->edges, e);
			e->next = v->edge;
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
	}
}

static void inc_edge_core_pog(POG *g, pog_node_t *v, pog_edge_t *e){
	pog_edge_t *f, *p;
	e->cov ++;
	p = ref_pogedgev(g->edges, v->edge);
	if(p == e){
		return;
	}
	f = NULL;
	while(1){
		if(p->cov >= e->cov){
			f = p;
		}
		if(p->next == offset_pogedgev(g->edges, e)){
			break;
		}
		p = ref_pogedgev(g->edges, p->next);
	}
	if(f == p){
		return;
	}
	p->next = e->next;
	if(f){
		e->next = f->next;
		f->next = offset_pogedgev(g->edges, e);
	} else {
		e->next = v->edge;
		v->edge = offset_pogedgev(g->edges, e);
	}
}

void check_dp_cell_btxs(POG *g, u4i nidx, u8i coff, u2i seqlex){
	pog_node_t *u, *v;
	pog_edge_t *e;
	u4i i, btx, beg, end, eidx, pass;
	v = ref_pognodev(g->nodes, nidx);
	btx = MAX_U4;
	beg = g->rows->buffer[coff + seqlex + 1];
	end = g->rows->buffer[coff + seqlex + 2];
	for(i=beg;i<end;i++){
		if(g->rows->buffer[coff + i] != POG_SCORE_MIN && g->btxs->buffer[coff + i] != btx){
			btx = g->btxs->buffer[coff + i];
			u = ref_pognodev(g->nodes, btx);
			eidx = u->edge;
			pass = 0;
			while(eidx){
				e = ref_pogedgev(g->edges, eidx);
				if(e->node == nidx){
					pass = 1;
					break;
				}
				eidx = e->next;
			}
			if(pass == 0){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
}

static inline int align_rd_pog(POG *g, u2i rid){
	pog_node_t *u, *v;
	pog_edge_t *e, *f;
	u4i seqoff, seqlen, seqlex, slen, seqinc;
	b4i score, x, xb, xe;
	b2i *qp, *row, *btds;
	u4i *btxs;
	u8i rmax;
	u4i nidx, xidx, eidx, coff, btx, bt, mnode;
	u4i i, j, flag, bb, bl;
	seqoff = g->seqs->rdoffs->buffer[rid];
	seqlen = g->seqs->rdlens->buffer[rid];
	seqlex = roundup_times(seqlen, 8); // 128 = 8 * sizeof(b2i)
	slen = seqlex >> 3;
	seqinc = seqlex + 8; // seqlex, max_idx, beg, end, and paddings
	// set query prof
	// fit buffer to 16 bytes aligned
	head_sl_b2v(g->qprof, g->qprof->n_head);
	clear_and_encap_b2v(g->qprof, seqlex * 5 + 8);
	head_sr_b2v(g->qprof, (16 - (((u8i)g->qprof->buffer) & 0xF)) >> 1);
	for(i=0;i<4;i++){
		qp = g->qprof->buffer + i * seqlex;
		bl = 4;
		for(j=0;j<seqlen;j++){
			bb = get_basebank(g->seqs->rdseqs, seqoff + j);
			if(0){
				if(bb == bl){
					if(bb == i){
						qp[j] = g->M - 1;
					} else {
						qp[j] = g->X + 1;
					}
				} else {
					if(bb == i){
						qp[j] = g->M;
					} else {
						qp[j] = g->X;
					}
				}
			} else {
				if(bb == i){
					qp[j] = g->M;
				} else {
					qp[j] = g->X;
				}
			}
			bl = bb;
		}
		for(;j<seqlex;j++){
			qp[j] = POG_SCORE_MIN;
		}
	}
	{
		__m128i X;
		qp = g->qprof->buffer + 4 * seqlex;
		if(1){
			X = _mm_set1_epi16(g->X);
			for(j=0;j<slen;j++){
				_mm_stream_si128(((__m128i*)qp) + j, X);
			}
		} else {
			for(j=0;j<seqlex;j++){
				qp[j] = g->X;
			}
		}
	}
	// init graph
	encap_pognodev(g->nodes, seqlen);
	encap_pogedgev(g->edges, seqlen);
	rmax = g->seqs->nseq + 1;
	bb = 1;
	for(i=0;i<g->nodes->size;i++){
		v = ref_pognodev(g->nodes, i);
		v->coff = 0;
		v->vst = 0;
		if(v->nin || v->edge){
			rmax ++; // estimate max rows
		}
	}
	// init first row
	// fit to 16 bytes aligned
	head_sl_b2v(g->rows, g->rows->n_head);
	clear_and_encap_b2v(g->rows, rmax * seqinc + 8);
	head_sr_b2v(g->rows, (16 - (((u8i)g->rows->buffer) & 0xF)) >> 1);
	head_sl_b2v(g->btds, g->btds->n_head);
	clear_and_encap_b2v(g->btds, rmax * seqinc + 8);
	head_sr_b2v(g->btds, (16 - (((u8i)g->btds->buffer) & 0xF)) >> 1);
	head_sl_u4v(g->btxs, g->btxs->n_head);
	clear_and_encap_u4v(g->btxs, rmax * seqinc + 4);
	head_sr_u4v(g->btxs, (16 - (((u8i)g->btxs->buffer) & 0xF)) >> 2);
	u = ref_pognodev(g->nodes, POG_HEAD_NODE);
	u->coff = g->rows->size;
	inc_b2v(g->rows, seqinc);
	inc_b2v(g->btds, seqinc);
	inc_u4v(g->btxs, seqinc);
	row = ref_b2v(g->rows, u->coff);
	btds = ref_b2v(g->btds, u->coff);
	btxs = ref_u4v(g->btxs, u->coff);
	if(g->aln_mode == POG_ALNMODE_OVERLAP){ // overlap alignment
		memset(row,  0, seqinc * sizeof(b2i));
		row[0] = g->T;
	} else { // global
		__m128i MIN;
		MIN = _mm_set1_epi16(POG_SCORE_MIN);
		for(j=0;j<slen;j++){
			_mm_store_si128(((__m128i*)row) + j, MIN);
		}
		row[0] = g->T;
	}
	memset(btds, 0, seqinc * sizeof(b2i));
	memset(btxs, 0, seqinc * sizeof(u4i));
	if(g->near_dialog && g->W_score <= 0){
		row[seqlex] = g->W;
		row[seqlex + 1] = 0;
		row[seqlex + 2] = num_min(Int(seqlen), 2 * g->W + 1);
	} else {
		row[seqlex] = 0;
		row[seqlex + 1] = 0;
		row[seqlex + 2] = seqlen;
	}
	clear_u4v(g->stack);
	push_u4v(g->stack, POG_HEAD_NODE);
	score = POG_SCORE_MIN;
	mnode = POG_TAIL_NODE; // point to POG_TAIL
	x = seqlen - 1;
	int my_print = (cns_debug > 2);
	while(pop_u4v(g->stack, &nidx)){
		u = ref_pognodev(g->nodes, nidx);
		if(my_print){
			fprintf(stderr, "NODEALN[%u:R%u:%u] %d-%d:%d:%d %d\n", nidx, u->rid, u->pos, g->rows->buffer[u->coff + seqlex + 1], g->rows->buffer[u->coff + seqlex + 2], g->rows->buffer[u->coff + seqlex], g->rows->buffer[u->coff + g->rows->buffer[u->coff + seqlex]], g->btds->buffer[u->coff + g->rows->buffer[u->coff + seqlex]]);
		}
		if(g->aln_mode == POG_ALNMODE_OVERLAP){
			if((g->rows->buffer[u->coff + seqlex + 2]) >= Int(seqlen) &&  g->rows->buffer[u->coff + x] > score){ // overlap alignment
				score = g->rows->buffer[u->coff + x];
				mnode = nidx;
			}
		}
		eidx = u->edge;
		while(eidx){
			e = ref_pogedgev(g->edges, eidx);
			eidx = e->next;
			v = ref_pognodev(g->nodes, e->node);
			qp = g->qprof->buffer + v->base * seqlex;
			coff = g->rows->size;
			inc_b2v(g->rows, seqinc);
			inc_b2v(g->btds, seqinc);
			inc_u4v(g->btxs, seqinc);
			if(g->sse == 1){
				v1_sse_band_row_rdaln_pog(g, nidx, e->node, seqlen, u->coff, coff, qp);
			} else {
				sse_band_row_rdaln_pog(g, nidx, e->node, seqlen, u->coff, coff, qp);
			}
			if(v->vst){
				merge_row_rdaln_pog(g, seqlen, coff, v->coff);
				//check_dp_cell_btxs(g, e->node, v->coff, seqlex);
				trunc_b2v(g->rows, seqinc);
				trunc_b2v(g->btds, seqinc);
				trunc_u4v(g->btxs, seqinc);
			} else {
				v->coff = coff;
			}
			v->vst ++;
			if(v->vst >= v->nin){
				push_u4v(g->stack, e->node);
			}
		}
	}
	v = ref_pognodev(g->nodes, POG_TAIL_NODE);
	if(g->W == 0){
		row = ref_b2v(g->rows, v->coff);
		j = 0;
		for(i=1;i<seqlen;i++){
			if(row[i] > row[j]){
				j = i;
			}
		}
		row[seqlex] = j;
		row[seqlex + 1] = 0;
		row[seqlex + 2] = seqlen;
	}
	if(g->aln_mode == POG_ALNMODE_OVERLAP){
		if(g->rows->buffer[v->coff + (g->rows->buffer[v->coff + seqlex])] + g->T > score){
			score = g->rows->buffer[v->coff + (g->rows->buffer[v->coff + seqlex])];
			mnode = POG_TAIL_NODE;
			x = (g->rows->buffer[v->coff + seqlex]);
		}
	} else {
		score = g->rows->buffer[v->coff + seqlen - 1];
		mnode = POG_TAIL_NODE;
		x = seqlen - 1;
	}
	xe = x;
	xb = x;
	nidx = POG_TAIL_NODE;
	if(x + 1 < (int)seqlen){
		v = ref_pognodev(g->nodes, nidx);
		for(i=seqlen-1;(int)i>x;i--){
			if(0){
				fprintf(stderr, "BT[%u] y=N%u_R%u_%u_%c x=%d:%c z=%d\n", rid, nidx, v->rid, v->pos, bit_base_table[v->base], i, bit_base_table[get_basebank(g->seqs->rdseqs, seqoff + i)], -1);
				fflush(stderr);
			}
			e = next_ref_pogedgev(g->edges);
			e->node = offset_pognodev(g->nodes, v);
			e->cov  = 1;
			v->nin ++;
			u = next_ref_pognodev(g->nodes);
			ZEROS(u);
			u->rid  = rid;
			u->pos  = i;
			u->base = get_basebank(g->seqs->rdseqs, seqoff + u->pos);
			u->aligned = offset_pognodev(g->nodes, u);
			u->edge = offset_pogedgev(g->edges, e);
			e->next = 0;
			v = u;
		}
		v->coff = ref_pognodev(g->nodes, POG_TAIL_NODE)->coff;
	} else {
		v = ref_pognodev(g->nodes, nidx);
		if(mnode != POG_TAIL_NODE){
			u = ref_pognodev(g->nodes, mnode);
			g->btds->buffer[v->coff + x] = 0;
			g->btxs->buffer[v->coff + x] = mnode;
		}
	}
	btx = v->coff + x;
	while(1){
		bt = get_b2v(g->btds, btx);
		if(my_print){
			pog_node_t *w;
			w = ref_pognodev(g->nodes, nidx);
			fprintf(stderr, "BT[%u] y=N%u_R%u_%u_%c x=%d:%c z=%d %d\n", rid, nidx, w->rid, w->pos, bit_base_table[w->base], x, x>=0? bit_base_table[get_basebank(g->seqs->rdseqs, seqoff + x)] : '*', bt, g->rows->buffer[btx]);
			fflush(stderr);
		}
		if(bt == POG_DP_BT_M){
			e = next_ref_pogedgev(g->edges);
			e->node = offset_pognodev(g->nodes, v);
			e->cov  = 1;
			v->nin ++;
			if(btx){
				u = next_ref_pognodev(g->nodes);
				ZEROS(u);
				u->rid  = rid;
				u->pos  = x;
				u->base = get_basebank(g->seqs->rdseqs, seqoff + u->pos);
				u->aligned = offset_pognodev(g->nodes, u);
				u->edge = offset_pogedgev(g->edges, e);
				e->next = 0;
			} else {
				u = ref_pognodev(g->nodes, POG_HEAD_NODE);
				eidx = u->edge;
				flag = 0;
				while(eidx){
					f = ref_pogedgev(g->edges, eidx);
					eidx = f->next;
					if(f->node == e->node){
						inc_edge_core_pog(g, u, f);
						//f->cov ++;
						v->nin --;
						flag = 1;
						break;
					}
				}
				if(flag == 0){
					add_edge_core_pog(g, u, e);
					//e->next = u->edge;
					//u->edge = offset_pogedgev(g->edges, e);
				} else {
					trunc_pogedgev(g->edges, 1);
				}
				xb = x;
				break;
			}
			x --;
			if(nidx > POG_TAIL_NODE){
				flag = 0;
				xidx = nidx;
				do { // TODO: found BUG, endless loop, xidx == 0
					v = ref_pognodev(g->nodes, xidx);
					if(v->nin && v->base == u->base){
						u->edge = 0;
						eidx = v->edge;
						while(eidx){
							f = ref_pogedgev(g->edges, eidx);
							eidx = f->next;
							if(f->node == e->node){
								f->cov ++;
								ref_pognodev(g->nodes, e->node)->nin --;
								flag = 1;
								break;
							}
						}
						if(flag == 0){
							flag = 1;
							add_edge_core_pog(g, v, e);
							//e->next = v->edge;
							//v->edge = offset_pogedgev(g->edges, e);
						}
						break;
					}
					xidx = v->aligned;
				} while(xidx != nidx);
				{
					u->aligned = v->aligned;
					v->aligned = offset_pognodev(g->nodes, u);
				}
				if(flag){
					u = v;
				}
			}
			v = u;
		} else if(bt & POG_DP_BT_I){
			u = next_ref_pognodev(g->nodes);
			ZEROS(u);
			u->rid  = rid;
			u->pos  = x;
			u->base = get_basebank(g->seqs->rdseqs, seqoff + u->pos);
			u->aligned = offset_pognodev(g->nodes, u);
			u->edge = g->edges->size;
			e = next_ref_pogedgev(g->edges);
			e->node = offset_pognodev(g->nodes, v);
			e->cov  = 1;
			e->next = 0;
			v->nin ++;
			x --;
			v = u;
		} else {
		}
		if(x < 0){
			nidx = 0;
			btx = 0;
			continue;
		}
		nidx = (bt & POG_DP_BT_I)? nidx : get_u4v(g->btxs, btx);
		if(nidx >= g->nodes->size){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		u = ref_pognodev(g->nodes, nidx);
		btx = u->coff + x;
		if(x < g->rows->buffer[u->coff + seqlex + 1] || x >= g->rows->buffer[u->coff + seqlex + 2]){
			xb = x;
			while(x >= 0){
				e = next_ref_pogedgev(g->edges);
				e->node = offset_pognodev(g->nodes, v);
				e->cov  = 1;
				v->nin ++;
				u = next_ref_pognodev(g->nodes);
				ZEROS(u);
				u->rid  = rid;
				u->pos  = x --;
				u->base = get_basebank(g->seqs->rdseqs, seqoff + u->pos);
				u->aligned = offset_pognodev(g->nodes, u);
				u->edge = offset_pogedgev(g->edges, e);
				e->next = 0;
				v = u;
			}
			{
				e = next_ref_pogedgev(g->edges);
				e->node = offset_pognodev(g->nodes, v);
				e->cov  = 1;
				v->nin ++;
				u = ref_pognodev(g->nodes, POG_HEAD_NODE);
				add_edge_core_pog(g, u, e);
				//e->next = u->edge;
				//u->edge = offset_pogedgev(g->edges, e);
				v = u;
			}
			break;
		}
	}
	if(cns_debug > 1){
		fprintf(stderr, "ALIGN[%03d] len=%u aligned=%d,%d score=%d\n", rid, g->seqs->rdlens->buffer[rid], xb, xe + 1, score);
	}
	return score;
}

static inline void push_pog(POG *g, char *seq, u2i len){
	push_seqbank(g->seqs, NULL, 0, seq, len);
}

static inline void fwdbitpush_pog(POG *g, u8i *bits, u8i off, u2i len){
	fwdbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
}

static inline void revbitpush_pog(POG *g, u8i *bits, u8i off, u2i len){
	revbitpush_seqbank(g->seqs, NULL, 0, bits, off, len);
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
		row2[end] = 0;
		bs[end] = 0;
		if(cns_debug > 2){
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
		if(cns_debug > 2 && bt){
			fprintf(stderr, " -- x = %d y = %d bt = %d in %s -- %s:%d --\n", x, y, bt, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
	}
	while(y >= 0){
		r[y] = 4;
		y --;
	}
	if(cns_debug > 1){
		fprintf(stderr, " -- max = %d x = %d in %s -- %s:%d --\n", max, x, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	}
	return max;
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
			str[j-b] = "ACGT-"[c];
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
	free(str);
}

static inline void realign_msa_pog(POG *g){
	u4i ridx;
	if(g->rW <= 0 || g->seqs->nseq < 3) return;
	if(cns_debug > 1){
		fprintf(stderr, "RAW MSA\n");
		print_msa_pog(g, stderr);
	}
	for(ridx=0;ridx<g->seqs->nseq;ridx++){
		if(cns_debug > 1){
			fprintf(stderr, "REALIGN[%u]\n", ridx);
		}
		realign_msa_pog_core(g, ridx, g->rW);
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

static inline void gen_cns_pog(POG *g){
	u1i *s, *r;
	u4i idx, lst, lsv, lsx, rid, max, i, mcnt, min_cnt, max_cnt, runlen, cnl, beg, end;
	u4i freqs[5][11], fmax1, fmax2;
	u2i cnts[5], *bls, cl, bx, bl, bm, dl, dm, *hcovs, corr;
	//u2i *vsts;
	float fl, fx1, fx2, norm;
	mcnt = num_min(g->msa_min_cnt, g->seqs->nseq);
	min_cnt = num_max(mcnt, UInt(g->msa_min_freq * g->seqs->nseq));
	max_cnt = g->seqs->nseq - min_cnt;
	s = ref_u1v(g->msa, g->msa_len * g->seqs->nseq);
	memset(s, 4, g->msa_len);
	clear_basebank(g->cns);
	//vsts = calloc(g->seqs->nseq, sizeof(u2i));
	bls  = calloc(g->seqs->nseq, sizeof(u2i));
	hcovs = g->hcovs->buffer;
	if(g->rW){
		memset(g->hcovs->buffer, 0, g->msa_len * sizeof(u2i));
		for(rid=0;rid<g->seqs->nseq;rid++){
			r = ref_u1v(g->msa, g->msa_len * rid);
			beg = 0;
			while(beg < g->msa_len && r[beg] == 4) beg ++;
			end = g->msa_len;
			while(end && r[end - 1] == 4) end --;
			for(i=beg;i<end;i++) g->hcovs->buffer[i] ++;
		}
	}
	fmax1 = 5;
	fmax2 = 10;
	for(i=0;i<fmax1;i++){
		memset(freqs[i], 0, (fmax2 + 1) * sizeof(u4i));
	}
	lst = lsv = MAX_U4;
	lsx = 0;
	runlen = 0;
	cl = 0;
	cnl = 0;
	end = 0;
	for(idx=0;idx<=g->msa_len;idx++){
		if(idx <= g->msa_len){
			memset(cnts, 0, 5 * 2);
			for(rid=0;rid<g->seqs->nseq;rid++){
				cnts[g->msa->buffer[g->msa_len * rid + idx]] ++;
			}
			max = 0;
			for(i=1;i<4;i++){
				if(cnts[i] > cnts[max]){
					max = i;
				}
			}
		}
		if(idx == g->msa_len || (hcovs[idx] >= min_cnt && cnts[4] + hcovs[idx] - g->seqs->nseq <= max_cnt && cnts[max] >= mcnt)){
			if(lst == MAX_U4){ lst = idx; lsv = idx; }
			if(idx < g->msa_len){
				s[idx] = max;
				end = idx;
			} else if(end + 1 == idx){
				end = idx;
			}
			if(idx == g->msa_len || s[lst] != max){
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
						if(cl <= fmax1 && bls[rid - 1] <= fmax2){
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
				//if(bl && cl != bl){
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
	//free(vsts);
	free(bls);
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
	u4i i, ridx, nidx, eidx, xidx, moff, ready, beg, end;
	int score;
	for(ridx=0;ridx<g->seqs->nseq;ridx++){
		score = align_rd_pog(g, ridx);
	}
	for(nidx=0;nidx<g->nodes->size;nidx++){
		u = ref_pognodev(g->nodes, nidx);
		u->vst = 0;
		u->btx = 0;
		u->coff = 0;
	}
	//fprintf(stderr, " -- check dup edges in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	//check_dup_edges_pog(g);
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
			if(v->btx) continue; // already pushed
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
					v->btx  = 1;
					push_u4v(g->stack, e->node);
					xidx = v->aligned;
					while(xidx != e->node){
						x = ref_pognodev(g->nodes, xidx);
						x->coff = moff;
						if(x->edge){
							push_u4v(g->stack, xidx);
							x->btx = 1;
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
	u = ref_pognodev(g->nodes, POG_TAIL_NODE);
	g->msa_len = u->coff;
	// generate MSA
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
							set_u1v(g->msa, x->rid * g->msa_len + x->coff, x->base);
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
	clear_and_encap_u2v(g->hcovs, g->msa_len);
	memset(g->hcovs->buffer, 0, g->msa_len * sizeof(u2i));
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
	gen_cns_pog(g);
	if(cns_debug > 1){
		print_msa_pog(g, stderr);
	}
	if(0){
		fprintf(stderr, " -- seqs\t%llu --\n", (u8i)g->seqs->rdseqs->cap); fflush(stderr);
		fprintf(stderr, " -- nodes\t%llu --\n", (u8i)g->nodes->cap); fflush(stderr);
		fprintf(stderr, " -- edges\t%llu --\n", (u8i)g->edges->cap); fflush(stderr);
		fprintf(stderr, " -- rows\t%llu --\n", (u8i)g->rows->cap); fflush(stderr);
		fprintf(stderr, " -- btds\t%llu --\n", (u8i)g->btds->cap); fflush(stderr);
		fprintf(stderr, " -- btxs\t%llu --\n", (u8i)g->btxs->cap); fflush(stderr);
		fprintf(stderr, " -- cns\t%llu --\n", (u8i)g->cns->cap); fflush(stderr);
	}
}

#endif
