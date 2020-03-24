/*
 * 
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

#ifndef __WTDBG_H_RJ
#define __WTDBG_H_RJ

#include "kbm.h"
#include "kbmpoa.h"
#include "filewriter.h"
#include "pgzf.h"
#include <getopt.h>
#include <regex.h>

#define WT_MAX_RD			0x3FFFFFFF // 1 G
#define WT_MAX_RDLEN		0x00FFFFFF // 16 Mb
#define WT_MAX_RDBIN		0x0000FFFF // 64 K bins
#define WT_MAX_NODE			0x000000FFFFFFFFFFLLU
#define WT_MAX_EDGE			0x000000FFFFFFFFFFLLU
#define WT_MAX_NODE_EDGES	0xFFFF
#define WT_MAX_EDGE_COV		0xFFF

typedef struct {
	u8i rid:30, dir:1, beg:16, end:16, closed:1;
} rd_frg_t;
define_list(rdfrgv, rd_frg_t);

typedef struct {
	u8i node;
	u8i rid:30, dir:1, beg:16, end:16, closed:1;
} rd_reg_t;
define_list(rdregv, rd_reg_t);

typedef struct {
	u8i rid:30, dir:1, beg:16, end:16, closed:1;
	u8i read_link;
} rd_rep_t;
define_list(rdrepv, rd_rep_t);

typedef struct {
	u8i node;
	u8i rid:30, dir:1, beg:16, end:16, closed:1;
	u8i read_link;
} reg_t;
define_list(regv, reg_t);

typedef struct { uint32_t x, y; } xy_t;
define_list(xyv, xy_t);

typedef struct { uint64_t idx:46, cnt:18; } vec_ref_t;
typedef struct { uint64_t idx:46, cnt:18; } ptr_ref_t;
typedef struct { uint64_t idx:46, cnt:18, fix:1, rank:45, score:18; } rnk_ref_t;
static const vec_ref_t VEC_REF_NULL = (vec_ref_t){0, 0};
static const ptr_ref_t PTR_REF_NULL = (ptr_ref_t){0, 0};
static const rnk_ref_t RNK_REF_NULL = (rnk_ref_t){0, 0, 0, 0, 0};
define_list(vecrefv, vec_ref_t);
define_list(ptrrefv, ptr_ref_t);
define_list(rnkrefv, rnk_ref_t);
#define ptrref_hashcode(E) u64hashcode((E).idx)
#define ptrref_hashequals(E1, E2) ((E1).idx == (E2).idx)
define_hashset(ptrrefhash, ptr_ref_t, ptrref_hashcode, ptrref_hashequals);

#define WT_EDGE_CLOSED_NULL	0
#define WT_EDGE_CLOSED_MASK	1
#define WT_EDGE_CLOSED_LESS	2
#define WT_EDGE_CLOSED_HARD	3
typedef struct {
	uint64_t node1:45, dir1:1, dir2:1, status:1, closed:2, flag:2, cov:12;
	uint64_t node2:45; int64_t off:19;
} edge_t;
define_list(edgev, edge_t);

static inline uint64_t _edge_hashcode(edge_t e){
	const uint64_t m = 0xc6a4a7935bd1e995LLU;
	const int r = 47;
	uint64_t h = 1023 ^ (16 * m);
	uint64_t k = (e.node1 << 1) | e.dir1;
	k *= m;
	k ^= k >> r;
	k *= m;
	h ^= k;
	h *= m;
	k = (e.node2 << 1) | e.dir2;
	k *= m;
	k ^= k >> r;
	k *= m;
	h ^= k;
	h *= m;
	h ^= h >> r;
	h *= m;
	h ^= h >> r;
	return h;
}
#define EDGEHASH(idx) ((edgev *)set->userdata)->buffer[idx]
#define edge_hashcode(E) _edge_hashcode(EDGEHASH(E))
#define edge_hashequals(E1, E2) (EDGEHASH(E1).node1 == EDGEHASH(E2).node1 && EDGEHASH(E1).node2 == EDGEHASH(E2).node2 \
			&& EDGEHASH(E1).dir1 == EDGEHASH(E2).dir1 && EDGEHASH(E1).dir2 == EDGEHASH(E2).dir2)
define_hashset(edgehash, uint64_t, edge_hashcode, edge_hashequals);

typedef struct { uint64_t idx:63, flg:1; uint64_t next; } edge_ref_t;
static const edge_ref_t EDGE_REF_NULL = (edge_ref_t){0x7FFFFFFFFFFFFFFLLU, 1, 0};
define_list(edgerefv, edge_ref_t);

#define MAX_REP_IDX	MAX_U4

typedef struct {
	u8i rep_idx:32;
	u4i unvisit:16, cov:16;
	u8i closed:1, single_in:1, bt_visit:45, rep_dir:1, bt_idx:16, init_end:1;
	vec_ref_t regs;
	ptr_ref_t edges[2];
} node_t;
define_list(nodev, node_t);

typedef struct {
	u8i idx:39, flg:1, cnt:24;
} hit_lnk_t;
define_list(hitlnkv, hit_lnk_t);

typedef struct {
	rd_frg_t  frgs[2];
	hit_lnk_t lnks[2]; // idx+flg: link, cnt[0]: mat, cnt[1]:cglen
} rd_hit_t;
define_list(rdhitv, rd_hit_t);

typedef struct {
	u8i visit:63, flag:1;
	hit_lnk_t hits; // point to the g->rdhits
	int clps[2];
	ptr_ref_t regs;
	//u2i corr_bincnt;
} read_t;
define_list(readv, read_t);

typedef struct {
	u8i node;
	//node_t *n;
	edge_ref_t edges[2];
	u4i dir:2, cov:30;
	int off;
} trace_t;
define_list(tracev, trace_t);
#define WT_TRACE_MSG_ZERO	0
#define WT_TRACE_MSG_ONE	1
#define WT_TRACE_MSG_MORE	2
#define WT_TRACE_MSG_VISITED	3
#define WT_TRACE_MSG_UNDEF	4

typedef struct {
	union {
		struct {
			u4i frg1:31, dir1:1;
			u4i frg2:31, dir2:1;
		};
		u8i key;
	};
	u4i cov:13, flag:2, tidx1:8, tidx2:8, weak:1, closed:2;
	b4i off;
} lnk_t;
define_list(lnkv, lnk_t);
#define lnk_hashcode(E) (E).key
#define lnk_hashequals(E1, E2) (E1).key == (E2).key
define_hashset(lnkhash, lnk_t, lnk_hashcode, lnk_hashequals);

typedef struct {
	u8i toff:46, tcnt:18; // extended traces and core traces
	u4i tx, ty; // core traces
	ptr_ref_t lnks[2];
	u4i len, length;
	u8i rep_idx:48, unvisit:16;
	u8i closed:1, single_in:1, bt_visit:46, rep_dir:1, bt_idx:16;
} frg_t;
define_list(frgv, frg_t);

typedef struct {
	edge_ref_t lnks[2];
	u4i frg;
	u4i dir:2, cov:30;
	int off;
	u4i tx, ty;
} path_t;
define_list(pathv, path_t);
#define WT_PATH_MSG_ZERO	0
#define WT_PATH_MSG_ONE	1
#define WT_PATH_MSG_MORE	2
#define WT_PATH_MSG_VISITED	3
#define WT_PATH_MSG_UNDEF	4

typedef struct {
	u8i node1:63, dir1:1;
	u8i node2:63, dir2:1;
	b8i off:42, len:22;
} seqlet_t;
define_list(seqletv, seqlet_t);

typedef struct {
	KBM      *kbm;
	KBMPar   *par, *rpar;
	regv     *regs;
	readv    *reads;

	rdhitv   *rdhits;
	BitsVec  *cigars;

	nodev    *nodes;
	edgev    *edges;
	edgehash *ehash;
	edgerefv *erefs;

	frgv     *frgs;
	lnkv     *lnks;
	edgerefv *lrefs;
	tracev   *traces;

	u8i      genome_size;
	u4i      num_index;
	//u4i      corr_mode, corr_min, corr_max;
	//u4i      corr_bsize, corr_bstep;
	//float    corr_cov, corr_gcov;
	int      node_order, mem_stingy;
	u4i      n_fix, only_fix; // first n sequences are accurate contigs; only_fix means whether to include other pacbio sequenes
	u4i      reglen, regovl, bestn;
	int      min_node_mats;
	int      max_overhang, chainning_hits, uniq_hit;
	float    node_max_conflict; // 0.25
	float    node_merge_cutoff;
	uint32_t max_node_cov, min_node_cov, exp_node_cov, min_edge_cov, max_edge_span;
	u4i      max_node_cov_sg, max_sg_end;
	int      cut_tip;
	int      store_low_cov_edge;
	int      rep_filter, rep_detach;
	uint32_t bub_step, tip_step, rep_step;
	int min_ctg_len, min_ctg_nds, minimal_output;

	vplist *utgs;
	u4i    major_nctg;
	vplist *ctgs;
} Graph;

static const char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};

static inline Graph* init_graph(KBM *kbm){
	Graph *g;
	u4i rid;
	g = malloc(sizeof(Graph));
	g->kbm = kbm;
	g->par = kbm->par;
	g->rpar = NULL;
	g->regs = init_regv(32);
	g->reads = init_readv(kbm->reads->size);
	g->reads->size = kbm->reads->size;
	for(rid=0;rid<g->reads->size;rid++){
		g->reads->buffer[rid].clps[0] = 0;
		g->reads->buffer[rid].clps[1] = g->kbm->reads->buffer[rid].bincnt;
	}
	g->nodes = init_nodev(32);
	g->rdhits = init_rdhitv(1024);
	g->cigars = init_bitsvec(1024, 3);
	g->edges = init_edgev(32);
	g->ehash = init_edgehash(1023);
	set_userdata_edgehash(g->ehash, g->edges);
	g->erefs = init_edgerefv(32);
	g->frgs = init_frgv(32);
	g->lnks = init_lnkv(32);
	g->lrefs = init_edgerefv(32);
	g->traces = init_tracev(32);
	g->genome_size = 1024 * 1024 * 1024LLU;
	g->num_index = 1;
	//g->corr_mode = 0;
	//g->corr_min  = 5;
	//g->corr_max  = 10;
	//g->corr_cov  = 0.75;
	//g->corr_gcov = 5.0;
	//g->corr_bsize = 2048;
	//g->corr_bstep = 2048 - 512;
	g->node_order = 0;
	g->n_fix = 0;
	g->only_fix = 0;
	g->reglen = 1024;
	g->regovl = 256;
	g->node_max_conflict = 0.25;
	g->node_merge_cutoff = 0.8;
	g->max_overhang = -1;
	g->bestn = 0;
	g->min_node_mats = 1;
	g->chainning_hits = 0;
	g->uniq_hit = 0;
	g->min_node_cov = 4;
	g->max_node_cov = 60;
	g->exp_node_cov = 40; // expected node cov
	g->min_edge_cov = 4;
	g->max_edge_span = 1024; // 1024 * 256 = 256 kb
	g->max_node_cov_sg = 2;
	g->max_sg_end = 5;
	g->store_low_cov_edge = 1;
	g->rep_filter = 1;
	g->rep_detach = 1;
	g->bub_step = 10;
	g->tip_step = 5;
	g->rep_step = 20;
	g->cut_tip = 1;
	g->min_ctg_len = 10000;
	g->min_ctg_nds = 5;
	g->minimal_output = 0;
	g->utgs = init_vplist(32);
	g->ctgs = init_vplist(32);
	g->major_nctg = 0;
	return g;
}

static inline void free_graph(Graph *g){
	u8i i;
	free_regv(g->regs);
	free_readv(g->reads);
	free_nodev(g->nodes);
	free_rdhitv(g->rdhits);
	free_bitsvec(g->cigars);
	free_edgev(g->edges);
	free_edgehash(g->ehash);
	free_edgerefv(g->erefs);
	free_frgv(g->frgs);
	free_lnkv(g->lnks);
	free_edgerefv(g->lrefs);
	free_tracev(g->traces);
	for(i=0;i<g->utgs->size;i++) free_tracev(g->utgs->buffer[i]);
	free_vplist(g->utgs);
	for(i=0;i<g->ctgs->size;i++) free_tracev(g->ctgs->buffer[i]);
	free_vplist(g->ctgs);
	free(g);
}

#define KBM_MAP2RDHIT_QUICK

static inline int map2rdhits_graph(Graph *g, kbm_map_t *hit){
	u4i k, f, d, p, n, add;
	rd_hit_t *rh, *hp, *hn;
	read_t *rd;
	if(hit->mat == 0) return 0;
	rh = next_ref_rdhitv(g->rdhits);
	rh->frgs[0] = (rd_frg_t){hit->qidx, hit->qdir, hit->qb, hit->qe, 0};
	rh->frgs[1] = (rd_frg_t){hit->tidx, hit->tdir, hit->tb, hit->te, 0};
	rh->lnks[0].cnt = hit->mat;
	rh->lnks[1].cnt = hit->cglen;
	add = 0;
	for(k=0;k<2;k++){
		rd = ref_readv(g->reads, rh->frgs[k].rid);
#ifdef KBM_MAP2RDHIT_QUICK
		{ // Just add it
			rh->lnks[k].idx = rd->hits.idx;
			rh->lnks[k].flg = rd->hits.flg;
			rd->hits.idx = offset_rdhitv(g->rdhits, rh);
			rd->hits.flg = k;
			rd->hits.cnt ++;
			add ++;
			continue;
		}
#endif
		hn = ref_rdhitv(g->rdhits, rd->hits.idx);
		f  = rd->hits.flg;
		hp = NULL;
		p = 0;
		n = 1;
		while(1){
			if(hn->lnks[0].cnt <= rh->lnks[0].cnt){ // hit->mat
				break;
			}
			p = f;
			hp = hn;
			if(hn->lnks[f].idx == 0) break;
			f = hn->lnks[p].flg;
			hn = ref_rdhitv(g->rdhits, hn->lnks[p].idx);
			n ++;
		}
		if(g->bestn && n > g->bestn){
			continue;
		}
		if(hp == NULL){
			rh->lnks[k].idx = rd->hits.idx;
			rh->lnks[k].flg = rd->hits.flg;
			rd->hits.idx = offset_rdhitv(g->rdhits, rh);
			rd->hits.flg = k;
		} else {
			rh->lnks[k].idx = hp->lnks[p].idx;
			rh->lnks[k].flg = hp->lnks[p].flg;
			hp->lnks[p].idx = offset_rdhitv(g->rdhits, rh);
			hp->lnks[p].flg = k;
		}
		rd->hits.cnt ++;
		if(g->bestn && rd->hits.cnt > g->bestn){
			hn = rh;
			p  = k;
			while(n < g->bestn){
				f  = hn->lnks[p].flg;
				hn = ref_rdhitv(g->rdhits, hn->lnks[p].idx);
				p = f;
				n ++;
			}
			hp = hn;
			f = p;
			while(hn->lnks[f].idx){
				d  = hn->lnks[f].flg;
				hn = ref_rdhitv(g->rdhits, hn->lnks[f].idx);
				hn->frgs[f].closed = 1;
				f = d;
			}
			hp->lnks[p].idx = 0;
			hp->lnks[p].flg = 0;
			rd->hits.cnt = n;
		}
		add ++;
	}
	if(add == 0){
		trunc_rdhitv(g->rdhits, 1);
		g->cigars->size -= hit->cglen;
	}
	return add;
}

static inline int is_dovetail_overlap(Graph *g, kbm_map_t *hit){
	read_t *q, *t;
	int overhangs[2][2];
	q = ref_readv(g->reads, hit->qidx);
	t = ref_readv(g->reads, hit->tidx);
	overhangs[1][0] = hit->tb - t->clps[0];
	overhangs[1][1] = (t->clps[1] - t->clps[0]) - hit->te;
	if(hit->qdir){
		overhangs[0][0] = (q->clps[1] - q->clps[0]) - hit->qe;
		overhangs[0][1] = hit->qb;
	} else {
		overhangs[0][0] = hit->qb;
		overhangs[0][1] = (q->clps[1] - q->clps[0]) - hit->qe;
	}
	if(overhangs[0][0] > g->max_overhang && overhangs[1][0] > g->max_overhang){
		return 0;
	}
	if(overhangs[0][1] > g->max_overhang && overhangs[1][1] > g->max_overhang){
		return 0;
	}
	return 1;
}

// qlen is measured by KBM_BIN_SIZE
static inline int hit2rdregs_graph(Graph *g, rdregv *regs, int qlen, kbm_map_t *hit, BitsVec *cigars, u4v *maps[3]){
	KBM *kbm;
	u8i ndoff;
	u4i bpos[2][2], npos[2][2], clen, ndbeg, qn, j, qbincnt;
	int tmp, bt, tlen, x, y, mat, beg, end, min_node_len, max_node_len;
	int mask, closed;
	kbm = g->kbm;
	mask = 0;
	if(g->max_overhang >= 0){
		if(!is_dovetail_overlap(g, hit)){
			mask = 1;
		}
	}
	qn = g->reglen;
	min_node_len = (qn - 1);
	max_node_len = (qn + 1);
	// translate into reverse sequence order
	if(qlen == 0){
		qbincnt = kbm->reads->buffer[hit->qidx].bincnt;
		qlen = qbincnt;
	} else if(qlen > kbm->reads->buffer[hit->qidx].bincnt){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	} else {
		qbincnt = qlen;
	}
	tlen = kbm->reads->buffer[hit->tidx].bincnt;
	if(hit->qdir){
		tmp = qlen - hit->qb;
		hit->qb = qlen - hit->qe;
		hit->qe = tmp;
	}
	{
		clear_u4v(maps[0]);
		clear_u4v(maps[1]);
		clear_u4v(maps[2]);
		clen = hit->cglen;
		x = -1; y = -1;
		while(clen){
			bt = get_bitsvec(cigars, hit->cgoff + clen - 1);
			push_u4v(maps[2], !(bt >> 2));
			bt = bt & 0x03;
			switch(bt){
				case 0: x ++; y ++; break;
				case 1: x ++; break;
				case 2: y ++; break;
			}
			push_u4v(maps[0], (x < 0)? 0 : x);
			push_u4v(maps[1], (y < 0)? 0 : y);
			clen --;
		}
		push_u4v(maps[0], x + 1);
		push_u4v(maps[1], y + 1);
		push_u4v(maps[2], 0);
#if DEBUG
		if(x + 1 + hit->qb != hit->qe || y + 1 + hit->tb != hit->te){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			print_hit_kbm(g->kbm, g->kbm->reads->buffer[hit->qidx].tag, g->kbm->reads->buffer[hit->qidx].bincnt * KBM_BIN_SIZE, hit, cigars, NULL, stderr);
			abort();
		}
#endif
	}
	bpos[0][0] = hit->qb;
	bpos[0][1] = hit->qe;
	bpos[1][0] = hit->tb;
	bpos[1][1] = hit->te;
#if 0
	fprintf(stdout, "BPOS\t%d\t%d\t%d\t%d\n", bpos[0][0], bpos[0][1], bpos[1][0], bpos[1][1]);
	for(j=0;j<maps[0]->size;j++){
		fprintf(stdout, "%d,%d\t", maps[0]->buffer[j], maps[1]->buffer[j]);
		if((j % 10) == 9) fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
#endif
	{
		ndoff = kbm->reads->buffer[hit->qidx].binoff;
		if(hit->qdir){
			ndbeg = qbincnt - bpos[0][0];
			ndbeg = ndbeg % qn;
			ndoff = ndoff + ((qbincnt - (ndbeg + bpos[0][0])) / qn) - 1;
		} else {
			ndbeg = (bpos[0][0] % qn)? qn - (bpos[0][0] % qn) : 0;
			ndoff = ndoff + ((ndbeg + bpos[0][0]) / qn);
		}
		x = 0;
		for(j=ndbeg+bpos[0][0];j+qn<=bpos[0][1];j+=qn){
			mat = 0;
			while(maps[0]->buffer[x] < j - bpos[0][0]){
				x ++;
			}
			npos[0][0] = maps[1]->buffer[x];
			while(maps[0]->buffer[x] == j - bpos[0][0]){
				if(maps[2]->buffer[x]) mat ++;
				x ++;
			}
			npos[0][1] = maps[1]->buffer[x];
			while(maps[0]->buffer[x] < j + qn - 1 - bpos[0][0]){
				if(maps[2]->buffer[x]) mat ++;
				x ++;
			}
			npos[1][0] = maps[1]->buffer[x];
			while(maps[0]->buffer[x + 1] == j + qn - 1 - bpos[0][0]){
				if(maps[2]->buffer[x]) mat ++;
				x ++;
			}
			npos[1][1] = maps[1]->buffer[x];
			// TODO: arbitrarily use outer boundary as matched region, need to test its effect
			// TODO: whether to remove regs outside clps
			if(mat >= g->min_node_mats){
				beg = (npos[0][0] + bpos[1][0]);
				end = (npos[1][1] + bpos[1][0] + 1);
				if(end - beg >= min_node_len && end - beg <= max_node_len){
					closed = 0;
				} else {
					closed = 1;
				}
#if DEBUG
				if(beg >= end || beg < 0){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
#endif
				push_rdregv(regs, (rd_reg_t){ndoff, hit->tidx, hit->qdir, beg, end, mask | closed});
			}
			if(hit->qdir){
				ndoff --;
			} else {
				ndoff ++;
			}
#if 0
			rd_reg_t *rg = ref_rdregv(regs, regs->size - 1);
			fprintf(stdout, "NPOS\tX=%d,%d\tY=%d,%d\t%d\t%d\t%d\t%d\n", j, j - bpos[0][0], j + qn - 1, j + qn - 1 - bpos[0][0], npos[0][0], npos[0][1], npos[1][0], npos[1][1]);
			fprintf(stdout, "REG\t%llu\t%s\t%c\t%d\t%d\t%d\n", rg->node, kbm->reads->buffer[rg->rid].tag, "+-"[rg->dir], rg->beg, rg->end,  rg->end - rg->beg);
#endif
		}
	}
	//if(!g->corr_mode){
	{
		ndoff = kbm->reads->buffer[hit->tidx].binoff;
		ndbeg = (bpos[1][0] % qn)? qn - (bpos[1][0] % qn): 0;
		ndoff = ndoff + ((ndbeg + bpos[1][0]) / qn);
		x = 0;
		for(j=ndbeg+bpos[1][0];j+qn<=bpos[1][1];j+=qn){
			mat = 0;
			while(maps[1]->buffer[x] < j - bpos[1][0]){
				x ++;
			}
			npos[0][0] = maps[0]->buffer[x];
			while(maps[1]->buffer[x] == j - bpos[1][0]){
				if(maps[2]->buffer[x]) mat ++;
				x ++;
			}
			npos[0][1] = maps[0]->buffer[x];
			while(maps[1]->buffer[x] < j + qn - 1 - bpos[1][0]){
				if(maps[2]->buffer[x]) mat ++;
				x ++;
			}
			npos[1][0] = maps[0]->buffer[x];
			while(maps[1]->buffer[x + 1] == j + qn - 1 - bpos[1][0]){
				if(maps[2]->buffer[x]) mat ++;
				x ++;
			}
			npos[1][1] = maps[0]->buffer[x];
			// TODO: arbitrarily use outer boundary as matched region, need to test its effect
			if(mat >= g->min_node_mats){
				if(hit->qdir){
					beg = qlen - (npos[1][1] + bpos[0][0] + 1);
					end = qlen - (npos[0][0] + bpos[0][0]);
				} else {
					beg = (npos[0][0] + bpos[0][0]);
					end = (npos[1][1] + bpos[0][0] + 1);
				}
				if(end - beg >= min_node_len && end - beg <= max_node_len){
					closed = 0;
				} else {
					closed = 1;
				}
#if DEBUG
				if(beg >= end || beg < 0){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
#endif
				push_rdregv(regs, (rd_reg_t){ndoff, hit->qidx, hit->qdir, beg, end, mask | closed});
#if 0
				fprintf(stdout, "NPOS\tX=%d,%d\tY=%d,%d\t%d\t%d\t%d\t%d\n", j, j - bpos[1][0], j + qn - 1, j + qn - 1 - bpos[1][0], npos[0][0], npos[0][1], npos[1][0], npos[1][1]);
				rd_reg_t *rg = ref_rdregv(regs, regs->size - 1);
				fprintf(stdout, "REG\t%llu\t%s\t%c\t%d\t%d\t%d\n", rg->node, kbm->reads->buffer[rg->rid].tag, "+-"[rg->dir], rg->beg, rg->end,  rg->end - rg->beg);
#endif
			}
			ndoff ++;
		}
	}
	if(hit->qdir){
		tmp = qlen - hit->qb;
		hit->qb = qlen - hit->qe;
		hit->qe = tmp;
	}
	return !mask;
}

thread_beg_def(mdbg);
Graph *g;
KBMAux *aux, *raux;
CTGCNS *cc;
reg_t reg;
rdregv *regs;
BitVec *rdflags;
u4i beg, end;
int raw;
FILE *alno;
int task;
thread_end_def(mdbg);

thread_beg_func(mdbg);
Graph *g;
KBM *kbm;
KBMAux *aux, *raux;
rdregv *regs;
BitVec *rdflags;
u4v *maps[3], *tidxs;
volatile reg_t *reg;
g = mdbg->g;
kbm = g->kbm;
reg = (reg_t*)&mdbg->reg;
aux = mdbg->aux;
raux = mdbg->raux;
regs = mdbg->regs;
rdflags = mdbg->rdflags;
maps[0] = init_u4v(32);
maps[1] = init_u4v(32);
maps[2] = init_u4v(32);
tidxs = init_u4v(16);
thread_beg_loop(mdbg);
if(mdbg->task == 1){
	if(reg->closed) continue;
	//if(g->corr_mode){
	if(0){
	} else {
		query_index_kbm(aux, NULL, reg->rid, kbm->rdseqs, (kbm->reads->buffer[reg->rid].seqoff + UInt(reg->beg)) * KBM_BIN_SIZE, UInt(reg->end - reg->beg) * KBM_BIN_SIZE);
		map_kbm(aux);
	}
	if(raux && aux->hits->size){ // refine
		kbm_read_t *rd;
		u4i i, j, tidx;
		clear_kbm(raux->kbm);
		bitpush_kbm(raux->kbm, NULL, 0, kbm->rdseqs->bits, (kbm->reads->buffer[reg->rid].seqoff + UInt(reg->beg)) * KBM_BIN_SIZE, UInt(reg->end - reg->beg) * KBM_BIN_SIZE);
		ready_kbm(raux->kbm);
		simple_index_kbm(raux->kbm, 0, raux->kbm->bins->size);
		clear_u4v(tidxs);
		sort_array(aux->hits->buffer, aux->hits->size, kbm_map_t, num_cmpgt(a.tidx, b.tidx));
		for(i=0;i<aux->hits->size;i++){
			if(i && tidxs->buffer[tidxs->size - 1] == aux->hits->buffer[i].tidx) continue;
			push_u4v(tidxs, aux->hits->buffer[i].tidx);
		}
		clear_kbmmapv(aux->hits);
		clear_bitsvec(aux->cigars);
		for(i=0;i<tidxs->size;i++){
			tidx = get_u4v(tidxs, i);
			rd = ref_kbmreadv(aux->kbm->reads, tidx);
			query_index_kbm(raux, rd->tag, tidx, aux->kbm->rdseqs, rd->seqoff * KBM_BIN_SIZE, rd->bincnt * KBM_BIN_SIZE);
			map_kbm(raux);
			for(j=0;j<raux->hits->size;j++){
				flip_hit_kbmaux(aux, raux, j);
			}
		}
	}
	sort_array(aux->hits->buffer, aux->hits->size, kbm_map_t, num_cmpgt(b.mat, a.mat));
}
thread_end_loop(mdbg);
free_u4v(maps[0]);
free_u4v(maps[1]);
free_u4v(maps[2]);
free_u4v(tidxs);
thread_end_func(mdbg);

typedef struct {
	int pos:19;
	u4i dir:1, spur:1, dep:11;
} rd_clp_t;
define_list(rdclpv, rd_clp_t);

static inline void clip_read_algo(int clps[2], rdclpv *brks, rdclpv *chis, int avg_dep, int min_dep){
	rd_clp_t *c;
	u4i i;
	int dep, x, y, max, mx, my;
	sort_array(brks->buffer, brks->size, rd_clp_t, num_cmpgt((a.pos << 1) | a.dir, (b.pos << 1) | b.dir));
	x = y = 0; max = 0; dep = 0; mx = my = 0;
	for(i=0;i<brks->size;i++){
		c = ref_rdclpv(brks, i);
		if(dep >= min_dep){
			y = c->pos;
			if(y - x > max){
				max = y - x;
				mx = x;
				my = y;
			}
		}
		if(c->dir){ // end of overlap
			c->dep = dep;
			dep --;
		} else {
			dep ++;
			c->dep = dep;
			if(dep == min_dep){
				x = c->pos;
			}
		}
		if(c->spur){
			push_rdclpv(chis, (rd_clp_t){c->pos - KBM_BIN_SIZE, 0, 0, c->dep});
			push_rdclpv(chis, (rd_clp_t){c->pos - 1,            1, 0, c->dep});
			push_rdclpv(chis, (rd_clp_t){c->pos,                0, 1, c->dep});
			push_rdclpv(chis, (rd_clp_t){c->pos + KBM_BIN_SIZE, 1, 0, c->dep});
		}
	}
	clps[0] = mx / KBM_BIN_SIZE;
	clps[1] = my / KBM_BIN_SIZE;
	if((int)chis->size < avg_dep){
		return;
	}
	sort_array(chis->buffer, chis->size, rd_clp_t, num_cmpgtx(a.pos, b.pos, b.dir, a.dir));
	dep = 0; max = 0; mx = 0;
	for(i=0;i<chis->size;i++){
		c = ref_rdclpv(chis, i);
		if(c->dir){
			if(c->spur){
				if(dep >= max){
					max = dep;
					mx = i;
				}
			}
			dep --;
		} else {
			dep ++;
			if(c->spur){
				if(dep >= max){
					max = dep;
					mx = i;
				}
			}
		}
	}
	if(max * 2 < avg_dep){
		return;
	}
	c = ref_rdclpv(chis, mx);
	if(c->dep >= avg_dep){
		return;
	}
	if(2 * c->dep > max + 1){
		return;
	}
	if(c->pos <= clps[0] * KBM_BIN_SIZE || c->pos >= clps[1] * KBM_BIN_SIZE){
		return;
	}
	if(c->pos - clps[0] > clps[1] - c->pos){
		clps[1] = (c->pos + 1) / KBM_BIN_SIZE;
	} else {
		clps[0] = (c->pos + 1) / KBM_BIN_SIZE;
	}
}

static inline void clip_read_core(Graph *g, u4i rid, hitlnkv *lnks, rdclpv *brks, rdclpv *chis){
	read_t *rd;
	hit_lnk_t *lnk;
	rd_hit_t *hit;
	rd_frg_t *f1, *f2;
	u4i i;
	int rlen, dlen, margin, min_dep, dep, avg, x, y;
	rd = ref_readv(g->reads, rid);
	rlen = g->kbm->reads->buffer[rid].bincnt;
	if(rlen == 0){
		return;
	}
	lnk = &rd->hits;
	clear_rdclpv(brks);
	margin = g->max_overhang > 0? g->max_overhang : 4;
	min_dep = 2;
	dep = 0;
	clear_hitlnkv(lnks);
	while(lnk->idx){
		hit = ref_rdhitv(g->rdhits, lnk->idx);
		f1 = hit->frgs + lnk->flg;
		if(f1->closed == 0) push_hitlnkv(lnks, *lnk);
		lnk = hit->lnks + lnk->flg;
	}
	for(i=0;i<lnks->size;i++){
		lnk = ref_hitlnkv(lnks, i);
		hit = ref_rdhitv(g->rdhits, lnk->idx);
		f1 = hit->frgs + lnk->flg;
		if(f1->closed == 0){
			// whether be contained
			if(f1->beg <= margin && f1->end + margin >= rlen){
				rd->clps[0] = 0;
				rd->clps[1] = rlen;
				return;
			}
			f2 = hit->frgs + !lnk->flg;
			dlen = g->kbm->reads->buffer[f2->rid].bincnt;
			if(f1->dir ^ f2->dir){
				x = (f1->beg > margin && f2->end + margin < dlen);
				y = (f1->end + margin < rlen && f2->beg > margin);
			} else {
				x = (f1->beg > margin && f2->beg > margin);
				y = (f1->end + margin < rlen && f2->end + margin < dlen);
			}
			if(x && y){
				//push_rdclpv(brks, (rd_clp_t){f1->beg, 0, 1, 0});
				//push_rdclpv(brks, (rd_clp_t){f1->end, 1, 1, 0});
			} else if(x){
				push_rdclpv(brks, (rd_clp_t){f1->beg * KBM_BIN_SIZE, 0, 1, 0});
				push_rdclpv(brks, (rd_clp_t){f1->end * KBM_BIN_SIZE, 1, 0, 0});
			} else if(y){
				push_rdclpv(brks, (rd_clp_t){f1->beg * KBM_BIN_SIZE, 0, 0, 0});
				push_rdclpv(brks, (rd_clp_t){f1->end * KBM_BIN_SIZE, 1, 1, 0});
			} else {
				dep += f1->end - f1->beg;
				push_rdclpv(brks, (rd_clp_t){f1->beg * KBM_BIN_SIZE, 0, 0, 0});
				push_rdclpv(brks, (rd_clp_t){f1->end * KBM_BIN_SIZE, 1, 0, 0});
			}
		}
	}
	clear_rdclpv(chis);
	avg = (dep + rlen - 1) / rlen;
	clip_read_algo(rd->clps, brks, chis, avg, min_dep);
}

thread_beg_def(mclp);
Graph *g;
int task;
thread_end_def(mclp);

thread_beg_func(mclp);
Graph *g;
rdclpv *brks, *chis;
hitlnkv *lnks;
u4i rid;
g = mclp->g;
brks = init_rdclpv(32);
chis = init_rdclpv(32);
lnks = init_hitlnkv(32);
thread_beg_loop(mclp);
if(mclp->task == 1){
	for(rid=mclp->t_idx;rid<mclp->g->reads->size;rid+=mclp->n_cpu){
		clip_read_core(mclp->g, rid, lnks, brks, chis);
	}
} else if(mclp->task == 2){
	hitlnkv *lnks, *chks;
	read_t *rd;
	hit_lnk_t *lnk;
	u4i i;
	lnks = init_hitlnkv(1024);
	chks = init_hitlnkv(1024);
	for(rid=mclp->t_idx;rid<mclp->g->reads->size;rid+=mclp->n_cpu){
		rd = ref_readv(g->reads, rid);
		if(rd->hits.idx == 0) continue;
		clear_hitlnkv(lnks);
		lnk = &rd->hits;
		while(1){
			if(lnk->idx == 0) break;
			push_hitlnkv(lnks, *lnk);
			lnk = g->rdhits->buffer[lnk->idx].lnks + lnk->flg;
		}
		if(lnks->size < 2) continue;
		sort_array(lnks->buffer, lnks->size, hit_lnk_t, num_cmpgt(g->rdhits->buffer[b.idx].lnks[0].cnt, g->rdhits->buffer[a.idx].lnks[0].cnt));
		lnk = &rd->hits;
		for(i=0;i<lnks->size;i++){
			lnk->idx = lnks->buffer[i].idx;
			lnk->flg = lnks->buffer[i].flg;
			lnk = g->rdhits->buffer[lnks->buffer[i].idx].lnks + lnks->buffer[i].flg;
		}
		lnk->idx = 0;
		lnk->flg = 0;
		{
			clear_hitlnkv(chks);
			lnk = &rd->hits;
			while(1){
				if(lnk->idx == 0) break;
				push_hitlnkv(chks, *lnk);
				lnk = g->rdhits->buffer[lnk->idx].lnks + lnk->flg;
			}
			if(lnks->size != chks->size){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			for(i=0;i<lnks->size;i++){
				if(lnks->buffer[i].idx != chks->buffer[i].idx || lnks->buffer[i].flg != chks->buffer[i].flg){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
			}
		}
		if(g->bestn && lnks->size > mclp->g->bestn){
			lnk = lnks->buffer + g->bestn - 1;
			lnk = g->rdhits->buffer[lnk->idx].lnks + lnk->flg;
			lnk->idx = 0;
			lnk->flg = 0;
			for(i=g->bestn;i<lnks->size;i++){
				lnk = lnks->buffer + i;
				g->rdhits->buffer[lnk->idx].frgs[lnk->flg].closed = 1;
			}
			rd->hits.cnt = g->bestn;
		}
	}
	free_hitlnkv(lnks);
	free_hitlnkv(chks);
}
thread_end_loop(mclp);
free_rdclpv(brks);
free_rdclpv(chis);
free_hitlnkv(lnks);
thread_end_func(mclp);

static inline u8i check_read_reg_conflict_core(Graph *g, rd_reg_t *hit, int *conflict){
	read_t *rd;
	reg_t *reg;
	u8i idx, pre;
	int x, y;
	rd = ref_readv(g->reads, hit->rid);
	idx = rd->regs.idx;
	pre = 0;
	*conflict = 0;
	while(idx){
		reg = ref_regv(g->regs, idx);
		if(hit->beg >= reg->beg) pre = idx;
		idx = reg->read_link;
		if(hit->end <= reg->beg) break;
		if(hit->beg >= reg->end) continue;
		if(reg->closed) continue;
		if(hit->beg <= reg->beg){
			x = reg->beg;
			if(hit->end >= reg->end){
				y = reg->end;
			} else y = hit->end;
		} else {
			x = hit->beg;
			if(hit->end <= reg->end){
				y = hit->end;
			} else y = reg->end;
		}
		if(x < y){
			if(x + (int)g->regovl < y) *conflict = 1;
		}
	}
	return pre;
}

thread_beg_def(mupd);
Graph *g;
rdregv *regs;
rnk_ref_t *nd;
u8v *ins;
u8i vt;
thread_end_def(mupd);

thread_beg_func(mupd);
Graph *g;
rdregv *regs;
rd_reg_t *hit;
u8i pre;
u4i i;
int conflict;
g = mupd->g;
regs = mupd->regs;
thread_beg_loop(mupd);
if(mupd->nd == NULL) continue;
clear_u8v(mupd->ins);
if(1){
	for(i=1;i<mupd->nd->cnt;i++){
		if(regs->buffer[mupd->nd->idx + i].rid == regs->buffer[mupd->nd->idx + i - 1].rid){
			regs->buffer[mupd->nd->idx + i].closed = 1;
			regs->buffer[mupd->nd->idx + i - 1].closed = 1;
		}
	}
}
for(i=0;i<mupd->nd->cnt;i++){
	hit = ref_rdregv(regs, mupd->nd->idx + i);
	conflict = 0;
	pre = check_read_reg_conflict_core(g, hit, &conflict);
	if(conflict){
		hit->closed = 1;
	}
	push_u8v(mupd->ins, pre);
}
thread_end_loop(mupd);
thread_end_func(mupd);

//TODO: reg_t to store regs, sort them by rank, and sort them by rid in another u4i/u8i array
//TODO: fast generate nodes and read_link, then select important intervals
static inline void mul_update_regs_graph(Graph *g, rdregv *regs, rnkrefv *nds, u4i ncpu, u8i upds[3]){
	u8i idx, vt, inc;
	u4i i, j, pass;
	read_t *rd;
	node_t *n;
	reg_t *reg, *r;
	rd_reg_t *hit;
	int conflict, closed;
	thread_prepare(mupd);
	for(i=0;i<g->reads->size;i++){
		g->reads->buffer[i].visit = 0; // in which round did the read be touched
	}
	if(0){
		encap_regv(g->regs, roundup_times(regs->size + 1, 8)); // make sure next_ref_regv in this function is thread-safe
	}
	if(0){
		renew_regv(g->regs, roundup_times(regs->size + 1, 8)); // cannot allocate memory when reg->size is very big, weird
		ZEROS(next_ref_regv(g->regs));
	}
	thread_beg_init(mupd, ncpu);
	mupd->g = g;
	mupd->regs = regs;
	mupd->nd = NULL;
	mupd->ins = init_u8v(64);
	mupd->vt = 0;
	thread_end_init(mupd);
	vt = 1;
	upds[0] = upds[1] = upds[2] = 0;
	inc = 0;
	for(i=0;i<1024&&i<nds->size;i++) inc += nds->buffer[i].cnt;
	for(idx=0;idx<nds->size+ncpu;idx++){
		if(g->regs->size + inc > g->regs->cap){
			thread_wait_all(mupd);
			encap_regv(g->regs, inc);
		}
		thread_wait_next(mupd);
		if(mupd->nd){
			pass = 0;
			for(j=0;j<mupd->nd->cnt;j++){
				hit = ref_rdregv(regs, mupd->nd->idx + j);
				rd = ref_readv(g->reads, hit->rid);
				if(rd->visit >= mupd->vt){ // previous node had changed the layout of this read, need to check conflict again
					mupd->ins->buffer[j] = check_read_reg_conflict_core(g, hit, &conflict);
					if(conflict){
						hit->closed = 1;
					}
				}
				rd->visit = idx;
				if(hit->closed == 0) pass ++;
			}
			do {
				closed = 0;
				if(0 && mupd->nd->cnt > g->max_node_cov){
					// Repetitive nodes should be in graph to avoid mis-assembling
					upds[1] ++;
					closed = 1;
					continue;
				}
				if(mupd->nd->fix){
					if(pass == 0){
						upds[1] ++;
						closed = 1;
						continue;
					} else {
						upds[0] ++;
					}
				} else {
					if(pass < g->min_node_cov || pass < (u4i)(mupd->nd->cnt * (1 - g->node_max_conflict))){
						upds[1] ++;
						closed = 1;
						continue;
					} else {
						upds[0] ++;
					}
				}
				n = next_ref_nodev(g->nodes);
				n->rep_idx = MAX_REP_IDX;
				n->unvisit = 0;
				n->closed = 0;
				n->single_in = 0;
				n->bt_visit = 0;
				n->bt_idx = 0;
				n->init_end = 0;
				n->regs.idx = g->regs->size;
				n->regs.cnt = 0;
				n->cov = 0;
				n->edges[0] = PTR_REF_NULL;
				n->edges[1] = PTR_REF_NULL;
				for(j=0;j<mupd->nd->cnt;j++){
					hit = ref_rdregv(regs, mupd->nd->idx + j);
					rd = ref_readv(g->reads, hit->rid);
					n->regs.cnt ++;
					reg = next_ref_regv(g->regs); // Now, it is thread-safe
					reg->node = g->nodes->size - 1;
					reg->rid  = hit->rid;
					reg->dir  = hit->dir;
					reg->beg  = hit->beg;
					reg->end  = hit->end;
					reg->closed = hit->closed | closed;
					if(reg->closed == 0) n->cov ++;
					reg->read_link = 0;
					if(mupd->ins->buffer[j]){
						r = ref_regv(g->regs, mupd->ins->buffer[j]);
						reg->read_link = r->read_link;
						r->read_link = g->regs->size - 1;
					} else {
						reg->read_link = rd->regs.idx;
						rd->regs.idx = g->regs->size - 1;
					}
					rd->regs.cnt ++;
				}
			} while(0);
		}
		if(idx < nds->size){
			mupd->nd = ref_rnkrefv(nds, idx);
			mupd->vt = idx;
			thread_wake(mupd);
		} else {
			mupd->nd = NULL;
		}
	}
	thread_beg_close(mupd);
	free_u8v(mupd->ins);
	thread_end_close(mupd);
}

static inline u8i chainning_hits_core(kbmmapv *hits, BitsVec *cigars, int uniq_hit, float aln_var){
	kbm_map_t HIT;
	u8i idx, bst, ite, lst, lsx, nrm;
	u4i state;
	sort_array(hits->buffer, hits->size, kbm_map_t, num_cmpgtxx(a.qidx, b.qidx, a.tidx, b.tidx, a.qdir, b.qdir));
	nrm = 0;
	for(idx=lst=lsx=0;idx<=hits->size;idx++){
		state = 0;
		if(idx == hits->size){
			state = 1;
		} else if(hits->buffer[lst].qidx == hits->buffer[idx].qidx && hits->buffer[lst].tidx == hits->buffer[idx].tidx){
			if(hits->buffer[lst].qdir == hits->buffer[idx].qdir){
				state = 0;
			} else {
				state = 2;
			}
		} else {
			state = 1;
		}
		if(state){
			if(idx > lst + 1){
				if(simple_chain_all_maps_kbm(hits->buffer + lst, idx - lst, cigars, &HIT, cigars, aln_var)){
					hits->buffer[lst++] = HIT;
					while(lst < idx){
						hits->buffer[lst++].mat = 0; // closed = 1
						nrm ++;
					}
				}
			}
			if(state == 1){
				// Choose the best hit
				if(uniq_hit && idx > lsx + 1){
					bst = lsx;
					for(bst=lsx;bst<idx;bst++){
						if(hits->buffer[bst].mat) break;
					}
					for(ite=bst+1;ite<idx;ite++){
						if(hits->buffer[ite].mat == 0) continue;
						if((hits->buffer[ite].aln > hits->buffer[bst].aln) || (hits->buffer[ite].aln == hits->buffer[bst].aln && hits->buffer[ite].mat > hits->buffer[bst].mat)){
							bst = ite;
						}
					}
					for(ite=lsx;ite<bst;ite++){
						hits->buffer[ite].mat = 0;
					}
					for(ite=bst+1;ite<idx;ite++){
						hits->buffer[ite].mat = 0;
					}
					nrm += idx - lsx - 1;
				}
				lsx = idx;
			}
			lst = idx;
		}
	}
	return nrm;
}

static inline u8i load_alignments_core(Graph *g, FileReader *pws, int raw, rdregv *regs, u4v *maps[3]){
	kbmmapv *hits;
	BitsVec *cigars;
	u4v *cgs;
	kbm_map_t *hit, HIT, *h;
	cuhash_t *cu;
	u8i nhit;
	u4i i, rid;
	int qlen, val, flg, nwarn, mwarn, ncol;
	char *cgstr, *qtag;
	mwarn = 20;
	nwarn = 0;
	cgs = init_u4v(4);
	hits = init_kbmmapv(32);
	cigars = g->chainning_hits? init_bitsvec(1024, 3) : g->cigars;
	memset(&HIT, 0, sizeof(kbm_map_t));
	hit = &HIT;
	nhit = 0;
	while((ncol = readtable_filereader(pws)) != -1){
		if((pws->n_line % 100000) == 0){
			fprintf(KBM_LOGF, "\r%llu", pws->n_line); fflush(KBM_LOGF);
		}
		if(pws->line->buffer[0] == '#'){
			//if(strncasecmp(pws->line->buffer, "#corr_mode=1", 12) == 0){
				//g->corr_mode = 1;
			//}
			continue;
		}
		if(ncol < 15) continue;
		if((cu = get_cuhash(g->kbm->tag2idx, get_col_str(pws, 0))) == NULL){
			if(nwarn < mwarn){
				fprintf(stderr, " -- Cannot find read \"%s\" in LINE:%llu in %s -- %s:%d --\n", get_col_str(pws, 0), pws->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				nwarn ++;
			}
			continue;
		} else rid = cu->val;
		hit->qidx = rid;
		qtag = get_col_str(pws, 0);
		hit->qdir = (get_col_str(pws, 1)[0] == '-');
		qlen = atoi(get_col_str(pws, 2));
		//if(g->corr_mode){
		if(0){
			//g->reads->buffer[hit->qidx].corr_bincnt = qlen / KBM_BIN_SIZE;
		} else if(qlen != (int)g->kbm->reads->buffer[hit->qidx].bincnt * KBM_BIN_SIZE){
			if(nwarn < mwarn){
				fprintf(stderr, " -- inconsisitent read length \"%s\" %d != %d in %s -- %s:%d --\n", qtag, qlen, g->kbm->reads->buffer[hit->qidx].bincnt * KBM_BIN_SIZE, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				nwarn ++;
			}
		}
		hit->qb = atoi(get_col_str(pws, 3)) / KBM_BIN_SIZE;
		hit->qe = atoi(get_col_str(pws, 4)) / KBM_BIN_SIZE;
		if((cu = get_cuhash(g->kbm->tag2idx, get_col_str(pws, 5))) == NULL){
			if(nwarn < mwarn){
				fprintf(stderr, " -- Cannot find read \"%s\" in LINE:%llu in %s -- %s:%d --\n", get_col_str(pws, 5), pws->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				nwarn ++;
			}
			continue;
		} else rid = cu->val;
		if(hit->qidx >= rid) continue;
		hit->tidx = rid;
		hit->tdir = 0; // col 6 always eq '+'
		// skil col 7
		hit->tb = atoi(get_col_str(pws, 8)) / KBM_BIN_SIZE;
		hit->te = atoi(get_col_str(pws, 9)) / KBM_BIN_SIZE;
		// skip col 10-13
		hit->mat = atoi(get_col_str(pws, 10));
		if(hit->mat < g->par->min_mat) continue;
		hit->aln = atoi(get_col_str(pws, 11)) / KBM_BIN_SIZE;
		if(hit->aln < g->par->min_aln) continue;
		if(hit->mat < hit->aln * KBM_BIN_SIZE * g->par->min_sim) continue;
		if(num_diff(hit->qe - hit->qb, hit->te - hit->tb) > (int)num_max(g->par->aln_var * hit->aln, 1.0)) continue;
		hit->cnt = atoi(get_col_str(pws, 12));
		hit->gap = atoi(get_col_str(pws, 13));
		if(g->chainning_hits){
			if(hits->size && peer_kbmmapv(hits)->qidx != hit->qidx){
				chainning_hits_core(hits, cigars, g->uniq_hit, g->kbm->par->aln_var);
				for(i=0;i<hits->size;i++){
					h = ref_kbmmapv(hits, i);
					if(h->mat == 0) continue;
					append_bitsvec(g->cigars, cigars, h->cgoff, h->cglen);
					h->cgoff = g->cigars->size - h->cglen;
					nhit ++;
					map2rdhits_graph(g, h);
				}
				clear_kbmmapv(hits);
				clear_bitsvec(cigars);
			}
		}
		hit->cgoff = cigars->size;
		clear_u4v(cgs);
		cgstr = get_col_str(pws, 14);
		if(cgstr[0] == '*'){ // no cigar
			continue;
		}
		val = 0;
		while(cgstr[0]){
			if(cgstr[0] >= '0' && cgstr[0] <= '9'){
				val = val * 10 + (cgstr[0] - '0');
			} else {
				flg = -1;
				switch(cgstr[0]){
					case 'M': flg = 0; break;
					case 'I': flg = 1; break;
					case 'D': flg = 2; break;
					case 'm': flg = 4; break;
					case 'i': flg = 5; break;
					case 'd': flg = 6; break;
					default:
						fprintf(stderr, " -- Bad cigar '%c' \"%s\" in LINE:%llu in %s -- %s:%d --\n", cgstr[0], get_col_str(pws, 14), pws->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						exit(1);
				}
				if(val == 0) val = 1;
				push_u4v(cgs, (val << 8) | flg);
				val = 0;
			}
			cgstr ++;
		}
		for(i=cgs->size;i>0;i--){
			val = cgs->buffer[i - 1] >> 8;
			flg = cgs->buffer[i - 1] & 0xFF;
			pushs_bitsvec(cigars, flg, val);
		}
		hit->cglen = cigars->size - hit->cgoff;
		if(raw){
			hit2rdregs_graph(g, regs, qlen / KBM_BIN_SIZE, hit, cigars, maps);
			clear_bitsvec(cigars);
		} else if(g->chainning_hits){
			push_kbmmapv(hits, *hit);
		} else {
			map2rdhits_graph(g, hit);
		}
	}
	if(g->chainning_hits){
			if(hits->size){
				chainning_hits_core(hits, cigars, g->uniq_hit, g->kbm->par->aln_var);
				for(i=0;i<hits->size;i++){
					h = ref_kbmmapv(hits, i);
					if(h->mat == 0) continue;
					append_bitsvec(g->cigars, cigars, h->cgoff, h->cglen);
					h->cgoff = g->cigars->size - h->cglen;
					nhit ++;
					map2rdhits_graph(g, h);
				}
				clear_kbmmapv(hits);
				clear_bitsvec(cigars);
			}
	}
	fprintf(KBM_LOGF, "\r%llu lines, %llu hits\n", pws->n_line, nhit);
	free_kbmmapv(hits);
	if(g->chainning_hits) free_bitsvec(cigars);
	free_u4v(cgs);
	return nhit;
}

static inline u8i proc_alignments_core(Graph *g, int ncpu, int raw, rdregv *regs, u4v *maps[3], char *prefix, char *dump_kbm){
	kbm_map_t *hit;
	kbm_read_t *pb;
	BitVec *rdflags;
	BufferedWriter *bw;
	FILE *alno, *kmlog;
	u8i nbp, mbp, nhit;
	u8i i, ib, ie, ic;
	u4i rid, qb, qe, ii, in;
	int reset_kbm, n_cpu;
	thread_prepare(mdbg);
	if(KBM_LOG) n_cpu = 1;
	else n_cpu = ncpu;
	ic = (g->kbm->bins->size + g->num_index - 1) / g->num_index;
	ie = 0;
	UNUSED(mbp);
#if 0
	if(g->corr_mode){
		mbp = g->genome_size * g->corr_gcov;
		qb = qe = g->kbm->reads->size / 2;
		nbp = g->kbm->reads->buffer[qb].bincnt * KBM_BIN_SIZE;
		while(nbp < mbp && qb && qe + 1 < g->kbm->reads->size){
			qb --;
			qe ++;
			nbp += g->kbm->reads->buffer[qb].bincnt * KBM_BIN_SIZE;
			nbp += g->kbm->reads->buffer[qe].bincnt * KBM_BIN_SIZE;
		}
		if(qe < g->kbm->reads->size) qe ++;
		fprintf(KBM_LOGF, "[%s] turn correct-mode on, reads[%u ~ %u = %u] (%llu bp), genome-size=%llu, corr-gcov=%0.2f, corr-dep=[%d,%d,%0.2f]\n", date(), qb, qe, qe - qb, nbp, g->genome_size, g->corr_gcov, g->corr_min, g->corr_max, g->corr_cov); fflush(KBM_LOGF);
	} else {
#else
	{
#endif
		qb = 0;
		qe = g->reads->size;
	}
	alno = open_file_for_write(prefix, ".alignments.gz", 1);
	bw = zopen_bufferedwriter(alno, 1024 * 1024, ncpu, 0);
#if 0
	if(g->corr_mode){
		beg_bufferedwriter(bw);
		fprintf(bw->out, "#corr_mode=1\n");
		end_bufferedwriter(bw);
	}
	rdflags = (!g->corr_mode && g->par->skip_contained)? init_bitvec(g->kbm->reads->size) : NULL;
	in = g->corr_mode? 1 : g->num_index;
#else
	rdflags = (g->par->skip_contained)? init_bitvec(g->kbm->reads->size) : NULL;
	in = g->num_index;
#endif
	if(g->kbm->seeds->size){
		reset_kbm = 0;
		if(in > 1){
			fprintf(KBM_LOGF, " ** WARNNING: change number of kbm index to 1 **\n"); fflush(KBM_LOGF);
			in = 1;
		}
	} else {
		reset_kbm = 1;
	}
	//fix_node = 0;
	nhit = 0;
	for(ii=0;ii<in;ii++){
		ib = ie;
		ie = num_min(ib + ic, g->kbm->bins->size);
		while(ie > 0 && ie < g->kbm->bins->size && g->kbm->bins->buffer[ie - 1].ridx == g->kbm->bins->buffer[ie].ridx) ie ++;
#if 0
		if(g->corr_mode == 0){
#else
		{
#endif
			qb = 0;
			qe = ie? g->kbm->bins->buffer[ie - 1].ridx : 0;
		}
		nbp = ((u8i)(ie - ib)) * KBM_BIN_SIZE;
		if(reset_kbm){
			reset_index_kbm(g->kbm);
			fprintf(KBM_LOGF, "[%s] indexing bins[(%llu,%llu)/%llu] (%llu/%llu bp), %d threads\n", date(), ib, ie, (u8i)g->kbm->bins->size, nbp, (u8i)g->kbm->rdseqs->size, ncpu); fflush(KBM_LOGF);
			kmlog = (in > 1)? NULL : open_file_for_write(prefix, ".kmerdep", 1);
			index_kbm(g->kbm, ib, ie, ncpu, kmlog);
			if(kmlog){
				fclose(kmlog);
				kmlog = NULL;
			}
			fprintf(KBM_LOGF, "[%s] Done\n", date()); fflush(KBM_LOGF);
			if(in == 1 && dump_kbm){
				FILE *dump;
				fprintf(KBM_LOGF, "[%s] dump kbm-index to %s ...", date(), dump_kbm); fflush(KBM_LOGF);
				dump = open_file_for_write(dump_kbm, NULL, 1);
				mem_dump_obj_file(g->kbm, 1, &kbm_obj_desc, 1, 0, dump);
				fclose(dump);
				fprintf(KBM_LOGF, " Done\n"); fflush(KBM_LOGF);
			}
			if(in == 1){
				u4i *deps;
				u8i hidx;
				kmlog = open_file_for_write(prefix, ".binkmer", 1);
				deps = calloc(KBM_BIN_SIZE + 1, 4);
				for(hidx=0;hidx<g->kbm->bins->size;hidx++){
					deps[g->kbm->bins->buffer[hidx].degree] ++;
				}
				for(hidx=0;hidx<KBM_BIN_SIZE;hidx++){
					fprintf(kmlog, "%u\n", deps[hidx]);
				}
				fclose(kmlog);
				free(deps);
				if(!g->minimal_output && 0){
					kbm_bin_t *bn;
					kmlog = open_file_for_write(prefix, ".closed_bins", 1);
					for(hidx=0;hidx<g->kbm->bins->size;hidx++){
						bn = ref_kbmbinv(g->kbm->bins, hidx);
						if(bn->closed){
							fprintf(kmlog, "%s_F_%d_%d\t%d\n", g->kbm->reads->buffer[bn->ridx].tag, bn->off * KBM_BIN_SIZE, KBM_BIN_SIZE, bn->degree);
						}
					}
					fclose(kmlog);
				}
			}
		}
		{
			{
				thread_beg_init(mdbg, n_cpu);
				mdbg->g = g;
				memset((void*)&mdbg->reg, 0, sizeof(reg_t));
				mdbg->reg.closed = 1;
				mdbg->aux = init_kbmaux(g->kbm);
				if(g->rpar){
					mdbg->raux = init_kbmaux(init_kbm(g->rpar));
				} else {
					mdbg->raux = NULL;
				}
#if 0
				if(g->corr_mode){
					KBMBlock *kb;
					POGPar par;
					kb = init_kbmblock(g->corr_bsize, g->corr_bstep);
					//mdbg->cc = init_ctgcns(kb, iter_kbmblock, info_kbmblock, 1, 1, 1, g->corr_max, 200, 100, 1, 96, 2, -5, -2, -4, -1, 16, 3, 0.5, g->corr_bsize - g->corr_bstep + KBM_BIN_SIZE);
					par = DEFAULT_POG_PAR;
					par.refmode = 1;
					mdbg->cc = init_ctgcns(kb, iter_kbmblock, info_kbmblock, 1, 1, g->corr_max, 200, 100, 1, g->corr_bsize - g->corr_bstep + KBM_BIN_SIZE, &par);
				} else {
#else
				{
#endif
					mdbg->cc = NULL;
				}
				mdbg->aux->par = (KBMPar*)malloc(sizeof(KBMPar));
				memcpy(mdbg->aux->par, g->par, sizeof(KBMPar));
				mdbg->regs = regs;
				mdbg->rdflags = rdflags;
				mdbg->beg = 0;
				mdbg->end = 0;
				mdbg->raw = raw;
				mdbg->alno = alno;
				thread_end_init(mdbg);
			}
			thread_beg_iter(mdbg);
			mdbg->task = 1;
			thread_end_iter(mdbg);
			for(rid=qb;rid<=qe+ncpu;rid++){
				if(rid < qe){
					if(!KBM_LOG && ((rid - qb) % 2000) == 0){ fprintf(KBM_LOGF, "\r%u|%llu", rid - qb, nhit); fflush(KBM_LOGF); }
					thread_wait_one(mdbg);
				} else {
					thread_wait_next(mdbg);
					pb = NULL;
				}
				if(mdbg->reg.closed == 0){
					KBMAux *aux = mdbg->aux;
					//if(g->corr_mode && mdbg->cc->cns->size){
						//g->reads->buffer[mdbg->reg.rid].corr_bincnt = mdbg->cc->cns->size / KBM_BIN_SIZE;
					//}
					if(alno){
						beg_bufferedwriter(bw);
						//if(g->corr_mode && mdbg->aux->cns->size){
							//fprintf(bw->out, "#corrected\t%s\t%u\t", mdbg->cc->tag->string, (u4i)mdbg->cc->cns->size);
							//println_fwdseq_basebank(mdbg->cc->cns, 0, mdbg->cc->cns->size, bw->out);
						//}
						for(i=0;i<mdbg->aux->hits->size;i++){
							hit = ref_kbmmapv(mdbg->aux->hits, i);
							fprint_hit_kbm(mdbg->aux, i, bw->out);
						}
						end_bufferedwriter(bw);
					}
					for(i=0;i<mdbg->aux->hits->size;i++){
						hit = ref_kbmmapv(mdbg->aux->hits, i);
						if(hit->mat == 0) continue;
						if(rdflags
							&& g->kbm->reads->buffer[hit->tidx].bincnt < g->kbm->reads->buffer[hit->qidx].bincnt
							&& (hit->tb <= 1 && hit->te + 1 >= (int)(g->kbm->reads->buffer[hit->tidx].bincnt))
							&& (hit->qb > 1 || hit->qe + 1 < (int)(g->kbm->reads->buffer[hit->qidx].bincnt))
							){
							one_bitvec(rdflags, hit->tidx);
						}
					}
					if(g->chainning_hits){
						chainning_hits_core(aux->hits, aux->cigars, g->uniq_hit, g->kbm->par->aln_var);
					}
					for(i=0;i<aux->hits->size;i++){
						hit = ref_kbmmapv(aux->hits, i);
						if(hit->mat == 0) continue;
						//hit->qb  /= KBM_BIN_SIZE;
						//hit->qe  /= KBM_BIN_SIZE;
						//hit->tb  /= KBM_BIN_SIZE;
						//hit->te  /= KBM_BIN_SIZE;
						//hit->aln /= KBM_BIN_SIZE;
						nhit ++;
						append_bitsvec(g->cigars, aux->cigars, hit->cgoff, hit->cglen);
						hit->cgoff = g->cigars->size - hit->cglen;
						if(raw){
							//hit2rdregs_graph(g, regs, g->corr_mode? mdbg->cc->cns->size / KBM_BIN_SIZE : 0, hit, mdbg->aux->cigars, maps);
							hit2rdregs_graph(g, regs, 0, hit, mdbg->aux->cigars, maps);
						} else {
							map2rdhits_graph(g, hit);
						}
					}
					if(KBM_LOG){
						fprintf(KBM_LOGF, "QUERY: %s\t+\t%d\t%d\n", g->kbm->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end);
						for(i=0;i<mdbg->aux->hits->size;i++){
							hit = ref_kbmmapv(mdbg->aux->hits, i);
								fprintf(KBM_LOGF, "\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->kbm->reads->buffer[hit->tidx].tag, "+-"[hit->qdir], g->kbm->reads->buffer[hit->tidx].bincnt * KBM_BIN_SIZE, hit->tb * KBM_BIN_SIZE, hit->te * KBM_BIN_SIZE, hit->aln * KBM_BIN_SIZE, hit->mat);
						}
					}
					mdbg->reg.closed = 1;
				}
				if(rid < qe && (rdflags == NULL || get_bitvec(rdflags, rid) == 0)){
					pb = ref_kbmreadv(g->kbm->reads, rid);
					mdbg->reg = (reg_t){0, rid, 0, 0, pb->bincnt, 0, 0};
					thread_wake(mdbg);
				}
			}
			{
				thread_beg_close(mdbg);
				free(mdbg->aux->par);
				free_kbmaux(mdbg->aux);
				//if(g->corr_mode){
					//free_kbmblock((KBMBlock*)mdbg->cc->obj);
					//free_ctgcns(mdbg->cc);
				//}
				if(mdbg->raux){
					free_kbm(mdbg->raux->kbm);
					free_kbmaux(mdbg->raux);
				}
				thread_end_close(mdbg);
			}
		}
		if(!KBM_LOG) fprintf(KBM_LOGF, "\r%u reads|total hits %llu\n", qe - qb, nhit);
		if(reset_kbm){
			reset_index_kbm(g->kbm);
		}
	}
	if(bw) close_bufferedwriter(bw);
	if(alno) fclose(alno);
	if(rdflags) free_bitvec(rdflags);
	return nhit;
}

static inline void build_nodes_graph(Graph *g, u8i maxbp, int ncpu, FileReader *pws, int rdclip, char *prefix, char *dump_kbm){
	kbm_map_t *hit, HIT;
	BitVec *rks;
	u4v *maps[3];
	FILE *clplog;
	u8i idx, rank, kcnts[256], upds[3], fix_node, nhit;
	u4i rid, i, dep;
	int raw;
	rdregv *regs;
	rnkrefv *nds;
	rnk_ref_t *nd;
	thread_preprocess(mclp);
	regs = init_rdregv(1024);
	nds = init_rnkrefv(1024);
	renew_rdhitv(g->rdhits, 1024);
	maps[0] = init_u4v(4);
	maps[1] = init_u4v(4);
	maps[2] = init_u4v(4);
	clear_regv(g->regs);
	next_ref_regv(g->regs);
	clear_rdhitv(g->rdhits);
	ZEROS(next_ref_rdhitv(g->rdhits));
	//clear_kbmmapv(g->pwalns);
	clear_bitsvec(g->cigars);
	raw = !((g->bestn > 0) || rdclip);
	//if(g->corr_mode){
		//raw = 1;
	//}
	fix_node = 0; // TODO: new code hasn't coped with contigs+longreads mode
	if(pws){
		nhit = load_alignments_core(g, pws, raw, regs, maps);
	} else {
		UNUSED(maxbp);
		nhit = proc_alignments_core(g, ncpu, raw, regs, maps, prefix, dump_kbm);
	}
	print_proc_stat_info(0);
	if(raw){
		fprintf(KBM_LOGF, "[%s] generated %llu regs\n", date(), (u8i)regs->size);
	} else {
#ifdef KBM_MAP2RDHIT_QUICK
		// sort rdhits and pick bestn
		{
			fprintf(KBM_LOGF, "[%s] sorting rdhits ... ", date()); fflush(KBM_LOGF);
			thread_beg_init(mclp, ncpu);
			mclp->g = g;
			mclp->task = 2;
			thread_end_init(mclp);
			thread_wake_all(mclp);
			thread_wait_all(mclp);
			thread_beg_close(mclp);
			thread_end_close(mclp);
			fprintf(KBM_LOGF, "Done\n"); fflush(KBM_LOGF);
		}
#endif
		// clip reads
		if(rdclip){
			fprintf(KBM_LOGF, "[%s] clipping ... ", date()); fflush(KBM_LOGF);
			clplog = open_file_for_write(prefix, ".clps", 1);
			thread_beg_init(mclp, ncpu);
			mclp->g = g;
			mclp->task = 1;
			thread_end_init(mclp);
			thread_wake_all(mclp);
			thread_wait_all(mclp);
			thread_beg_close(mclp);
			thread_end_close(mclp);
			u8i tot, clp;
			tot = clp = 0;
			for(rid=0;rid<g->reads->size;rid++){
				tot += g->kbm->reads->buffer[rid].bincnt * KBM_BIN_SIZE;
				clp += (g->reads->buffer[rid].clps[1] - g->reads->buffer[rid].clps[0]) * KBM_BIN_SIZE;
				fprintf(clplog, "%s\t%d\t%d\t%d\n", g->kbm->reads->buffer[rid].tag, g->kbm->reads->buffer[rid].bincnt * KBM_BIN_SIZE, g->reads->buffer[rid].clps[0] * KBM_BIN_SIZE, g->reads->buffer[rid].clps[1] * KBM_BIN_SIZE);
			}
			fclose(clplog);
			fprintf(KBM_LOGF, "%.2f%% bases\n", ((tot - clp) * 100.0) / tot); fflush(KBM_LOGF);
		}
		fprintf(KBM_LOGF, "[%s] generating regs ... ", date()); fflush(KBM_LOGF);
		rd_hit_t *rh;
		u8i cgoff;
		cgoff = 0;
		hit = &HIT;
		hit->cnt = 0;
		hit->gap = 0;
		for(idx=0;idx<g->rdhits->size;idx++){
			rh = ref_rdhitv(g->rdhits, idx);
			hit->cgoff = cgoff;
			hit->cglen = rh->lnks[1].cnt;
			cgoff += hit->cglen;
			if(rh->frgs[0].closed && rh->frgs[1].closed){
				continue;
			}
			hit->qidx = rh->frgs[0].rid;
			hit->tidx = rh->frgs[1].rid;
			hit->qdir = rh->frgs[0].dir;
			hit->tdir = rh->frgs[1].dir;
			hit->qb = rh->frgs[0].beg;
			hit->qe = rh->frgs[0].end;
			hit->tb = rh->frgs[1].beg;
			hit->te = rh->frgs[1].end;
			hit->mat = rh->lnks[0].cnt;
			hit->aln = num_min(hit->qe - hit->qb, hit->te - hit->tb);
			//hit2rdregs_graph(g, regs, g->reads->buffer[hit->qidx].corr_bincnt, hit, g->cigars, maps);
			hit2rdregs_graph(g, regs, g->kbm->reads->buffer[hit->qidx].bincnt, hit, g->cigars, maps);
		}
		free_bitsvec(g->cigars); g->cigars = init_bitsvec(1024, 3);
		fprintf(KBM_LOGF, "%llu\n", (u8i)regs->size); fflush(KBM_LOGF);
	}
	free_u4v(maps[0]);
	free_u4v(maps[1]);
	free_u4v(maps[2]);
	// add node itself
	//if(!g->corr_mode){
	{
		rks = init_bitvec(g->kbm->bins->size);
		for(idx=0;idx<regs->size;idx++){
			one_bitvec(rks, regs->buffer[idx].node);
		}
		for(idx=0;idx<g->kbm->bins->size;idx++){
			if(get_bitvec(rks, idx)){
				rid = g->kbm->bins->buffer[idx].ridx;
				kbm_read_t *rd = ref_kbmreadv(g->kbm->reads, rid);
				push_rdregv(regs, (rd_reg_t){idx, rid, 0, (idx - rd->binoff) * g->reglen, (idx + 1 - rd->binoff) * g->reglen, 0});
			}
		}
		free_bitvec(rks);
	}
	// generating nodes
	fprintf(KBM_LOGF, "[%s] sorting regs ... ", date()); fflush(KBM_LOGF);
	psort_array(regs->buffer, regs->size, rd_reg_t, ncpu, num_cmpgtxx((a.node << 30) | a.rid, (b.node << 30) | b.rid, a.beg, b.beg, b.end, a.end));
	fprintf(KBM_LOGF, " Done\n"); fflush(KBM_LOGF);
	u4v *kbcnts;
	kbcnts = init_u4v(1024);
	rank = 0xFFFFFFFFFFFFFFFFLLU;
	nd = NULL;
	fprintf(KBM_LOGF, "[%s] generating intervals ... ", date()); fflush(KBM_LOGF);
	for(idx=0;idx<regs->size;idx++){
		if(0){
			kbm_read_t *rd1, *rd2;
			rd_reg_t *r;
			r = ref_rdregv(regs, idx);
			rd1 = ref_kbmreadv(g->kbm->reads, g->kbm->bins->buffer[r->node].ridx);
			rd2 = ref_kbmreadv(g->kbm->reads, r->rid);
			fprintf(stdout, "%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\n", rd1->tag, '+', (int)((r->node) - rd1->binoff) * g->reglen,
					(int)(r->node - rd1->binoff + 1) * g->reglen,
					rd2->tag, "+-"[r->dir], r->beg, r->end);
		}
		if(regs->buffer[idx].node != rank){
			if(nd){
				while(nd->cnt >= kbcnts->size) push_u4v(kbcnts, 0);
				kbcnts->buffer[nd->cnt] ++;
			}
			nd = next_ref_rnkrefv(nds);
			nd->idx = idx;
			nd->rank = regs->buffer[idx].node;
			nd->fix = regs->buffer[idx].node < fix_node;
			nd->cnt = 1;
			//nd->score = (regs->buffer[idx].beg >= g->reads->buffer[regs->buffer[idx].rid].clps[0] && regs->buffer[idx].end <= g->reads->buffer[regs->buffer[idx].rid].clps[1]);
			rank = regs->buffer[idx].node;
		} else {
			nd->cnt ++;
			//nd->score += (regs->buffer[idx].beg >= g->reads->buffer[regs->buffer[idx].rid].clps[0] && regs->buffer[idx].end <= g->reads->buffer[regs->buffer[idx].rid].clps[1]);
		}
	}
	if(nd){
		while(nd->cnt >= kbcnts->size) push_u4v(kbcnts, 0);
		kbcnts->buffer[nd->cnt] ++;
	}
	// find medean k-bin depth
	u4i gcov;
	for(gcov=0,rank=0;gcov<kbcnts->size&&rank<nds->size;gcov++){
		rank += kbcnts->buffer[gcov];
	}
	free_u4v(kbcnts);
	//psort_array(nds->buffer, nds->size, rnk_ref_t, ncpu, num_cmpgtx(b.cnt, a.cnt, a.rank, b.rank));
	psort_array(nds->buffer, nds->size, rnk_ref_t, ncpu, num_cmpgtx(num_diff(a.cnt, gcov), num_diff(b.cnt, gcov), a.rank, b.rank));
	fprintf(KBM_LOGF, " %llu intervals\n", (u8i)nds->size); fflush(KBM_LOGF);
	//psort_array(nds->buffer, nds->size, rnk_ref_t, ncpu, num_cmpgtx(b.score, a.score, a.rank, b.rank));
	if(0){
		memset(kcnts, 0, sizeof(uint64_t) * 256);
		for(idx=0;idx<nds->size;idx++){
			dep = num_min(nds->buffer[idx].cnt, 255);
			kcnts[dep] ++;
		}
		for(i=1;i<51;i++){
			fprintf(KBM_LOGF, "%10llu", (long long unsigned int)kcnts[i]);
			if(((i - 1) % 10) == 9) fprintf(KBM_LOGF, "\n");
		}
		if(((i - 1) % 10) != 0) fprintf(KBM_LOGF, "\n");
	}
	fprintf(KBM_LOGF, "[%s] selecting important intervals from %llu intervals\n", date(), (u8i)nds->size);
	mul_update_regs_graph(g, regs, nds, ncpu, upds);
	free_rdregv(regs);
	free_rnkrefv(nds);
	encap_regv(g->regs, 1);
	memset(g->regs->buffer + g->regs->size, 0xFF, sizeof(reg_t));
	if(!KBM_LOG) fprintf(KBM_LOGF, "[%s] Intervals: kept %llu, discarded %llu\n", date(), upds[0], upds[1]);
	print_proc_stat_info(0);
}

uint32_t estimate_genome_size(Graph *g, unsigned long long tot_bp, FILE *out){
	uint64_t kcnts[256];
	node_t *n;
	uint64_t nid, sum, ncnt, pmax;
	uint32_t i, dep, peak, mid;
	float avg;
	ncnt = g->nodes->size;
	memset(kcnts, 0, sizeof(uint64_t) * 256);
	sum = 0;
	for(nid=0;nid<ncnt;nid++){
		n = ref_nodev(g->nodes, nid);
		dep = num_min(n->regs.cnt, 255);
		sum += n->regs.cnt;
		kcnts[dep] ++;
	}
	mid = pmax = 0;
	while(mid < 255){
		pmax += kcnts[mid];
		if(pmax >= ncnt / 2) break;
		mid ++;
	}
	avg = 1.0 * sum / (ncnt + 1);
	fprintf(out, "[%s] median node depth = %d\n", date(), mid);
	return mid;
	//TODO: calculate the genome coverage
	for(i=1;i<51;i++){
		fprintf(out, "%10llu", (long long unsigned int)kcnts[i]);
		if(((i - 1) % 10) == 9) fprintf(out, "\n");
	}
	if(((i - 1) % 10) != 0) fprintf(out, "\n");
	return avg;
	pmax = 0; peak = avg;
	for(i=g->min_node_cov+1;i<256;i++){
		if(kcnts[i] > pmax){ peak = i; pmax = kcnts[i]; }
		else if(i > avg && kcnts[i] < 0.8 * pmax) break;
	}
	fprintf(out, "[%s] depth peak = %d\n", date(), peak);
	fprintf(out, "[%s] genome size = %llu\n", date(), tot_bp / peak);
	return peak;
}

// MUST be called before build_edges
static u8i mask_nodes_by_cov_graph(Graph *g, FILE *out){
	node_t *n;
	u8i ret, i;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->regs.cnt > g->max_node_cov || n->regs.cnt < g->min_node_cov){
			n->closed = 1;
			ret ++;
			if(out) fprintf(out, "MASK_COV\tN%llu\t%u\t%u\n", (u8i)i, (u4i)n->regs.cnt, n->cov);
		}
	}
	return ret;
}

static inline void remove_all_edges_graph(Graph *g){
	node_t *n;
	uint64_t nid;
	free_edgev(g->edges);
	g->edges = init_edgev(32);
	free_edgehash(g->ehash);
	g->ehash = init_edgehash(1023);
	set_userdata_edgehash(g->ehash, g->edges);
	free_edgerefv(g->erefs);
	g->erefs = init_edgerefv(32);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		n->edges[0] = n->edges[1] = PTR_REF_NULL;
	}
}

typedef struct {
	uint64_t idx:46; int off:18;
} edge_off_t;
define_list(edgeoffv, edge_off_t);

static inline int estimate_edge_length(edge_off_t *ps, uint32_t size, uint32_t idxs[2]){
	int64_t tot, var;
	uint32_t i, b, e, mi;
	int v, max, len, avg, std;
	 idxs[0] = 0; idxs[1] = size;
	if(size == 0){ return 0; }
	if(size <= 2){ return ps[size/2].off; }
	b = 0; e = size;
	for(i=b,tot=0;i<e;i++) tot += ps[i].off;
	len = tot / (e - b);
	//fprintf(stdout, "b=%d\te=%d\tlen=%d\n", b, e, len);
	while(b + 2 < e){
		max = 0; mi = 0;
		for(i=b+1;i<e;i++){
			if(ps[i].off - ps[i-1].off > max){
				max = ps[i].off - ps[i-1].off;
				mi = i;
			}
		}
		if(max < len * 0.5) break;
		else if(max < 100) break;
		if(mi - b > e - mi) e = mi;
		else b = mi;
		for(i=b,tot=0;i<e;i++) tot += ps[i].off;
		avg = tot / (e - b);
		if(num_diff(avg, len) < num_max(avg * 0.2, 50)) break;
		len = avg;
	}
	//fprintf(stdout, "b=%d\te=%d\tlen=%d\n", b, e, len);
	if(0){
		if(b + 1 < e){
			for(i=b,var=0;i<e;i++){
				v = ((int)ps[i].off) - ((int)len);
				var += v * v;
			}
			std = sqrt(var / (e - b));
			//fprintf(stdout, "std=%d\n", std);
			b = 0; e = size;
			while(b < e && num_diff(ps[b].off, len) > 3 * std) b ++;
			while(b < e && num_diff(ps[e - 1].off, len) > 3 * std) e --;
		}
		idxs[0] = b;
		idxs[1] = e;
	}
	return len;
}

static inline int calculate_edge_cov_off_graph(Graph *g, edge_t *e, edgeoffv *offs){
	node_t *n1, *n2;
	reg_t *r1, *r2;
	uint32_t dir1, dir2, cov, idxs[2];
	int off;
	n1 = ref_nodev(g->nodes, e->node1);
	n2 = ref_nodev(g->nodes, e->node2);
	r1 = ref_regv(g->regs, n1->regs.idx);
	r2 = ref_regv(g->regs, n2->regs.idx);
	cov = e->cov;
	clear_edgeoffv(offs);
	while(1){
		if(r1->rid > r2->rid){
			r2 ++;
			if(r2->node != e->node2) break;
		} else if(r1->rid < r2->rid){
			r1 ++;
			if(r1->node != e->node1) break;
		} else {
			if(r1->beg < r2->beg){
				if(r1->node < r2->node){
					dir1 = r1->dir;
					dir2 = r2->dir;
				} else {
					dir1 = !r2->dir;
					dir2 = !r1->dir;
				}
				off = r2->beg - r1->end;
			} else {
				if(r2->node < r1->node){
					dir1 = r2->dir;
					dir2 = r1->dir;
				} else {
					dir1 = !r1->dir;
					dir2 = !r2->dir;
				}
				off = ((int)r1->beg) - r2->end;
			}
			if(dir1 == e->dir1 && dir2 == e->dir2){
				push_edgeoffv(offs, (edge_off_t){e - g->edges->buffer, off});
			}
			r1 ++;
			if(r1->node != e->node1) break;
			r2 ++;
			if(r2->node != e->node2) break;
		}
	}
	e->off = estimate_edge_length(offs->buffer, offs->size, idxs);
	e->cov = idxs[1] - idxs[0];
	if(cov != e->cov) return 1;
	else return 0;
}

static inline void build_edges_graph(Graph *g, int ncpu, FILE *log){
	read_t *rd;
	node_t *n;
	edge_t *E, *e;
	vplist *regs;
	reg_t *r1, *r2;
	edge_ref_t f1, f2;
	edgeoffv *offs;
	uint64_t idx, lst, cnt, x, *u;
	uint32_t rid, i, idxs[2];
	int exists;
	UNUSED(log);
	clear_edgev(g->edges);
	E = next_ref_edgev(g->edges);
	memset(E, 0, sizeof(edge_t));
	E->closed = WT_EDGE_CLOSED_MASK;
	E->cov = 1;
	//E->status = 0;
	clear_edgehash(g->ehash);
	offs = init_edgeoffv(32);
	regs = init_vplist(32);
	for(rid=0;rid<g->kbm->reads->size;rid++){
		rd = ref_readv(g->reads, rid);
		if(rd->regs.cnt < 2) continue;
		clear_vplist(regs);
		idx = rd->regs.idx;
		while(idx){
			r2 = ref_regv(g->regs, idx);
			idx = r2->read_link;
			if(g->nodes->buffer[r2->node].closed) continue;
			if(r2->closed) continue;
			if(!(r2->beg >= rd->clps[0] && r2->end <= rd->clps[1])) continue;
			for(i=regs->size;i>0;i--){
				r1 = (reg_t*)get_vplist(regs, i - 1);
				if(r1->end + g->max_edge_span < r2->beg) break;
				E = ref_edgev(g->edges, 0);
				if(r1->node < r2->node){
					E->node1 = r1->node;
					E->node2 = r2->node;
					E->dir1  = r1->dir;
					E->dir2  = r2->dir;
				} else {
					E->node1 = r2->node;
					E->node2 = r1->node;
					E->dir1  = !r2->dir;
					E->dir2  = !r1->dir;
				}
				E->off = ((int)r2->beg) - r1->end;
				u = prepare_edgehash(g->ehash, 0, &exists);
				if(exists){
					//if(g->edges->buffer[*u].cov < WT_MAX_EDGE_COV) g->edges->buffer[*u].cov ++;
				} else {
					*u = g->edges->size;
					push_edgev(g->edges, *E);
				}
				if(i == regs->size){
					e = ref_edgev(g->edges, *u);
					//if(e->cov < WT_MAX_EDGE_COV) e->cov ++;
					//if(e->cov >= 1){
						e->closed = 0;
					//}
				}
				push_edgeoffv(offs, (edge_off_t){*u, ((int)r2->beg) - r1->end});
			}
			push_vplist(regs, r2);
		}
	}
	free_vplist(regs);
	if(g->edges->size == 0){
		free_edgeoffv(offs);
		return;
	}
	psort_array(offs->buffer, offs->size, edge_off_t, ncpu, num_cmpgtx(a.idx, b.idx, a.off, b.off));
	lst = 0;
	for(idx=1;idx<=offs->size;idx++){
		if(idx < offs->size && offs->buffer[idx].idx == offs->buffer[lst].idx) continue;
		if(1){
			g->edges->buffer[offs->buffer[lst].idx].off = offs->buffer[(lst+idx)/2].off;
			g->edges->buffer[offs->buffer[lst].idx].cov = idx - lst;
		} else {
			g->edges->buffer[offs->buffer[lst].idx].off = estimate_edge_length(offs->buffer + lst, idx - lst, idxs);
			g->edges->buffer[offs->buffer[lst].idx].cov = idxs[1] - idxs[0];
		}
		if(0){
			uint64_t m;
			for(m=lst;m<idx;m++) fprintf(stdout, "%u\t", offs->buffer[m].off);
			fprintf(stdout, "\n");
			for(m=lst+idxs[0];m<lst+idxs[1];m++) fprintf(stdout, "%u\t", offs->buffer[m].off);
			fprintf(stdout, "\n");
		}
		lst = idx;
	}
	free_edgeoffv(offs);
	clear_edgerefv(g->erefs);
	push_edgerefv(g->erefs, EDGE_REF_NULL);
	for(idx=1;idx<g->edges->size;idx++){
		if(g->edges->buffer[idx].cov < g->min_edge_cov){
			if(g->store_low_cov_edge) g->edges->buffer[idx].closed = WT_EDGE_CLOSED_LESS;
			else continue;
		}
		if(g->nodes->buffer[g->edges->buffer[idx].node1].closed || g->nodes->buffer[g->edges->buffer[idx].node2].closed){
			g->edges->buffer[idx].closed = WT_EDGE_CLOSED_HARD;
		} else if(g->edges->buffer[idx].closed == WT_EDGE_CLOSED_MASK){
			g->edges->buffer[idx].closed = WT_EDGE_CLOSED_LESS;
		//} else {
			//g->edges->buffer[idx].closed = WT_EDGE_CLOSED_NULL;
		}
		push_edgerefv(g->erefs, (edge_ref_t){idx, 0, 0});
		push_edgerefv(g->erefs, (edge_ref_t){idx, 1, 0});
	}
	psort_array(g->erefs->buffer + 1, g->erefs->size - 1, edge_ref_t, ncpu, num_cmpgtx(
		(a.flg? ((g->edges->buffer[a.idx].node2 << 1) | !g->edges->buffer[a.idx].dir2) : ((g->edges->buffer[a.idx].node1 << 1) | g->edges->buffer[a.idx].dir1)),
		(b.flg? ((g->edges->buffer[b.idx].node2 << 1) | !g->edges->buffer[b.idx].dir2) : ((g->edges->buffer[b.idx].node1 << 1) | g->edges->buffer[b.idx].dir1)),
		g->edges->buffer[a.idx].off, g->edges->buffer[b.idx].off));
	f1.idx = g->nodes->size; f1.flg = 0;
	cnt = 0;
	for(lst=idx=1;idx<g->erefs->size;idx++){
		if(g->erefs->buffer[idx].flg){
			f2.idx =  g->edges->buffer[g->erefs->buffer[idx].idx].node2;
			f2.flg = !g->edges->buffer[g->erefs->buffer[idx].idx].dir2;
		} else {
			f2.idx =  g->edges->buffer[g->erefs->buffer[idx].idx].node1;
			f2.flg =  g->edges->buffer[g->erefs->buffer[idx].idx].dir1;
		}
		if(f1.idx == f2.idx && f1.flg == f2.flg) continue;
		if(lst < idx){
			n = ref_nodev(g->nodes, f1.idx);
			n->edges[f1.flg].idx = lst;
			n->edges[f1.flg].cnt = 0;
			for(x=lst;x+1<idx;x++){
				g->erefs->buffer[x].next = x + 1;
				if(g->edges->buffer[g->erefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) n->edges[f1.flg].cnt ++;
			}
			if(g->edges->buffer[g->erefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) n->edges[f1.flg].cnt ++;
		}
		lst = idx;
		f1 = f2;
	}
	if(lst < idx){
		n = ref_nodev(g->nodes, f1.idx);
		n->edges[f1.flg].idx = lst;
		n->edges[f1.flg].cnt = 0;
		for(x=lst;x+1<idx;x++){
			g->erefs->buffer[x].next = x + 1;
			if(g->edges->buffer[g->erefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) n->edges[f1.flg].cnt ++;
		}
		if(g->edges->buffer[g->erefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) n->edges[f1.flg].cnt ++;
	}
}

static inline void load_nodes_graph(Graph *g, FileReader *clp, FileReader *nds){
	read_t *rd;
	node_t *n;
	reg_t *reg, *r;
	uint64_t nid;
	uint32_t i, nreg, rid;
	char *str, *tok;
	int ncol, closed, nwarn;
	clear_nodev(g->nodes);
	clear_regv(g->regs);
	next_ref_regv(g->regs);
	nwarn = 0;
	while((ncol = readtable_filereader(clp)) != -1){
		if(clp->line->string[0] == '#') continue;
		if(ncol < 4) continue;
		if((rid = getval_cuhash(g->kbm->tag2idx, get_col_str(clp, 0))) == MAX_U4){
			if(nwarn < 10){
				fprintf(stderr, " -- Cannot find read \"%s\" in LINE:%llu in %s -- %s:%d --\n", get_col_str(clp, 0), clp->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				nwarn ++;
			}
			continue;
		}
		rd = ref_readv(g->reads, rid);
		rd->clps[0] = atoi(get_col_str(clp, 2)) / KBM_BIN_SIZE;
		rd->clps[1] = atoi(get_col_str(clp, 3)) / KBM_BIN_SIZE;
	}
	nwarn = 0;
	while((ncol = readtable_filereader(nds)) != -1){
		if(nds->line->string[0] == '#') continue;
		if(ncol < 2) continue;
		nreg = atoi(get_col_str(nds, 1));
		//if(nreg == 0) continue;
		//if(nreg < g->min_node_cov) continue;
		//node_id ~ N(\d+)(\*?)
		nid = atol(get_col_str(nds, 0) + 1);
		if(get_col_str(nds, 0)[get_col_len(nds, 0) - 1] == '*'){ // assert get_col_len(nds, 0) > 0
			closed = 1;
		} else {
			closed = 0;
		}
		while(g->nodes->size < nid){
			n = next_ref_nodev(g->nodes);
			ZEROS(n);
			n->closed = 1;
		}
		n = next_ref_nodev(g->nodes);
		ZEROS(n);
		n->closed = closed;
		n->regs.idx = g->regs->size;
		n->regs.cnt = nreg;
		for(i=0;i<nreg;i++){
			reg = next_ref_regv(g->regs);
			reg->node = nid;
			reg->closed = 0;
			reg->read_link = 0;
			str = get_col_str(nds, 2 + i);
			if(0){
				if(str[get_col_len(nds, 2 + i)-1] == '*'){
					reg->closed = 1;
				}
				tok = index(str, '_'); *tok = '\0';
				reg->rid = getval_cuhash(g->kbm->tag2idx, str);
				str = tok + 1; *tok = '_'; tok = index(str, '_'); *tok = '\0';
				reg->dir = (str[0] == 'R');
				str = tok + 1; *tok = '_'; tok = index(str, '_'); *tok = '\0';
				reg->beg = atoi(str);
				str = tok + 1;
				reg->end = atoi(str);
				reg->end += reg->beg;
			} else {
				if(str[get_col_len(nds, 2 + i)-1] == '*'){
					reg->closed = 1;
				}
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->end = atoi(tok);
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->beg = atoi(tok);
				reg->end += reg->beg;
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->dir = (tok[0] == 'R');
				if((reg->rid = getval_cuhash(g->kbm->tag2idx, str)) == MAX_U4){
					g->regs->size --;
					if(nwarn < 10){
						fprintf(stderr, " -- Cannot find read \"%s\" in LINE:%llu in %s -- %s:%d --\n", str, nds->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						nwarn ++;
					}
					continue;
				}
			}
			rd = ref_readv(g->reads, reg->rid);
			if(rd->regs.idx){
				r = ref_regv(g->regs, rd->regs.idx);
				if(r->beg > reg->beg){
					reg->read_link = rd->regs.idx;
					rd->regs.idx = g->regs->size - 1;
				} else {
					while(1){
						if(r->read_link == 0) break;
						if(g->regs->buffer[r->read_link].beg > reg->beg) break;
						r = ref_regv(g->regs, r->read_link);
					}
					reg->read_link = r->read_link;
					r->read_link = g->regs->size - 1;
				}
			} else {
				rd->regs.idx = g->regs->size - 1;
			}
			rd->regs.cnt ++;
		}
	}
	encap_regv(g->regs, 1);
	g->regs->buffer[g->regs->size].node = WT_MAX_NODE;
}

#endif
