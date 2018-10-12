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

#include "kbm.h"
#include "filewriter.h"
#include <getopt.h>
#include <regex.h>

#define WT_MAX_RD			0x03FFFFFF
#define WT_MAX_RDLEN		0x0003FFFF
#define WT_MAX_NODE			0x000000FFFFFFFFFFLLU
#define WT_MAX_EDGE			0x000000FFFFFFFFFFLLU
#define WT_MAX_NODE_EDGES	0xFFFF
#define WT_MAX_EDGE_COV		0x7FFF

typedef struct {
	u8i rid:26, dir:1, beg:18, end:18, closed:1;
} rd_frg_t;
define_list(rdfrgv, rd_frg_t);

typedef struct {
	u8i node;
	u8i rid:26, dir:1, beg:18, end:18, closed:1;
} rd_reg_t;
define_list(rdregv, rd_reg_t);

typedef struct {
	u8i rid:26, dir:1, beg:18, end:18, closed:1;
	u8i read_link;
} rd_rep_t;
define_list(rdrepv, rd_rep_t);

typedef struct {
	u8i node;
	u8i rid:26, dir:1, beg:18, end:18, closed:1;
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
	u8i idx:63, flg:1;
} hit_lnk_t;
define_list(hitlnkv, hit_lnk_t);

typedef struct {
	rd_frg_t  frgs[2];
	hit_lnk_t lnks[2]; // link to next hit of the read
} rd_hit_t;
define_list(rdhitv, rd_hit_t);

typedef struct {
	u8i visit:63, flag:1;
	hit_lnk_t hits; // point to the g->rdhits
	int clps[2];
	ptr_ref_t regs;
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
	KBMPar   *par;
	regv     *regs;
	readv    *reads;

	rdhitv   *rdhits;
	kbmmapv  *pwalns;
	BitsVec  *cigars;

	nodev    *nodes;
	edgev    *edges;
	edgehash *ehash;
	edgerefv *erefs;

	frgv     *frgs;
	lnkv     *lnks;
	edgerefv *lrefs;
	tracev   *traces;

	int      node_order;
	uint32_t n_fix, only_fix; // first n sequences are accurate contigs; only_fix means whether to include other pacbio sequenes
	uint32_t reglen, regovl, bestn;
	int      min_node_mats;
	int      max_overhang, chainning_hits;
	float    node_max_conflict; // 0.25
	float    node_merge_cutoff;
	uint32_t max_node_cov, min_node_cov, exp_node_cov, min_edge_cov;
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

Graph* init_graph(KBM *kbm){
	Graph *g;
	u4i rid;
	g = malloc(sizeof(Graph));
	g->kbm = kbm;
	g->par = kbm->par;
	g->regs = init_regv(32);
	g->reads = init_readv(kbm->reads->size);
	g->reads->size = kbm->reads->size;
	for(rid=0;rid<g->reads->size;rid++){
		g->reads->buffer[rid].clps[0] = 0;
		g->reads->buffer[rid].clps[1] = g->kbm->reads->buffer[rid].rdlen;
	}
	g->nodes = init_nodev(32);
	g->rdhits = init_rdhitv(1024);
	g->pwalns = init_kbmmapv(1024);
	g->cigars = init_bitsvec(1024, 3);
	g->edges = init_edgev(32);
	g->ehash = init_edgehash(1023);
	set_userdata_edgehash(g->ehash, g->edges);
	g->erefs = init_edgerefv(32);
	g->frgs = init_frgv(32);
	g->lnks = init_lnkv(32);
	g->lrefs = init_edgerefv(32);
	g->traces = init_tracev(32);
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
	g->min_node_cov = 4;
	g->max_node_cov = 60;
	g->exp_node_cov = 40; // expected node cov
	g->min_edge_cov = 4;
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

void free_graph(Graph *g){
	uint64_t i;
	free_regv(g->regs);
	free_readv(g->reads);
	free_nodev(g->nodes);
	free_rdhitv(g->rdhits);
	free_kbmmapv(g->pwalns);
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

int map2rdhits_graph(Graph *g, kbm_map_t *hit){
	rd_hit_t *rh;
	rh = next_ref_rdhitv(g->rdhits);
	rh->frgs[0] = (rd_frg_t){hit->qidx, hit->qdir, hit->qb, hit->qe, 0};
	rh->frgs[1] = (rd_frg_t){hit->tidx, hit->tdir, hit->tb, hit->te, 0};
	rh->lnks[0] = g->reads->buffer[hit->qidx].hits;
	g->reads->buffer[hit->qidx].hits = (hit_lnk_t){g->rdhits->size - 1, 0};
	rh->lnks[1] = g->reads->buffer[hit->tidx].hits;
	g->reads->buffer[hit->tidx].hits = (hit_lnk_t){g->rdhits->size - 1, 1};
	return 1;
}

int is_dovetail_overlap(Graph *g, kbm_map_t *hit){
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

int hit2rdregs_graph(Graph *g, rdregv *regs, kbm_map_t *hit, BitsVec *cigars, u4v *maps[3]){
	KBM *kbm;
	u8i ndoff;
	u4i bpos[2][2], npos[2][2], clen, ndbeg, qn, j;
	int tmp, bt, qlen, tlen, x, y, mat, beg, end, min_node_len, max_node_len;
	int mask, closed;
	kbm = g->kbm;
	mask = 0;
	if(g->max_overhang >= 0){
		if(!is_dovetail_overlap(g, hit)){
			mask = 1;
		}
	}
	qn = g->reglen / KBM_BIN_SIZE;
	min_node_len = (qn - 1) * KBM_BIN_SIZE;
	max_node_len = (qn + 1) * KBM_BIN_SIZE;
	// translate into reverse sequence order
	qlen = kbm->reads->buffer[hit->qidx].bincnt * KBM_BIN_SIZE;
	tlen = kbm->reads->buffer[hit->tidx].bincnt * KBM_BIN_SIZE;
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
			push_u4v(maps[0], x);
			push_u4v(maps[1], y);
			clen --;
		}
		push_u4v(maps[0], x + 1);
		push_u4v(maps[1], y + 1);
		push_u4v(maps[2], 0);
#if 0
		if(x + 1 + (hit->qb / KBM_BIN_SIZE) != (hit->qe / KBM_BIN_SIZE) || y + 1 + (hit->tb / KBM_BIN_SIZE) != (hit->te / KBM_BIN_SIZE)){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			print_hit_kbm(g->kbm, hit, cigars, stderr);
			abort();
		}
#endif
	}
	bpos[0][0] = hit->qb / KBM_BIN_SIZE;
	bpos[0][1] = hit->qe / KBM_BIN_SIZE;
	bpos[1][0] = hit->tb / KBM_BIN_SIZE;
	bpos[1][1] = hit->te / KBM_BIN_SIZE;
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
			ndbeg = kbm->reads->buffer[hit->qidx].bincnt - bpos[0][0];
			ndbeg = ndbeg % qn;
			ndoff = ndoff + ((kbm->reads->buffer[hit->qidx].bincnt - (ndbeg + bpos[0][0])) / qn) - 1;
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
				beg = (npos[0][0] + bpos[1][0]) * KBM_BIN_SIZE;
				end = (npos[1][1] + bpos[1][0] + 1) * KBM_BIN_SIZE;
				if(end - beg >= min_node_len && end - beg <= max_node_len){
					closed = 0;
				} else {
					closed = 1;
				}
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
					beg = qlen - (npos[1][1] + bpos[0][0] + 1) * KBM_BIN_SIZE;
					end = qlen - (npos[0][0] + bpos[0][0]) * KBM_BIN_SIZE;
				} else {
					beg = (npos[0][0] + bpos[0][0]) * KBM_BIN_SIZE;
					end = (npos[1][1] + bpos[0][0] + 1) * KBM_BIN_SIZE;
				}
				if(end - beg >= min_node_len && end - beg <= max_node_len){
					closed = 0;
				} else {
					closed = 1;
				}
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
KBMAux *aux;
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
KBMAux *aux;
rdregv *regs;
BitVec *rdflags;
u4v *maps[3];
kmeroffv *kmers[2];
kbm_map_t *hit;
u4i ridx, i, todo;
volatile reg_t *reg;
g = mdbg->g;
kbm = g->kbm;
reg = (reg_t*)&mdbg->reg;
aux = mdbg->aux;
regs = mdbg->regs;
rdflags = mdbg->rdflags;
maps[0] = init_u4v(32);
maps[1] = init_u4v(32);
maps[2] = init_u4v(32);
kmers[0] = adv_init_kmeroffv(64, 0, 1);
kmers[1] = adv_init_kmeroffv(64, 0, 1);
thread_beg_loop(mdbg);
if(mdbg->task == 1){
	if(reg->closed) continue;
	query_index_kbm(aux, NULL, reg->rid, kbm->rdseqs, kbm->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, kmers);
	map_kbm(aux);
	sort_array(aux->hits->buffer, aux->hits->size, kbm_map_t, num_cmpgt(b.mat, a.mat));
} else {
	for(ridx=mdbg->beg+mdbg->t_idx;ridx<mdbg->end;ridx+=mdbg->n_cpu){
		if(rdflags){
			todo = 1;
			thread_beg_syn(mdbg);
			todo = !get_bitvec(rdflags, ridx);
			thread_end_syn(mdbg);
			if(!todo) continue;
		}
		query_index_kbm(aux, NULL, ridx, kbm->rdseqs, kbm->reads->buffer[ridx].rdoff, kbm->reads->buffer[ridx].rdlen, kmers);
		map_kbm(aux);
		sort_array(aux->hits->buffer, aux->hits->size, kbm_map_t, num_cmpgt(b.mat, a.mat));
		{
			thread_beg_syn(mdbg);
			if(!KBM_LOG && ((ridx - mdbg->beg) % 1000) == 0){ fprintf(KBM_LOGF, "\r%u|%llu", ridx - mdbg->beg, (u8i)g->pwalns->size); fflush(KBM_LOGF); }
			u8i cgoff = g->cigars->size;
			if(mdbg->raw){
			} else {
				append_bitsvec(g->cigars, mdbg->aux->cigars, 0, mdbg->aux->cigars->size);
			}
			for(i=0;i<mdbg->aux->hits->size;i++){
				hit = ref_kbmmapv(mdbg->aux->hits, i);
				if(mdbg->alno){
					fprint_hit_kbm(mdbg->aux, i, mdbg->alno);
				}
				if(rdflags
					&& g->kbm->reads->buffer[hit->tidx].bincnt < g->kbm->reads->buffer[hit->qidx].bincnt
					&& (hit->tb <= KBM_BSIZE && hit->te + KBM_BSIZE >= (int)(g->kbm->reads->buffer[hit->tidx].bincnt * KBM_BSIZE))
					&& (hit->qb > KBM_BSIZE || hit->qe + KBM_BSIZE < (int)(g->kbm->reads->buffer[hit->qidx].bincnt * KBM_BSIZE))
					){
					one_bitvec(rdflags, hit->tidx);
				}
				if(mdbg->raw){
					hit2rdregs_graph(g, regs, hit, mdbg->aux->cigars, maps);
				} else {
					push_kbmmapv(g->pwalns, *hit);
					peer_kbmmapv(g->pwalns)->cgoff += cgoff;
				}
			}
			if(KBM_LOG){
				fprintf(KBM_LOGF, "QUERY: %s\t+\t%d\t%d\n", g->kbm->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end);
				for(i=0;i<mdbg->aux->hits->size;i++){
					hit = ref_kbmmapv(mdbg->aux->hits, i);
					fprintf(KBM_LOGF, "\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->kbm->reads->buffer[hit->tidx].tag, "+-"[hit->qdir], g->kbm->reads->buffer[hit->tidx].rdlen, hit->tb, hit->te, hit->aln, hit->mat);
				}
			}
			thread_end_syn(mdbg);
		}
	}
}
thread_end_loop(mdbg);
free_u4v(maps[0]);
free_u4v(maps[1]);
free_u4v(maps[2]);
free_kmeroffv(kmers[0]);
free_kmeroffv(kmers[1]);
thread_end_func(mdbg);

thread_beg_def(mbio);
int     bidx;
FILE  *bios[2], *out;
size_t  buf_size;
thread_end_def(mbio);

thread_beg_func(mbio);
char   *buffs[2];
size_t  blens[2], bsize[2];
int bidx;
mbio->once = 0;
mbio->bidx = 0;
buffs[0] = NULL;
buffs[1] = NULL;
blens[0] = 0;
blens[1] = 0;
mbio->bios[0] = open_memstream(buffs + 0, blens + 0);
mbio->bios[1] = open_memstream(buffs + 1, blens + 1);
mbio->bidx = 0;
thread_beg_loop(mbio);
bidx = mbio->bidx;
bsize[0] = ftell(mbio->bios[0]);
bsize[1] = ftell(mbio->bios[1]);
if(bsize[bidx] >= mbio->buf_size){
	thread_beg_syn(mbio);
}
if(bsize[!bidx]){
	fwrite(buffs[!bidx], 1, bsize[!bidx], mbio->out);
	fseek(mbio->bios[!bidx], 0, SEEK_SET);
}
if(bsize[bidx] >= mbio->buf_size){
	mbio->bidx = !bidx;
	thread_end_syn(mbio);
} else if(bsize[bidx]){
	thread_beg_syn(mbio);
	mbio->bidx = !bidx;
	thread_end_syn(mbio);
}
thread_end_loop(mbio);
{
	bsize[0] = ftell(mbio->bios[0]);
	bsize[1] = ftell(mbio->bios[1]);
	bidx = mbio->bidx;
	if(bsize[!bidx]){
		fwrite(buffs[!bidx], 1, bsize[!bidx], mbio->out);
	}
	if(bsize[bidx]){
		fwrite(buffs[bidx], 1, bsize[bidx], mbio->out);
	}
}
fclose(mbio->bios[0]);
fclose(mbio->bios[1]);
if(buffs[0]) free(buffs[0]);
if(buffs[1]) free(buffs[1]);
thread_end_func(mbio);

typedef struct {
	int pos:19;
	u4i dir:1, spur:1, dep:11;
} rd_clp_t;
define_list(rdclpv, rd_clp_t);

void clip_read_algo(int clps[2], rdclpv *brks, rdclpv *chis, int avg_dep, int min_dep){
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
	clps[0] = mx;
	clps[1] = my;
	if((int)chis->size < avg_dep){
		return;
	}
	sort_array(chis->buffer, chis->size, rd_clp_t, a.pos > b.pos);
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
	if(c->pos <= clps[0] || c->pos >= clps[1]){
		return;
	}
	if(c->pos - clps[0] > clps[1] - c->pos){
		clps[1] = c->pos;
	} else {
		clps[0] = c->pos;
	}
}

void clip_read_core(Graph *g, u4i rid, hitlnkv *lnks, rdclpv *brks, rdclpv *chis){
	read_t *rd;
	hit_lnk_t *lnk, *l1, *l2, *l3;
	rd_hit_t *hit, *hit1, *hit2, *hit3;
	rd_frg_t *f1, *f2;
	u4i i, j, h;
	int rlen, dlen, margin, min_dep, dep, avg, x, y;
	rd = ref_readv(g->reads, rid);
	rlen = g->kbm->reads->buffer[rid].bincnt * KBM_BIN_SIZE;
	if(rlen == 0){
		return;
	}
	lnk = &rd->hits;
	clear_rdclpv(brks);
	margin = g->max_overhang > 0? g->max_overhang : KBM_BIN_SIZE * 4;
	min_dep = 2;
	dep = 0;
	clear_hitlnkv(lnks);
	while(lnk->idx){
		hit = ref_rdhitv(g->rdhits, lnk->idx);
		f1 = hit->frgs + lnk->flg;
		if(f1->closed == 0) push_hitlnkv(lnks, *lnk);
		lnk = hit->lnks + lnk->flg;
	}
	if(0 && g->chainning_hits){
		sort_array(
			lnks->buffer, lnks->size, hit_lnk_t,
			num_cmpgtxx(
				g->rdhits->buffer[a.idx].frgs[!a.flg].rid, g->rdhits->buffer[b.idx].frgs[!b.flg].rid,
				g->rdhits->buffer[a.idx].frgs[0].dir, g->rdhits->buffer[b.idx].frgs[0].dir,
				g->rdhits->buffer[a.idx].frgs[!a.flg].beg, g->rdhits->buffer[b.idx].frgs[!b.flg].beg
			)
		);
		push_hitlnkv(lnks, (hit_lnk_t){0, 0});
		l1 = ref_hitlnkv(lnks, 0);
		hit1 = ref_rdhitv(g->rdhits, l1->idx);
		for(i=1,j=0;i<lnks->size;i++){
			l2 = ref_hitlnkv(lnks, i);
			hit2 = ref_rdhitv(g->rdhits, l2->idx);
			if(hit1->frgs[!l1->flg].rid != hit2->frgs[!l2->flg].rid || hit1->frgs[0].dir != hit2->frgs[0].dir){
				for(h=j+1;h<i;h++){
					l3 = ref_hitlnkv(lnks, h);
					hit3 = ref_rdhitv(g->rdhits, l3->idx);
					hit3->frgs[0].closed = 1;
					hit3->frgs[1].closed = 1;
					if(hit1->frgs[0].beg > hit3->frgs[l1->flg ^ l3->flg].beg) hit1->frgs[0].beg = hit3->frgs[l1->flg ^ l3->flg].beg;
					if(hit1->frgs[0].end < hit3->frgs[l1->flg ^ l3->flg].end) hit1->frgs[0].end = hit3->frgs[l1->flg ^ l3->flg].end;
					if(hit1->frgs[1].beg > hit3->frgs[(l1->flg == l3->flg)].beg) hit1->frgs[0].beg = hit3->frgs[(l1->flg == l3->flg)].beg;
					if(hit1->frgs[1].end < hit3->frgs[(l1->flg == l3->flg)].end) hit1->frgs[0].end = hit3->frgs[(l1->flg == l3->flg)].end;
				}
				j = i;
				l1 = l2;
				hit1 = hit2;
			}
		}
		lnks->size --;
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
			dlen = g->kbm->reads->buffer[f2->rid].bincnt * KBM_BIN_SIZE;
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
				push_rdclpv(brks, (rd_clp_t){f1->beg, 0, 1, 0});
				push_rdclpv(brks, (rd_clp_t){f1->end, 1, 0, 0});
			} else if(y){
				push_rdclpv(brks, (rd_clp_t){f1->beg, 0, 0, 0});
				push_rdclpv(brks, (rd_clp_t){f1->end, 1, 1, 0});
			} else {
				dep += f1->end - f1->beg;
				push_rdclpv(brks, (rd_clp_t){f1->beg, 0, 0, 0});
				push_rdclpv(brks, (rd_clp_t){f1->end, 1, 0, 0});
			}
		}
	}
	clear_rdclpv(chis);
	avg = (dep + rlen - 1) / rlen;
	clip_read_algo(rd->clps, brks, chis, avg, min_dep);
}

thread_beg_def(mclp);
Graph *g;
thread_end_def(mclp);

thread_beg_func(mclp);
rdclpv *brks, *chis;
hitlnkv *lnks;
u4i rid;
brks = init_rdclpv(32);
chis = init_rdclpv(32);
lnks = init_hitlnkv(32);
thread_beg_loop(mclp);
for(rid=mclp->t_idx;rid<mclp->g->reads->size;rid+=mclp->n_cpu){
	clip_read_core(mclp->g, rid, lnks, brks, chis);
}
thread_end_loop(mclp);
free_rdclpv(brks);
free_rdclpv(chis);
free_hitlnkv(lnks);
thread_end_func(mclp);

u8i check_read_reg_conflict_core(Graph *g, rd_reg_t *hit, int *conflict){
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
void mul_update_regs_graph(Graph *g, rdregv *regs, rnkrefv *nds, u4i ncpu, u8i upds[3]){
	u8i idx, vt;
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
	encap_regv(g->regs, regs->size + 1); // make sure next_ref_regv in this function is thread-safe
	thread_beg_init(mupd, ncpu);
	mupd->g = g;
	mupd->regs = regs;
	mupd->nd = NULL;
	mupd->ins = init_u8v(64);
	mupd->vt = 0;
	thread_end_init(mupd);
	vt = 1;
	upds[0] = upds[1] = upds[2] = 0;
	for(idx=0;idx<nds->size+ncpu;idx++){
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

void build_nodes_graph(Graph *g, u8i maxbp, int ncpu, FileReader *pws, int rdclip, char *prefix, char *dump_kbm){
	kbm_map_t *hit, HIT;
	BitVec *rks, *rdflags;
	u4v *maps[3];
	BufferedWriter *bw;
	FILE *alno, *clplog, *kmlog;
	cuhash_t *cu;
	u8i idx, rank, kcnts[256], upds[3], nbp, fix_node;
	u4i nqry, rid, i, dep, ib, ie, mi, max_node_cov, round;
	int n_cpu, ncol, reset_kbm, raw;
	kbm_read_t *pb;
	rdregv *regs;
	rnkrefv *nds;
	rnk_ref_t *nd;
	thread_preprocess(mdbg);
	thread_preprocess(mclp);
	if(KBM_LOG) n_cpu = 1;
	else n_cpu = ncpu;
	regs = init_rdregv(1024);
	nds = init_rnkrefv(1024);
	free_rdhitv(g->rdhits); g->rdhits = init_rdhitv(1024);
	maps[0] = init_u4v(4);
	maps[1] = init_u4v(4);
	maps[2] = init_u4v(4);
	clear_regv(g->regs);
	next_ref_regv(g->regs);
	clear_rdhitv(g->rdhits);
	ZEROS(next_ref_rdhitv(g->rdhits));
	clear_kbmmapv(g->pwalns);
	clear_bitsvec(g->cigars);
	raw = !(g->chainning_hits || (g->bestn > 0) || rdclip);
	fix_node = 0; // TODO: new code hasn't coped with contigs+longreads mode
	if(pws){
		u4v *cgs;
		u8i nhit;
		int qlen, val, flg, nwarn, mwarn;
		char *cgstr, *qtag;
		mwarn = 20;
		nwarn = 0;
		cgs = init_u4v(4);
		memset(&HIT, 0, sizeof(kbm_map_t));
		hit = &HIT;
		nhit = 0;
		while((ncol = readtable_filereader(pws)) != -1){
			if((pws->n_line % 100000) == 0){
				fprintf(KBM_LOGF, "\r%llu", pws->n_line); fflush(KBM_LOGF);
			}
			if(pws->line->buffer[0] == '#') continue;
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
			if(qlen != (int)g->kbm->reads->buffer[hit->qidx].rdlen){
				if(nwarn < mwarn){
					fprintf(stderr, " -- inconsisitent read length \"%s\" %d != %d in %s -- %s:%d --\n", qtag, qlen, g->kbm->reads->buffer[hit->qidx].rdlen, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					nwarn ++;
				}
				continue;
			}
			hit->qb = atoi(get_col_str(pws, 3));
			hit->qe = atoi(get_col_str(pws, 4));
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
			hit->tb = atoi(get_col_str(pws, 8));
			hit->te = atoi(get_col_str(pws, 9));
			// skip col 10-13
			hit->mat = atoi(get_col_str(pws, 10));
			if(hit->mat < g->par->min_mat) continue;
			hit->aln = atoi(get_col_str(pws, 11));
			if(hit->aln < g->par->min_aln) continue;
			hit->cnt = atoi(get_col_str(pws, 12));
			hit->gap = atoi(get_col_str(pws, 13));
			hit->cgoff = g->cigars->size;
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
				while(val > 0){
					push_bitsvec(g->cigars, flg);
					val --;
				}
			}
			hit->cglen = g->cigars->size - hit->cgoff;
			nhit ++;
			if(raw){
				hit2rdregs_graph(g, regs, hit, g->cigars, maps);
				clear_bitsvec(g->cigars);
			} else {
				push_kbmmapv(g->pwalns, *hit);
			}
		}
		fprintf(KBM_LOGF, "\r%llu lines, %llu hits\n", pws->n_line, nhit);
		free_u4v(cgs);
	} else {
		nqry = g->only_fix? g->n_fix : g->reads->size;
		ib = g->only_fix? g->kbm->reads->buffer[g->n_fix].binoff : 0;
		mi = g->kbm->bins->size;
		ie = 0;
		if(maxbp < g->kbm->rdseqs->size){
			max_node_cov = 1.0 * g->max_node_cov * maxbp / g->kbm->rdseqs->size;
			if(max_node_cov < 10) max_node_cov = 10;
		} else max_node_cov = g->max_node_cov;
		round = 0;
		alno   = open_file_for_write(prefix, ".alignments", 1);
		bw = open_bufferedwriter(alno, 1024 * 1024);
		rdflags = g->par->skip_contained? init_bitvec(g->kbm->reads->size) : NULL;
		thread_beg_init(mdbg, n_cpu);
		mdbg->g = g;
		memset((void*)&mdbg->reg, 0, sizeof(reg_t));
		mdbg->reg.closed = 1;
		mdbg->aux = init_kbmaux(g->kbm);
		mdbg->aux->par = g->par;
		mdbg->regs = regs;
		mdbg->rdflags = rdflags;
		mdbg->beg = 0;
		mdbg->end = 0;
		mdbg->raw = raw;
		mdbg->alno = alno;
		thread_end_init(mdbg);
		reset_kbm = 0;
		while(ib < mi){
			nbp = 0;
			round ++;
			ie = num_min(ib + maxbp / KBM_BIN_SIZE, mi);
			nbp = ((u8i)(ie - ib)) * KBM_BSIZE;
			if(round == 1 && g->kbm->seeds->size){
				// already indexed, maybe loaded from kbm-index file
				reset_kbm = 0;
			} else {
				reset_index_kbm(g->kbm);
				fprintf(KBM_LOGF, "[%s] indexing bins[%u,%u] (%llu bp), %d threads\n", date(), ib, ie, nbp, ncpu); fflush(KBM_LOGF);
				index_kbm(g->kbm, ib, ie, ncpu);
				reset_kbm = 1;
				fprintf(KBM_LOGF, "[%s] Done\n", date()); fflush(KBM_LOGF);
				if(round == 1 && dump_kbm){
					FILE *dump;
					fprintf(KBM_LOGF, "[%s] dump kbm index to %s ...", date(), dump_kbm); fflush(KBM_LOGF);
					dump = open_file_for_write(dump_kbm, NULL, 1);
					mem_dump_obj_file(g->kbm, 1, &kbm_obj_desc, 1, 0, dump);
					fclose(dump);
					fprintf(KBM_LOGF, " Done\n"); fflush(KBM_LOGF);
				}
			}
			ib = ie;
			if(round == 1) fix_node = 0;
			if(round == 1){
				kbm_kmer_t *u;
				size_t iter_ptr;
				u4i *deps, hidx;
				deps = calloc(KBM_MAX_KCNT + 1, 4);
				kmlog = open_file_for_write(prefix, ".kmerdep", 1);
				for(hidx=0;hidx<KBM_N_HASH;hidx++){
					iter_ptr = 0;
					while((u = ref_iter2_kbmhash(g->kbm->hashs[hidx], &iter_ptr))){
						deps[u->tot] ++;
					}
				}
				for(hidx=0;hidx<=KBM_MAX_KCNT;hidx++){
					fprintf(kmlog, "%u\t%u\n", hidx, deps[hidx]);
				}
				fclose(kmlog);
				free(deps);
				kmlog = open_file_for_write(prefix, ".binkmer", 1);
				deps = calloc(KBM_BIN_SIZE + 1, 4);
				for(hidx=0;hidx<g->kbm->bins->size;hidx++){
					deps[g->kbm->bins->buffer[hidx].degree] ++;
				}
				for(hidx=0;hidx<256;hidx++){
					fprintf(kmlog, "%u\n", deps[hidx]);
				}
				fclose(kmlog);
				free(deps);
				if(!g->minimal_output){
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
			if(0){
				thread_beg_iter(mdbg);
				mdbg->beg = 0;
				mdbg->end = nqry;
				mdbg->task = 2;
				thread_wake(mdbg);
				thread_end_iter(mdbg);
			} else {
				thread_beg_iter(mdbg);
				mdbg->task = 1;
				thread_end_iter(mdbg);
				for(rid=0;rid<=nqry+ncpu;rid++){
					if(rid < nqry){
						if(!KBM_LOG && (rid % 2000) == 0){ fprintf(KBM_LOGF, "\r%u|%llu", rid, (u8i)g->pwalns->size); fflush(KBM_LOGF); }
						thread_wait_one(mdbg);
					} else {
						thread_wait_next(mdbg);
						pb = NULL;
					}
					if(mdbg->reg.closed == 0){
						u8i cgoff = g->cigars->size;
						if(raw){
						} else {
							if(1){
								append_bitsvec(g->cigars, mdbg->aux->cigars, 0, mdbg->aux->cigars->size);
							} else {
								for(i=0;i<mdbg->aux->cigars->size;i++){
									push_bitsvec(g->cigars, get_bitsvec(mdbg->aux->cigars, i));
								}
							}
						}
						if(alno){
							beg_bufferedwriter(bw);
							for(i=0;i<mdbg->aux->hits->size;i++){
								hit = ref_kbmmapv(mdbg->aux->hits, i);
								fprint_hit_kbm(mdbg->aux, i, bw->out);
							}
							end_bufferedwriter(bw);
						}
						for(i=0;i<mdbg->aux->hits->size;i++){
							hit = ref_kbmmapv(mdbg->aux->hits, i);
							if(rdflags
								&& g->kbm->reads->buffer[hit->tidx].bincnt < g->kbm->reads->buffer[hit->qidx].bincnt
								&& (hit->tb <= KBM_BSIZE && hit->te + KBM_BSIZE >= (int)(g->kbm->reads->buffer[hit->tidx].bincnt * KBM_BSIZE))
								&& (hit->qb > KBM_BSIZE || hit->qe + KBM_BSIZE < (int)(g->kbm->reads->buffer[hit->qidx].bincnt * KBM_BSIZE))
								){
								one_bitvec(rdflags, hit->tidx);
							}
							if(raw){
								hit2rdregs_graph(g, regs, hit, mdbg->aux->cigars, maps);
							} else {
								push_kbmmapv(g->pwalns, *hit);
								peer_kbmmapv(g->pwalns)->cgoff += cgoff;
							}
						}
						if(KBM_LOG){
							fprintf(KBM_LOGF, "QUERY: %s\t+\t%d\t%d\n", g->kbm->reads->buffer[mdbg->reg.rid].tag, mdbg->reg.beg, mdbg->reg.end);
							for(i=0;i<mdbg->aux->hits->size;i++){
								hit = ref_kbmmapv(mdbg->aux->hits, i);
								fprintf(KBM_LOGF, "\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->kbm->reads->buffer[hit->tidx].tag, "+-"[hit->qdir], g->kbm->reads->buffer[hit->tidx].rdlen, hit->tb, hit->te, hit->aln, hit->mat);
							}
						}
						mdbg->reg.closed = 1;
					}
					if(rid < nqry && (rdflags == NULL || get_bitvec(rdflags, rid) == 0)){
						pb = ref_kbmreadv(g->kbm->reads, rid);
						mdbg->reg = (reg_t){0, rid, 0, 0, pb->rdlen, 0, 0};
						thread_wake(mdbg);
					}
				}
			}
			if(ib < mi && !KBM_LOG) fprintf(KBM_LOGF, "\r");
		}
		thread_beg_close(mdbg);
		free_kbmaux(mdbg->aux);
		thread_end_close(mdbg);
		if(!KBM_LOG) fprintf(KBM_LOGF, "\r%u reads|total hits %llu\n", nqry, (u8i)g->pwalns->size);
		if(bw) close_bufferedwriter(bw);
		if(alno) fclose(alno);
		if(rdflags) free_bitvec(rdflags);
		if(reset_kbm){
			reset_index_kbm(g->kbm);
		}
	}
	print_proc_stat_info(0);
	if(raw){
	} else {
		if(g->chainning_hits){
			u8i lst, nch, nmg;
			fprintf(KBM_LOGF, "[%s] chainning ... ", date()); fflush(KBM_LOGF);
			psort_array(g->pwalns->buffer, g->pwalns->size, kbm_map_t, ncpu, num_cmpgtxx(a.qidx, b.qidx, a.tidx, b.tidx, a.qdir, b.qdir));
			nch = nmg = 0;
			for(idx=lst=0;idx<=g->pwalns->size;idx++){
				if(idx == g->pwalns->size || g->pwalns->buffer[lst].qidx != g->pwalns->buffer[idx].qidx ||
					g->pwalns->buffer[lst].tidx != g->pwalns->buffer[idx].tidx || g->pwalns->buffer[lst].qdir != g->pwalns->buffer[idx].qdir){
					if(idx > lst + 1){
						if(0){ // TODO: remove this block after debug
							u8i i;
							fprintf(KBM_LOGF, "++\n");
							for(i=lst;i<idx;i++){
								print_hit_kbm(g->kbm, g->pwalns->buffer + i, g->cigars, KBM_LOGF);
							}
							fprintf(KBM_LOGF, "--\n");
						}
						if(simple_chain_all_maps_kbm(g->pwalns->buffer + lst, idx - lst, g->cigars, &HIT, g->cigars, g->kbm->par->aln_var)){
							nch += idx - lst;
							nmg ++;
							g->pwalns->buffer[lst++] = HIT;
							if(0){
								fprintf(KBM_LOGF, "MERGED\t");
								print_hit_kbm(g->kbm, &HIT, g->cigars, KBM_LOGF);
							}
							while(lst < idx){
								g->pwalns->buffer[lst++].mat = 0; // closed = 1
							}
						}
					}
					lst = idx;
				}
			}
			fprintf(KBM_LOGF, " %llu hits into %llu\n", nch, nmg); fflush(KBM_LOGF);
		}
		if(g->bestn > 0){
			fprintf(KBM_LOGF, "[%s] picking best %d hits for each read ... ", date(), g->bestn); fflush(KBM_LOGF);
			psort_array(g->pwalns->buffer, g->pwalns->size, kbm_map_t, ncpu, num_cmpgtx(b.aln, a.aln, b.mat, a.mat));
			u4i *cnts;
			u8i size;
			cnts = calloc(g->reads->size, sizeof(u4i)); // all zeros
			size = 0;
			for(idx=0;idx<g->pwalns->size;idx++){
				hit = ref_kbmmapv(g->pwalns, idx);
				if(hit->mat == 0) continue;
				//if(cnts[hit->qidx] < g->bestn || cnts[hit->tidx] < g->bestn){
				if(cnts[hit->qidx] < g->bestn && cnts[hit->tidx] < g->bestn){
					size ++;
					cnts[hit->qidx] ++;
					cnts[hit->tidx] ++;
				} else {
					hit->mat = 0;
				}
			}
			free(cnts);
			fprintf(KBM_LOGF, "%llu hits\n", size);
		}
		// clip reads
		if(rdclip){
			fprintf(KBM_LOGF, "[%s] clipping ... ", date()); fflush(KBM_LOGF);
			for(idx=0;idx<g->pwalns->size;idx++){
				hit = ref_kbmmapv(g->pwalns, idx);
				if(hit->mat == 0) continue;
				map2rdhits_graph(g, hit);
			}
			clplog = open_file_for_write(prefix, ".clps", 1);
			thread_beg_init(mclp, ncpu);
			mclp->g = g;
			thread_end_init(mclp);
			thread_wake_all(mclp);
			thread_beg_close(mclp);
			thread_end_close(mclp);
			u8i tot, clp;
			tot = clp = 0;
			for(rid=0;rid<g->reads->size;rid++){
				tot += g->kbm->reads->buffer[rid].rdlen;
				clp += g->reads->buffer[rid].clps[1] - g->reads->buffer[rid].clps[0];
				fprintf(clplog, "%s\t%d\t%d\t%d\n", g->kbm->reads->buffer[rid].tag, g->kbm->reads->buffer[rid].rdlen, g->reads->buffer[rid].clps[0], g->reads->buffer[rid].clps[1]);
			}
			fclose(clplog);
			fprintf(KBM_LOGF, "%.2f%% bases\n", ((tot - clp) * 100.0) / tot); fflush(KBM_LOGF);
		}
		if(0){
			for(idx=0;idx<g->pwalns->size;idx++){
				hit = ref_kbmmapv(g->pwalns, idx);
				if(hit->mat == 0) continue;
				hit2rdregs_graph(g, regs, hit, g->cigars, maps);
			}
		} else {
			thread_fast_run(mhit, ncpu, EXPR(
				u8i i;
				rdregv *rs;
				u4v *ms[3];
				rs = init_rdregv(1024);
				ms[0] = init_u4v(1024);
				ms[1] = init_u4v(1024);
				ms[2] = init_u4v(1024);
				for(i=TIDX;i<g->pwalns->size;i+=NCPU){
					clear_rdregv(rs);
					hit2rdregs_graph(g, rs, ref_kbmmapv(g->pwalns, i), g->cigars, ms);
					if(rs->size){
						thread_beg_syn(mhit);
						append_rdregv(regs, rs);
						thread_end_syn(mhit);
					}
				}
				free_rdregv(rs);
				free_u4v(ms[0]);
				free_u4v(ms[1]);
				free_u4v(ms[2]);
			));
		}
		free_kbmmapv(g->pwalns); g->pwalns = init_kbmmapv(1024);
		free_bitsvec(g->cigars); g->cigars = init_bitsvec(1024, 3);
	}
	free_u4v(maps[0]);
	free_u4v(maps[1]);
	free_u4v(maps[2]);
	fprintf(KBM_LOGF, "[%s] generated %llu regs\n", date(), (u8i)regs->size);
	// add node itself
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
	// generating nodes
	fprintf(KBM_LOGF, "[%s] sorting regs ... ", date()); fflush(KBM_LOGF);
	psort_array(regs->buffer, regs->size, rd_reg_t, ncpu, num_cmpgtxx((a.node << 26) | a.rid, (b.node << 26) | b.rid, a.beg, b.beg, b.end, a.end));
	fprintf(KBM_LOGF, " Done\n"); fflush(KBM_LOGF);
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
	psort_array(nds->buffer, nds->size, rnk_ref_t, ncpu, num_cmpgtx(b.cnt, a.cnt, a.rank, b.rank));
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

void remove_all_edges_graph(Graph *g){
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

int estimate_edge_length(edge_off_t *ps, uint32_t size, uint32_t idxs[2]){
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
		for(i=b,tot=0;i<e;i++) tot += ps[i].off; avg = tot / (e - b);
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

int calculate_edge_cov_off_graph(Graph *g, edge_t *e, edgeoffv *offs){
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

void build_edges_graph(Graph *g, int ncpu, FILE *log){
	read_t *rd;
	node_t *n;
	edge_t *E;
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
	//E->cov = 1;
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
			for(i=0;i<regs->size;i++){
				r1 = (reg_t*)get_vplist(regs, i);
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
					if(g->edges->buffer[*u].cov < WT_MAX_EDGE_COV) g->edges->buffer[*u].cov ++;
				} else {
					*u = g->edges->size;
					push_edgev(g->edges, *E);
				}
				if(i + 1 == regs->size){
					g->edges->buffer[*u].closed = 0;
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

void load_nodes_graph(Graph *g, FileReader *clp, FileReader *nds){
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
		rd->clps[0] = atoi(get_col_str(clp, 2));
		rd->clps[1] = atoi(get_col_str(clp, 3));
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

void print_node_edges_cov_graph(Graph *g, FILE *out){
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	u8i idx, nid;
	u4i k, i;
	u4v *covs;
	covs = init_u4v(32);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		clear_u4v(covs);
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				e = ref_edgev(g->edges, f->idx);
				push_u4v(covs, e->cov);
				idx = f->next;
			}
		}
		if(covs->size == 0) continue;
		sort_array(covs->buffer, covs->size, u4i, num_cmpgt(b, a));
		fprintf(out, "NODE_COV\tN%llu\t%u\t%u", nid, n->cov, n->regs.cnt);
		for(i=0;i<covs->size;i++){
			fprintf(out, "\t%u", covs->buffer[i]);
		}
		fprintf(out, "\n");
	}
	free_u4v(covs);
}

// MUST be called before build_edges
uint64_t mask_nodes_by_cov_graph(Graph *g, FILE *out){
	node_t *n;
	uint64_t ret, i;
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

static inline void cut_edge_core_graph(Graph *g, edge_t *e, int closed_val){
	//if(e->closed == closed_val) return;
	if(e->closed) return;
	e->closed = closed_val;
	ref_nodev(g->nodes, e->node1)->edges[e->dir1].cnt --;
	ref_nodev(g->nodes, e->node2)->edges[!e->dir2].cnt --;
}

#define cut_edge_graph(g, e) cut_edge_core_graph(g, e, 1)

static inline void cut_lnk_core_graph(Graph *g, lnk_t *e, int closed_val){
	if(e->closed) return;
	e->closed = closed_val;
	ref_frgv(g->frgs, e->frg1)->lnks[e->dir1].cnt --;
	ref_frgv(g->frgs, e->frg2)->lnks[!e->dir2].cnt --;
}

#define cut_lnk_graph(g, e) cut_lnk_core_graph(g, e, 1)

static inline void revive_edge_graph(Graph *g, edge_t *e){
	if(e->closed == WT_EDGE_CLOSED_NULL) return;
	e->closed = WT_EDGE_CLOSED_NULL;
	ref_nodev(g->nodes, e->node1)->edges[e->dir1].cnt ++;
	ref_nodev(g->nodes, e->node2)->edges[!e->dir2].cnt ++;
}

static inline void revive_lnk_graph(Graph *g, lnk_t *e){
	if(e->closed == WT_EDGE_CLOSED_NULL) return;
	e->closed = WT_EDGE_CLOSED_NULL;
	ref_frgv(g->frgs, e->frg1)->lnks[e->dir1].cnt ++;
	ref_frgv(g->frgs, e->frg2)->lnks[!e->dir2].cnt ++;
}

void del_node_edges_graph(Graph *g, node_t *n){
	edge_ref_t *f;
	edge_t *e;
	uint64_t idx;
	uint32_t k;
	for(k=0;k<2;k++){
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = g->edges->buffer + f->idx;
			cut_edge_core_graph(g, e, WT_EDGE_CLOSED_HARD);
		}
	}
}

void del_frg_lnks_graph(Graph *g, frg_t *n){
	edge_ref_t *f;
	lnk_t *e;
	uint64_t idx;
	uint32_t k;
	for(k=0;k<2;k++){
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = g->lnks->buffer + f->idx;
			cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
		}
	}
}

void del_node_graph(Graph *g, node_t *n){
	del_node_edges_graph(g, n);
	n->closed = 1;
}

u8i mask_nodes_by_edge_cov_graph(Graph *g, u4i min_node_cov, float min_edge_cov_ratio, FILE *out){
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	u8i idx, nid, ret;
	u4i max, k;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		//if(n->regs.cnt < min_node_cov) continue;
		if(n->cov < min_node_cov) continue;
		max = 0;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				e = ref_edgev(g->edges, f->idx);
				if(e->cov > max) max = e->cov;
				idx = f->next;
			}
		}
		if(max < (u4i)(n->regs.cnt * min_edge_cov_ratio)){
			del_node_graph(g, n);
			if(out){
				fprintf(out, "NODE_EDGE_COV\tN%llu\t%u\t%u\n", nid, n->regs.cnt, max);
			}
			ret ++;
		}
	}
	return ret;
}

#define MAX_BT_NIDX	0xFFFFFFU
typedef struct {
	u8i node:46, visit:1, closed:1, cov:11, fix:1, sub_dir:1;
	u8i flag;
	u8i bt_dir:1, bt_open:10, bt_nidx:24, bt_score:20, bt_step:8, bt_hit:1;
	ptr_ref_t edges[2];
} subnode_t;
#define subnode_hashcode(E) u64hashcode((E).node)
#define subnode_hashequals(E1, E2) (E1).node == (E2).node
define_hashset(subnodehash, subnode_t, subnode_hashcode, subnode_hashequals);

typedef struct {
	subnode_t *node;
	uint32_t cov:28, visit:1, fwd:1, dir:1, closed:1;
	uint32_t next;
} subedge_t;
define_list(subedgev, subedge_t);

int evaluate_node_connectivity_graph(Graph *g, uint64_t nid, u4v *rds, subnodehash *nodes, subedgev *edges, ptrrefv *stack){
	node_t *nd;
	read_t *rd;
	reg_t  *rg;
	subnode_t N, *n, *n1, *n2;
	subedge_t *e;
	ptr_ref_t *p;
	u8i idx, edx, aim;
	uint32_t i, k, k1, k2, cnt;
	int exists;
	// collect reads containing nid
	clear_u4v(rds);
	nd = ref_nodev(g->nodes, nid);
	for(i=0;i<nd->regs.cnt;i++){
		rg = ref_regv(g->regs, nd->regs.idx + i);
		push_u4v(rds, rg->rid);
	}
	// prepare nodes in subgraph
	clear_subnodehash(nodes);
	clear_subedgev(edges);
	next_ref_subedgev(edges);
	memset(&N, 0, sizeof(subnode_t));
	N.cov = 1;
	for(i=0;i<rds->size;i++){
		rd = ref_readv(g->reads, rds->buffer[i]);
		rg = NULL;
		idx = rd->regs.idx;
		while(idx){
			rg = ref_regv(g->regs, idx);
			idx = rg->read_link;
			if(rg->closed) continue;
			N.node = rg->node;
			n = prepare_subnodehash(nodes, N, &exists);
			if(exists){
				n->cov ++;
			} else {
				*n = N;
			}
		}
	}
	// mask low cov nodes
	reset_iter_subnodehash(nodes);
	while((n = ref_iter_subnodehash(nodes))){
		if(n->cov < g->min_node_cov) n->closed = 1;
	}
	// build edges
	for(i=0;i<rds->size;i++){
		rd = ref_readv(g->reads, rds->buffer[i]);
		n1 = NULL;
		k1 = 0;
		idx = rd->regs.idx;
		while(idx){
			rg = ref_regv(g->regs, idx);
			idx = rg->read_link;
			if(rg->closed) continue;
			N.node = rg->node;
			n2 = get_subnodehash(nodes, N);
			k2 = rg->dir;
			if(n2->closed) continue;
			if(n1){
				// link n1 to n2
				edx = n1->edges[k1].idx;
				while(edx){
					e = ref_subedgev(edges, edx);
					if(e->node == n2 && e->dir == k2){
						e->cov ++;
						break;
					}
					edx = e->next;
				}
				if(edx == 0){
					edx = edges->size;
					e = next_ref_subedgev(edges);
					e->node = n2;
					e->dir = k2;
					e->cov = 1;
					e->next = n1->edges[k1].idx;
					n1->edges[k1].idx = edx;
					n1->edges[k1].cnt ++;
				}
				// link rev n2 to rev n1
				edx = n2->edges[!k2].idx;
				while(edx){
					e = ref_subedgev(edges, edx);
					if(e->node == n1 && e->dir == !k1){
						e->cov ++;
						break;
					}
					edx = e->next;
				}
				if(edx == 0){
					edx = edges->size;
					e = next_ref_subedgev(edges);
					e->node = n1;
					e->dir = !k1;
					e->cov = 1;
					e->next = n2->edges[!k2].idx;
					n2->edges[!k2].idx = edx;
					n2->edges[!k2].cnt ++;
				}
			}
			n1 = n2;
			k1 = k2;
		}
	}
	// find the nid node
	N.node = nid;
	n = get_subnodehash(nodes, N);
	n->visit = 1;
	// checking whether its out-edges collapse into one node
	for(k=0;k<2;k++){
		if(n->edges[k].cnt > 64) return 0;
		if(n->edges[k].cnt < 2) continue;
		idx = n->edges[k].idx;
		cnt = 0;
		while(idx){
			e = ref_subedgev(edges, idx);
			idx = e->next;
			if(e->cov == 1) continue; // don't track low cov out-edges
			cnt ++;
		}
		aim = 0xFFFFFFFFFFFFFFFFLLU >> (64 - cnt);
		cnt = 0;
		exists = 0;
		if(k){
			reset_iter_subnodehash(nodes);
			while((n1 = ref_iter_subnodehash(nodes))){ n1->flag = 0; }
		}
		idx = n->edges[k].idx;
		while(idx){
			e = ref_subedgev(edges, idx);
			idx = e->next;
			if(e->cov == 1) continue; // don't track low cov out-edges
			e->node->flag |= 1LLU << cnt;
			cnt ++;
			reset_iter_subnodehash(nodes);
			while((n1 = ref_iter_subnodehash(nodes))){ n1->visit = 0; }
			n->visit = 1;
			clear_ptrrefv(stack);
			push_ptrrefv(stack, (ptr_ref_t){offset_subnodehash(nodes, e->node), e->dir});
			while(stack->size){
				p = peer_ptrrefv(stack);
				n1 = nodes->array + p->idx;
				k1 = p->cnt;
				stack->size --;
				if(n1->flag == aim){ exists = 1; break; }
				if(n1->visit) continue;
				n1->visit = 1;
				edx = n1->edges[k1].idx;
				while(edx){
					e = ref_subedgev(edges, edx);
					edx = e->next;
					if(e->node->visit) continue;
					e->node->flag |= n1->flag;
					push_ptrrefv(stack, (ptr_ref_t){offset_subnodehash(nodes, e->node), e->dir});
				}
				if(exists) break;
			}
			if(exists) break;
		}
		if(exists == 0) return 0;
	}
	return 1;
}

void print_subgraph_dot(Graph *g, u8i id, subnodehash *nodes, subedgev *edges, FILE *out){
	subnode_t *n;
	subedge_t *e;
	u8i idx;
	int k;
	fprintf(out, "digraph N%llu {\n", id);
	fprintf(out, " N%llu [style=filled fillcolor=yellow]\n", id);
	reset_iter_subnodehash(nodes);
	while((n = ref_iter_subnodehash(nodes))){
		if(n->closed) continue;
		fprintf(out, "N%llu [label=\"N%llu(%llu)\"]\n", (u8i)n->node, (u8i)n->node, (u8i)g->nodes->buffer[n->node].rep_idx);
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				e = ref_subedgev(edges, idx);
				idx = e->next;
				fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d\"]\n", (u8i)n->node, (u8i)e->node->node, "+-"[k], "+-"[e->dir], e->cov);
			}
		}
	}
	fprintf(out, "}\n");
}

thread_beg_def(mrep);
Graph *g;
u8i ret;
u8v *reps;
thread_end_def(mrep);

thread_beg_func(mrep);
subnodehash *nodes;
subedgev *edges;
u4v *rds;
ptrrefv *stack;
u8i nid, tidx, ncpu;
nodes = init_subnodehash(1023);
edges = init_subedgev(32);
rds = init_u4v(32);
stack = init_ptrrefv(32);
tidx = mrep->t_idx;
ncpu = mrep->n_cpu;
thread_beg_loop(mrep);
for(nid=tidx;nid<mrep->g->nodes->size;nid+=ncpu){
	if(mrep->g->nodes->buffer[nid].closed) continue;
	if(evaluate_node_connectivity_graph(mrep->g, nid, rds, nodes, edges, stack) == 0){
		if(0){
			print_subgraph_dot(mrep->g, nid, nodes, edges, stdout);
		}
		mrep->g->nodes->buffer[nid].closed = 1;
		push_u8v(mrep->reps, nid);
		mrep->ret ++;
	}
}
thread_end_loop(mrep);
free_subnodehash(nodes);
free_subedgev(edges);
free_u4v(rds);
free_ptrrefv(stack);
thread_end_func(mrep);

u8i mask_nodes_by_connectivity_graph(Graph *g, int ncpu, FILE *out){
	node_t *n;
	u8i ret, i;
	thread_preprocess(mrep);
	ret = 0;
	thread_beg_init(mrep, ncpu);
	mrep->g   = g;
	mrep->ret = 0;
	mrep->reps = init_u8v(32);
	thread_end_init(mrep);
	thread_wake_all(mrep);
	thread_wait_all(mrep);
	thread_beg_close(mrep);
	if(out){
		for(i=0;i<mrep->reps->size;i++){
			n = ref_nodev(g->nodes, mrep->reps->buffer[i]);
			fprintf(out, "N%llu\t%u\tconn\n", (u8i)mrep->reps->buffer[i], (u4i)n->regs.cnt);
		}
	}
	ret += mrep->ret;
	free_u8v(mrep->reps);
	thread_end_close(mrep);
	return ret;
}

// remove regs which have no high cov edges with other regs in one read
thread_beg_def(mrdk);
Graph *g;
u4i ridx;
u8v *masks;
thread_end_def(mrdk);

thread_beg_func(mrdk);
Graph *g;
u4i ridx;
u8v *regs;
u4v *gidxs, *gcnts;
UUhash *hash;
UUhash_t *UU;
read_t *rd;
reg_t *reg;
node_t *n;
edge_ref_t *f;
edge_t *e;
u8i idx, fidx, nidx, hidx;
u4i i, k, gid, max;
g = mrdk->g;
regs = init_u8v(32);
gidxs = init_u4v(32);
gcnts = init_u4v(32);
hash = init_UUhash(1023);
thread_beg_loop(mrdk);
clear_u8v(mrdk->masks);
ridx = mrdk->ridx;
if(ridx == MAX_U4) continue;
clear_u8v(regs);
clear_u4v(gidxs);
clear_u4v(gcnts);
clear_UUhash(hash);
rd = ref_readv(g->reads, ridx);
idx = rd->regs.idx;
while(idx){
	push_u4v(gidxs, regs->size);
	push_u4v(gcnts, 0);
	reg = ref_regv(g->regs, idx);
	if(reg->closed == 0){
		put_UUhash(hash, (UUhash_t){reg->node, regs->size});
		push_u8v(regs, idx);
	}
	idx = reg->read_link;
}
for(i=0;i<regs->size;i++){
	idx = regs->buffer[i];
	gid = gidxs->buffer[i];
	reg = ref_regv(g->regs, idx);
	n = ref_nodev(g->nodes, reg->node);
	for(k=0;k<2;k++){
		fidx = n->edges[k].idx;
		while(fidx){
			f = ref_edgerefv(g->erefs, fidx);
			fidx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(e->cov < g->min_edge_cov) continue;
			nidx = f->flg? e->node1 : e->node2;
			if((UU = get_UUhash(hash, nidx)) == NULL) continue;
			if((hidx = UU->val) <= i) continue;
			gidxs->buffer[hidx] = gid;
		}
	}
}
for(i=0;i<gidxs->size;i++) gcnts->buffer[gidxs->buffer[i]] ++;
max = 0;
gid = 0;
for(i=0;i<gcnts->size;i++){
	if(gcnts->buffer[i] > max){
		max = gcnts->buffer[i];
		gid = i;
	}
}
for(i=0;i<gidxs->size;i++){
	if(gidxs->buffer[i] != gid){
		push_u8v(mrdk->masks, regs->buffer[i]);
	}
}
thread_end_loop(mrdk);
free_u8v(regs);
free_u4v(gidxs);
free_u4v(gcnts);
free_UUhash(hash);
thread_end_func(mrdk);

u8i mask_read_weak_regs_graph(Graph *g, int ncpu){
	u8i i, ret;
	u4i j;
	thread_preprocess(mrdk);
	thread_beg_init(mrdk, ncpu);
	mrdk->g = g;
	mrdk->ridx = MAX_U4;
	mrdk->masks = init_u8v(32);
	thread_end_init(mrdk);
	ret = 0;
	for(i=0;i<g->reads->size+ncpu;i++){
		if(i < g->reads->size){
			thread_wait_one(mrdk);
		} else {
			thread_wait_next(mrdk);
		}
		if(mrdk->masks->size){
			for(j=0;j<mrdk->masks->size;j++){
				ref_regv(g->regs, mrdk->masks->buffer[j])->closed = 1;
			}
			ret += mrdk->masks->size;
			clear_u8v(mrdk->masks);
		}
		if(i < g->reads->size){
			mrdk->ridx = i;
			thread_wake(mrdk);
		}
	}
	thread_beg_close(mrdk);
	free_u8v(mrdk->masks);
	thread_end_close(mrdk);
	return ret;
}

uint64_t mask_possible_tip_nodes_graph(Graph *g){
	node_t *n;
	reg_t *r;
	uint64_t ret, i;
	uint32_t j, cnt;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		cnt = 0;
		for(j=0;j<n->regs.cnt;j++){
			r = ref_regv(g->regs, n->regs.idx + j);
			if(r->read_link == 0) continue;
			if(n->regs.idx + j == g->reads->buffer[r->rid].regs.idx) continue;
			cnt ++;
		}
		if(cnt < g->min_edge_cov){
			n->closed = 1;
			ret ++;
		}
	}
	return ret;
}

void print_node_edges_graph(Graph *g, u8i nid, int dir, FILE *out){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	uint64_t idx;
	n = ref_nodev(g->nodes, nid);
	idx = n->edges[dir].idx;
	while(idx){
		f = ref_edgerefv(g->erefs, idx);
		idx = f->next;
		e = ref_edgev(g->edges, f->idx);
		if(f->flg){
			fprintf(out, "N%llu\t%c\tN%llu\t%c\t%d\t%d\n", (u8i)e->node2, "+-"[!e->dir1], (u8i)e->node1, "+-"[e->dir1], e->cov, e->off);
		} else {
			fprintf(out, "N%llu\t%c\tN%llu\t%c\t%d\t%d\n", (u8i)e->node1, "+-"[e->dir1], (u8i)e->node2, "+-"[e->dir2], e->cov, e->off);
		}
	}
}

edge_ref_t* first_living_edge_graph(Graph *g, node_t *n, int dir, int *info){
	edge_ref_t *f, *ret;
	uint64_t idx;
	ret = NULL;
	if(info){
		*info = WT_TRACE_MSG_ZERO;
		if(n->edges[dir].cnt == 0) return NULL;
		idx = n->edges[dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			if(g->edges->buffer[f->idx].closed) continue;
			if(ret){ *info = WT_TRACE_MSG_MORE; return NULL; }
			else { *info = WT_TRACE_MSG_ONE; ret = f; }
		}
	} else {
		if(n->edges[dir].cnt == 0) return NULL;
		idx = n->edges[dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			if(g->edges->buffer[f->idx].closed) continue;
			if(ret){ return NULL; }
			else { ret = f; }
		}
	}
	return ret;
}

edge_ref_t* first_living_lnk_graph(Graph *g, frg_t *n, int dir, int *info){
	edge_ref_t *f, *ret;
	uint64_t idx;
	ret = NULL;
	if(info){
		*info = WT_TRACE_MSG_ZERO;
		if(n->lnks[dir].cnt == 0) return NULL;
		idx = n->lnks[dir].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			if(g->lnks->buffer[f->idx].closed) continue;
			if(ret){ *info = WT_TRACE_MSG_MORE; return NULL; }
			else { *info = WT_TRACE_MSG_ONE; ret = f; }
		}
	} else {
		if(n->lnks[dir].cnt == 0) return NULL;
		idx = n->lnks[dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			if(g->lnks->buffer[f->idx].closed) continue;
			if(ret){ return NULL; }
			else { ret = f; }
		}
	}
	return ret;
}

#define count_living_edges_graph(g, n, dir) (n)->edges[dir].cnt

#define count_living_lnks_graph(g, n, dir) (n)->lnks[dir].cnt

// dir = 2 means either strand
edge_ref_t* edge_node2node_graph(Graph *g, u8i node1, int dir1, u8i node2, int dir2){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	uint64_t idx;
	int dire;
	n = ref_nodev(g->nodes, node1);
	if(dir1 > 1){
		dir1 = 0; dire = 2;
	} else {
		dire = dir1 + 1;
	}
	while(dir1 < dire){
		idx = n->edges[dir1].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(f->flg){
				if(e->node1 == node2 && (dir2 > 1? 1 : (dir2 == (!e->dir1)))) return f;
			} else {
				if(e->node2 == node2 && (dir2 > 1? 1 : (dir2 == e->dir2))) return f;
			}
		}
		dir1 ++;
	}
	return NULL;
}

uint64_t linear_trace_graph(Graph *g, tracev *path, uint64_t max_step, int *msg){
	trace_t *t;
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	uint64_t step;
	int dir, info;
	if(path->size == 0){
		if(msg) *msg = 3;
		return 0;
	}
	t = ref_tracev(path, path->size - 1);
	step = 0;
	while(step < max_step){
		f = first_living_edge_graph(g, ref_nodev(g->nodes, t->node), t->dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = info; break; }
		e = g->edges->buffer + f->idx;
		n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
		dir = f->flg? !e->dir1 : e->dir2;
		t->edges[t->dir] = *f;
		t = next_ref_tracev(path);
		t->node = n - g->nodes->buffer;
		t->dir = dir;
		t->edges[!t->dir] = (edge_ref_t){f->idx, !f->flg, 0};
		t->edges[t->dir] = EDGE_REF_NULL;
		f = first_living_edge_graph(g, n, !dir, &info);
		if(info != WT_TRACE_MSG_ONE){  path->size --; if(msg) *msg = -1 - info; break; }
		step ++;
	}
	return step;
}

uint64_t linear_path_graph(Graph *g, pathv *path, int max_len, int *msg){
	path_t *t;
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	int len;
	int dir, info;
	if(path->size == 0){
		if(msg) *msg = 3;
		return 0;
	}
	t = ref_pathv(path, path->size - 1);
	len = ref_frgv(g->frgs, t->frg)->len;
	while(len < max_len){
		f = first_living_lnk_graph(g, ref_frgv(g->frgs, t->frg), t->dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = info; break; }
		e = g->lnks->buffer + f->idx;
		len += e->off;
		n = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
		dir = f->flg? !e->dir1 : e->dir2;
		t->lnks[t->dir] = *f;
		t = next_ref_pathv(path);
		t->frg = n - g->frgs->buffer;
		t->dir = dir;
		t->lnks[!t->dir] = (edge_ref_t){f->idx, !f->flg, 0};
		t->lnks[t->dir] = EDGE_REF_NULL;
		f = first_living_lnk_graph(g, n, !dir, &info);
		if(info != WT_TRACE_MSG_ONE){  path->size --; if(msg) *msg = -1 - info; break; }
		len += n->len;
	}
	return len;
}

int cal_offset_traces_graph(Graph *g, tracev *path, u8i beg, u8i end, int offset){
	trace_t *t;
	node_t *n;
	reg_t *r;
	edge_t *e;
	uint64_t i;
	int off;
	off = offset;
	for(i=beg;i<end;i++){
		t = ref_tracev(path, i);
		t->off = off;
		n = ref_nodev(g->nodes, t->node);
		r = ref_regv(g->regs, n->regs.idx);
		if(t->edges[t->dir].idx == EDGE_REF_NULL.idx){
			off += r->end - r->beg;
		} else {
			e = ref_edgev(g->edges, t->edges[t->dir].idx);
			off += r->end - r->beg + e->off;
		}
	}
	return off;
}

int cal_offset_paths_graph(Graph *g, pathv *path, u8i beg, u8i end){
	path_t *t;
	frg_t *n;
	lnk_t *e;
	uint64_t i;
	int off;
	off = 0;
	for(i=beg;i<end;i++){
		t = ref_pathv(path, i);
		t->off = off;
		n = ref_frgv(g->frgs, t->frg);
		if(t->lnks[t->dir].idx == EDGE_REF_NULL.idx){
			off += n->length;
		} else {
			e = ref_lnkv(g->lnks, t->lnks[t->dir].idx);
			off += ((i == beg)? n->length : n->len) + e->off;
		}
	}
	return off;
}

uint64_t true_linear_unique_trace_graph(Graph *g, tracev *path, uint64_t max_step, uint64_t visit, int *msg){
	trace_t *t;
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	uint64_t step;
	int dir, info;
	if(path->size == 0){
		if(msg) *msg = WT_TRACE_MSG_ZERO;
		return 0;
	}
	step = 0;
	if(msg) *msg = WT_TRACE_MSG_ONE;
	while(step < max_step){
		t = ref_tracev(path, path->size - 1);
		f = first_living_edge_graph(g, ref_nodev(g->nodes, t->node), !t->dir, &info);
		if(info == WT_TRACE_MSG_MORE){ if(path->size > 1) path->size --; if(msg) *msg = -1 - info; break; }
		n = ref_nodev(g->nodes, t->node);
		n->bt_visit = visit;
		f = first_living_edge_graph(g, ref_nodev(g->nodes, t->node),  t->dir, &info);
		if(info == WT_TRACE_MSG_ZERO){ if(msg) *msg = info; break; }
		else if(info == WT_TRACE_MSG_MORE){ if(path->size > 1) path->size --; if(msg) *msg = info; break; }
		e = g->edges->buffer + f->idx;
		n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
		if(n->bt_visit == visit){ if(msg) *msg = WT_TRACE_MSG_VISITED; break; }
		dir = f->flg? !e->dir1 : e->dir2;
		t->edges[t->dir] = *f;
		t = next_ref_tracev(path);
		t->node = n - g->nodes->buffer;
		t->dir = dir;
		t->edges[!t->dir] = (edge_ref_t){f->idx, !f->flg, 0};
		t->edges[t->dir] = EDGE_REF_NULL;
		step ++;
	}
	return step;
}

uint64_t true_linear_unique_path_graph(Graph *g, pathv *path, uint64_t max_step, uint64_t visit, int *msg){
	path_t *t;
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	uint64_t step;
	int dir, info;
	if(path->size == 0){
		if(msg) *msg = WT_TRACE_MSG_ZERO;
		return 0;
	}
	if(msg) *msg = WT_TRACE_MSG_ONE;
	step = 0;
	while(step < max_step){
		t = ref_pathv(path, path->size - 1);
		f = first_living_lnk_graph(g, ref_frgv(g->frgs, t->frg), !t->dir, &info);
		if(info == WT_TRACE_MSG_MORE){ if(path->size > 1) path->size --; if(msg) *msg = -1 - info; break; }
		n = ref_frgv(g->frgs, t->frg);
		n->bt_visit = visit;
		f = first_living_lnk_graph(g, ref_frgv(g->frgs, t->frg),  t->dir, &info);
		if(info == WT_TRACE_MSG_ZERO){ if(msg) *msg = info; break; }
		else if(info == WT_TRACE_MSG_MORE){ if(path->size > 1) path->size --; if(msg) *msg = info; break; }
		e = g->lnks->buffer + f->idx;
		n = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
		if(n->bt_visit == visit){ if(msg) *msg = WT_TRACE_MSG_VISITED; break; }
		dir = f->flg? !e->dir1 : e->dir2;
		t->lnks[t->dir] = *f;
		t = next_ref_pathv(path);
		t->frg = n - g->frgs->buffer;
		t->dir = dir;
		t->lnks[!t->dir] = (edge_ref_t){f->idx, !f->flg, 0};
		t->lnks[t->dir] = EDGE_REF_NULL;
		step ++;
	}
	return step;
}

uint64_t count_linear_trace_graph(Graph *g, trace_t *t, uint64_t max_step, int *msg){
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	uint64_t step;
	int dir, info;
	step = 0;
	if(msg) *msg = WT_TRACE_MSG_ONE;
	while(step < max_step){
		f = first_living_edge_graph(g, ref_nodev(g->nodes, t->node), t->dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = info; break; }
		e = g->edges->buffer + f->idx;
		n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
		dir = f->flg? !e->dir1 : e->dir2;
		t->node = n - g->nodes->buffer;
		t->dir = dir;
		f = first_living_edge_graph(g, n, !dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = -1 - info; break; }
		step ++;
	}
	return step;
}

int count_linear_path_graph(Graph *g, path_t *t, int max_len, int *msg){
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	int len;
	int dir, info;
	if(msg) *msg = WT_TRACE_MSG_ONE;
	len = ref_frgv(g->frgs, t->frg)->len;
	while(len < max_len){
		f = first_living_lnk_graph(g, ref_frgv(g->frgs, t->frg), t->dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = info; break; }
		e = g->lnks->buffer + f->idx;
		len += e->off;
		n = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
		dir = f->flg? !e->dir1 : e->dir2;
		t->frg = n - g->frgs->buffer;
		t->dir = dir;
		f = first_living_lnk_graph(g, n, !dir, &info);
		if(info != WT_TRACE_MSG_ONE){ if(msg) *msg = -1 - info; break; }
		len += n->len;
	}
	return len;
}

u4i del_isolated_nodes_graph(Graph *g, FILE *log){
	node_t *n;
	u4i ret, i;
	int f, r;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		first_living_edge_graph(g, n, 0, &f);
		if(f != WT_TRACE_MSG_ZERO) continue;
		first_living_edge_graph(g, n, 1, &r);
		if(r != WT_TRACE_MSG_ZERO) continue;
		n->closed = 1;
		if(log) fprintf(log, "DEL_ISO\tN%u\n", i);
		ret ++;
	}
	return ret;
}

uint64_t cut_binary_edges_graph(Graph *g){
	UUhash *hash;
	UUhash_t *u;
	node_t *n;
	edge_ref_t *f;
	edge_t *e, *p;
	uint64_t idx, nid, ret;
	ret = 0;
	hash = init_UUhash(15);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		clear_UUhash(hash);
		idx = n->edges[0].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(e->closed < WT_EDGE_CLOSED_LESS){
				put_UUhash(hash, (UUhash_t){f->flg? e->node1 : e->node2, f->idx});
			}
		}
		idx = n->edges[1].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(e->closed >= WT_EDGE_CLOSED_LESS) continue;
			if((u = get_UUhash(hash, f->flg? e->node1 : e->node2)) == NULL) continue;
			p = ref_edgev(g->edges, u->val);
			if(0){
				if(p->cov > e->cov) cut_edge_core_graph(g, e, WT_EDGE_CLOSED_HARD);
				else cut_edge_core_graph(g, p, WT_EDGE_CLOSED_HARD);
				ret ++;
			} else {
				cut_edge_core_graph(g, e, WT_EDGE_CLOSED_HARD);
				cut_edge_core_graph(g, p, WT_EDGE_CLOSED_HARD);
				ret += 2;
			}
		}
	}
	free_UUhash(hash);
	return ret;
}

uint64_t _old_cut_binary_lnks_graph(Graph *g, FILE *info){
	UUhash *hash;
	UUhash_t *u;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e, *p;
	uint64_t idx, nid, ret;
	ret = 0;
	hash = init_UUhash(15);
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		clear_UUhash(hash);
		idx = n->lnks[0].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = ref_lnkv(g->lnks, f->idx);
			if(e->closed < WT_EDGE_CLOSED_LESS){
				put_UUhash(hash, (UUhash_t){f->flg? e->frg1 : e->frg2, f->idx});
			}
		}
		idx = n->lnks[1].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = ref_lnkv(g->lnks, f->idx);
			if(e->closed >= WT_EDGE_CLOSED_LESS) continue;
			if((u = get_UUhash(hash, f->flg? e->frg1 : e->frg2)) == NULL) continue;
			p = ref_lnkv(g->lnks, u->val);
			if(info){
				fprintf(info, "BINARY_LINK\tF%d\t%c\tF%d\t%c\t%d\t%d\n", e->frg1, "+-"[e->dir1], e->frg2, "+-"[e->dir2], e->cov, e->off);
				fprintf(info, "BINARY_LINK\tF%d\t%c\tF%d\t%c\t%d\t%d\n", p->frg1, "+-"[p->dir1], p->frg2, "+-"[p->dir2], p->cov, p->off);
			}
			if(1){
				if(p->cov > e->cov) cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
				else if(p->cov < e->cov) cut_lnk_core_graph(g, p, WT_EDGE_CLOSED_HARD);
				else {
					cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
					cut_lnk_core_graph(g, p, WT_EDGE_CLOSED_HARD);
					ret ++;
				}
				ret ++;
			} else {
				cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
				cut_lnk_core_graph(g, p, WT_EDGE_CLOSED_HARD);
				ret += 2;
			}
		}
	}
	free_UUhash(hash);
	return ret;
}

uint64_t cut_binary_lnks_graph(Graph *g, FILE *info){
	UUhash *hash;
	UUhash_t *u;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e, *p;
	u8i idx, nid, wid, ret;
	u4i k;
	int exists;
	ret = 0;
	hash = init_UUhash(15);
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		clear_UUhash(hash);
		for(k=0;k<2;k++){
			idx = n->lnks[k].idx;
			while(idx){
				f = ref_edgerefv(g->lrefs, idx);
				idx = f->next;
				e = ref_lnkv(g->lnks, f->idx);
				if(e->closed != WT_EDGE_CLOSED_HARD){
					wid = f->flg? e->frg1 : e->frg2;
					u = prepare_UUhash(hash, wid, &exists);
					if(exists){
						p = ref_lnkv(g->lnks, u->val);
						if(info){
							fprintf(info, "BINARY_LINK\tF%d\t%c\tF%d\t%c\t%d\t%d\n", e->frg1, "+-"[e->dir1], e->frg2, "+-"[e->dir2], e->cov, e->off);
							fprintf(info, "BINARY_LINK\tF%d\t%c\tF%d\t%c\t%d\t%d\n", p->frg1, "+-"[p->dir1], p->frg2, "+-"[p->dir2], p->cov, p->off);
						}
						if(p->cov < e->cov){
							cut_lnk_core_graph(g, p, WT_EDGE_CLOSED_HARD);
							u->val = f->idx;
						} else if(p->cov > e->cov){
							cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
						} else {
							cut_lnk_core_graph(g, p, WT_EDGE_CLOSED_HARD);
							cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
							ret ++;
						}
						ret ++;
					} else {
						u->key = wid;
						u->val = f->idx;
					}
				}
			}
		}
	}
	free_UUhash(hash);
	return ret;
}

u8i cut_low_cov_lnks_graph(Graph *g, int low_cov){
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e;
	u8i idx, nid, ret;
	u4i k;
	int max_cov;
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		for(k=0;k<2;k++){
			max_cov = 0;
			idx = n->lnks[k].idx;
			while(idx){
				f = ref_edgerefv(g->lrefs, idx);
				idx = f->next;
				e = ref_lnkv(g->lnks, f->idx);
				if(e->cov > max_cov) max_cov = e->cov;
			}
			if(max_cov <= low_cov || max_cov < (int)g->min_edge_cov) continue;
			idx = n->lnks[k].idx;
			while(idx){
				f = ref_edgerefv(g->lrefs, idx);
				idx = f->next;
				e = ref_lnkv(g->lnks, f->idx);
				if(e->cov <= low_cov){
					cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint32_t rescue_low_cov_transitive_edges_graph(Graph *g, uint64_t nid, u8v *edges, UUhash *hash){
	node_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	edge_t *e, *e1, *e2;
	reg_t *r;
	UUhash_t *u;
	uint64_t idx, nid2, nid3;
	uint32_t i, k, k2, k3, k4, ret;
	int off1, off2, yes;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		clear_UUhash(hash);
		clear_u8v(edges);
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			e = ref_edgev(g->edges, f->idx);
			//if(e->closed < WT_EDGE_CLOSED_LESS){
			if(e->closed == 0){
				put_UUhash(hash, (UUhash_t){f->flg? e->node1: e->node2, idx});
			} else if(e->closed == WT_EDGE_CLOSED_LESS){
				push_u8v(edges, idx);
			}
			idx = f->next;
		}
		if(edges->size <= 1) continue;
		for(i=0;i<edges->size;i++){
			f = ref_edgerefv(g->erefs, edges->buffer[i]);
			e = ref_edgev(g->edges, f->idx);
			//if(e->status) continue;
			if(e->closed != WT_EDGE_CLOSED_LESS) continue;
			if(f->flg){ nid2 = e->node1; k2 = !e->dir1; }
			else      { nid2 = e->node2; k2 =  e->dir2; }
			yes = 0;
			w = ref_nodev(g->nodes, nid2);
			idx = w->edges[k2].idx;
			while(idx){
				f2 = ref_edgerefv(g->erefs, idx);
				idx = f2->next;
				if(f2->flg){ nid3 = g->edges->buffer[f2->idx].node1; k3 = !g->edges->buffer[f2->idx].dir1; }
				else       { nid3 = g->edges->buffer[f2->idx].node2; k3 =  g->edges->buffer[f2->idx].dir2; }
				e2 = ref_edgev(g->edges, f2->idx);
				//if(e2->closed >= WT_EDGE_CLOSED_LESS) continue;
				if(e2->closed) continue;
				if((u = get_UUhash(hash, nid3)) == NULL) continue;
				f3 = ref_edgerefv(g->erefs, u->val);
				e1 = ref_edgev(g->edges, f3->idx);
				k4 = f3->flg? !e1->dir1 : e1->dir2;
				if(k3 != k4) continue;
				v = ref_nodev(g->nodes, nid3);
				off1 = off2 = (g->regs->buffer[n->regs.idx].end - g->regs->buffer[n->regs.idx].beg) + (g->regs->buffer[v->regs.idx].end - g->regs->buffer[v->regs.idx].beg);
				off1 += e1->off;
				off2 += e->off + e2->off + (g->regs->buffer[w->regs.idx].end - g->regs->buffer[w->regs.idx].beg);
				r = ref_regv(g->regs, w->regs.idx);
				//if(e->off + e2->off + (r->end - r->beg) >= longest) continue;
				if(num_diff(off1, off2) >= num_min(off1, off2)) continue;
				yes = 1;
				//revive_edge_graph(g, e);
				e->flag |= (1 << f->flg);
				ret ++;
				break;
			}
			if(0){
				if(yes) continue;
				idx = w->edges[!k2].idx;
				while(idx){
					f2 = ref_edgerefv(g->erefs, idx);
					idx = f2->next;
					if(f2->flg){ nid3 = g->edges->buffer[f2->idx].node1; k3 = !g->edges->buffer[f2->idx].dir1; }
					else       { nid3 = g->edges->buffer[f2->idx].node2; k3 =  g->edges->buffer[f2->idx].dir2; }
					e2 = ref_edgev(g->edges, f2->idx);
					if(e2->closed >= WT_EDGE_CLOSED_LESS) continue;
					if((u = get_UUhash(hash, nid3)) == NULL) continue;
					f3 = ref_edgerefv(g->erefs, u->val);
					e1 = ref_edgev(g->edges, f3->idx);
					k4 = f3->flg? !e1->dir1 : e1->dir2;
					if(k3 != k4) continue;
					v = ref_nodev(g->nodes, nid3);
					off1 = off2 = (g->regs->buffer[n->regs.idx].end - g->regs->buffer[n->regs.idx].beg) + (g->regs->buffer[w->regs.idx].end - g->regs->buffer[w->regs.idx].beg);
					off1 += e->off;
					off2 += e1->off + e2->off + (g->regs->buffer[v->regs.idx].end - g->regs->buffer[v->regs.idx].beg);
					r = ref_regv(g->regs, w->regs.idx);
					if(num_diff(off1, off2) >= num_min(off1, off2)) continue;
					yes = 1;
					//revive_edge_graph(g, e);
					e->flag |= (1 << f->flg);
					ret ++;
					break;
				}
			}
		}
	}
	return ret;
}

uint64_t rescue_low_cov_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	uint64_t i, nid, ret;
	ret = 0;
	edges = init_u8v(32);
	hash = init_UUhash(13);
	for(i=0;i<g->edges->size;i++) g->edges->buffer[i].flag = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		rescue_low_cov_transitive_edges_graph(g, nid, edges, hash);
	}
	for(i=0;i<g->edges->size;i++){
		//if(g->edges->buffer[i].flag == 3){
		if(g->edges->buffer[i].flag == 3){
			revive_edge_graph(g, g->edges->buffer + i);
			ret ++;
		}
		g->edges->buffer[i].flag = 0;
	}
	free_u8v(edges);
	free_UUhash(hash);
	return ret;
}

uint32_t rescue_low_cov_tip_edges_core(Graph *g, uint64_t nid){
	node_t *n, *w, *ww;
	edge_t *e, *ee;
	edge_ref_t *f;
	uint64_t idx, wid;
	uint32_t k, dir, ret;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	if(n->edges[0].cnt == 0 && n->edges[1].cnt == 0) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		if(n->edges[k].cnt) continue;
		idx = n->edges[k].idx;
		ee = NULL;
		ww = NULL;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			if(f->flg){
				wid = e->node1; dir = !e->dir1;
			} else {
				wid = e->node2; dir = e->dir2;
			}
			w = ref_nodev(g->nodes, wid);
			//if(w->edges[!dir].cnt) continue;
			if(w->edges[dir].cnt == 0) continue;
			if(ee == NULL || e->cov > ee->cov || (e->cov == ee->cov && w->regs.cnt > ww->regs.cnt)){ ee = e; ww = w; }
		}
		if(ee){ revive_edge_graph(g, ee); ret ++; }
	}
	return ret;
}

uint64_t rescue_low_cov_tip_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	uint64_t nid, ret;
	ret = 0;
	edges = init_u8v(32);
	hash = init_UUhash(13);
	for(nid=0;nid<g->nodes->size;nid++){
		ret += rescue_low_cov_tip_edges_core(g, nid);
	}
	free_u8v(edges);
	free_UUhash(hash);
	return ret;
}

int rescue_mercy_edge_core_graph(Graph *g, u4i rid, BitVec *tips[2]){
	read_t *rd;
	edge_t *e;
	edge_ref_t *f;
	reg_t *r, *s;
	u8i idx;
	int c, d, ret;
	rd = ref_readv(g->reads, rid);
	if(rd->regs.cnt < 2) return 0;
	idx = rd->regs.idx;
	s = NULL;
	c = 0;
	ret = 0;
	while(idx){
		r = ref_regv(g->regs, idx);
		idx = r->read_link;
		if(r->closed) continue;
		if(!(r->beg >= rd->clps[0] && r->end <= rd->clps[1])) continue;
		if(g->nodes->buffer[r->node].closed) continue;
		d = 0;
		if(get_bitvec(tips[0], r->node)){
			d |= 1;
		}
		if(get_bitvec(tips[1], r->node)){
			d |= 2;
		}
		if(d == 3) continue;
		if(s){
			if(d == (1 << (!r->dir))){
				f = edge_node2node_graph(g, s->node, (0b100 >> c) & 0x01, r->node, !((0b100 >> d) & 0x01));
				if(f){
					e = ref_edgev(g->edges, f->idx);
					revive_edge_graph(g, e);
					zero_bitvec(tips[e->dir1], e->node1);
					zero_bitvec(tips[!e->dir2], e->node1);
					ret ++;
#ifdef __DEBUG__
				} else {
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
#endif
				}
				s = NULL;
				continue;
			} else {
				s = NULL;
			}
		}
		if(d == (1 << (r->dir))){
			s = r;
			c = d;
		}
	}
	return ret;
}

u8i rescue_mercy_edges_graph(Graph *g){
	BitVec *tips[2];
	node_t *n;
	u8i i, ret;
	ret = 0;
	tips[0] = init_bitvec(g->nodes->size);
	tips[1] = init_bitvec(g->nodes->size);
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		if(n->edges[0].cnt == 0) one_bitvec(tips[0], i);
		if(n->edges[1].cnt == 0) one_bitvec(tips[1], i);
	}
	for(i=0;i<g->reads->size;i++){
		ret += rescue_mercy_edge_core_graph(g, i, tips);
	}
	free_bitvec(tips[0]);
	free_bitvec(tips[1]);
	return ret;
}

uint32_t rescue_weak_tip_lnks_core(Graph *g, uint64_t nid){
	frg_t *n, *w, *ww;
	lnk_t *e, *ee;
	edge_ref_t *f;
	uint64_t idx, wid;
	uint32_t k, dir, ret;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	if(n->lnks[0].cnt == 0 && n->lnks[1].cnt == 0) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		if(n->lnks[k].cnt) continue;
		idx = n->lnks[k].idx;
		ee = NULL;
		ww = NULL;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = ref_lnkv(g->lnks, f->idx);
			if(e->weak == 0) continue;
			if(f->flg){
				wid = e->frg1; dir = !e->dir1;
			} else {
				wid = e->frg2; dir = e->dir2;
			}
			w = ref_frgv(g->frgs, wid);
			if(w->lnks[!dir].cnt) continue;
			if(ee == NULL){ ee = e; ww = w; }
			else { ee = NULL; break; }
		}
		if(ee){ revive_lnk_graph(g, ee); ret ++; }
	}
	return ret;
}

uint64_t rescue_weak_tip_lnks2_graph(Graph *g){
	uint64_t nid, ret;
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		ret += rescue_weak_tip_lnks_core(g, nid);
	}
	return ret;
}

u8i rescue_weak_tip_lnks_graph(Graph *g){
	u8v *weaks[2];
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	u8i nid, i, ret;
	u4i k, eidx, idx;
	weaks[0] = init_u8v(g->frgs->size);
	weaks[1] = init_u8v(g->frgs->size);
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		for(k=0;k<2;k++){
			if(n->lnks[k].cnt){
				eidx = 0;
			} else {
				idx = n->lnks[k].idx;
				eidx = 0;
				while(idx){
					f = ref_edgerefv(g->lrefs, idx);
					e = ref_lnkv(g->lnks, f->idx);
					if(e->weak){
						if(eidx == 0) eidx = f->idx;
						else eidx = MAX_VALUE_U4;
					}
					idx = f->next;
				}
			}
			push_u8v(weaks[k], eidx);
		}
	}
	for(k=0;k<2;k++){
		for(i=0;i<weaks[k]->size;i++){
			if(weaks[k]->buffer[i] == 0 || weaks[k]->buffer[i] == MAX_VALUE_U4) continue;
			e = ref_lnkv(g->lnks, weaks[k]->buffer[i]);
			if(i != e->frg1) continue;
			if(weaks[0]->buffer[e->frg2] == weaks[k]->buffer[i] || weaks[1]->buffer[e->frg2] == weaks[k]->buffer[i]){
				ret ++;
				revive_lnk_graph(g, e);
			}
		}
	}
	free_u8v(weaks[0]);
	free_u8v(weaks[1]);
	return ret;
}

static inline int _scoring_edge_orders(Graph *g, uint64_t fidx){
	edge_ref_t *f;
	edge_t *e;
	node_t *n;
	int score;
	f = ref_edgerefv(g->erefs, fidx);
	e = ref_edgev(g->edges, f->idx);
	n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
	score = e->off + (n->regs.cnt * -5) + (e->cov * -5);
	return score;
}


uint32_t reduce_transitive_edges_core_graph(Graph *g, uint64_t nid, u8v *edges, UUhash *hash, uint32_t closed_val){
	node_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	edge_t *e, *e1, *e2;
	UUhash_t *u;
	uint64_t idx, nid2, nid3;
	uint32_t i, k, k2, k3, k4, ret;
	int off1, off2;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		clear_UUhash(hash);
		clear_u8v(edges);
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			e = ref_edgev(g->edges, f->idx);
			if(e->closed < closed_val){
				put_UUhash(hash, (UUhash_t){f->flg? e->node1: e->node2, (idx << 1) | e->closed});
				push_u8v(edges, idx);
			}
			idx = f->next;
		}
		if(edges->size <= 1) continue;
		// sort the edges by composition of e->off and e->cov
		//sort_array(edges->buffer, edges->size, uint64_t, num_cmpgt(_scoring_edge_orders(g, a), _scoring_edge_orders(g, b)));
		for(i=0;i<edges->size;i++){
			f = ref_edgerefv(g->erefs, edges->buffer[i]);
			e = ref_edgev(g->edges, f->idx);
			if(f->flg){ nid2 = e->node1; k2 = !e->dir1; }
			else      { nid2 = e->node2; k2 =  e->dir2; }
			w = ref_nodev(g->nodes, nid2);
			idx = w->edges[k2].idx;
			while(idx){
				f2 = ref_edgerefv(g->erefs, idx);
				idx = f2->next;
				if(f2->flg){ nid3 = g->edges->buffer[f2->idx].node1; k3 = !g->edges->buffer[f2->idx].dir1; }
				else       { nid3 = g->edges->buffer[f2->idx].node2; k3 =  g->edges->buffer[f2->idx].dir2; }
				e2 = ref_edgev(g->edges, f2->idx);
				//if(e2->closed) continue;
				if(e2->closed >= closed_val) continue;
				if((u = get_UUhash(hash, nid3)) == NULL) continue;
				if(u->val & 0x01) continue; // already deleted
				f3 = ref_edgerefv(g->erefs, u->val >> 1);
				e1 = ref_edgev(g->edges, f3->idx);
				k4 = f3->flg? !e1->dir1 : e1->dir2;
				if(k3 != k4) continue;
				v = ref_nodev(g->nodes, nid3);
				off1 = off2 = (g->regs->buffer[n->regs.idx].end - g->regs->buffer[n->regs.idx].beg) + (g->regs->buffer[v->regs.idx].end - g->regs->buffer[v->regs.idx].beg);
				off1 += e1->off;
				off2 += e->off + e2->off + (g->regs->buffer[w->regs.idx].end - g->regs->buffer[w->regs.idx].beg);
				// check whether off1 and off2 diff too much
				if(num_diff(off1, off2) >= num_min(off1, off2)) continue;
				u->val |= 1;
			}
		}
		reset_iter_UUhash(hash);
		while((u = ref_iter_UUhash(hash))){
			if(u->val & 0x01){
				e = ref_edgev(g->edges, ref_edgerefv(g->erefs, u->val >> 1)->idx);
				if(e->closed == WT_EDGE_CLOSED_NULL){
					cut_edge_graph(g, e);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint64_t reduce_transitive_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	uint64_t nid, ret;
	ret = 0;
	edges = init_u8v(32);
	hash = init_UUhash(13);
	for(nid=0;nid<g->nodes->size;nid++){
		ret += reduce_transitive_edges_core_graph(g, nid, edges, hash, 2);
	}
	free_u8v(edges);
	free_UUhash(hash);
	return ret;
}

void set_init_ends_graph(Graph *g){
	node_t *n;
	u8i nid;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->edges[0].cnt == 0 || n->edges[1].cnt == 0){
			n->init_end = 1;
		}
	}
}

// node_t->unvisit is used to indicate vacant, inplay and eliminated

u4i myers_transitive_reduction_core_graph(Graph *g, u8i nid, float _fuzz, int closed){
	node_t *n, *w, *x;
	edge_ref_t *f, *f2;
	edge_t *e, *e2;
	u8i idx, idx2;
	u4i k, k2, ret;
	int longest, len, fuzz;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		longest = 0;
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			e = ref_edgev(g->edges, f->idx);
			//if(e->closed == 0 && !(e->flag & (1 << f->flg))){
			if(e->closed <= closed){
				w = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
				w->unvisit = 1; // inplay
				len = e->off + g->reglen;
				if(len > longest) longest = len;
			}
			idx = f->next;
		}
		fuzz = (int)_fuzz;
		if(longest * (_fuzz - fuzz) > fuzz) fuzz = longest * (_fuzz - fuzz);
		longest += fuzz;
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			e = ref_edgev(g->edges, f->idx);
			//if(e->closed == 0 && !(e->flag & (1 << f->flg))){
			if(e->closed <= closed){
				w = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
				if(w->unvisit == 1){
					k2 = f->flg? !e->dir1 : e->dir2;
					idx2 = w->edges[k2].idx;
					while(idx2){
						f2 = ref_edgerefv(g->erefs, idx2);
						e2 = ref_edgev(g->edges, f2->idx);
						// TODO: check the strand
						//if(e2->closed == 0 && !(e2->flag & (1 << f2->flg))){
						if(e2->closed <= closed){
							x = ref_nodev(g->nodes, f2->flg? e2->node1 : e2->node2);
							if(x->unvisit == 1){
								len = e->off + e2->off + g->reglen;
								if(len <= longest){
									x->unvisit = 2; // eliminated
								}
							}
						}
						idx2 = f2->next;
					}
				}
			}
			idx = f->next;
		}
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			e = ref_edgev(g->edges, f->idx);
			//if(e->closed == 0 && !(e->flag & (1 << f->flg))){
			if(e->closed <= closed){
				w = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
				//if(w->unvisit == 1){
					k2 = f->flg? !e->dir1 : e->dir2;
					idx2 = w->edges[k2].idx;
					while(idx2){
						f2 = ref_edgerefv(g->erefs, idx2);
						e2 = ref_edgev(g->edges, f2->idx);
						//if(e2->closed == 0 && !(e2->flag & (1 << f2->flg))){
						if(e2->closed <= closed){
							x = ref_nodev(g->nodes, f2->flg? e2->node1 : e2->node2);
							if(x->unvisit == 1){
								if(e2->off <= fuzz || idx2 == w->edges[k2].idx){
									x->unvisit = 2; // eliminated
								}
							}
						}
						idx2 = f2->next;
					}
				//}
			}
			idx = f->next;
		}
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			e = ref_edgev(g->edges, f->idx);
			//if(e->closed == 0 && !(e->flag & (1 << f->flg))){
			if(e->closed <= closed){
				w = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
				if(w->unvisit == 2){
					e->flag |= 1 << f->flg;
					ret ++;
				}
				w->unvisit = 0;
			}
			idx = f->next;
		}
	}
	return ret;
}

u8i myers_transitive_reduction_graph(Graph *g, float fuzz){
	node_t *n;
	edge_t *e;
	u8i nid, eid, ret;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		n->unvisit = 0; // vacant
	}
	for(eid=0;eid<g->edges->size;eid++){
		e = ref_edgev(g->edges, eid);
		e->flag = 0;
	}
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		myers_transitive_reduction_core_graph(g, nid, fuzz, WT_EDGE_CLOSED_NULL);
	}
	ret = 0;
	for(eid=0;eid<g->edges->size;eid++){
		e = ref_edgev(g->edges, eid);
		if(e->flag == 3){
			cut_edge_graph(g, e);
			ret ++;
		}
		e->flag = 0;
	}
	return ret;
}

u4i myers_transitive_reduction_core_frg_graph(Graph *g, u8i nid, float _fuzz){
	frg_t *n, *w, *x;
	edge_ref_t *f, *f2;
	lnk_t *e, *e2;
	u8i idx, idx2;
	u4i k, k2, ret;
	int longest, len, fuzz;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		longest = 0;
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			//if(e->closed == 0 && !(e->flag & (1 << f->flg))){
			if(e->closed == 0){
				w = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
				w->unvisit = 1; // inplay
				len = e->off + g->reglen;
				if(len > longest) longest = len;
			}
			idx = f->next;
		}
		fuzz = (int)_fuzz;
		if(longest * (_fuzz - fuzz) > fuzz) fuzz = longest * (_fuzz - fuzz);
		longest += fuzz;
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			//if(e->closed == 0 && !(e->flag & (1 << f->flg))){
			if(e->closed == 0){
				w = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
				if(w->unvisit == 1){
					k2 = f->flg? !e->dir1 : e->dir2;
					idx2 = w->lnks[k2].idx;
					while(idx2){
						f2 = ref_edgerefv(g->lrefs, idx2);
						e2 = ref_lnkv(g->lnks, f2->idx);
						// TODO: check the strand
						//if(e2->closed == 0 && !(e2->flag & (1 << f2->flg))){
						if(e2->closed == 0){
							x = ref_frgv(g->frgs, f2->flg? e2->frg1 : e2->frg2);
							if(x->unvisit == 1){
								len = e->off + e2->off + g->reglen;
								if(len <= longest){
									x->unvisit = 2; // eliminated
								}
							}
						}
						idx2 = f2->next;
					}
				}
			}
			idx = f->next;
		}
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			//if(e->closed == 0 && !(e->flag & (1 << f->flg))){
			if(e->closed == 0){
				w = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
				//if(w->unvisit == 1){
					k2 = f->flg? !e->dir1 : e->dir2;
					idx2 = w->lnks[k2].idx;
					while(idx2){
						f2 = ref_edgerefv(g->lrefs, idx2);
						e2 = ref_lnkv(g->lnks, f2->idx);
						//if(e2->closed == 0 && !(e2->flag & (1 << f2->flg))){
						if(e2->closed == 0){
							x = ref_frgv(g->frgs, f2->flg? e2->frg1 : e2->frg2);
							if(x->unvisit == 1){
								if(e2->off <= fuzz || idx2 == w->lnks[k2].idx){
									x->unvisit = 2; // eliminated
								}
							}
						}
						idx2 = f2->next;
					}
				//}
			}
			idx = f->next;
		}
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			//if(e->closed == 0 && !(e->flag & (1 << f->flg))){
			if(e->closed == 0){
				w = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
				if(w->unvisit == 2){
					e->flag |= 1 << f->flg;
					ret ++;
				}
				w->unvisit = 0;
			}
			idx = f->next;
		}
	}
	return ret;
}

u8i myers_transitive_reduction_frg_graph(Graph *g, float fuzz){
	frg_t *n;
	lnk_t *e;
	u8i nid, eid, ret;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		n->unvisit = 0; // vacant
	}
	for(eid=0;eid<g->lnks->size;eid++){
		e = ref_lnkv(g->lnks, eid);
		e->flag = 0;
	}
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		myers_transitive_reduction_core_frg_graph(g, nid, fuzz);
	}
	ret = 0;
	for(eid=0;eid<g->lnks->size;eid++){
		e = ref_lnkv(g->lnks, eid);
		if(e->flag == 3){
			cut_lnk_graph(g, e);
			ret ++;
		}
		e->flag = 0;
	}
	return ret;
}

u4i detach_repetitive_frg_core_graph(Graph *g, u8i nid, u4i max_dist, u8i visit, u8v *nds, u8v *heap){
	frg_t *n, *w, *x;
	edge_ref_t *f;
	lnk_t *e;
	u8i idx, wid, xid, eid;
	u4i k, i, j, bidx, bcnt[2];
	int off;
	n = ref_frgv(g->frgs, nid);
	n->bt_visit = visit;
	clear_u8v(nds);
	clear_u8v(heap);
	for(k=0;k<2;k++){
		bcnt[k] = 0;
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			idx = f->next;
			if(e->closed) continue;
			wid = f->flg? e->frg1 : e->frg2;
			w = ref_frgv(g->frgs, wid);
			if(w->bt_visit == visit){
				// TODO: remove the fprintf after debug
				fprintf(stderr, " -- F%llu:%c -> F%llu:%c already visited in %s -- %s:%d --\n", nid, "+-"[k], wid, "+-"[w->rep_dir], __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				return 0;
			}
			push_u8v(nds, wid);
			w->bt_visit = visit;
			w->rep_dir  = f->flg? !e->dir1 : e->dir2;
			w->bt_idx = (bcnt[k] << 1) | k;
			bcnt[k] ++;
			off = n->len / 2 + e->off;
			if(off < 0) off = 0;
			w->rep_idx  = off + w->len; // rep_idx is used as dist
			if(w->rep_idx < max_dist){
				array_heap_push(heap->buffer, heap->size, heap->cap, u8i, wid, num_cmp(g->frgs->buffer[a].rep_idx, g->frgs->buffer[b].rep_idx));
			}
		}
	}
	if(bcnt[0] == 0 || bcnt[1] == 0) return 0;
	if(bcnt[0] < 2 && bcnt[1] < 2) return 0;
	// extending branches to max_dist in length
	while(heap->size){
		wid = heap->buffer[0];
		array_heap_remove(heap->buffer, heap->size, heap->cap, u8i, 0, num_cmp(g->frgs->buffer[a].rep_idx, g->frgs->buffer[b].rep_idx));
		w = ref_frgv(g->frgs, wid);
		idx = w->lnks[w->rep_dir].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			idx = f->next;
			if(e->closed) continue;
			xid = f->flg? e->frg1 : e->frg2;
			x = ref_frgv(g->frgs, xid);
			if(x->bt_visit == visit) continue;
			off = w->rep_idx + e->off;
			if(off < 0) off = 0;
			if(off > (int)max_dist) continue;
			push_u8v(nds, xid);
			x->bt_visit = visit;
			x->rep_dir = f->flg? !e->dir1 : e->dir2;
			x->bt_idx = w->bt_idx;
			x->rep_idx = off + x->len;
			if(x->rep_idx < max_dist){
				array_heap_push(heap->buffer, heap->size, heap->cap, u8i, xid, num_cmp(g->frgs->buffer[a].rep_idx, g->frgs->buffer[b].rep_idx));
			}
		}
	}
	// find cross links
	encap_and_zeros_u8v(heap, bcnt[0] * bcnt[1]); // reuse heap as matrix
	for(i=0;i<nds->size;i++){
		wid = nds->buffer[i];
		w = ref_frgv(g->frgs, wid);
		idx = w->lnks[!w->rep_dir].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			idx = f->next;
			if(e->closed == WT_EDGE_CLOSED_HARD) continue;
			xid = f->flg? e->frg1 : e->frg2;
			if(xid == nid) continue;
			x = ref_frgv(g->frgs, xid);
			if(x->bt_visit != visit) continue;
			if(((w->bt_idx ^ x->bt_idx) & 0x01) == 0) continue;
			if(w->bt_idx & 0x01){
				if(x->bt_idx & 0x01){
					continue;
				} else {
					bidx = (x->bt_idx >> 1) * bcnt[1] + (w->bt_idx >> 1);
				}
			} else {
				if(x->bt_idx & 0x01){
					bidx = (w->bt_idx >> 1) * bcnt[1] + (x->bt_idx >> 1);
				} else {
					continue;
				}
			}
			eid = heap->buffer[bidx];
			if(eid == 0 || e->cov > g->lnks->buffer[eid].cov){
				heap->buffer[bidx] = eid = f->idx;
			}
		}
	}
	// find one-left to one-right links
	clear_u8v(nds); // will reuse nds to store bidx
	for(i=0;i<bcnt[0];i++){
		k = bcnt[1];
		for(j=0;j<bcnt[1];j++){
			if(heap->buffer[i * bcnt[1] + j]){
				if(k < bcnt[1]){
					k = bcnt[1];
					break;
				} else {
					k = j;
				}
			}
		}
		for(j=0;k<bcnt[1]&&j<bcnt[0];j++){
			if(j == i) continue;
			if(heap->buffer[j * bcnt[1] + k]){
				k = bcnt[1];
				break;
			}
		}
		if(k == bcnt[1]) continue;
		push_u8v(nds, i * bcnt[1] + k);
	}
	// detach repeats
	for(i=0;i<nds->size;i++){
		//bid[0] = nds->buffer[i] / bcnt[0];
		//bid[1] = nds->buffer[i] % bcnt[0];
		for(k=0;k<2;k++){
			idx = n->lnks[k].idx;
			bidx = k? (nds->buffer[i] % bcnt[0]) : (nds->buffer[i] / bcnt[0]);
			while(idx){
				f = ref_edgerefv(g->lrefs, idx);
				e = ref_lnkv(g->lnks, f->idx);
				idx = f->next;
				if(e->closed) continue;
				wid = f->flg? e->frg1 : e->frg2;
				w = ref_frgv(g->frgs, wid);
				if((w->bt_idx >> 1) == bidx){
					//cut_lnk_graph(g, e);
					cut_lnk_core_graph(g, e, WT_EDGE_CLOSED_HARD); // in case of loop cut and revive
				}
			}
		}
	}
	for(i=0;i<nds->size;i++){
		eid = heap->buffer[nds->buffer[i]];
		e = ref_lnkv(g->lnks, eid);
		revive_lnk_graph(g, e);
	}
	return nds->size;
}

u4i detach_repetitive_frg_graph(Graph *g, u4i max_dist){
	u8v *nds;
	u8v *heap;
	frg_t *n;
	u8i nid, visit;
	u4i ret;
	nds = init_u8v(32);
	heap = init_u8v(32);
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		n->bt_visit = 0;
	}
	visit = 0;
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->lnks[0].cnt == 0 || n->lnks[1].cnt == 0) continue;
		if(n->lnks[0].cnt < 2 && n->lnks[1].cnt < 2) continue;
		ret += detach_repetitive_frg_core_graph(g, nid, max_dist, ++visit, nds, heap);
	}
	free_u8v(nds);
	free_u8v(heap);
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		n->bt_visit = 0;
	}
	return ret;
}

uint32_t reduce_transitive_lnks_core_graph(Graph *g, uint64_t nid, u8v *lnks, UUhash *hash, uint32_t closed_val){
	frg_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	lnk_t *e, *e1, *e2;
	UUhash_t *u;
	uint64_t idx, nid2, nid3;
	uint32_t i, k, k2, k3, k4, ret;
	int off1, off2;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	ret = 0;
	for(k=0;k<2;k++){
		clear_UUhash(hash);
		clear_u8v(lnks);
		idx = n->lnks[k].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			e = ref_lnkv(g->lnks, f->idx);
			if(e->closed < closed_val){
				put_UUhash(hash, (UUhash_t){f->flg? e->frg1: e->frg2, (idx << 1) | e->closed});
				push_u8v(lnks, idx);
			}
			idx = f->next;
		}
		if(lnks->size <= 1) continue;
		// sort the lnks by composition of e->off and e->cov
		//sort_array(lnks->buffer, lnks->size, uint64_t, num_cmpgt(_scoring_lnk_orders(g, a), _scoring_lnk_orders(g, b)));
		for(i=0;i<lnks->size;i++){
			f = ref_edgerefv(g->lrefs, lnks->buffer[i]);
			e = ref_lnkv(g->lnks, f->idx);
			if(f->flg){ nid2 = e->frg1; k2 = !e->dir1; }
			else      { nid2 = e->frg2; k2 =  e->dir2; }
			w = ref_frgv(g->frgs, nid2);
			idx = w->lnks[k2].idx;
			while(idx){
				f2 = ref_edgerefv(g->lrefs, idx);
				idx = f2->next;
				if(f2->flg){ nid3 = g->lnks->buffer[f2->idx].frg1; k3 = !g->lnks->buffer[f2->idx].dir1; }
				else       { nid3 = g->lnks->buffer[f2->idx].frg2; k3 =  g->lnks->buffer[f2->idx].dir2; }
				e2 = ref_lnkv(g->lnks, f2->idx);
				//if(e2->closed) continue;
				if(e2->closed >= closed_val) continue;
				if((u = get_UUhash(hash, nid3)) == NULL) continue;
				if(u->val & 0x01) continue; // already deleted
				f3 = ref_edgerefv(g->lrefs, u->val >> 1);
				e1 = ref_lnkv(g->lnks, f3->idx);
				k4 = f3->flg? !e1->dir1 : e1->dir2;
				if(k3 != k4) continue;
				v = ref_frgv(g->frgs, nid3);
				off1 = off2 = n->len + v->len;
				off1 += e1->off;
				off2 += e->off + e2->off + w->len;
				// check whether off1 and off2 diff too much
				if(num_diff(off1, off2) >= num_min(off1, off2)) continue;
				u->val |= 1;
			}
		}
		reset_iter_UUhash(hash);
		while((u = ref_iter_UUhash(hash))){
			if(u->val & 0x01){
				e = ref_lnkv(g->lnks, ref_edgerefv(g->lrefs, u->val >> 1)->idx);
				if(e->closed == WT_EDGE_CLOSED_NULL){
					cut_lnk_graph(g, e);
					ret ++;
				}
			}
		}
	}
	return ret;
}

uint64_t reduce_transitive_lnks_graph(Graph *g){
	u8v *lnks;
	UUhash *hash;
	uint64_t nid, ret;
	ret = 0;
	lnks = init_u8v(32);
	hash = init_UUhash(13);
	for(nid=0;nid<g->frgs->size;nid++){
		ret += reduce_transitive_lnks_core_graph(g, nid, lnks, hash, 2);
	}
	free_u8v(lnks);
	free_UUhash(hash);
	return ret;
}

uint64_t trim_tip_core_graph(Graph *g, uint16_t max_step, tracev *path, uint64_t nid, int hard_trim){
	trace_t *t, T;
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	uint64_t ret, idx;
	uint32_t i, dir, step, step2, found, n_in;
	int msg1, msg2;
	if(g->cut_tip == 0) return 0;
	ret = 0;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	if(n->edges[0].cnt == 1 && n->edges[1].cnt == 0){
		dir = 0;
	} else if(n->edges[0].cnt == 0 && n->edges[1].cnt == 1){
		dir = 1;
	} else {
		return 0;
	}
	clear_tracev(path);
	t = next_ref_tracev(path);
	t->node = nid;
	t->edges[0] = EDGE_REF_NULL;
	t->edges[1] = EDGE_REF_NULL;
	t->dir = dir;
	step = linear_trace_graph(g, path, max_step, &msg1) + 1;
	if(step > max_step) return 0;
	//if(msg1 != -1 - WT_TRACE_MSG_MORE && msg1 != WT_TRACE_MSG_MORE) return 0;
	if(msg1 == WT_TRACE_MSG_MORE){
		if(!hard_trim) return 0;
		path->size --;
	} else if(msg1 == -1 - WT_TRACE_MSG_MORE){
		t = ref_tracev(path, path->size); // please see linear_trace_graph
		n = ref_nodev(g->nodes, t->node);
		dir = !t->dir;
		n_in = 0;
		idx = n->edges[dir].idx;
		found = 0;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			if(f->idx == t->edges[dir].idx) continue;
			e = g->edges->buffer + f->idx;
			if(e->closed) continue;
			n_in ++;
			T.node = f->flg? e->node1 : e->node2;
			T.dir  = f->flg? !e->dir1 : e->dir2;
			step2 = count_linear_trace_graph(g, &T, step + 1, &msg2) + 1;
			if(msg2 != WT_TRACE_MSG_ZERO) step2 ++;
			if(step2 > step){
				found = 1;
				break; 
			} else if(step2 == step && (e->cov > g->edges->buffer[t->edges[dir].idx].cov || (step == 1 && e->cov == g->edges->buffer[t->edges[dir].idx].cov && e->off >= g->edges->buffer[t->edges[dir].idx].off))){
				found = 1;
				break; 
			}
			//if(step2 + 1 >= step && msg2 != WT_TRACE_MSG_ZERO){ found = 1; break; }
		}
		if(!found) return 0;
	} else return 0;
	for(i=0;i<path->size;i++){
		del_node_graph(g, ref_nodev(g->nodes, path->buffer[i].node));
		//del_node_edges_graph(g, ref_nodev(g->nodes, path->buffer[i].node));
		ret ++;
	}
	return ret;
}

uint64_t trim_tips_graph(Graph *g, uint16_t max_step, int hard_trim){
	tracev *path;
	node_t *n;
	uint64_t ret, nid;
	ret = 0;
	path = init_tracev(32);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(trim_tip_core_graph(g, max_step, path, nid, hard_trim)) ret ++;
	}
	free_tracev(path);
	return ret;
}

// careful sharp_tip -> blunt_tip -> sharp_tip and so on
u4i trim_blunt_tip_core_graph(Graph *g, u8i nid){
	node_t *n;
	int k;
	if(g->cut_tip == 0) return 0;
	n = ref_nodev(g->nodes, nid);
	if(n->edges[0].cnt && n->edges[1].cnt) return 0;
	k = (n->edges[0].cnt == 0);
	if(n->edges[k].cnt < 2) return 0;
	del_node_graph(g, n);
	return 1;
}

uint64_t trim_blunt_tips_graph(Graph *g){
	node_t *n;
	uint64_t ret, nid;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(trim_blunt_tip_core_graph(g, nid)) ret ++;
	}
	return ret;
}

uint64_t trim_frgtip_core_graph(Graph *g, int max_len, pathv *path, uint64_t nid){
	path_t *t, T;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e;
	uint64_t ret, idx;
	uint32_t i, dir, found, n_in;
	int len, len2;
	int msg1, msg2;
	if(g->cut_tip == 0) return 0;
	ret = 0;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	first_living_lnk_graph(g, n, 0, &msg1);
	first_living_lnk_graph(g, n, 1, &msg2);
	if(msg1 != WT_TRACE_MSG_ZERO){
		if(msg2 != WT_TRACE_MSG_ZERO) return 0;
		dir = 0;
	} else if(msg2 != WT_TRACE_MSG_ZERO){
		dir = 1;
	} else return 0;
	clear_pathv(path);
	t = next_ref_pathv(path);
	t->frg = nid;
	t->lnks[0] = EDGE_REF_NULL;
	t->lnks[1] = EDGE_REF_NULL;
	t->dir = dir;
	len = linear_path_graph(g, path, max_len, &msg1) + 1;
	if(len > max_len) return 0;
	//if(msg1 != -1 - WT_TRACE_MSG_MORE && msg1 != WT_TRACE_MSG_MORE) return 0;
	if(msg1 == WT_TRACE_MSG_MORE){
		path->size --;
	} else if(msg1 == -1 - WT_TRACE_MSG_MORE){
		t = ref_pathv(path, path->size); // please see linear_path_graph
		n = ref_frgv(g->frgs, t->frg);
		dir = !t->dir;
		n_in = 0;
		idx = n->lnks[dir].idx;
		found = 0;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			if(f->idx == t->lnks[dir].idx) continue;
			e = g->lnks->buffer + f->idx;
			if(e->closed) continue;
			n_in ++;
			T.frg = f->flg? e->frg1 : e->frg2;
			T.dir  = f->flg? !e->dir1 : e->dir2;
			len2 = count_linear_path_graph(g, &T, len + 1, &msg2) + 1;
			if(msg2 != WT_TRACE_MSG_ZERO) len2 ++;
			if(len2 >= len){ found = 1; break; }
		}
		if(!found) return 0;
	} else return 0;
	for(i=0;i<path->size;i++){
		del_frg_lnks_graph(g, ref_frgv(g->frgs, path->buffer[i].frg));
		ret ++;
	}
	return ret;
}

uint64_t trim_frgtips_graph(Graph *g, int max_len){
	pathv *path;
	frg_t *n;
	uint64_t ret, nid;
	ret = 0;
	path = init_pathv(32);
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		if(trim_frgtip_core_graph(g, max_len, path, nid)) ret ++;
	}
	free_pathv(path);
	return ret;
}

typedef struct {
	node_t *n;
	edge_t *e; // incoming edge
	uint64_t dir:1, ind:1, step:8, bt:16, ending:16, score:20, keep:2;
} bt_t;
define_list(btv, bt_t);
#define WT_MAX_BTIDX	0xFFFF

uint32_t pop_bubble_backtrace_graph(Graph *g, btv *bts, uint32_t idx){
	bt_t *bt;
	uint32_t ret;
	while(idx){
		bt = ref_btv(bts, idx);
		bt->step = 0;
		idx = bt->bt;
	}
	ret = 0;
	for(idx=1;idx<bts->size;idx++){
		bt = ref_btv(bts, idx);
		if(bt->step == 0) continue;
		cut_edge_graph(g, bt->e);
		ret ++;
	}
	return ret;
}

typedef struct {
	frg_t *n;
	lnk_t *e; // incoming edge
	uint64_t dir:1, ind:1, step:8, bt:16, ending:16, score:22;
} frg_bt_t;
define_list(frgbtv, frg_bt_t);
uint32_t pop_frg_bubble_backtrace_graph(Graph *g, frgbtv *bts, uint32_t idx){
	frg_bt_t *bt;
	uint32_t ret;
	while(idx){
		bt = ref_frgbtv(bts, idx);
		bt->step = 0;
		idx = bt->bt;
	}
	ret = 0;
	for(idx=1;idx<bts->size;idx++){
		bt = ref_frgbtv(bts, idx);
		if(bt->step == 0) continue;
		cut_lnk_graph(g, bt->e);
		ret ++;
	}
	return ret;
}

uint32_t pop_bubble2_backtrace_graph(Graph *g, btv *bts, uint32_t _idx){
	bt_t *bt;
	uint32_t ret, i, idx;
	for(i=1;i<bts->size;i++){
		bt = ref_btv(bts, i);
		bt->keep = 2;
	}
	idx = _idx;
	while(idx){
		bt = ref_btv(bts, idx);
		bt->keep = 1;
		idx = bt->bt;
	}
	for(i=1;i<bts->size;i++){
		bt = ref_btv(bts, i);
		if(bt->keep == 1) continue;
		if(bt->ending == 0) continue;
		if(bts->buffer[bt->ending].keep == 1){
			idx = i;
			while(idx){
				bt = ref_btv(bts, idx);
				if(bt->keep != 2) break;
				idx = bt->bt;
				bt->keep = 0;
			}
		}
	}
	ret = 0;
	for(idx=1;idx<bts->size;idx++){
		bt = ref_btv(bts, idx);
		if(bt->keep) continue;
		cut_edge_graph(g, bt->e);
		ret ++;
	}
	return ret;
}

uint32_t pop_frg_bubble2_backtrace_graph(Graph *g, frgbtv *bts, uint32_t _idx){
	frg_bt_t *bt;
	uint32_t ret, i, idx;
	for(i=1;i<bts->size;i++){
		bt = ref_frgbtv(bts, i);
		if(bt->ending == _idx){
			idx = i;
			while(idx){
				bt = ref_frgbtv(bts, idx);
				bt->step = 0;
				idx = bt->bt;
			}
		}
	}
	idx = _idx;
	while(idx){
		bt = ref_frgbtv(bts, idx);
		bt->step = 1;
		idx = bt->bt;
	}
	ret = 0;
	for(idx=1;idx<bts->size;idx++){
		bt = ref_frgbtv(bts, idx);
		if(bt->step != 0) continue;
		cut_lnk_graph(g, bt->e);
		ret ++;
	}
	return ret;
}

uint32_t safe_cut_redundant_edges_graph(Graph *g, btv *bts, bt_t *b1, bt_t *b2){
	uint32_t ret;
	ret = 0;
	if(0){
		ret = 1;
		cut_edge_graph(g, b1->e);
		b1 = ref_btv(bts, b1->bt);
	}
	while(1){
		if(b1->step >= b2->step){
			if(b1 == b2) break;
			ret ++;
			cut_edge_graph(g, b1->e);
			b1 = ref_btv(bts, b1->bt);
		} else {
			b2 = ref_btv(bts, b2->bt);
		}
	}
	return ret;
}

uint32_t safe_cut_redundant_lnks_graph(Graph *g, frgbtv *bts, frg_bt_t *b1, frg_bt_t *b2){
	uint32_t ret;
	ret = 0;
	if(0){
		ret = 1;
		cut_lnk_graph(g, b1->e);
		b1 = ref_frgbtv(bts, b1->bt);
	}
	while(1){
		if(b1->step >= b2->step){
			if(b1 == b2) break;
			ret ++;
			cut_lnk_graph(g, b1->e);
			b1 = ref_frgbtv(bts, b1->bt);
		} else {
			b2 = ref_frgbtv(bts, b2->bt);
		}
	}
	return ret;
}

uint32_t pop_bubble_core_graph(Graph *g, uint16_t max_step, btv *bts, u4v *heap, uint64_t nid, uint32_t dir, uint64_t visit, int safe){
	bt_t *bt, *tb;
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	uint64_t ret, idx;
	uint32_t bidx, i, lst, unclosed;
	ret = 0;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	if(count_living_edges_graph(g, n, dir) < 2) return 0;
	clear_btv(bts);
	next_ref_btv(bts);
	bt = next_ref_btv(bts);
	bt->n = n;
	bt->e = NULL;
	bt->dir = dir;
	bt->ind = safe? 0 : 1;
	bt->step = 0;
	bt->bt = 0;
	bt->score = 0;
	bt->ending = 0;
	clear_u4v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	n->single_in = 1;
	unclosed = 0;
	while(heap->size && heap->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, uint32_t, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
		encap_btv(bts, bts->buffer[bidx].n->edges[bts->buffer[bidx].dir].cnt);
		bt = ref_btv(bts, bidx);
		if(bt->step >= max_step) return 0;
		if(bt->ind && bt->n->single_in == 0) bt->ind = 0;
		lst = bts->size;
		idx = bt->n->edges[bt->dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = g->edges->buffer + f->idx;
			if(e->closed) continue;
			tb = next_ref_btv(bts);
			tb->n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
			if(tb->n == bts->buffer[1].n) return 0;
			tb->e = e;
			tb->dir = f->flg? !e->dir1 : e->dir2;
			tb->step = bt->step + 1;
			tb->bt   = bidx;
			tb->ind = 0;
			tb->score = bt->score + e->cov;
			tb->ending = 0;
		}
		if(bt->ind && (bt->bt == 0 || lst + 1 == bts->size)){
			for(i=lst;i<bts->size;i++) bts->buffer[i].ind = 1;
		}
		for(i=lst;i<bts->size;i++){
			tb = ref_btv(bts, i);
			if(tb->n->bt_visit != visit){
				tb->n->bt_visit = visit;
				tb->n->unvisit = count_living_edges_graph(g, tb->n, !tb->dir);
				if(tb->n->unvisit == 1) tb->n->single_in = 1;
				else tb->n->single_in = 0;
				tb->n->bt_idx = i;
				unclosed ++;
			} else {
				tb->ending = tb->n->bt_idx;
				if(tb->dir == bts->buffer[tb->n->bt_idx].dir){
					if(tb->ind && bts->buffer[tb->n->bt_idx].ind){
						if(tb->step == bts->buffer[tb->n->bt_idx].step){
							return safe_cut_redundant_edges_graph(g, bts, tb, ref_btv(bts, tb->n->bt_idx));
						} else {
							return safe_cut_redundant_edges_graph(g, bts, ref_btv(bts, tb->n->bt_idx), tb);
						}
					} else if(tb->ind){
						return safe_cut_redundant_edges_graph(g, bts, tb, ref_btv(bts, tb->n->bt_idx));
					} else if(bts->buffer[tb->n->bt_idx].ind){
						return safe_cut_redundant_edges_graph(g, bts, ref_btv(bts, tb->n->bt_idx), tb);
					}
				} else {
					// circle
					return 0;
				}
			}
			tb->n->unvisit --;
			if(tb->n->unvisit == 0){
				if(tb->step > bts->buffer[tb->n->bt_idx].step){
					bts->buffer[tb->n->bt_idx].ending = i;
					tb->n->bt_idx = i;
					tb->ending = 0;
				}
				if(count_living_edges_graph(g, tb->n, tb->dir)){
					array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
				}
				unclosed --;
			}
		}
		if(heap->size == 1 && unclosed == 0){
			return pop_bubble2_backtrace_graph(g, bts, heap->buffer[0]);
		}
	}
	return 0;
}

uint64_t pop_bubbles_graph(Graph *g, uint16_t max_step, int safe){
	btv *bts;
	u4v *heap;
	node_t *n;
	uint64_t nid, visit, ret, _ret;
	int dir;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++) g->nodes->buffer[nid].bt_visit = 0;
	bts = init_btv(32);
	heap = init_u4v(32);
	visit = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		for(dir=0;dir<2;dir++){
			_ret = pop_bubble_core_graph(g, max_step, bts, heap, nid, dir, ++visit, safe);
			if(_ret) ret ++;
		}
	}
	free_btv(bts);
	free_u4v(heap);
	return ret;
}

uint32_t pop_frg_bubble_core_graph(Graph *g, uint16_t max_step, frgbtv *bts, u4v *heap, uint64_t nid, uint32_t dir, uint64_t visit){
	frg_bt_t *bt, *tb;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e;
	uint64_t ret, idx;
	uint32_t bidx, i, lst, unclosed;
	ret = 0;
	n = ref_frgv(g->frgs, nid);
	if(n->closed) return 0;
	if(count_living_lnks_graph(g, n, dir) < 2) return 0;
	clear_frgbtv(bts);
	next_ref_frgbtv(bts);
	bt = next_ref_frgbtv(bts);
	bt->n = n;
	bt->e = NULL;
	bt->dir = dir;
	bt->ind = 1;
	//bt->ind = 0;
	bt->step = 0;
	bt->bt = 0;
	bt->score = 0;
	bt->ending = 0;
	clear_u4v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	n->single_in = 1;
	unclosed = 0;
	while(heap->size && heap->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, uint32_t, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
		encap_frgbtv(bts, bts->buffer[bidx].n->lnks[bts->buffer[bidx].dir].cnt);
		bt = ref_frgbtv(bts, bidx);
		if(bt->step >= max_step) return 0;
		if(bt->ind && bt->n->single_in == 0) bt->ind = 0;
		lst = bts->size;
		idx = bt->n->lnks[bt->dir].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = g->lnks->buffer + f->idx;
			if(e->closed) continue;
			tb = next_ref_frgbtv(bts);
			tb->n = ref_frgv(g->frgs, f->flg? e->frg1 : e->frg2);
			if(tb->n == bts->buffer[1].n) return 0;
			tb->e = e;
			tb->dir = f->flg? !e->dir1 : e->dir2;
			tb->step = bt->step + 1;
			tb->bt   = bidx;
			tb->ind = 0;
			tb->score = bt->score + e->cov;
			tb->ending = 0;
		}
		if(bt->ind && (bt->bt == 0 || lst + 1 == bts->size)){
			for(i=lst;i<bts->size;i++) bts->buffer[i].ind = 1;
		}
		for(i=lst;i<bts->size;i++){
			tb = ref_frgbtv(bts, i);
			if(tb->n->bt_visit != visit){
				tb->n->bt_visit = visit;
				tb->n->unvisit = count_living_lnks_graph(g, tb->n, !tb->dir);
				if(tb->n->unvisit == 1) tb->n->single_in = 1;
				else tb->n->single_in = 0;
				tb->n->bt_idx = i;
				unclosed ++;
			} else {
				tb->ending = tb->n->bt_idx;
				if(tb->dir == bts->buffer[tb->n->bt_idx].dir){
					if(tb->ind && bts->buffer[tb->n->bt_idx].ind){
						if(tb->step == bts->buffer[tb->n->bt_idx].step){
							return safe_cut_redundant_lnks_graph(g, bts, tb, ref_frgbtv(bts, tb->n->bt_idx));
						} else {
							return safe_cut_redundant_lnks_graph(g, bts, ref_frgbtv(bts, tb->n->bt_idx), tb);
						}
					} else if(tb->ind){
						return safe_cut_redundant_lnks_graph(g, bts, tb, ref_frgbtv(bts, tb->n->bt_idx));
					} else if(bts->buffer[tb->n->bt_idx].ind){
						return safe_cut_redundant_lnks_graph(g, bts, ref_frgbtv(bts, tb->n->bt_idx), tb);
					}
				}
			}
			tb->n->unvisit --;
			if(tb->n->unvisit == 0){
				if(tb->step > bts->buffer[tb->n->bt_idx].step) tb->n->bt_idx = i;
				if(count_living_lnks_graph(g, tb->n, tb->dir)){
					array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
				}
				unclosed --;
			}
		}
		if(heap->size == 1 && unclosed == 0){
			return pop_frg_bubble2_backtrace_graph(g, bts, heap->buffer[0]);
		}
	}
	return 0;
}

uint64_t pop_frg_bubbles_graph(Graph *g, uint16_t max_step){
	frgbtv *bts;
	u4v *heap;
	frg_t *n;
	uint64_t nid, visit, ret, _ret;
	int dir;
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++) g->frgs->buffer[nid].bt_visit = 0;
	bts = init_frgbtv(32);
	heap = init_u4v(32);
	visit = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		for(dir=0;dir<2;dir++){
			_ret = pop_frg_bubble_core_graph(g, max_step, bts, heap, nid, dir, ++visit);
			if(_ret) ret ++;
		}
	}
	free_frgbtv(bts);
	free_u4v(heap);
	return ret;
}

u8i remove_boomerangs_frg_graph(Graph *g, u4i max_frg_len){
	frg_t *frg;
	u8i i, ret;
	ret = 0;
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		if(frg->closed) continue;
		if(frg->len > max_frg_len) continue;
		if(frg->lnks[0].cnt == 0 && frg->lnks[1].cnt > 1){
		} else if(frg->lnks[1].cnt == 0 && frg->lnks[0].cnt > 1){
		} else continue;
		ret ++;
		del_frg_lnks_graph(g, frg);
	}
	return ret;
}

u8i cut_weak_branches_frg_graph(Graph *g){
	frg_t *frg1, *frg2;
	lnk_t *lnk;
	u8v *cuts;
	u8i i, ret;
	cuts = init_u8v(32);
	for(i=0;i<g->lnks->size;i++){
		lnk = ref_lnkv(g->lnks, i);
		if(lnk->weak == 0) continue;
		frg1 = ref_frgv(g->frgs, lnk->frg1);
		frg2 = ref_frgv(g->frgs, lnk->frg2);
		if(frg1->lnks[lnk->dir1].cnt > 1){
			push_u8v(cuts, i);
		} else if(frg2->lnks[!lnk->dir2].cnt > 1){
			push_u8v(cuts, i);
		}
	}
	ret = cuts->size;
	for(i=0;i<cuts->size;i++){
		cut_lnk_graph(g, ref_lnkv(g->lnks, cuts->buffer[i]));
	}
	free_u8v(cuts);
	return ret;
}

u4i resolve_yarn_core_graph(Graph *g, u4i max_step, btv *bts, u4v *heap, u8i nid, u4i dir, u8i visit){
	bt_t *bt, *tb;
	node_t *n, *m;
	edge_ref_t *f;
	edge_t *e;
	uint64_t ret, idx, tip_idx;
	uint32_t bidx, i, lst, tip;
	int n_in;
	ret = 0;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	if(count_living_edges_graph(g, n, dir) < 2) return 0;
	clear_btv(bts);
	next_ref_btv(bts);
	bt = next_ref_btv(bts);
	bt->n = n;
	bt->e = NULL;
	bt->dir = dir;
	bt->step = 0;
	bt->bt = 0;
	bt->score = 0;
	bt->ending = 0;
	clear_u4v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	tip = 0; tip_idx = WT_MAX_BTIDX;
	n_in = 1;
	while(heap->size && bts->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, uint32_t, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
		encap_btv(bts, bts->buffer[bidx].n->edges[bts->buffer[bidx].dir].cnt);
		bt = ref_btv(bts, bidx);
		if(bt->step >= max_step){
			return 0;
		}
		lst = bts->size;
		idx = bt->n->edges[bt->dir].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = g->edges->buffer + f->idx;
			if(e->closed) continue;
			tb = next_ref_btv(bts);
			tb->n = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
			//if(tb->n == bts->buffer[1].n) return 0;
			tb->e = e;
			tb->dir = f->flg? !e->dir1 : e->dir2;
			tb->step = bt->step + 1;
			tb->bt   = bidx;
			tb->score = bt->score + e->cov;
			tb->ending = 0;
		}
		for(i=lst;i<bts->size;i++){
			tb = ref_btv(bts, i);
			if(tb->n->bt_visit != visit){
				tb->n->bt_visit = visit;
				tb->n->unvisit = count_living_edges_graph(g, tb->n, !tb->dir);
				n_in += tb->n->unvisit - 1;
				tb->n->bt_idx = i;
				if(count_living_edges_graph(g, tb->n, tb->dir)){
					array_heap_push(heap->buffer, heap->size, heap->cap, uint32_t, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
				} else if(tip == 0){
					tip = 1;
					tip_idx = i;
				} else {
					return 0;
				}
			} else {
				tb->ending = tb->n->bt_idx;
				n_in --;
			}
		}
		if(n_in == 1){
			if(heap->size == 0 && tip == 1){
				return pop_bubble_backtrace_graph(g, bts, tip? tip_idx : heap->buffer[0]);
			} else if(heap->size == 1 && tip == 0){
				tip_idx = heap->buffer[0];
				tb = ref_btv(bts, heap->buffer[0]);
				idx = tb->n->edges[tb->dir].idx;
				while(idx){
					f = ref_edgerefv(g->erefs, idx);
					idx = f->next;
					e = g->edges->buffer + f->idx;
					if(e->closed) continue;
					m = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
					if(m->bt_visit == visit) continue;
					else tip ++;
				}
				if(tip == 0) return 0;
				return pop_bubble_backtrace_graph(g, bts, tip_idx);
			}
		}
	}
	return 0;
}

// very complicated local region, like yarn, but with Single In edge and Single Out edges
u8i resolve_yarns_graph(Graph *g, u4i max_step){
	btv *bts;
	u4v *heap;
	node_t *n;
	uint64_t nid, visit, ret, _ret;
	int dir;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++) g->nodes->buffer[nid].bt_visit = 0;
	bts = init_btv(32);
	heap = init_u4v(32);
	visit = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(n->edges[0].cnt <= 1 && n->edges[1].cnt > 1){
			dir = 1;
		} else if(n->edges[1].cnt <= 1 && n->edges[0].cnt > 1){
			dir = 0;
		} else continue;
		_ret = resolve_yarn_core_graph(g, max_step, bts, heap, nid, dir, ++visit);
		if(_ret) ret ++;
	}
	free_btv(bts);
	free_u4v(heap);
	return ret;
}

u8i mask_all_branching_nodes_graph(Graph *g){
	node_t *n;
	u8i node, ret;
	ret = 0;
	for(node=0;node<g->nodes->size;node++){
		n = ref_nodev(g->nodes, node);
		if(n->closed) continue;
		if(n->edges[0].cnt > 1 || n->edges[1].cnt > 1){
			n->rep_idx = 1;
			ret ++;
		} else {
			n->rep_idx = 0;
		}
	}
	for(node=0;node<g->nodes->size;node++){
		n = ref_nodev(g->nodes, node);
		if(n->closed) continue;
		if(n->rep_idx == 0) continue;
		del_node_graph(g, n);
	}
	return ret;
}

uint64_t gen_unitigs_graph(Graph *g){
	tracev *path;
	u4v *lens;
	trace_t *t;
	node_t *n;
	uint64_t nid, nutg, i;
	for(i=0;i<g->utgs->size;i++) free_tracev(g->utgs->buffer[i]);
	clear_vplist(g->utgs);
	lens = init_u4v(1024);
	nutg = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		g->nodes->buffer[nid].bt_visit = 0;
		g->nodes->buffer[nid].rep_idx  = MAX_REP_IDX;
	}
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(n->bt_visit) continue;
		path = init_tracev(4);
		nutg ++;
		t = next_ref_tracev(path);
		t->node = nid;
		t->edges[0] = EDGE_REF_NULL;
		t->edges[1] = EDGE_REF_NULL;
		t->dir = 0;
		true_linear_unique_trace_graph(g, path, 0xFFFFFFFFFFFFFFFFLLU, nutg, NULL);
		reverse_tracev(path);
		for(i=0;i<path->size;i++) path->buffer[i].dir = !path->buffer[i].dir;
		true_linear_unique_trace_graph(g, path, 0xFFFFFFFFFFFFFFFFLLU, nutg, NULL);
		push_u4v(lens, cal_offset_traces_graph(g, path, 0, path->size, 0));
		for(i=0;i<path->size;i++){
			ref_nodev(g->nodes, path->buffer[i].node)->rep_idx = g->utgs->size;
		}
		push_vplist(g->utgs, path);
	}
	fprintf(KBM_LOGF, "[%s] ", date()); num_n50(lens, KBM_LOGF); fprintf(KBM_LOGF, "\n");
	free_u4v(lens);
	return nutg;
}

seqletv* path2seqlets_graph(Graph *g, pathv *path){
	seqletv *qs;
	path_t *p1, *p2;
	frg_t *frg1, *frg2;
	edge_ref_t *f;
	lnk_t *l;
	trace_t *t1, *t2;
	b8i off;
	int len;
	u4i i, j;
	qs = init_seqletv(4);
	cal_offset_paths_graph(g, path, 0, path->size);
	p1   = NULL;
	frg1 = NULL;
	for(i=0;i<path->size;i++){
		p2 = ref_pathv(path, i);
		frg2 = ref_frgv(g->frgs, p2->frg);
		cal_offset_traces_graph(g, g->traces, frg2->toff, frg2->toff + frg2->tcnt, 0);
		p2->tx = 0;
		p2->ty = frg2->tcnt - 1;
		if(p1){
			f = p1->lnks + p1->dir;
			l = ref_lnkv(g->lnks, f->idx);
			if(f->flg){
				if(p1->dir){
					p1->tx = l->tidx2;
				} else {
					p1->ty = frg1->tcnt - 1 - l->tidx2;
				}
				if(p2->dir){
					p2->ty = frg2->tcnt - 1 - l->tidx1;
				} else {
					p2->tx = l->tidx1;
				}
			} else {
				if(p1->dir){
					p1->tx = l->tidx1;
				} else {
					p1->ty = frg1->tcnt - 1 - l->tidx1;
				}
				if(p2->dir){
					p2->ty = frg2->tcnt - 1 - l->tidx2;
				} else {
					p2->tx = l->tidx2;
				}
			}
			if(p1->ty >= frg1->tcnt || p2->ty >= frg2->tcnt){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
			off += l->off;
		}
		p1 = p2;
		frg1 = frg2;
	}
	p1   = NULL;
	frg1 = NULL;
	for(i=0;i<path->size;i++){
		p2 = ref_pathv(path, i);
		frg2 = ref_frgv(g->frgs, p2->frg);
		if(p1){
			t1 = ref_tracev(g->traces, frg1->toff + (p1->dir? p1->tx : p1->ty));
			t2 = ref_tracev(g->traces, frg2->toff + (p2->dir? p2->ty : p2->tx));
			off = p1->off + (p1->dir? (int)(frg1->len - (t1->off + g->reglen)) : t1->off);
			f = p1->lnks + p1->dir;
			l = ref_lnkv(g->lnks, f->idx);
			len = l->off + (p1->dir? t1->off + g->reglen : frg1->len - t1->off) + (p2->dir? frg2->len - t2->off : t2->off + g->reglen);
			push_seqletv(qs, (seqlet_t){t1->node, p1->dir ^ t1->dir, t2->node, p2->dir ^ t2->dir, off, len});
		}
		if(p2->dir){
			for(j=p2->tx+1;j<=p2->ty;j++){
				t1 = ref_tracev(g->traces, frg2->toff + frg2->tcnt - 1 - j);
				t2 = ref_tracev(g->traces, frg2->toff + frg2->tcnt - 0 - j);
				off = p2->off + frg2->len - (t2->off + g->reglen);
				len = t2->off - t1->off + g->reglen;
				push_seqletv(qs, (seqlet_t){t2->node, !t2->dir, t1->node, !t1->dir, off, len});
			}
		} else {
			for(j=p2->tx+1;j<=p2->ty;j++){
				t1 = ref_tracev(g->traces, frg2->toff + j - 1);
				t2 = ref_tracev(g->traces, frg2->toff + j);
				off = p2->off + t1->off;
				len = t2->off - t1->off + g->reglen;
				push_seqletv(qs, (seqlet_t){t1->node, t1->dir, t2->node, t2->dir, off, len});
			}
		}
		p1 = p2;
		frg1 = frg2;
	}
	return qs;
}

uint64_t gen_contigs_graph(Graph *g, FILE *out){
	pathv *path;
	seqletv *qs;
	seqlet_t *q;
	path_t *t;
	frg_t *n;
	uint64_t nid, nctg, i;
	for(i=0;i<g->ctgs->size;i++) free_tracev(g->ctgs->buffer[i]);
	clear_vplist(g->ctgs);
	nctg = 0;
	for(nid=0;nid<g->frgs->size;nid++) g->frgs->buffer[nid].bt_visit = 0;
	path = init_pathv(4);
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		if(n->closed) continue;
		if(n->bt_visit) continue;
		nctg ++;
		clear_pathv(path);
		t = next_ref_pathv(path);
		t->frg = nid;
		t->lnks[0] = EDGE_REF_NULL;
		t->lnks[1] = EDGE_REF_NULL;
		t->dir = 0;
		true_linear_unique_path_graph(g, path, 0xFFFFFFFFFFFFFFFFLLU, nctg, NULL);
		reverse_pathv(path);
		for(i=0;i<path->size;i++) path->buffer[i].dir = !path->buffer[i].dir;
		true_linear_unique_path_graph(g, path, 0xFFFFFFFFFFFFFFFFLLU, nctg, NULL);
		qs = path2seqlets_graph(g, path);
		if(qs->size + 1 < (u4i)g->min_ctg_nds){
			free_seqletv(qs);
			continue;
		}
		q = ref_seqletv(qs, qs->size - 1);
		if((int)q->off + (int)q->len < (int)g->min_ctg_len){
			free_seqletv(qs);
			continue;
		}
		if(out){
			for(i=0;i<path->size;i++){
				t = ref_pathv(path, i);
				fprintf(out, "ctg%d\tF%d\t%c\t%d\n", (int)g->ctgs->size, t->frg, "+-*@"[t->dir], t->off);
			}
		}
		push_vplist(g->ctgs, qs);
	}
	free_pathv(path);
	g->major_nctg = g->ctgs->size;
	return nctg;
}

u8i gen_complex_contigs_graph(Graph *g){
	tracev *ts;
	trace_t *t;
	seqletv *qs;
	seqlet_t *q;
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	u8i i, idx, mi, cnt;
	u4i j, k, mk, mc;
	int mf;
	for(i=0;i<g->nodes->size;i++) g->nodes->buffer[i].unvisit = 1;
	cnt = 0;
	for(i=0;i<g->utgs->size;i++){
		ts = (tracev*)get_vplist(g->utgs, i);
		for(j=0;j<ts->size;j++){
			t = ref_tracev(ts, j);
			n = ref_nodev(g->nodes, t->node);
			n->unvisit = 0;
		}
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->unvisit == 0) continue;
		if(n->regs.cnt < g->min_node_cov){
			n->unvisit = 0;
			continue;
		}
	}
	cnt = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->unvisit == 0) continue;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			mi = mk = mc = 0; mf = MAX_VALUE_B4;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				if(e->node1 == e->node2) continue;
				if(e->cov < g->min_edge_cov) continue;
				if(e->off < mf){
					mi = f - g->erefs->buffer; mk = k; mc = e->cov; mf = e->off;
				}
			}
			if(mf == MAX_VALUE_B4) continue;
			f = ref_edgerefv(g->erefs, mi);
			if(f->flg) continue;
			e = ref_edgev(g->edges, f->idx);
			if(g->nodes->buffer[e->node2].unvisit == 0) continue;
			qs = init_seqletv(1);
			q = next_ref_seqletv(qs);
			q->node1 = i;
			q->dir1  = k;
			q->node2 = e->node2;
			q->dir2  = e->dir2;
			q->off   = 0;
			q->len   = g->reglen * 2 + e->off;
			push_vplist(g->ctgs, qs);
			cnt ++;
		}
	}
	return cnt;
}

void n50_stat_contigs_graph(Graph *g){
	seqletv *qs;
	seqlet_t *q;
	u4v *lens;
	int len;
	u8i i;
	lens = init_u4v(g->major_nctg + 1);
	for(i=0;i<g->major_nctg;i++){
		qs = (seqletv*)get_vplist(g->ctgs, i);
		q = ref_seqletv(qs, qs->size - 1);
		len = (int)q->off + (int)q->len;
		if(len <= 0){
			fprintf(stderr, " -- seqlet[ctg%llu off=%d, len=%d, sum=%d]  in %s -- %s:%d --\n", i, (int)q->off, (int)q->len, len, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			continue;
		}
		push_u4v(lens, len);
	}
	fprintf(KBM_LOGF, "[%s] Estimated: ", date()); num_n50(lens, KBM_LOGF); fprintf(KBM_LOGF, "\n");
	free_u4v(lens);
}

// after gen_contigs
u4i count_isolated_reads_graph(Graph *g){
	seqletv *ts;
	seqlet_t *t;
	node_t *n;
	read_t *rd;
	reg_t *r;
	u8i i, j, cnt, idx;
	int fnd;
	for(i=0;i<g->nodes->size;i++) g->nodes->buffer[i].unvisit = 1;
	cnt = 0;
	for(i=0;i<g->ctgs->size;i++){
		ts = (seqletv*)get_vplist(g->ctgs, i);
		for(j=0;j<ts->size;j++){
			t = ref_seqletv(ts, j);
			n = ref_nodev(g->nodes, t->node1);
			n->unvisit = 0;
			n = ref_nodev(g->nodes, t->node2);
			n->unvisit = 0;
		}
	}
	for(i=0;i<g->reads->size;i++){
		rd = ref_readv(g->reads, i);
		fnd = 0;
		idx = rd->regs.idx;
		while(idx){
			r = ref_regv(g->regs, idx);
			idx = r->read_link;
			n = ref_nodev(g->nodes, r->node);
			if(n->unvisit == 0){ fnd = 1; break; }
		}
		if(fnd == 0) cnt ++;
	}
	return cnt;
}

void print_dot_subgraph(Graph *g, subnodehash *nodes, subedgev *edges, FILE *out){
	subnode_t *n1, *n2;
	subedge_t *e;
	u4i k, idx;
	fprintf(out, "digraph {\n");
	fprintf(out, " rankdir=LR\n");
	reset_iter_subnodehash(nodes);
	while((n1 = ref_iter_subnodehash(nodes))){
		if(n1->closed) continue;
		fprintf(out, "N%llu [label=\"N%llu(%llu) %d:%d:%d\" style=filled fillcolor=\"%s\" color=\"%s\"]\n", (u8i)n1->node, (u8i)n1->node, (u8i)g->nodes->buffer[n1->node].rep_idx, n1->cov, n1->visit, n1->bt_open, n1->fix? "yellow" : "white", n1->visit? "green" : (n1->cov > 2? "blue" : "black"));
	}
	reset_iter_subnodehash(nodes);
	while((n1 = ref_iter_subnodehash(nodes))){
		if(n1->closed) continue;
		for(k=0;k<2;k++){
			idx = n1->edges[k].idx;
			while(idx){
				e = ref_subedgev(edges, idx);
				idx = e->next;
				if(e->fwd == 0) continue;
				if(e->closed) continue;
				n2 = e->node;
				fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d:%d\" color=\"%s\" %s]\n", (u8i)n1->node, (u8i)n2->node, "+-"[k], "+-"[e->dir], e->cov, e->visit, e->cov > 1? "blue" : "black", e->visit? "style=dashed":"");
			}
		}
	}
	fprintf(out, "}\n");
	fflush(out);
}

void fprintf_dot_subgraph(Graph *g, subnodehash *nodes, subedgev *edges, char *name_prefix, char *name_suffix){
	FILE *out;
	out = open_file_for_write(name_prefix, name_suffix, 1);
	print_dot_subgraph(g, nodes, edges, out);
	fclose(out);
}

typedef struct {
	u4i node:31, dir:1;
	u4i flag;
	u4i prev;
	u4i step;
	int score;
} sg_heap_t;
define_list(sgheapv, sg_heap_t);

typedef struct {
	u4i node:31, dir:1;
	u4i group:30, solid:1, closed:1;
} sg_tip_t;
define_list(sgtipv, sg_tip_t);

subedge_t* find_edge_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
	subnode_t *n;
	subedge_t *e;
	u8i idx;
	n = ref_subnodehash(nodes, node1);
	idx = n->edges[dir1].idx;
	e = ref_subedgev(edges, idx);
	if(offset_subnodehash(nodes, e->node) == node2 && e->dir == dir2){
		return e;
	}
	while((idx = e->next)){
		e = ref_subedgev(edges, idx);
		if(offset_subnodehash(nodes, e->node) == node2 && e->dir == dir2){
			return e;
		}
	}
	return NULL;
}

int cut_edge_core_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
	subnode_t *n;
	subedge_t *e, *p;
	u8i idx;
	n = ref_subnodehash(nodes, node1);
	idx = n->edges[dir1].idx;
	e = ref_subedgev(edges, idx);
	if(offset_subnodehash(nodes, e->node) == node2 && e->dir == dir2){
		e->closed = 1;
		n->edges[dir1].idx = e->next;
		n->edges[dir1].cnt --;
		return 1;
	}
	while((idx = e->next)){
		p = e;
		e = ref_subedgev(edges, idx);
		if(offset_subnodehash(nodes, e->node) == node2 && e->dir == dir2){
			e->closed = 1;
			p->next = e->next;
			n->edges[dir1].cnt --;
			return 1;
		}
	}
	return 0;
}

int cut_edge_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
	return cut_edge_core_subgraph(nodes, edges, node1, dir1, node2, dir2)
		+ cut_edge_core_subgraph(nodes, edges, node2, !dir2, node1, !dir1);
}

u4i unitigs2frgs_graph(Graph *g, int ncpu){
	frg_t *frg;
	node_t *n;
	tracev *ts;
	trace_t *t;
	u4i i, j, tid, ret;
	{
		clear_frgv(g->frgs);
		clear_lnkv(g->lnks);
		clear_edgerefv(g->lrefs);
		clear_tracev(g->traces);
	}
	ret = 0;
	for(tid=0;tid<g->utgs->size;tid++){
		ret ++;
		frg = next_ref_frgv(g->frgs);
		frg->toff = g->traces->size;
		frg->lnks[0] = PTR_REF_NULL;
		frg->lnks[1] = PTR_REF_NULL;
		frg->closed = 0;
		ts = (tracev*)get_vplist(g->utgs, tid);
		frg->tx = 0;
		frg->ty = ts->size;
		for(i=frg->tx;i<frg->ty;i++) ts->buffer[i].cov = 0;
		append_array_tracev(g->traces, ts->buffer + frg->tx, frg->ty - frg->tx);
		frg->tcnt = frg->ty - frg->tx;
		frg->length = frg->len  = cal_offset_traces_graph(g, g->traces, frg->toff + frg->tx, frg->toff + frg->ty, 0);
	}
	psort_array(g->frgs->buffer, g->frgs->size, frg_t, ncpu, num_cmpgt(b.length, a.length));
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		n->unvisit = 1;
		n->rep_idx = MAX_REP_IDX;
		n->bt_visit = 0;
	}
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		for(j=frg->tx;j<frg->ty;j++){
			t = ref_tracev(g->traces, frg->toff + j);
			n = ref_nodev(g->nodes, t->node);
			n->rep_idx = i;
			n->rep_dir = t->dir;
			n->bt_visit = j;
		}
	}
	return ret;
}

int scan_rd_lnk_core(Graph *g, u4i rid, lnk_t *lnk, u8v *regids){
	read_t *rd;
	reg_t  *r, *rl, *rr;
	node_t *n, *nl, *nr;
	frg_t  *fl, *fr;
	trace_t *tl, *tr;
	u8i idx;
	u4i i, tmp;
	rd = ref_readv(g->reads, rid);
	idx = rd->regs.idx;
	clear_u8v(regids);
	while(idx){
		push_u8v(regids, idx);
		idx = ref_regv(g->regs, idx)->read_link;
	}
	if(regids->size < 2) return 0;
	rl = NULL; nl = NULL;
	for(i=0;i<regids->size;i++){
		r = ref_regv(g->regs, regids->buffer[i]);
		n = ref_nodev(g->nodes, r->node);
		if(n->rep_idx == MAX_REP_IDX) continue;
		if(rl && n->rep_idx != nl->rep_idx) break;
		rl = r;
		nl = n;
	}
	if(i == regids->size) return 0;
	rr = NULL; nr = NULL;
	for(i=regids->size;i>0;i--){
		r = ref_regv(g->regs, regids->buffer[i - 1]);
		n = ref_nodev(g->nodes, r->node);
		if(n->rep_idx == MAX_REP_IDX) continue;
		if(rr && n->rep_idx != nr->rep_idx) break;
		rr = r;
		nr = n;
	}
	fl = ref_frgv(g->frgs, nl->rep_idx);
	tl = ref_tracev(g->traces, fl->toff + nl->bt_visit);
	fr = ref_frgv(g->frgs, nr->rep_idx);
	tr = ref_tracev(g->traces, fr->toff + nr->bt_visit);
	lnk->frg1 = nl->rep_idx;
	lnk->frg2 = nr->rep_idx;
	if(lnk->frg1 == lnk->frg2) return 0;
	lnk->dir1 = rl->dir ^ nl->rep_dir;
	lnk->dir2 = rr->dir ^ nr->rep_dir;
	lnk->tidx1 = lnk->dir1? nl->bt_visit : fl->tcnt - 1 - nl->bt_visit;
	lnk->tidx2 = lnk->dir2? fr->tcnt - 1 - nr->bt_visit : nr->bt_visit;
	lnk->cov = 1; // directed link
	//ln->weak = 0;
	//lnk->closed = 0;
	lnk->off = rr->beg - rl->end;
	lnk->off -= lnk->dir1? tl->off : ((b4i)fl->len) - (tl->off + rl->end - rl->beg);
	lnk->off -= lnk->dir2? ((b4i)fr->len) - (tr->off + rr->end - rr->beg) : tr->off;
	if(lnk->off + (int)g->reglen < 0) return 0;
	if(lnk->frg1 > lnk->frg2){
		swap_tmp(lnk->frg1, lnk->frg2, tmp);
		swap_tmp(lnk->dir1, lnk->dir2, tmp);
		lnk->dir1 = !lnk->dir1;
		lnk->dir2 = !lnk->dir2;
		swap_tmp(lnk->tidx1, lnk->tidx2, tmp);
	}
	return 1;
}

int scan_nd_lnk_core(Graph *g, u8i nid, lnk_t *lnk){
	node_t *n, *w;
	edge_ref_t *f;
	edge_t *e;
	frg_t *fl, *fr;
	trace_t *tl, *tr;
	u8i idx, wid;
	u4i fids[2], dirs[2], dir, covs[2], tidx[2], k;
	int offs[2];
	n = ref_nodev(g->nodes, nid);
	if(n->rep_idx != MAX_REP_IDX) return 0;
	for(k=0;k<2;k++){
		fids[k] = 0;
		dirs[k] = 0;
		covs[k] = 0;
		tidx[k] = 0;
		offs[k] = 0;
		idx = n->edges[k].idx;
		while(idx){
			f = ref_edgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_edgev(g->edges, f->idx);
			wid = f->flg? e->node1 : e->node2;
			w = ref_nodev(g->nodes, wid);
			if(w->rep_idx == MAX_REP_IDX) continue;
			dir = f->flg? !e->dir1 : e->dir2;
			dir = w->rep_dir ^ dir;
			if(fids[k] == 0){
				fids[k] = w->rep_idx;
				dirs[k] = dir;
				covs[k] = 1;
				tidx[k] = w->bt_visit; // bt_visit is used to store nodes's tidx on frg
				offs[k] = e->off;
			} else if(fids[k] == w->rep_idx && dirs[k] == dir){
				covs[k] ++;
			} else {
				fids[k] = MAX_U4;
				break;
			}
		}
		if(fids[k] == 0) return 0;
		if(fids[k] == MAX_U4) return 0;
		if(covs[k] < g->min_edge_cov) return 0;
	}
	if(fids[0] == fids[1]) return 0;
	lnk->cov = 0; // indirected link
	if(fids[0] < fids[1]){
		lnk->frg1 = fids[0];
		lnk->frg2 = fids[1];
		lnk->dir1 = !dirs[0];
		lnk->dir2 = dirs[1];
		lnk->tidx1 = tidx[0];
		lnk->tidx2 = tidx[1];
	} else {
		lnk->frg1 = fids[1];
		lnk->frg2 = fids[0];
		lnk->dir1 = !dirs[1];
		lnk->dir2 = dirs[0];
		lnk->tidx1 = tidx[1];
		lnk->tidx2 = tidx[0];
	}
	lnk->off = offs[0] + offs[1] + g->reglen;
	fl = ref_frgv(g->frgs, lnk->frg1);
	fr = ref_frgv(g->frgs, lnk->frg2);
	tl = ref_tracev(g->traces, fl->toff + lnk->tidx1);
	tr = ref_tracev(g->traces, fr->toff + lnk->tidx2);
	lnk->off -= lnk->dir1? tl->off : ((int)fl->len) - (int)(tl->off + g->reglen);
	lnk->off -= lnk->dir2? ((int)fr->len) - (int)(tr->off + g->reglen) : tr->off;
	if(lnk->off + (int)g->reglen < 0) return 0;
	return 1;
}

thread_beg_def(mlnk);
Graph *g;
lnkv *lnks;
int task;
thread_end_def(mlnk);

thread_beg_func(mlnk);
u8v *regids;
u8i nid;
u4i rid;
lnk_t LNK;
memset(&LNK, 0, sizeof(lnk_t));
LNK.cov = 1;
LNK.weak = 0;
LNK.closed = 0;
regids = init_u8v(32);
thread_beg_loop(mlnk);
if(mlnk->task == 1){
	for(rid=mlnk->t_idx;rid<mlnk->g->reads->size;rid+=mlnk->n_cpu){
		if(scan_rd_lnk_core(mlnk->g, rid, &LNK, regids)){
			push_lnkv(mlnk->lnks, LNK);
		}
	}
} else if(mlnk->task == 2){
	for(nid=mlnk->t_idx;nid<mlnk->g->nodes->size;nid+=mlnk->n_cpu){
		if(scan_nd_lnk_core(mlnk->g, nid, &LNK)){
			push_lnkv(mlnk->lnks, LNK);
		}
	}
}
thread_end_loop(mlnk);
free_u8v(regids);
thread_end_func(mlnk);

u4i gen_lnks_graph(Graph *g, int ncpu, FILE *log){
	frg_t *frg;
	trace_t *t;
	node_t *n;
	lnkv *lnks;
	lnk_t *l;
	edge_ref_t F1, F2;
	u8i lst, idx;
	u4i i, j, m, v1, v2, cnt, cov;
	int x;
	thread_preprocess(mlnk);
	clear_lnkv(g->lnks);
	memset(next_ref_lnkv(g->lnks), 0, sizeof(lnk_t));
	clear_edgerefv(g->lrefs);
	push_edgerefv(g->lrefs, EDGE_REF_NULL);
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		n->unvisit = 1;
		n->rep_idx = MAX_REP_IDX;
		n->bt_visit = 0;
	}
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		for(j=frg->tx;j<frg->ty;j++){
			t = ref_tracev(g->traces, frg->toff + j);
			n = ref_nodev(g->nodes, t->node);
			n->rep_idx = i;
			n->rep_dir = t->dir;
			n->bt_visit = j;
		}
	}
	thread_beg_init(mlnk, ncpu);
	mlnk->g = g;
	mlnk->lnks = init_lnkv(1024);
	mlnk->task = 0;
	thread_end_init(mlnk);
	thread_apply_all(mlnk, mlnk->task = 1);
	//thread_apply_all(mlnk, mlnk->task = 2);
	lnks = init_lnkv(1024);
	thread_beg_close(mlnk);
	append_lnkv(lnks, mlnk->lnks);
	free_lnkv(mlnk->lnks);
	thread_end_close(mlnk);
	psort_array(lnks->buffer, lnks->size, lnk_t, ncpu, num_cmpgtx(a.key, b.key, a.off, b.off));
	if(log){
		for(i=0;i<lnks->size;i++){
			l = ref_lnkv(lnks, i);
			fprintf(log, "F%d[%c:%d] -> F%d[%c:%d] = %d, cov=%d\n", l->frg1, "+-"[l->dir1], l->tidx1, l->frg2, "+-"[l->dir2], l->tidx2, l->off, l->cov);
		}
	}
	cov = 0;
	for(i=j=0;i<=lnks->size;i++){
		if(i == lnks->size || lnks->buffer[i].key != lnks->buffer[j].key){
			push_edgerefv(g->lrefs, (edge_ref_t){g->lnks->size, 0, 0});
			push_edgerefv(g->lrefs, (edge_ref_t){g->lnks->size, 1, 0});
			l = next_ref_lnkv(g->lnks);
			m = (((u8i)i) + j) / 2;
			if(cov && lnks->buffer[m].cov == 0){
				for(v1=1;v1+j<=m;v1++){
					if(lnks->buffer[m-v1].cov){
						v1 |= 0x80000000U;
						break;
					}
				}
				for(v2=1;m+v2<i;v2++){
					if(lnks->buffer[m+v2].cov){
						v2 |= 0x80000000U;
						break;
					}
				}
				if(v1 & 0x80000000U){
					if(v2 & 0x80000000U){
						if((v1 & 0x7FFFFFFFU) <= (v2 & 0x7FFFFFFFU)){
							m = m - (v1 & 0x7FFFFFFFU);
						} else {
							m = m + (v2 & 0x7FFFFFFFU);
						}
					} else {
						m = m - (v1 & 0x7FFFFFFFU);
					}
				} else {
					if(v2 & 0x80000000U){
						m = m + (v2 & 0x7FFFFFFFU);
					} else {
						fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					}
				}
			}
			*l = lnks->buffer[m];
			l->cov  = cov;
			l->weak = (l->cov < g->max_node_cov_sg);
			j = i;
			cov = lnks->buffer[i].cov;
		} else {
			cov += lnks->buffer[i].cov;
		}
	}
	free_lnkv(lnks);

	// sort lrefs
	psort_array(g->lrefs->buffer + 1, g->lrefs->size - 1, edge_ref_t, ncpu, num_cmpgt(
		(a.flg? ((g->lnks->buffer[a.idx].frg2 << 1) | !g->lnks->buffer[a.idx].dir2) : ((g->lnks->buffer[a.idx].frg1 << 1) | g->lnks->buffer[a.idx].dir1)),
		(b.flg? ((g->lnks->buffer[b.idx].frg2 << 1) | !g->lnks->buffer[b.idx].dir2) : ((g->lnks->buffer[b.idx].frg1 << 1) | g->lnks->buffer[b.idx].dir1))
		));
	push_edgerefv(g->lrefs, (edge_ref_t){g->lnks->size, 0, 0}); memset(next_ref_lnkv(g->lnks), 0, sizeof(lnk_t)); g->lnks->size --;
	g->lrefs->size --;
	F1.idx = g->lnks->size; F1.flg = 0;
	cnt = 0;
	// update frg->lnks
	for(lst=idx=1;idx<=g->lrefs->size;idx++){
		if(g->lrefs->buffer[idx].flg){
			F2.idx =  g->lnks->buffer[g->lrefs->buffer[idx].idx].frg2;
			F2.flg = !g->lnks->buffer[g->lrefs->buffer[idx].idx].dir2;
		} else {
			F2.idx =  g->lnks->buffer[g->lrefs->buffer[idx].idx].frg1;
			F2.flg =  g->lnks->buffer[g->lrefs->buffer[idx].idx].dir1;
		}
		if(F1.idx == F2.idx && F1.flg == F2.flg) continue;
		if(lst < idx){
			frg = ref_frgv(g->frgs, F1.idx);
			frg->lnks[F1.flg].idx = lst;
			for(x=lst;x+1<(int)idx;x++){
				g->lrefs->buffer[x].next = x + 1;
				if(g->lnks->buffer[g->lrefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) frg->lnks[F1.flg].cnt ++;
			}
			if(g->lnks->buffer[g->lrefs->buffer[x].idx].closed == WT_EDGE_CLOSED_NULL) frg->lnks[F1.flg].cnt ++;
		}
		lst = idx;
		F1 = F2;
	}
	return g->lnks->size - 1;
}

int gen_seq_traces_graph(Graph *g, tracev *path, String *seq){
	trace_t *t1, *t2;
	reg_t *reg, *r1, *r2;
	edge_t *e;
	uint32_t i;
	int inc, found;
	clear_string(seq);
	t1 = NULL;
	for(i=0;i<path->size;i++){
		t2 = ref_tracev(path, i);
		if(t1){
			inc = 0;
			r1 = ref_regv(g->regs, ref_nodev(g->nodes, t1->node)->regs.idx);
			r2 = ref_regv(g->regs, ref_nodev(g->nodes, t2->node)->regs.idx);
			e = ref_edgev(g->edges, t1->edges[t1->dir].idx);
			do {
				inc = 0;
				found = 0;
				while(r1->node == t1->node && r2->node == t2->node){
					if(r1->rid > r2->rid){
						r2 ++;
					} else if(r1->rid < r2->rid){
						r1 ++;
					} else {
						if(r1->beg < r2->beg){
							if(t1->dir ^ r1->dir){ r1++; r2++; continue; }
							inc = r2->beg - r1->end;
							if(inc <= 0) break;
							encap_string(seq, inc);
							seq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[r1->rid].rdoff + r1->end, inc, seq->string + seq->size);
							seq->size += inc;
						} else {
							if(!(t1->dir ^ r1->dir)){ r1++; r2++; continue; }
							inc = r1->beg - r2->end;
							if(inc <= 0) break;
							encap_string(seq, inc);
							revseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[r1->rid].rdoff + r2->end, inc, seq->string + seq->size);
							seq->size += inc;
						}
						inc = 0;
						found = 1; break;
					}
				}
				if(found == 0){ inc = e->off; break; }
			} while(0);
			if(inc > 0){ inc = 0; while(inc++ < e->off) add_char_string(seq, 'N'); }
			else if(inc < 0){
				if(seq->size + inc < 0) seq->size = 0;
				else seq->size += inc;
				seq->string[seq->size] = '\0';
			}
		}
		t1 = t2;
		reg = ref_regv(g->regs, ref_nodev(g->nodes, t1->node)->regs.idx);
		inc = reg->end - reg->beg;
		encap_string(seq, inc);
		if(t1->dir ^ reg->dir) revseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[reg->rid].rdoff + reg->beg, inc, seq->string + seq->size);
		else                   seq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[reg->rid].rdoff + reg->beg, inc, seq->string + seq->size);
		seq->size += inc;
	}
	return seq->size;
}

typedef struct {
	uint64_t rid:26, dir:1, beg:18, end:18, view:1;
} lay_reg_t;
define_list(layregv, lay_reg_t);

typedef struct {
	u4i tidx;
	u8i roff:48, rcnt:16;
} lay_t;
define_list(layv, lay_t);

void gen_lay_regs_core_graph(Graph *g, seqlet_t *q, layregv *regs){
	node_t *n1, *n2;
	reg_t *r1, *r2;
	uint32_t rid, beg, end;
	n1 = ref_nodev(g->nodes, q->node1);
	n2 = ref_nodev(g->nodes, q->node2);
	r1 = ref_regv(g->regs, n1->regs.idx);
	r2 = ref_regv(g->regs, n2->regs.idx);
	while(r1->node == q->node1 && r2->node == q->node2){
		if(r1->rid > r2->rid){
			r2 ++;
		} else if(r1->rid < r2->rid){
			r1 ++;
		} else {
			rid = r1->rid;
			if(r1->beg < r2->beg){
				if(q->dir1 ^ r1->dir){ r1 ++; r2 ++; continue; }
				beg = r1->beg; end = r2->end;
				push_layregv(regs, (lay_reg_t){rid, 0, beg, end, 0});
			} else {
				if(!(q->dir1 ^ r1->dir)){ r1 ++; r2 ++; continue; }
				beg = r2->beg; end = r1->end;
				push_layregv(regs, (lay_reg_t){rid, 1, beg, end, 0});
			}
			r1 ++; r2 ++;
		}
	}
	//if(regs->size == 0){
		//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	//}
}

thread_beg_def(mlay);
Graph *g;
seqletv *path;
layv    *lays;
layregv *regs;
FILE *log;
thread_end_def(mlay);

thread_beg_func(mlay);
seqlet_t *let;
lay_t *lay;
u8i i;
thread_beg_loop(mlay);
clear_layregv(mlay->regs);
for(i=mlay->t_idx;i<mlay->path->size;i+=mlay->n_cpu){
	let = ref_seqletv(mlay->path, i);
	lay = ref_layv(mlay->lays, i);
	lay->tidx = mlay->t_idx;
	lay->roff = mlay->regs->size;
	gen_lay_regs_core_graph(mlay->g, let, mlay->regs);
	lay->rcnt = mlay->regs->size - lay->roff;
	sort_array(mlay->regs->buffer + lay->roff, lay->rcnt, lay_reg_t, num_cmpgt(b.end - b.beg, a.end - a.beg));
	if(lay->rcnt == 0 && mlay->log){
		thread_beg_syn(mlay);
		fprintf(mlay->log, " -- N%llu(%c) -> N%llu(%c) has no read path --\n", (u8i)let->node1, "+-"[let->dir1], (u8i)let->node2, "+-"[let->dir2]); fflush(mlay->log);
		thread_end_syn(mlay);
	}
}
thread_end_loop(mlay);
thread_end_func(mlay);

u8i print_ctgs_graph(Graph *g, u8i uid, u8i beg, u8i end, char *prefix, char *lay_suffix, u4i ncpu, FILE *log){
	FILE *o_lay;
	BufferedWriter *bw;
	layv *lays;
	layregv *regs;
	seqletv *path;
	seqlet_t *t;
	lay_t *lay;
	lay_reg_t *reg;
	u8i i, ret;
	u4i j, c, len;
	thread_preprocess(mlay);
	o_lay = open_file_for_write(prefix, lay_suffix, 1);
	bw = open_bufferedwriter(o_lay, 8 * 1024 * 1024);
	lays = init_layv(32);
	thread_beg_init(mlay, ncpu);
	mlay->g = g;
	mlay->path = NULL;
	mlay->lays = lays;
	mlay->regs = init_layregv(32);
	mlay->log  = log;
	thread_end_init(mlay);
	ret = 0;
	for(i=beg;i<end;i++){
		path = (seqletv*)get_vplist(g->ctgs, i);
		if(path->size == 0) continue;
		clear_and_inc_layv(lays, path->size);
		thread_apply_all(mlay, EXPR(mlay->path = path));
		uid ++;
		ret ++;
		len = path->buffer[path->size - 1].off + path->buffer[path->size - 1].len;
		{
			beg_bufferedwriter(bw);
			fprintf(bw->out, ">ctg%llu nodes=%llu len=%u\n", uid, (u8i)path->size + 1, len);
			if(log) fprintf(log, "OUTPUT_CTG\tctg%d -> ctg%d nodes=%llu len=%u\n", (int)i, (int)uid, (u8i)path->size + 1, len);
			for(j=0;j<lays->size;j++){
				if((j % 100) == 0){
					end_bufferedwriter(bw);
					beg_bufferedwriter(bw);
				}
				lay = ref_layv(lays, j);
				if(lay->rcnt == 0) continue;
				t = ref_seqletv(path, j);
				fprintf(bw->out, "E\t%d\tN%llu\t%c\tN%llu\t%c\n", (int)t->off, (u8i)t->node1, "+-"[t->dir1], (u8i)t->node2, "+-"[t->dir2]);
				regs = thread_access(mlay, (j % mlay->n_cpu))->regs;
				for(c=0;c<lay->rcnt;c++){
					reg = ref_layregv(regs, lay->roff + c);
					fprintf(bw->out, "%c\t%s\t%c\t%d\t%d\t", "Ss"[reg->view], g->kbm->reads->buffer[reg->rid].tag, "+-"[reg->dir], reg->beg, reg->end - reg->beg);
					if(reg->dir){
						print_revseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, bw->out);
					} else {
						print_seq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[reg->rid].rdoff + reg->beg, reg->end - reg->beg, bw->out);
					}
					fprintf(bw->out, "\n");
				}
			}
			end_bufferedwriter(bw);
		}
	}
	close_bufferedwriter(bw);
	fclose(o_lay);
	thread_beg_close(mlay);
	free_layregv(mlay->regs);
	thread_end_close(mlay);
	fprintf(KBM_LOGF, "[%s] output %u contigs\n", date(), (u4i)ret);
	free_layv(lays);
	return uid;
}

uint32_t print_traces_graph(Graph *g, tracev *path, FILE *out){
	String *str;
	trace_t *t1, *t2;
	node_t *n1, *n2;
	reg_t *r1, *r2;
	edge_ref_t *f;
	edge_t *e;
	int offset, fst;
	uint64_t beg, end;
	uint32_t i, rid;
	if(path->size < 2) return 0;
	str = init_string(1024);
	offset = 0;
	t1 = ref_tracev(path, 0);
	for(i=1;i<path->size;i++){
		t2 = ref_tracev(path, i);
		{
			n1 = ref_nodev(g->nodes, t1->node);
			n2 = ref_nodev(g->nodes, t2->node);
			f  = t1->edges + t1->dir;
			e  = ref_edgev(g->edges, f->idx);
			r1 = ref_regv(g->regs, n1->regs.idx);
			r2 = ref_regv(g->regs, n2->regs.idx);
			fprintf(out, "E\t%d\tN%llu\t%c\tN%llu\t%c\n", offset, t1->node, "+-"[t1->dir], t2->node, "+-"[t2->dir]);
			fst = 1;
			while(r1->node == t1->node && r2->node == t2->node){
				if(r1->rid > r2->rid){
					r2 ++;
				} else if(r1->rid < r2->rid){
					r1 ++;
				} else {
					rid = r1->rid;
					if(r1->beg < r2->beg){
						if(t1->dir ^ r1->dir){ r1 ++; r2 ++; continue; }
						beg = r1->beg; end = r2->end;
						fprintf(out, "S\t%s\t", g->kbm->reads->buffer[rid].tag);
						fprintf(out, "+\t%d\t%d\t", (int)beg, (int)(end - beg));
						encap_string(str, end - beg);
						fwdseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[rid].rdoff + beg, end - beg, str->string);
						fputs(str->string, out);
						//beg += g->kbm->reads->buffer[rid].rdoff;
						//end += g->kbm->reads->buffer[rid].rdoff;
						//for(j=beg;j<end;j++){
							//fputc(bit_base_table[bits2bit(g->kbm->rdseqs->bits, j)], out);
						//}
					} else {
						if(!(t1->dir ^ r1->dir)){ r1 ++; r2 ++; continue; }
						beg = r2->beg; end = r1->end;
						fprintf(out, "S\t%s\t", g->kbm->reads->buffer[rid].tag);
						fprintf(out, "-\t%d\t%d\t", (int)beg, (int)(end - beg));
						encap_string(str, end - beg);
						revseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[rid].rdoff + beg, end - beg, str->string);
						fputs(str->string, out);
						//beg += g->kbm->reads->buffer[rid].rdoff;
						//end += g->kbm->reads->buffer[rid].rdoff;
						//for(j=end;j>beg;j--){
							//fputc(bit_base_table[bits2revbit(g->kbm->rdseqs->bits, (j-1))], out);
						//}
					}
					fputc('\n', out);
					if(fst){
						offset += end - beg; fst = 0;
					}
					r1 ++; r2 ++;
				}
			}
		}
		t1 = t2;
	}
	free_string(str);
	return offset;
}

uint64_t print_utgs_graph(Graph *g, char *prefix, char *utg, char *lay){
	FILE *o_seq, *o_lay, *files[4];
	tracev *path;
	String *seq;
	char *str;
	uint64_t i, uid, cnt, tot;
	char ch;
	int beg, end;
	files[0] = open_file_for_write(prefix, utg, 1);
	str = catstr(2, utg, ".filtered");
	files[1] = open_file_for_write(prefix, str, 1);
	free(str);
	files[2] = open_file_for_write(prefix, lay, 1);
	str = catstr(2, lay, ".filtered");
	files[3] = open_file_for_write(prefix, str, 1);
	free(str);
	seq = init_string(1024);
	tot = cnt = 0;
	for(i=uid=0;i<g->utgs->size;i++){
		path = (tracev*)get_vplist(g->utgs, i);
		if(gen_seq_traces_graph(g, path, seq) < g->min_ctg_len || (int)path->size < g->min_ctg_nds){
			o_seq = files[1];
			o_lay = files[3];
		} else {
			o_seq = files[0];
			o_lay = files[2];
			cnt ++;
			tot += seq->size;
		}
		uid ++;
		fprintf(o_seq, ">utg%llu len=%d nodes=%llu beg=N%llu end=N%llu\n", (unsigned long long)uid, seq->size, (unsigned long long)path->size,
			path->buffer[0].node, path->buffer[path->size - 1].node);
		for(beg=0;beg<seq->size;beg+=100){
			end = beg + 100;
			if(end > seq->size) end = seq->size;
			ch = seq->string[end];
			seq->string[end] = '\0';
			fprintf(o_seq, "%s\n", seq->string + beg);
			seq->string[end] = ch;
		}
		fprintf(o_lay, ">utg%llu len=%d nodes=%llu\n", (unsigned long long)uid, seq->size, (unsigned long long)path->size);
		print_traces_graph(g, path, o_lay);
	}
	free_string(seq);
	fprintf(KBM_LOGF, "[%s] %llu unitigs (>= %d bp), total %llu bp\n", date(), (unsigned long long)cnt, g->min_ctg_len, (unsigned long long)tot);
	fclose(files[0]);
	fclose(files[1]);
	fclose(files[2]);
	fclose(files[3]);
	return uid;
}

/*
 * For debug in  GDB
 * local_dot_node, local_dot_step, and print_local_dot_graph()
 */

static uint64_t local_dot_node = 1;
static uint32_t local_dot_step = 10;

void get_subgraph_nodes_graph(Graph *g, ptrrefhash *nodes, u8v *stack, uint16_t max_step, uint32_t closed_val){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	ptr_ref_t *p, *pp;
	uint64_t nid, idx;
	uint32_t k, cnt;
	int exists;
	clear_u8v(stack);
	reset_iter_ptrrefhash(nodes);
	while((p = ref_iter_ptrrefhash(nodes))){
		p->cnt = 0;
		push_u8v(stack, p->idx);
	}
	while(stack->size){
		p = get_ptrrefhash(nodes, (ptr_ref_t){stack->buffer[--stack->size], 0});
		if(p->cnt >> 16) continue;
		if((p->cnt & 0xFFFF) >= max_step) continue;
		n = ref_nodev(g->nodes, p->idx);
		cnt = p->cnt;
		p->cnt |= 1U << 16;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = ref_edgev(g->edges, f->idx);
				if(e->closed >= closed_val) continue;
				nid = f->flg? e->node1 : e->node2;
				pp = prepare_ptrrefhash(nodes, (ptr_ref_t){nid, 0}, &exists);
				if(exists) continue;
				pp->idx = nid; pp->cnt = cnt + 1;
				push_u8v(stack, nid);
			}
		}
	}
}

uint64_t print_local_dot_graph(Graph *g, FILE *out){
	ptrrefhash *hash;
	u8v *stack;
	ptr_ref_t *p;
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	uint32_t j, k, max;
	hash = init_ptrrefhash(1023);
	stack = init_u8v(32);
	put_ptrrefhash(hash, (ptr_ref_t){local_dot_node, 0});
	get_subgraph_nodes_graph(g, hash, stack, local_dot_step, 1);
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	reset_iter_ptrrefhash(hash);
	while((p = ref_iter_ptrrefhash(hash))){
		i = p->idx;
		n = ref_nodev(g->nodes, i);
		r = NULL; max = 0;
		for(j=0;j<n->regs.cnt;j++){
			rr = ref_regv(g->regs, n->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r == NULL) continue;
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	reset_iter_ptrrefhash(hash);
	while((p = ref_iter_ptrrefhash(hash))){
		i = p->idx;
		n = ref_nodev(g->nodes, i);
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				if(e->closed) continue;
				if(f->flg){
					//if(!exists_ptrrefhash(hash, (ptr_ref_t){e->node1, 0})) continue;
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1]);
				} else {
					//if(!exists_ptrrefhash(hash, (ptr_ref_t){e->node2, 0})) continue;
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2]);
				}
			}
		}
	}
	fprintf(out, "}\n");
	return 0;
}

uint64_t print_dot_full_graph(Graph *g, FILE *out){
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	uint32_t j, k, max;
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		//if(n->closed) continue;
		r = NULL; max = 0;
		for(j=0;j<n->regs.cnt;j++){
			rr = ref_regv(g->regs, n->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r == NULL) continue;
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		//if(n->closed) continue;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				//if(e->closed) continue;
				if(f->flg){
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1], e->closed? " style=dashed" : "");
				} else {
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2], e->closed? " style=dashed" : "");
				}
			}
		}
	}
	fprintf(out, "}\n");
	return 0;
}

uint64_t print_dot_graph(Graph *g, FILE *out){
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	uint32_t j, k, max;
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		r = NULL; max = 0;
		for(j=0;j<n->regs.cnt;j++){
			rr = ref_regv(g->regs, n->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r == NULL) continue;
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = g->edges->buffer + f->idx;
				if(e->closed) continue;
				if(f->flg){
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1]);
				} else {
					fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2]);
				}
			}
		}
	}
	fprintf(out, "}\n");
	return 0;
}

uint64_t print_nodes_graph(Graph *g, FILE *out){
	node_t *n;
	reg_t *r;
	unsigned long long i;
	uint32_t j;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed){
			fprintf(out, "N%llu*\t%u", i, n->regs.cnt);
		} else {
			fprintf(out, "N%llu\t%u", i, n->regs.cnt);
		}
		for(j=0;j<n->regs.cnt;j++){
			r = ref_regv(g->regs, n->regs.idx + j);
			fprintf(out, "\t%s_%c_%d_%d", g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
			if(r->closed) fputc('*', out);
		}
		fprintf(out, "\n");
	}
	return i;
}

uint64_t print_reads_graph(Graph *g, FILE *out){
	read_t *rd;
	reg_t  *r;
	uint64_t idx;
	uint32_t i;
	for(i=0;i<g->kbm->reads->size;i++){
		rd = ref_readv(g->reads, i);
		fprintf(out, "%s\t%d\t%u", g->kbm->reads->buffer[i].tag, g->kbm->reads->buffer[i].rdlen, rd->regs.cnt);
		idx = rd->regs.idx;
		while(idx){
			r = ref_regv(g->regs, idx);
			fprintf(out, "\tN%llu%s:%c_%d_%d", (unsigned long long)r->node, (r->closed? "*" : (g->nodes->buffer[r->node].closed? "!" : "")), "FR"[r->dir], r->beg, r->end - r->beg);
			idx = r->read_link;
		}
		fprintf(out, "\n");
	}
	return i;
}

uint64_t print_frgs_nodes_graph(Graph *g, FILE *out){
	frg_t *frg;
	trace_t *t;
	node_t *n;
	u4i i, j;
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		if(frg->closed) continue;
		fprintf(out, "F%u\t%d\t%d\t%u\t%u", i, frg->length, frg->len, frg->tcnt, frg->ty - frg->tx);
		for(j=0;j<frg->tcnt;j++){
			t = ref_tracev(g->traces, frg->toff + j);
			if(j < frg->tx || j >= frg->ty){
				n = ref_nodev(g->nodes, t->node);
				if(n->rep_idx == MAX_REP_IDX){
					fprintf(out, "\tn%llu:%c:%d:::%d", t->node, "+-"[t->dir], t->off, t->cov);
				} else {
					fprintf(out, "\tn%llu:%c:%d:F%llu:%c:%d", t->node, "+-"[t->dir], t->off, (u8i)n->rep_idx, "+-"[n->rep_dir], t->cov);
				}
			} else {
				fprintf(out, "\tN%llu:%c:%d", t->node, "+-"[t->dir], t->off);
			}
		}
		fprintf(out, "\n");
	}
	return i;
}

uint64_t print_frgs_dot_graph(Graph *g, FILE *out){
	frg_t *frg;
	trace_t *t1, *t2;
	node_t *n1, *n2;
	reg_t *r1, *r2, *rr;
	edge_ref_t *f;
	lnk_t *e;
	unsigned long long i, idx;
	uint32_t j, k, max;
	fprintf(out, "digraph {\n");
	fprintf(out, "node [shape=record]\n");
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		if(frg->closed) continue;
		//if(frg->ty - frg->tx < (u4i)g->min_ctg_nds) continue;
		t1 = ref_tracev(g->traces, frg->toff + frg->tx);
		t2 = ref_tracev(g->traces, frg->toff + frg->ty - 1);
		n1 = ref_nodev(g->nodes, t1->node);
		n2 = ref_nodev(g->nodes, t2->node);
		r1 = NULL; max = 0;
		for(j=0;j<n1->regs.cnt;j++){
			rr = ref_regv(g->regs, n1->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r1 = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r1 == NULL) continue;
		r2 = NULL; max = 0;
		for(j=0;j<n2->regs.cnt;j++){
			rr = ref_regv(g->regs, n2->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r2 = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r2 == NULL) continue;
		fprintf(out, "F%llu [label=\"{F%llu %u %u/%u | { {N%llu:%c | %s | %c_%d_%d} | {N%llu:%c | %s | %c_%d_%d}}}\"]\n", i, i, frg->ty - frg->tx, frg->len, frg->length,
			t1->node, "+-"[t1->dir], g->kbm->reads->buffer[r1->rid].tag, "FR"[r1->dir], r1->beg, r1->end - r1->beg,
			t2->node, "+-"[t2->dir], g->kbm->reads->buffer[r2->rid].tag, "FR"[r2->dir], r2->beg, r2->end - r2->beg);
	}
	for(i=0;i<g->frgs->size;i++){
		frg = ref_frgv(g->frgs, i);
		if(frg->closed) continue;
		//if(frg->ty - frg->tx < (u4i)g->min_ctg_nds) continue;
		for(k=0;k<2;k++){
			idx = frg->lnks[k].idx;
			while(idx){
				f = ref_edgerefv(g->lrefs, idx);
				idx = f->next;
				e = g->lnks->buffer + f->idx;
				if(e->closed) continue;
				if(f->flg){
					//if(g->frgs->buffer[e->frg1].ty - g->frgs->buffer[e->frg1].tx < (u4i)g->min_ctg_nds) continue;
					fprintf(out, "F%llu -> F%llu [label=\"%c%c:%d:%d\" color=%s style=%s]\n", i, (u8i)e->frg1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1], e->weak? "dashed" : "solid");
				} else {
					//if(g->frgs->buffer[e->frg2].ty - g->frgs->buffer[e->frg2].tx < (u4i)g->min_ctg_nds) continue;
					fprintf(out, "F%llu -> F%llu [label=\"%c%c:%d:%d\" color=%s style=%s]\n", i, (u8i)e->frg2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2], e->weak? "dashed" : "solid");
				}
			}
		}
	}
	fprintf(out, "}\n");
	return 0;
}

typedef uint64_t (*graph_print_func)(Graph *g, FILE *out);

uint64_t generic_print_graph(Graph *g, graph_print_func func, char *prefix, char *suffix){
	FILE *out;
	char *file;
	uint64_t cnt;
	{
		fprintf(KBM_LOGF, "[%s] output \"%s%s\".", date(), prefix, suffix? suffix : ""); fflush(KBM_LOGF);
		file = malloc(strlen(prefix) + (suffix? strlen(suffix) : 0) + 2);
		sprintf(file, "%s%s", prefix, suffix? suffix : "");
		out = fopen(file, "w");
		cnt = func(g, out);
		fclose(out);
		free(file);
		fprintf(KBM_LOGF, " Done.\n"); fflush(KBM_LOGF);
	}
	return cnt;
}

static struct option prog_opts[] = {
	{"cpu",                              1, 0, 't'},
	{"input",                            1, 0, 'i'},
	{"err-free-seq",                     1, 0, 'I'},
	{"force",                            0, 0, 'f'},
	{"prefix",                           1, 0, 'o'},
	{"kmer-fsize",                       1, 0, 'k'},
	{"kmer-psize",                       1, 0, 'p'},
	{"kmer-depth-max",                   1, 0, 'K'},
	{"kmer-depth-min",                   1, 0, 'E'},
	{"kmer-depth-min-filter",            0, 0, 'F'},
	{"kmer-subsampling",                 1, 0, 'S'},
	{"dp-max-gap",                       1, 0, 'X'},
	{"dp-max-var",                       1, 0, 'Y'},
	{"dp-penalty-gap",                   1, 0, 'x'},
	{"dp-penalty-var",                   1, 0, 'y'},
	{"aln-min-length",                   1, 0, 'l'},
	{"aln-min-match",                    1, 0, 'm'},
	{"aln-max-var",                      1, 0, 's'},
	{"verbose",                          0, 0, 'v'},
	{"quiet",                            0, 0, 'q'},
	{"help",                             0, 0, 1000}, // detailed document
	{"tidy-reads",                       1, 0, 'L'},
	{"keep-name",                        0, 0, 1001},
	{"err-free-nodes",                   0, 0, 1002},
	{"limit-input",                      1, 0, 1003},
	{"node-len",                         1, 0, 1004},
	{"node-ovl",                         1, 0, 1005},
	{"node-drop",                        1, 0, 1006},
	{"edge-min",                         1, 0, 'e'},
	{"node-min",                         1, 0, 1007},
	{"node-max",                         1, 0, 1008},
	{"ttr-cutoff-depth",                 1, 0, 1009},
	{"ttr-cutoff-ratio",                 1, 0, 1010},
	{"dump-kbm",                         1, 0, 1011},
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
	{0, 0, 0, 0}
};

int usage(int level){
	printf(
	"WTDBG: De novo assembler for long noisy sequences\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 2.1 (20181007)\n"
#ifdef TIMESTAMP
	"Compiled: %s\n"
#endif
	"Usage: wtdbg2 [options]\n"
	"Options:\n"
	" -i <string> Long reads sequences file (REQUIRED; can be multiple), []\n"
	" -I <string> Error-free sequences file (can be multiple), []\n"
	" -o <string> Prefix of output files (REQUIRED), []\n"
	" -t <int>    Number of threads, 0 for all cores, [4]\n"
	" -f          Force to overwrite output files\n"
	" -L <int>    Choose the longest subread and drop reads shorter than <int> (5000 recommended for PacBio) [0]\n"
	" -k <int>    Kmer fsize, 0 <= k <= 25, [0]\n"
	" -p <int>    Kmer psize, 0 <= p <= 25, [21]\n"
	"             k + p <= 25, seed is <k-mer>+<p-homopolymer-compressed>\n"
	" -K <float>  Filter high frequency kmers, maybe repetitive, [1000]\n"
	"             if K >= 1, take the integer value as cutoff, MUST <= 65535\n"
	"             else, mask the top fraction part high frequency kmers\n"
	" -E <int>    Min kmer frequency, [2]\n"
	" -F          Filter low frequency kmers by a 4G-bytes array (max_occ=3 2-bits). Here, -E must greater than 1\n"
	" -S <int>    Subsampling kmers, 1/(<-S>) kmers are indexed, [4]\n"
	"             -S is very useful in saving memeory and speeding up\n"
	"             please note that subsampling kmers will have less matched length\n"
	" -X <int>    Max number of bin(256bp) in one gap, [4]\n"
	" -Y <int>    Max number of bin(256bp) in one deviation, [4]\n"
	" -x <int>    penalty for BIN gap, [-7]\n"
	" -y <int>    penalty for BIN deviation, [-21]\n"
	" -l <float>  Min length of alignment, [2048]\n"
	" -m <float>  Min matched, [200]\n"
	" -A          Keep contained reads during alignment\n"
	" -s <float>  Max length variation of two aligned fragments, [0.2]\n"
	" -e <int>    Min read depth of a valid edge, [3]\n"
	" -q          Quiet\n"
	" -v          Verbose (can be multiple)\n"
	" --help      Show more options\n"
#ifdef TIMESTAMP
	, TOSTR(TIMESTAMP)
#endif
	);
	if(level > 0){
		printf(
	" ** more options **\n"
	" --cpu <int>\n"
	"   See -t 0, default: all cores\n"
	" --input <string> +\n"
	"   See -i\n"
	" --err-free-seq <string> +\n"
	"   See -I. Error-free sequences will be firstly token for nodes, if --err-free-nodes is specified, only select nodes from those sequences\n"
	" --force\n"
	"   See -f\n"
	" --prefix <string>\n"
	"   See -o\n"
	" --kmer-fsize <int>\n"
	"   See -k 0\n"
	" --kmer-psize <int>\n"
	"   See -p 21\n"
	" --kmer-depth-max <float>\n"
	"   See -K 1000\n"
	" --kmer-depth-min <int>\n"
	"   See -E\n"
	" --kmer-depth-min-filter\n"
	"   See -F\n"
	"   `wtdbg` uses a 4 Gbytes array to counting the occurence (0-3) of kmers in the way of counting-bloom-filter. It will reduce memory space largely\n"
	"    Orphaned kmers won't appear in building kbm-index\n"
	" --kmer-subampling <int>\n"
	"   See -S 1\n"
	" --aln-kmer-sampling <int>\n"
	"   Select no more than n seeds in a query bin, default: 256\n"
	" --dp-max-gap <int>\n"
	"   See -X 4\n"
	" --dp-max-var <int>\n"
	"   See -Y 4\n"
	" --dp-penalty-gap <int>\n"
	"   See -x -7\n"
	" --dp-penalty-var <int>\n"
	"   See -y -21\n"
	" --aln-min-length <int>\n"
	"   See -l 2048\n"
	" --aln-min-match <int>\n"
	"   See -m 200. Here the num of matches counting basepair of the matched kmer's regions\n"
	" --aln-max-var <float>\n"
	"   See -s 0.2\n"
	" --aln-dovetail <int>\n"
	"   Retain dovetail overlaps only, the max overhang size is <--aln-dovetail>, the value should be times of 256, -1 to disable filtering, default: 256\n"
	" --aln-strand <int>\n"
	"   1: forward, 2: reverse, 3: both. Please don't change the deault vaule 3, unless you exactly know what you are doing\n"
	" --aln-maxhit <int>\n"
	"   Max n hits for each read in build graph, default: 1000\n"
	" --aln-bestn <int>\n"
	"   Use best n hits for each read in build graph, 0: keep all, default: 500\n"
	"   <prefix>.alignments always store all alignments\n"
	" -A, --aln-noskip\n"
	"   Even a read was contained in previous alignment, still align it against other reads\n"
	" --verbose +\n"
	"   See -v. -vvvv will display the most detailed information\n"
	" --quiet\n"
	"   See -q\n"
	" -L <int>, --tidy-reads=<int>\n"
	"   Default: 0. Pick longest subreads if possible. Filter reads less than <--tidy-reads>. Rename reads into 'S%%010d' format. The first read is named as S0000000001\n"
	"   Set to 0 bp to disable tidy. Suggested vaule is 5000 for pacbio reads\n"
	" --keep-name\n"
	"   Keep orignal read names even with --tidy-reads\n"
	" --err-free-nodes\n"
	"   Select nodes from error-free-sequences only. E.g. you have contigs assembled from NGS-WGS reads, and long noisy reads.\n"
	"   You can type '--err-free-seq your_ctg.fa --input your_long_reads.fa --err-free-nodes' to perform assembly somehow act as long-reads scaffolding\n"
	" --limit-input <int>\n"
	"   Limit the input sequences to at most <int> M bp. Usually for test\n"
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
	" --drop-low-cov-edges\n"
	"   Don't attempt to rescue low coverage edges\n"
	" --node-min <int>\n"
	"   Min depth of a intreval to be selected as valid node. Defaultly, this value is automaticly the same with --edge-min.\n"
	" --node-max <int>\n"
	"   Nodes with too high depth will be regarded as repetitive, and be masked. Default: 200, more than 200 reads contain this node\n"
	" --ttr-cutoff-depth <int>, 0\n"
	" --ttr-cutoff-ratio <float>, 0.5\n"
	"   Tiny Tandom Repeat. A node located inside ttr will bring noisy in graph, should be masked. The pattern of such nodes is:\n"
	"   depth >= <--ttr-cutoff-depth>, and none of their edges have depth greater than depth * <--ttr-cutoff-ratio 0.5>\n"
	"   set --ttr-cutoff-depth 0 to disable ttr masking\n"
	" --dump-kbm <string>\n"
	"   Dump kbm index into file for loaded by `kbm` or `wtdbg`\n"
	" --load-kbm <string>\n"
	"   Instead of reading sequences and building kbm index, which is time-consumed, loading kbm-index from already dumped file.\n"
	"   Please note that, once kbm-index is mmaped by kbm -R <kbm-index> start, will just get the shared memory in minute time.\n"
	"   See `kbm` -R <your_seqs.kbmidx> [start | stop]\n"
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
	return 1;
}

int main(int argc, char **argv){
	Graph *g;
	KBMPar *par;
	KBM *kbm;
	FileReader *fr;
	BioSequence *seqs[2], *seq;
	cplist *pbs, *ngs, *pws;
	FILE *evtlog;
	char *prefix, *dump_kbm, *load_kbm, *load_nodes, *load_clips;
	char regtag[14];
	int len, tag_size, asyn_read;
	u8i tot_bp, cnt, bub, tip, rep, yarn, max_bp, max_idx_bp, nfix, opt_flags;
	uint32_t i, j, k;
	int c, opt_idx, ncpu, only_fix, node_cov, max_node_cov, exp_node_cov, min_bins, edge_cov, store_low_cov_edge, reglen, regovl, bub_step, tip_step, rep_step;
	int frgtip_len, ttr_n_cov;
	int quiet, tidy_reads, tidy_rdtag, less_out, tip_like, cut_tip, rep_filter, out_alns, cnn_filter, log_rep, rep_detach, del_iso, rdclip, chainning, bestn, rescue_low_edges;
	int min_ctg_len, min_ctg_nds, max_trace_end, max_overhang, overwrite, node_order, fast_mode;
	float node_drop, node_mrg, ttr_e_cov, fval;
	pbs = init_cplist(4);
	ngs = init_cplist(4);
	pws = init_cplist(4);
	asyn_read = 1;
	ncpu = 4;
	tidy_reads = 0;
	tidy_rdtag = -1;
	fast_mode = 0;
	max_bp = 0;
	max_idx_bp = 0LLU * 1000 * 1000 * 1000; // unlimited
	reglen = 1024;
	regovl = 256;
	node_drop = 0.25;
	node_mrg = 0.9;
	only_fix = 0;
	node_cov = 0;
	max_node_cov = 200;
	exp_node_cov = 40;
	min_bins = 1;
	edge_cov = 3;
	rdclip = 1;
	chainning = 1;
	bestn = 500;
	ttr_n_cov = 0;
	ttr_e_cov = 0.5;
	dump_kbm = NULL;
	load_kbm = NULL;
	load_clips = NULL;
	load_nodes = NULL;
	store_low_cov_edge = 1;
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
	par->min_aln = 1024 * 2;
	par->min_mat = 200;
	opt_flags = 0;
	while((c = getopt_long(argc, argv, "ht:i:I:fo:FE:k:p:K:S:X:Y:x:y:l:m:s:vqe:L:A", prog_opts, &opt_idx)) != -1){
		switch(c){
			case 't': ncpu = atoi(optarg); break;
			case 'i': push_cplist(pbs, optarg); break;
			case 'I': push_cplist(ngs, optarg); par->rd_len_order = 0; break;
			case 'f': overwrite = 1; break;
			case 'o': prefix = optarg; break;
			case 'k': par->ksize = atoi(optarg); opt_flags |= (1 << 1); break;
			case 'p': par->psize = atoi(optarg); opt_flags |= (1 << 0); break;
			case 'K': fval = atof(optarg);
					if(fval > 1) par->kmax = fval;
					else { par->kmax = 0; par->ktop = fval; }
					break;
			case 'E': par->kmin = atoi(optarg); break;
			case 'F': par->use_kf = 1; break;
			case 'S': par->kmer_mod = atoi(optarg) * KBM_N_HASH; opt_flags |= (1 << 2);break;
			case 'X': par->max_bgap = atoi(optarg); break;
			case 'Y': par->max_bvar = atoi(optarg); break;
			case 'x': par->pgap = atoi(optarg); break;
			case 'y': par->pvar = atoi(optarg); break;
			case 'l': par->min_aln = atoi(optarg); break;
			case 'm': par->min_mat = atoi(optarg); break;
			case 's': par->aln_var = atof(optarg); break;
			case 'v': KBM_LOG ++; break;
			case 'q': quiet = 1; break;
			case 'h': return usage(0);
			case 1000: return usage(1);
			case 'L':  tidy_reads = atoi(optarg); break;
			case 1001: tidy_rdtag = 0; break;
			case 1002: only_fix = 1; break;
			case 1003: max_bp = atol(optarg); break;
			case 1004: reglen = atoi(optarg); break;
			case 1005: regovl = atoi(optarg); break;
			case 1006: node_drop = atof(optarg); break;
			case 'e':  edge_cov = atoi(optarg); break;
			case 1007: node_cov = atoi(optarg); break;
			case 1008: max_node_cov = atoi(optarg); break;
			case 1009: ttr_n_cov = atoi(optarg); break;
			case 1010: ttr_e_cov = atof(optarg); break;
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
			case 'A':  par->skip_contained = 0; break;
			case 1031: min_bins = atoi(optarg); break;
			case 1032: rescue_low_edges = 1; break;
			case 1033: rescue_low_edges = 0; break;
			default: return usage(0);
		}
	}
	if (optind == 1) return usage(0);
	if(prefix == NULL) {
		fprintf(stderr, "ERROR: please specify the output prefix with -o\n");
		return 1;
	}
	if(load_kbm == NULL && pbs->size + ngs->size == 0) {
		fprintf(stderr, "ERROR: please specify the input with -i, -I or --load-kbm\n");
		return 1;
	}
	if((reglen % KBM_BIN_SIZE)){
		reglen = ((reglen + KBM_BIN_SIZE - 1) / KBM_BIN_SIZE) * KBM_BIN_SIZE;
		fprintf(stderr, " ** Adjust -j to %d\n", reglen);
	}
	if(node_cov == 0) node_cov = edge_cov;
	if(!overwrite && file_exists(prefix)){
		fprintf(stderr, "File exists! '%s'\n\n", prefix);
		return usage(0);
	}
	if(max_idx_bp == 0) max_idx_bp = 0xFFFFFFFFFFFFFFFFLLU;
	if(par->ksize + par->psize > 25){
		fprintf(stderr, " -- Invalid kmer size %d+%d=%d > 25 in %s -- %s:%d --\n", par->ksize, par->psize, par->ksize + par->psize,  __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(quiet){
		int devnull;
		devnull = open("/dev/null", O_WRONLY);
		dup2(devnull, STDERR_FILENO);
	}
	if(tidy_rdtag == -1){
		tidy_rdtag = tidy_reads;
	}
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
		fprintf(KBM_LOGF, "[%s] Done. %u sequences, %llu bp, parameter('-S %d')\n", date(), (u4i)kbm->reads->size, (u8i)kbm->rdseqs->size, kbm->par->kmer_mod / KBM_N_HASH);
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
		nfix = 0;
		tot_bp = kbm->rdseqs->size;
	} else {
		kbm = init_kbm(par);
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
		}
		regfree(&reg);
		free_biosequence(seqs[0]);
		free_biosequence(seqs[1]);
		if(!KBM_LOG){ fprintf(KBM_LOGF, "\r%u reads", (unsigned)kbm->reads->size); fflush(KBM_LOGF); }
		ready_kbm(kbm);
		fprintf(KBM_LOGF, "\n[%s] Done, %u reads, %llu bp, %u bins\n", date(), (unsigned)kbm->reads->size, tot_bp, (u4i)kbm->bins->size); fflush(KBM_LOGF);
	}
	print_proc_stat_info(0);
	g = init_graph(kbm);
	g->node_order = node_order;
	g->reglen = reglen;
	g->regovl = regovl;
	g->max_overhang = max_overhang;
	g->node_max_conflict = node_drop;
	g->node_merge_cutoff = node_mrg;
	g->min_node_cov = node_cov;
	g->max_node_cov_sg = node_cov;
	g->max_node_cov = max_node_cov;
	g->exp_node_cov = exp_node_cov;
	g->min_node_mats = min_bins;
	g->min_edge_cov = edge_cov;
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
	g->chainning_hits = (rdclip && chainning);
	g->bestn = bestn;
	g->minimal_output = less_out;
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
		////uint32_t mid = estimate_genome_size(g, tot_bp, KBM_LOGF);
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
	if(!less_out) generic_print_graph(g, print_dot_full_graph,   prefix, ".1.dot");
	fprintf(KBM_LOGF, "[%s] graph clean\n", date()); fflush(KBM_LOGF);
	if(0){
		cnt = mask_read_weak_regs_graph(g, ncpu);
		fprintf(KBM_LOGF, "[%s] masked %llu regions(%d bp) as unreliable, total regs %llu\n", date(), (unsigned long long)cnt, reglen, (u8i)g->regs->size);
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
	cnt = myers_transitive_reduction_graph(g, 100.2f);
	set_init_ends_graph(g);
	fprintf(KBM_LOGF, "[%s] cut %llu transitive edges\n", date(), (unsigned long long)cnt);
	if(del_iso){
		cnt = del_isolated_nodes_graph(g, evtlog);
		if(cnt){
			fprintf(KBM_LOGF, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
		}
	}
	if(!less_out) generic_print_graph(g, print_dot_graph,   prefix, ".2.dot");
	{
		bub = tip = rep = yarn = 0;
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
		if(bub + tip){ fprintf(KBM_LOGF, "[%s] %llu bubbles; %llu tips; %llu yarns;\n", date(), bub, tip, yarn); fflush(KBM_LOGF); }
	}
	if(del_iso){
		cnt = del_isolated_nodes_graph(g, evtlog);
		fprintf(KBM_LOGF, "[%s] deleted %llu isolated nodes\n", date(), (unsigned long long)cnt);
	}
	if(!less_out) generic_print_graph(g, print_dot_graph,   prefix, ".3.dot");
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
	if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.dot");
	if(1){
		cnt = rescue_weak_tip_lnks_graph(g);
		fprintf(KBM_LOGF, "[%s] rescue %llu weak links\n", date(), (unsigned long long)cnt);
	}
	cnt = cut_binary_lnks_graph(g, evtlog);
	fprintf(KBM_LOGF, "[%s] deleted %llu binary links\n", date(), (unsigned long long)cnt);
	//cnt = reduce_transitive_lnks_graph(g);
	cnt = myers_transitive_reduction_frg_graph(g, 10000.1f);
	fprintf(KBM_LOGF, "[%s] cut %llu transitive links\n", date(), (unsigned long long)cnt);
	cnt = remove_boomerangs_frg_graph(g, 30 * 1000);
	fprintf(KBM_LOGF, "[%s] remove %llu boomerangs\n", date(), (unsigned long long)cnt);
	cnt = 0;
	while(1){
		u8i c;
		if((c = detach_repetitive_frg_graph(g, 100 * 1000)) == 0){
			break;
		}
		cnt += c;
	}
	fprintf(KBM_LOGF, "[%s] detached %llu repeat-associated paths\n", date(), (unsigned long long)cnt);
	cnt = cut_weak_branches_frg_graph(g);
	fprintf(KBM_LOGF, "[%s] remove %llu weak branches\n", date(), (unsigned long long)cnt);
	//cnt = cut_low_cov_lnks_graph(g, 1);
	//fprintf(KBM_LOGF, "[%s] deleted %llu low cov links\n", date(), (unsigned long long)cnt);
	//if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.2.dot");
	cnt = trim_frgtips_graph(g, frgtip_len);
	fprintf(KBM_LOGF, "[%s] cut %llu tips\n", date(), (unsigned long long)cnt);
	//if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".frg.3.dot");
	cnt = pop_frg_bubbles_graph(g, bub_step);
	fprintf(KBM_LOGF, "[%s] pop %llu bubbles\n", date(), (unsigned long long)cnt);
	cnt = trim_frgtips_graph(g, frgtip_len);
	fprintf(KBM_LOGF, "[%s] cut %llu tips\n", date(), (unsigned long long)cnt);
	if(!less_out) generic_print_graph(g, print_frgs_dot_graph, prefix, ".ctg.dot");
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
	cnt = print_ctgs_graph(g, 0, 0, g->major_nctg, prefix, ".ctg.lay", ncpu, evtlog);
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
	free_graph(g);
	fprintf(KBM_LOGF, "[%s] Program Done\n", date());
	END_STAT_PROC_INFO(stderr);
	return 0;
}
