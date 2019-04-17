#ifndef __WTDBG_H_RJ
#define __WTDBG_H_RJ

#include "kbm.h"
#include "filewriter.h"
#include "pgzf.h"
#include <getopt.h>
#include <regex.h>

#define WT_MAX_RD			0x3FFFFFFF // 30 bits, 1G
#define WT_MAX_RDLEN		0x03FFFFFF // 24 bits, 16M
#define WT_MAX_RDBIN		0x0003FFFF // 16 bits, 64K
#define WT_MAX_BIDX			0x000000FFFFFFFFFFLLU // 40 bits, 1T
#define WT_MAX_KBIN			0xF // 4 bits, 16
#define WT_MAX_NODE			0x000000FFFFFFFFFFLLU
#define WT_MAX_NODE_COV		0xFFFF
#define WT_MAX_EDGE			0x000000FFFFFFFFFFLLU
#define WT_MAX_NODE_EDGES	0xFFFF
#define WT_MAX_EDGE_COV		0x3FF
#define WT_MAX_EDGE_LEN		0x3FF
#define MAX_UTG_IDX			0xFFFFFFFFU // 32 bits, 4G

typedef struct {
	u8i node:40, dir:1, off:16, len:4, closed:3;
} rd_reg_t;
define_list(rdregv, rd_reg_t);

typedef struct {
	u8i bidx:40, dir:1, off:16, len:4, closed:3;
	u8i rid;
} rb_reg_t;
define_list(rbregv, rb_reg_t);

#define WT_MAX_RDREG_AUX	0xFFF // rdreg_aux < WT_MAX_RDREG_AUX, can be directly located by g->rdregs->buffer[read_t->regs.idx + reg_t->rdreg_aux]
typedef union {
	u8i value;
	struct { u8i rid:30, dir:1, off:16, len:4, closed:1, rdreg_aux:12; };
} reg_t;
define_list(regv, reg_t);

typedef struct { uint64_t idx:40, cnt:24; } vec_ref_t;
typedef struct { uint64_t idx:40, cnt:24; } ptr_ref_t;
typedef struct { uint64_t idx:40, cnt:24, fix:1, rank:45, score:18; } rnk_ref_t;
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
	u8i node:40, dir:1, closed:2, cov:10;
	b8i off:11;
	u8i next:40, flg:1, status:1, flag:2, port:20;
} mono_edge_t;

typedef struct {
	mono_edge_t es[2];
} edge_t;
define_list(edgev, edge_t);

static inline uint64_t _edge_hashcode(edge_t e){
	const uint64_t m = 0xc6a4a7935bd1e995LLU;
	const int r = 47;
	uint64_t h = 1023 ^ (16 * m);
	uint64_t k = (e.es[0].node << 1) | e.es[0].dir;
	k *= m; k ^= k >> r; k *= m; h ^= k; h *= m; k = (e.es[1].node << 1) | e.es[1].dir; k *= m; k ^= k >> r; k *= m; h ^= k; h *= m; h ^= h >> r; h *= m; h ^= h >> r;
	return h;
}
#define EDGEHASH(idx) ((edgev *)set->userdata)->buffer[idx]
#define edge_hashcode(E) _edge_hashcode(EDGEHASH(E))
#define edge_hashequals(E1, E2) (EDGEHASH(E1).es[0].node == EDGEHASH(E2).es[0].node && EDGEHASH(E1).es[1].node == EDGEHASH(E2).es[1].node \
			&& EDGEHASH(E1).es[0].dir == EDGEHASH(E2).es[0].dir && EDGEHASH(E1).es[1].dir == EDGEHASH(E2).es[1].dir)
define_hashset(edgehash, u8i, edge_hashcode, edge_hashequals);

typedef struct { u8i idx:40, flg:1, flg2:1, cnt:22; } edge_ref_t;
static const edge_ref_t EDGE_REF_NULL = (edge_ref_t){0, 0, 0, 0};
define_list(edgerefv, edge_ref_t);

typedef struct { u8i idx:63, flg:1, next; } link_ref_t;
static const link_ref_t LINK_REF_NULL = (link_ref_t){0x7FFFFFFFFFFFFFFLLU, 1, 0};
define_list(linkrefv, link_ref_t);

typedef struct {
	u4i utg_idx;
	u4i utg_dir:1, bt_idx:31;
	u8i closed:1, bt_visit:45, bt_dir:1, cov:16, init_end:2;
	vec_ref_t regs;
	edge_ref_t erefs[2];
} node_t;
define_list(nodev, node_t);

typedef struct {
	vec_ref_t regs;
} read_t;
define_list(readv, read_t);

typedef struct {
	u2i refidx, refdir:1, mat:15;
	u2i qb, qe;
	u4i tb, te;
} read_map_t;
define_list(readmapv, read_map_t);

typedef struct {
	u8i node:40, dir:1, cov:23;
	edge_ref_t edges[2];
	u4i off;
} trace_t;
define_list(tracev, trace_t);
#define WT_TRACE_MSG_ZERO	0
#define WT_TRACE_MSG_ONE	1
#define WT_TRACE_MSG_MORE	2
#define WT_TRACE_MSG_VISITED	3
#define WT_TRACE_MSG_UNDEF	4

typedef struct {
	u8i node1:63, dir1:1;
	u8i node2:63, dir2:1;
	b8i off:42, len:22;
} seqlet_t;
define_list(seqletv, seqlet_t);

typedef struct { uint64_t node:63, dir:1; } utg_node_t;

typedef struct {
	utg_node_t nodes[2];
	u4i len, cnt:24, bt_nin:8;
	u8i bt_vst:40, bt_idx:20, bt_dir:1, sel:2, bt_use:1, closed:1;
	u4i bt_ptr, bt_cin:8, bt_cov:12, bt_sum:12;
} utg_t;
define_list(utgv, utg_t);

typedef struct {
	KBM      *kbm;
	KBMPar   *par, *rpar;
	rdregv   *rdregs;
	readv    *reads;
	readmapv *rdmaps;
	cplist   *reftags;
	cuhash   *ref2idx;

	regv     *regs;
	nodev    *nodes;
	edgev    *edges;
	edgehash *ehash;

	u8i      genome_size;
	u4i      num_index;
	//u4i      corr_mode, corr_min, corr_max;
	//u4i      corr_bsize, corr_bstep;
	//float    corr_cov, corr_gcov;
	int      node_order, mem_stingy;
	u4i      n_fix, only_fix; // first n sequences are accurate contigs; only_fix means whether to include other pacbio sequenes
	u4i      reglen, regovl, bestn;
	int      min_node_mats;
	//int      max_overhang, chainning_hits, uniq_hit;
	float    node_max_conflict; // 0.25
	//float    node_merge_cutoff;
	u4i      max_node_cov, min_node_cov, exp_node_cov, min_edge_cov;
	u4i      max_node_cov_sg, max_sg_end;
	int      cut_tip;
	int      store_low_cov_edge;
	int      rep_filter, rep_detach;
	uint32_t bub_step, tip_step, rep_step;
	u4i      min_ctg_len, min_ctg_nds, minimal_output;
	utgv *utgs;
	u4i major_nctg;
} Graph;

static const char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};

static inline Graph* init_graph(KBM *kbm){
	Graph *g;
	g = malloc(sizeof(Graph));
	g->kbm = kbm;
	g->par = kbm->par;
	g->rpar = NULL;
	g->rdregs = init_rdregv(32);
	g->reads = init_readv(kbm->reads->size);
	g->reads->size = kbm->reads->size;
	g->rdmaps = NULL;
	g->reftags = NULL;
	g->ref2idx = NULL;
	g->regs = init_regv(32);
	g->nodes = init_nodev(32);
	g->edges = init_edgev(32);
	g->ehash = init_edgehash(1023);
	set_userdata_edgehash(g->ehash, g->edges);
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
	//g->node_merge_cutoff = 0.8;
	//g->max_overhang = -1;
	//g->bestn = 0;
	g->min_node_mats = 1;
	//g->chainning_hits = 0;
	//g->uniq_hit = 0;
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
	g->utgs = init_utgv(32);
	g->major_nctg = 0;
	return g;
}

static inline void free_graph(Graph *g){
	free_rdregv(g->rdregs);
	free_readv(g->reads);
	if(g->rdmaps) free_readmapv(g->rdmaps);
	if(g->reftags){
		u4i i;
		for(i=0;i<g->reftags->size;i++){
			free(g->reftags->buffer[i]);
		}
		free_cplist(g->reftags);
	}
	if(g->ref2idx){
		free_cuhash(g->ref2idx);
	}
	free_regv(g->regs);
	free_nodev(g->nodes);
	free_edgev(g->edges);
	free_edgehash(g->ehash);
	free_utgv(g->utgs);
	free(g);
}

static inline void load_readmaps_paf(Graph *g, FileReader *fr){
	read_map_t *map;
	cuhash_t *u;
	u4i ridx, rdlen, mat;
	int ncol, exists;
	if(g->kbm->tag2idx == NULL) return;
	g->ref2idx = init_cuhash(13);
	g->rdmaps = init_readmapv(g->reads->size);
	g->reftags = init_cplist(32);
	while((ncol = readtable_filereader(fr)) >= 0){
		if(ncol < 12) continue;
		ridx = getval_cuhash(g->kbm->tag2idx, get_col_str(fr, 0));
		if(ridx == MAX_U4) continue;
		map = ref_readmapv(g->rdmaps, ridx);
		mat = atoi(get_col_str(fr, 9));
		if(mat <= map->mat) continue;
		u = prepare_cuhash(g->ref2idx, get_col_str(fr, 5), &exists);
		if(!exists){
			u->key = strdup(get_col_str(fr, 5));
			u->val = g->reftags->size;
			push_cplist(g->reftags, u->key);
		}
		map->refidx = u->val;
		map->refdir = (get_col_str(fr, 4)[0] == '-');
		map->mat = mat;
		rdlen = atoi(get_col_str(fr, 1));
		if(map->refdir){
			map->qb = (rdlen - atoi(get_col_str(fr, 3))) / KBM_BIN_SIZE;
			map->qe = (rdlen - atoi(get_col_str(fr, 2)) + KBM_BIN_SIZE - 1) / KBM_BIN_SIZE;
		} else {
			map->qb = atoi(get_col_str(fr, 2)) / KBM_BIN_SIZE;
			map->qe = (atoi(get_col_str(fr, 3)) + KBM_BIN_SIZE - 1) / KBM_BIN_SIZE;
		}
		map->tb = atol(get_col_str(fr, 7));
		map->te = atol(get_col_str(fr, 8));
	}
}

// qlen is measured by KBM_BIN_SIZE
static inline void kbmap2regs_graph(Graph *g, rbregv *rbregs, int qlen, kbm_map_t *hit, BitsVec *cigars, u4v *maps[3]){
	KBM *kbm;
	u8i ndoff;
	u4i bpos[2][2], npos[2][2], clen, ndbeg, qn, j, qbincnt;
	int tmp, bt, tlen, x, y, mat, beg, end, min_node_len, max_node_len, closed;
	kbm = g->kbm;
	qn = g->reglen;
	min_node_len = (qn - 1);
	max_node_len = (qn + 1);
	// translate into reverse sequence order
	if(qlen == 0){
		qbincnt = kbm->reads->buffer[hit->qidx].bincnt;
		qlen = qbincnt;
#ifdef DEBUG
	} else if(qlen > kbm->reads->buffer[hit->qidx].bincnt){
		fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
#endif
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
				if(closed){
				//} else if((ndoff % g->reglen) == 0){
				} else {
					push_rbregv(rbregs, (rb_reg_t){ndoff, hit->qdir, beg, end - beg, 0, hit->tidx});
				}
			}
			if(hit->qdir){
				ndoff --;
			} else {
				ndoff ++;
			}
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
				if(closed){
				//} else if((ndoff % g->reglen) == 0){
				} else {
					push_rbregv(rbregs, (rb_reg_t){ndoff, hit->qdir, beg, end - beg, 0, hit->qidx});
				}
			}
			ndoff ++;
		}
	}
	if(hit->qdir){
		tmp = qlen - hit->qb;
		hit->qb = qlen - hit->qe;
		hit->qe = tmp;
	}
}

thread_beg_def(mdbg);
Graph *g;
KBMAux *aux, *raux;
u4i rid;
u4i beg, end;
FILE *alno;
int task;
thread_end_def(mdbg);

thread_beg_func(mdbg);
Graph *g;
KBM *kbm;
KBMAux *aux, *raux;
kbm_read_t *rd, *rd2;
u4v *maps[3], *tidxs;
u4i rid;
g = mdbg->g;
kbm = g->kbm;
aux = mdbg->aux;
raux = mdbg->raux;
maps[0] = init_u4v(32);
maps[1] = init_u4v(32);
maps[2] = init_u4v(32);
tidxs = init_u4v(16);
thread_beg_loop(mdbg);
if(mdbg->task == 1){
	rid = mdbg->rid;
	if(rid == MAX_U4) continue;
	rd = ref_kbmreadv(kbm->reads, rid);
	if(rd->closed) continue;
	query_index_kbm(aux, NULL, rid, kbm->rdseqs, kbm_read_seqoff(rd), kbm_read_seqlen(rd));
	map_kbm(aux);
	if(raux && aux->hits->size){ // refine
		u4i i, j, tidx;
		clear_kbm(raux->kbm);
		bitpush_kbm(raux->kbm, NULL, 0, kbm->rdseqs->bits, kbm_read_seqoff(rd), kbm_read_seqlen(rd));
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
			rd2 = ref_kbmreadv(aux->kbm->reads, tidx);
			query_index_kbm(raux, rd2->tag, tidx, aux->kbm->rdseqs, kbm_read_seqoff(rd2), kbm_read_seqlen(rd2));
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

/*
// Roundly, Matthew effect, make sure the selected k-bins are less affected by each other
static inline u4i select_read_rbregs_core(Graph *g, u4i rid, rbregv *regs, u2i *kbcovs, u1i *kbrks, u4v *idxs, int finish){
	read_t *rd;
	rb_reg_t *r1, *r2;
	u4i i, j, idx, closed, ret;
	int ovl;
	rd = ref_readv(g->reads, rid);
	clear_and_encap_u4v(idxs, rd->regs.cnt);
	for(i=0;i<rd->regs.cnt;i++){
		idxs->buffer[i] = i;
		regs->buffer[rd->regs.idx + i].closed = 1;
	}
	idxs->size = rd->regs.cnt;
	if(g->kbm->reads->buffer[rid].bincnt * g->max_node_cov < rd->regs.cnt){
		// reject to cal highly repetitive read
		return 0;
	}
	sort_array(idxs->buffer, rd->regs.cnt, u4i,
		num_cmpgtxx(kbrks[regs->buffer[rd->regs.idx + b].bidx], kbrks[regs->buffer[rd->regs.idx + a].bidx],
			kbcovs[regs->buffer[rd->regs.idx + b].bidx], kbcovs[regs->buffer[rd->regs.idx + a].bidx],
			regs->buffer[rd->regs.idx + a].bidx, regs->buffer[rd->regs.idx + b].bidx));
	ret = 0;
	for(i=0;i<rd->regs.cnt;i++){
		idx = idxs->buffer[i];
		r1 = ref_rbregv(regs, rd->regs.idx + idx);
		closed = 0;
		for(j=idx+1;j<rd->regs.cnt;j++){
			r2 = ref_rbregv(regs, rd->regs.idx + j);
			ovl = Int(r1->off + r1->len) - Int(r2->off);
			if(ovl <= Int(g->regovl)) break;
			if(r2->closed) continue;
			closed = 1;
			break;
		}
		if(closed) continue;
		for(j=idx;j>0;j--){
			r2 = ref_rbregv(regs, rd->regs.idx + j - 1);
			ovl = Int(r2->off + r2->len) - Int(r1->off);
			if(ovl <= Int(g->regovl)) break;
			if(r2->closed) continue;
			closed = 1;
			break;
		}
		if(closed) continue;
		r1->closed = 0;
		ret ++;
		// TODO: if exceeds, extra codes need to locate rd_reg_t from reg_t
		if(ret >= WT_MAX_RDREG_AUX){
			break;
		}
	}
	static int print_rst = 0;
	if(print_rst){
		print_read_rbregs(stdout, g, rid, regs, kbcovs, 1);
	}
	if(finish){
		sort_array(regs->buffer + rd->regs.idx, rd->regs.cnt, rb_reg_t, num_cmpgtx(a.closed, b.closed, a.off, b.off));
		rd->regs.cnt = ret;
	}
	return ret;
}
*/

static inline u8i load_alignments_core(Graph *g, FileReader *pws, rbregv *regs, u4v *maps[3]){
	BitsVec *cigars;
	u4v *cgs;
	kbm_map_t *hit, HIT;
	cuhash_t *cu;
	u8i nhit;
	u4i i, rid, qlen;
	int val, flg, nwarn, mwarn, ncol;
	char *cgstr, *qtag;
	mwarn = 20;
	nwarn = 0;
	cgs = init_u4v(4);
	cigars = init_bitsvec(1024, 3);
	memset(&HIT, 0, sizeof(kbm_map_t));
	hit = &HIT;
	nhit = 0;
	while((ncol = readtable_filereader(pws)) != -1){
		if((pws->n_line % 100000) == 0){
			fprintf(KBM_LOGF, "\r%llu", pws->n_line); fflush(KBM_LOGF);
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
		qlen = atol(get_col_str(pws, 2));
		if(qlen != (u4i)(g->kbm->reads->buffer[hit->qidx].bincnt) * KBM_BIN_SIZE){
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
		if(hit->qidx >= rid){
			continue;
		}
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
		clear_bitsvec(cigars);
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
		kbmap2regs_graph(g, regs, qlen / KBM_BIN_SIZE, hit, cigars, maps);
		nhit ++;
	}
	fprintf(KBM_LOGF, "\r%llu lines, %llu hits\n", pws->n_line, nhit);
	free_bitsvec(cigars);
	free_u4v(cgs);
	return nhit;
}

static inline u8i load_alignments_sam_core(Graph *g, FileReader *sam, rbregv *regs, u4v *maps[3]){
	BitsVec *cigars;
	u4v *cgs, *bts;
	kbm_map_t *hit, HIT;
	cuhash_t *cu;
	u8i nhit;
	u4i i, rid, flag, op, ol, lens[2];
	u4i qlen, tlen, qdel;
	int nwarn, mwarn, ncol;
	char *cgstr, *qtag;
	mwarn = 20;
	nwarn = 0;
	cgs = init_u4v(4);
	bts = init_u4v(64);
	cigars = init_bitsvec(1024, 3);
	memset(&HIT, 0, sizeof(kbm_map_t));
	hit = &HIT;
	nhit = 0;
	while((ncol = readtable_filereader(sam)) != -1){
		if((sam->n_line % 100000) == 0){
			fprintf(KBM_LOGF, "\r%llu", sam->n_line); fflush(KBM_LOGF);
		}
		if(sam->line->buffer[0] == '@'){
			continue;
		}
		qtag = get_col_str(sam, 0);
		if((cu = get_cuhash(g->kbm->tag2idx, qtag)) == NULL){
			if(nwarn < mwarn){
				fprintf(stderr, " -- Cannot find read \"%s\" in LINE:%llu in %s -- %s:%d --\n", qtag, sam->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				nwarn ++;
			}
			continue;
		} else rid = cu->val;
		hit->qidx = rid;
		flag = atol(get_col_str(sam, 1));
		hit->qdir = (flag & 0x10) >> 4;
		qlen = g->kbm->reads->buffer[hit->qidx].bincnt * KBM_BIN_SIZE;
		hit->qb = 0;
		hit->qe = 0;
		if((cu = get_cuhash(g->kbm->tag2idx, get_col_str(sam, 2))) == NULL){
			if(nwarn < mwarn){
				fprintf(stderr, " -- Cannot find read \"%s\" in LINE:%llu in %s -- %s:%d --\n", get_col_str(sam, 2), sam->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				nwarn ++;
			}
			continue;
		} else rid = cu->val;
		if(hit->qidx >= rid) continue;
		hit->tidx = rid;
		tlen = g->kbm->reads->buffer[hit->tidx].bincnt * KBM_BIN_SIZE;
		hit->tdir = 0;
		hit->tb = atol(get_col_str(sam, 3));
		hit->te = hit->tb;
		hit->mat = 0;
		hit->aln = 0;
		hit->cnt = 0;
		hit->gap = 0;
		clear_u4v(cgs);
		cgstr = get_col_str(sam, 5);
		op = ol = 0;
		lens[0] = lens[1] = 0;
		while(cgstr[0]){
			if(cgstr[0] >= '0' && cgstr[0] <= '9'){
				ol = ol * 10 + (cgstr[0] - '0');
			} else {
				switch(cgstr[0]){
					case '=':
					case 'X':
					case 'M': op = 0b11; break;
					case 'S':
					case 'I': op = 0b10; break;
					case 'N':
					case 'H':
					case 'D': op = 0b01; break;
					case 'P': op = 0b00; break;
					default:
						fprintf(stderr, " -- Bad cigar '%c' \"%s\" in LINE:%llu in %s -- %s:%d --\n", cgstr[0], get_col_str(sam, 5), sam->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						exit(1);
				}
				if(ol == 0) ol = 1;
				{
					if(op & 0b01){
						lens[1] += ol;
					} else if(op & 0b10){
						lens[0] += ol;
					}
				}
				push_u4v(cgs, (ol << 8) | op);
				ol = 0;
			}
			cgstr ++;
		}
		// trim cigar right
		while(cgs->size){
			op = cgs->buffer[cgs->size - 1] >> 8;
			ol = cgs->buffer[cgs->size - 1] & 0xFF;
			if(op == 0b11) break;
			cgs->size --;
		}
		ol = op = i = 0;
		if(hit->qdir){ // KBM trim read into BIN size, but SAM not. Need to be careful in reverse strand
			qdel = lens[0] - qlen;
			while(i < cgs->size){
				ol = cgs->buffer[i] >> 8;
				op = cgs->buffer[i] & 0xFF;
				i ++;
				if(op & 0b10){
					if(ol >= qdel){
						ol -= qdel;
						//hit->qe += qdel;
						if(op & 0b01) hit->te += qdel;
					} else {
						qdel -= ol;
						//hit->qe += ol;
						if(op & 0b01) hit->te += ol;
					}
				} else if(op & 0b01){
					hit->te += ol;
				}
			}
		}
		// trim cigar left
		while(1){
			if(ol == 0){
				if(i >= cgs->size) break;
				ol = cgs->buffer[i] >> 8;
				op = cgs->buffer[i] & 0xFF;
				i ++;
			}
			if(op == 0b11) break;
			hit->te += ol * (op & 0b01);
			hit->qe += ol * ((op >> 1) & 0b01);
		}
		//call bin cigars
		hit->qb = hit->qe;
		hit->tb = hit->te;
		u4i qlst, tlst, mat;
		u4i bin_min_mat = 21;
		qlst = hit->qe / KBM_BIN_SIZE;
		tlst = hit->te / KBM_BIN_SIZE;
		mat = 0;
		clear_u4v(bts);
		while(1){
			if(ol == 0){
				if(i >= cgs->size) break;
				ol = cgs->buffer[i] >> 8;
				op = cgs->buffer[i] & 0xFF;
				i ++;
			}
			if(op == 0b11){
				while(ol){
					u4i l, d, qidx, ql, tidx, tl;
					qidx = hit->qe / KBM_BIN_SIZE;
					ql = ((~hit->qe) & 0xFFU) + 1;
					tidx = hit->te / KBM_BIN_SIZE;
					tl = ((~hit->te) & 0xFFU) + 1;
					l = num_min(ql, tl);
					l = num_min(l, ol);
					while(1){
						d = ((qlst < qidx)? 1 : 0) | ((tlst < tidx)? 2 : 0);
						if(d == 0)  break;
						if(mat < bin_min_mat){
							d |= 4;
						}
						qlst += d & 0x01;
						tlst += (d >> 1) & 0x01;
						if((d & 0x3) == 3) d &= 0x4;
						push_u4v(bts, d);
						mat = 0;
					}
					hit->qe += l;
					hit->te += l;
					ol -= l;
					mat += l;
					hit->mat += l;
				}
			} else {
				hit->te += ol * (op & 0b01);
				hit->qe += ol * ((op >> 1) & 0b01);
			}
			ol = 0;
		}
		hit->qb = hit->qb / KBM_BIN_SIZE;
		hit->tb = hit->qb / KBM_BIN_SIZE;
		// trim bin cigars right
		while(bts->size && (bts->buffer[bts->size - 1] & 0x04)) bts->size --;
		// trim bin cigars left
		for(i=0;i<bts->size;i++){
			op = bts->buffer[i];
			if(op & 0x04){
				if(op & 0x1) hit->qb ++;
				if(op & 0x2) hit->tb ++;
			} else {
				break;
			}
		}
		// tidy bin cigars
		hit->qe = hit->qb;
		hit->te = hit->tb;
		clear_bitsvec(cigars);
		hit->cgoff = cigars->size;
		op = 0;
		for(;i<bts->size;i++){
			if(bts->buffer[i] & 0x1) hit->qe ++;
			if(bts->buffer[i] & 0x2) hit->te ++;
			if(bts->buffer[i] & 0x4) hit->gap ++;
			if(bts->buffer[i] == 0 || bts->buffer[i] & 0x4){
				if(op) push_bitsvec(cigars, op);
				push_bitsvec(cigars, bts->buffer[i]);
				op = 0;
			} else if(op){
				if(op == bts->buffer[i]){
					push_bitsvec(cigars, op);
				} else {
					push_bitsvec(cigars, 0);
					op = 0;
				}
			} else {
				op = bts->buffer[i];
			}
		}
		if(op) push_bitsvec(cigars, op);
		hit->cglen = cigars->size - hit->cgoff;
		hit->aln = num_min(hit->qe - hit->qb, hit->te - hit->tb);
		if(hit->mat < g->par->min_mat) continue;
		if(hit->aln < g->par->min_aln) continue;
		if(hit->mat < hit->aln * KBM_BIN_SIZE * g->par->min_sim) continue;
		if(num_diff(hit->qe - hit->qb, hit->te - hit->tb) > (int)num_max(g->par->aln_var * hit->aln, 1.0)) continue;
		kbmap2regs_graph(g, regs, qlen / KBM_BIN_SIZE, hit, cigars, maps);
		nhit ++;
	}
	fprintf(KBM_LOGF, "\r%llu lines, %llu hits\n", sam->n_line, nhit);
	free_bitsvec(cigars);
	free_u4v(cgs);
	free_u4v(bts);
	return nhit;
}

static inline u8i proc_alignments_core(Graph *g, int ncpu, rbregv *regs, u4v *maps[3], char *prefix, char *dump_kbm){
	kbm_map_t *hit;
	kbm_read_t *pb;
	BitVec *rdflags;
	BufferedWriter *bw;
	FILE *alno, *kmlog;
	u8i nbp, nhit;
	u8i i, ib, ie, ic;
	u4i rid, qb, qe, ii, in;
	int reset_kbm, n_cpu;
	thread_prepare(mdbg);
	if(KBM_LOG) n_cpu = 1;
	else n_cpu = ncpu;
	ic = (g->kbm->bins->size + g->num_index - 1) / g->num_index;
	ie = 0;
	alno = open_file_for_write(prefix, ".alignments.gz", 1);
	bw = zopen_bufferedwriter(alno, 1024 * 1024, ncpu, 0);
	rdflags = (g->par->skip_contained)? init_bitvec(g->kbm->reads->size) : NULL;
	thread_beg_init(mdbg, n_cpu);
	mdbg->g = g;
	mdbg->aux = init_kbmaux(g->kbm);
	if(g->rpar){
		mdbg->raux = init_kbmaux(init_kbm(g->rpar));
	} else {
		mdbg->raux = NULL;
	}
	mdbg->aux->par = (KBMPar*)malloc(sizeof(KBMPar));
	memcpy(mdbg->aux->par, g->par, sizeof(KBMPar));
	mdbg->rid = MAX_U4;
	mdbg->beg = 0;
	mdbg->end = 0;
	mdbg->alno = alno;
	thread_end_init(mdbg);
	in = g->num_index;
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
		ie = ib + ic;
		while(ie > 0 && ie < g->kbm->bins->size && g->kbm->bins->buffer[ie - 1].ridx == g->kbm->bins->buffer[ie].ridx) ie ++;
		qb = 0;
		qe = ie? g->kbm->bins->buffer[ie - 1].ridx : 0;
		nbp = ((u8i)(ie - ib)) * KBM_BSIZE;
		if(reset_kbm){
			reset_index_kbm(g->kbm);
			fprintf(KBM_LOGF, "[%s] indexing bins[%llu,%llu] (%llu bp), %d threads\n", date(), ib, ie, nbp, ncpu); fflush(KBM_LOGF);
			kmlog = (in > 1)? NULL : open_file_for_write(prefix, ".kmerdep", 1);
			index_kbm(g->kbm, ib, ie, ncpu, kmlog);
			if(kmlog){
				fclose(kmlog);
				kmlog = NULL;
			}
			fprintf(KBM_LOGF, "[%s] Done\n", date()); fflush(KBM_LOGF);
			if(in == 1 && dump_kbm){
				FILE *dump;
				fprintf(KBM_LOGF, "[%s] dump kbm index to %s ...", date(), dump_kbm); fflush(KBM_LOGF);
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
			}
		}
		{
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
				if(mdbg->rid != MAX_U4){
					KBMAux *aux = mdbg->aux;
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
						if(hit->mat == 0) continue;
						if(rdflags
							&& g->kbm->reads->buffer[hit->tidx].bincnt < g->kbm->reads->buffer[hit->qidx].bincnt
							&& (hit->tb <= 1 && hit->te + 1 >= (int)(g->kbm->reads->buffer[hit->tidx].bincnt))
							&& (hit->qb > 1 || hit->qe + 1 < (int)(g->kbm->reads->buffer[hit->qidx].bincnt))
							){
							one_bitvec(rdflags, hit->tidx);
						}
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
						//append_bitsvec(g->cigars, aux->cigars, hit->cgoff, hit->cglen);
						//hit->cgoff = g->cigars->size - hit->cglen;
						kbmap2regs_graph(g, regs, 0, hit, mdbg->aux->cigars, maps);
					}
					if(KBM_LOG){
						fprintf(KBM_LOGF, "QUERY: %s\t+\t%d\n", g->kbm->reads->buffer[mdbg->rid].tag, g->kbm->reads->buffer[mdbg->rid].bincnt * KBM_BIN_SIZE);
						for(i=0;i<mdbg->aux->hits->size;i++){
							hit = ref_kbmmapv(mdbg->aux->hits, i);
								fprintf(KBM_LOGF, "\t%s\t%c\t%d\t%d\t%d\t%d\t%d\n", g->kbm->reads->buffer[hit->tidx].tag, "+-"[hit->qdir], g->kbm->reads->buffer[hit->tidx].bincnt * KBM_BIN_SIZE, hit->tb * KBM_BIN_SIZE, hit->te * KBM_BIN_SIZE, hit->aln * KBM_BIN_SIZE, hit->mat);
						}
					}
					mdbg->rid = MAX_U4;
				}
				if(rid < qe && (rdflags == NULL || get_bitvec(rdflags, rid) == 0)){
					mdbg->rid = rid;
					thread_wake(mdbg);
				}
			}
		}
		if(!KBM_LOG) fprintf(KBM_LOGF, "\r%u reads|total hits %llu\n", qe - qb, nhit);
		if(reset_kbm){
			reset_index_kbm(g->kbm);
		}
	}
	thread_beg_close(mdbg);
	free(mdbg->aux->par);
	free_kbmaux(mdbg->aux);
	if(mdbg->raux){
		free_kbm(mdbg->raux->kbm);
		free_kbmaux(mdbg->raux);
	}
	thread_end_close(mdbg);
	if(bw) close_bufferedwriter(bw);
	if(alno) fclose(alno);
	if(rdflags) free_bitvec(rdflags);
	return nhit;
}

static inline void kbcovs_read_rdregs(Graph *g, rbregv *regs, u2i *kbcovs, int all){
	u8i i;
	rb_reg_t *r;
	UNUSED(g);
	for(i=0;i<regs->size;i++){
		r = ref_rbregv(regs, i);
		if(!all && r->closed) continue;
		if(kbcovs[r->bidx] < MAX_U2) kbcovs[r->bidx] ++;
	}
}

static inline u8i mask_grouped_minor_kbins_by_kbcovs(Graph *g, rbregv *regs, u2i *kbcovs, BitVec *kbits){
	rb_reg_t *r;
	u8i rgidx, bidx, ret;
	ret = 0;
	for(rgidx=0;rgidx<regs->size;rgidx++){
		r = ref_rbregv(regs, rgidx);
		bidx = g->kbm->reads->buffer[r->rid].binoff + r->off;
		if(get_bitvec(kbits, bidx)){
		} else if(kbcovs[bidx] > kbcovs[r->bidx]){
		} else if(kbcovs[bidx] == kbcovs[r->bidx] && r->bidx > bidx){
		} else {
			one_bitvec(kbits, bidx);
			ret ++;
		}
	}
	return ret;
}

static inline void print_read_rbregs(FILE *out, Graph *g, u4i rid, rbregv *regs, u2i *kbcovs, int all){
	read_t *rd;
	rb_reg_t *r;
	u4i i;
	rd = ref_readv(g->reads, rid);
	fprintf(out, "%s", g->kbm->reads->buffer[rid].tag);
	for(i=0;i<rd->regs.cnt;i++){
		r = ref_rbregv(regs, rd->regs.idx + i);
		if(r->closed && !all) continue;
		if(r->closed == 0){
			fprintf(out, "\t\e[1;32mK%llu(%u)_%c_%u_%u\e[0m", (u8i)r->bidx, kbcovs[r->bidx], "FR"[r->dir], (u4i)r->off, (u4i)r->len);
		} else {
			fprintf(out, "\tK%llu(%u)_%c_%u_%u*", (u8i)r->bidx, kbcovs[r->bidx], "FR"[r->dir], (u4i)r->off, (u4i)r->len);
		}
	}
	fprintf(out, "\n");
}

static inline void print_reads_rbregs(FILE *out, Graph *g, u4i rbeg, u4i rend, rbregv *regs, u2i *kbcovs, int all){
	u4i rid;
	for(rid=rbeg;rid<rend;rid++){
		print_read_rbregs(out, g, rid, regs, kbcovs, all);
	}
}

void print_avg_kbcovs(Graph *g, u2i *kbcovs, u2i *kboris){
	u8i i, cnt, sum, ori;
	cnt = sum = ori = 0;
	for(i=0;i<g->kbm->bins->size;i++){
		if(kbcovs[i] >= g->min_node_cov){
			cnt ++;
			sum += kbcovs[i];
			ori += kboris[i];
		}
	}
	printf("KBCOV: (%llu, %llu) / %llu = (%0.3f, %0.3f)\n", sum, ori, cnt, 1.0 * sum / cnt, 1.0 * ori / cnt);
}

static inline u4i del_dup_bidx_read_rbregs(Graph *g, rbregv *regs, u4i rid, u4v *idxs){
	read_t *rd;
	rb_reg_t *r1, *r2;
	u4i i, ret;
	rd = ref_readv(g->reads, rid);
	if(rd->regs.cnt == 0) return 0;
	if(0){
		r1 = ref_rbregv(regs, rd->regs.idx + 0);
		for(i=1;i<rd->regs.cnt;i++){
			r2 = ref_rbregv(regs, rd->regs.idx + i);
			if(r1->off > r2->off){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
			if(r2->closed){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
	clear_u4v(idxs);
	ret = 0;
	for(i=0;i<rd->regs.cnt;i++){
		push_u4v(idxs, i);
	}
	sort_array(idxs->buffer, idxs->size, u4i, num_cmpgt(regs->buffer[rd->regs.idx + a].bidx, regs->buffer[rd->regs.idx + b].bidx));
	r1 = ref_rbregv(regs, rd->regs.idx + idxs->buffer[0]);
	for(i=1;i<idxs->size;i++){
		r2 = ref_rbregv(regs, rd->regs.idx + idxs->buffer[i]);
		if(r1->bidx == r2->bidx){
			if(r1->closed == 0){
				r1->closed = 1;
				ret ++;
			}
			r2->closed = 1;
			ret ++;
		}
		r1 = r2;
	}
	sort_array(regs->buffer + rd->regs.idx, rd->regs.cnt, rb_reg_t, num_cmpgtx(a.closed, b.closed, a.off, b.off));
	rd->regs.cnt -= ret;
	if(0){
		r1 = ref_rbregv(regs, rd->regs.idx + 0);
		for(i=1;i<rd->regs.cnt;i++){
			r2 = ref_rbregv(regs, rd->regs.idx + i);
			if(r1->off > r2->off){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
		}
	}
	return ret;
}

static inline void affact_array_rbregs_graph_core(Graph *g, rbregv *regs, ptr_ref_t *kb, ptr_ref_t *rs){
	read_t *rd;
	rb_reg_t *r1, *r2;
	u4i j, m, beg, end;
	int ovl;
	r1 = ref_rbregv(regs, kb->idx);
	m = kb->cnt;
	rs->idx = 0;
	rs->cnt = 0;
	// TODO: if exceeds WT_MAX_RDREG_AUX, extra codes need to locate rd_reg_t from reg_t
	if(m > WT_MAX_RDREG_AUX) return;
	rd = ref_readv(g->reads, r1->rid);
	for(j=m+1;j<rd->regs.cnt;j++){
		r2 = ref_rbregv(regs, rd->regs.idx + j);
		ovl = Int(r1->off + r1->len) - Int(r2->off);
		if(ovl < Int(g->regovl)) break;
	}
	end = j;
	for(j=m;j>0;j--){
		r2 = ref_rbregv(regs, rd->regs.idx + j - 1);
		ovl = Int(r2->off + g->reglen) - Int(r1->off);
		//if(r1->off - r2->off > Int(g->reglen + 1)) break;
		if(ovl < Int(g->regovl)) break;
	}
	beg = j;
	rs->idx = rd->regs.idx + beg;
	rs->cnt = end - beg;
}

static inline int check_conflict_rbreg_graph_core(Graph *g, rbregv *regs, ptr_ref_t *kb){
	read_t *rd;
	rb_reg_t *r1, *r2;
	u4i j, m;
	int ovl;
	r1 = ref_rbregv(regs, kb->idx);
	m = kb->cnt;
	// TODO: if exceeds WT_MAX_RDREG_AUX, extra codes need to locate rd_reg_t from reg_t
	if(m > WT_MAX_RDREG_AUX) return 1;
	rd = ref_readv(g->reads, r1->rid);
	for(j=m+1;j<rd->regs.cnt;j++){
		r2 = ref_rbregv(regs, rd->regs.idx + j);
		ovl = Int(r1->off + r1->len) - Int(r2->off);
		if(ovl <= Int(g->regovl)) break;
		if(r2->closed) continue;
		return 1;
	}
	for(j=m;j>0;j--){
		r2 = ref_rbregv(regs, rd->regs.idx + j - 1);
		ovl = Int(r2->off + g->reglen) - Int(r1->off);
		if(ovl <= Int(g->regovl)) break;
		//ovl = Int(r2->off + r2->len) - Int(r1->off);
		//if(r1->off - r2->off > Int(g->reglen + 1)) break;
		//if(ovl <= Int(g->regovl)) continue;
		if(r2->closed) continue;
		return 1;
	}
	return 0;
}

static inline u4i update_rbregs_graph(Graph *g, rbregv *regs, ptr_ref_t *kbins, u8i off, u4i len, u4v *pass){
	rb_reg_t *r1, *r2;
	u4i i;
	clear_u4v(pass);
	for(i=0;i<len;i++){
		r1 = ref_rbregv(regs, kbins[off + i].idx);
		if(check_conflict_rbreg_graph_core(g, regs, kbins + off + i)){
		} else {
			push_u4v(pass, i);
		}
	}
	if(pass->size < g->min_node_cov || pass->size < (u4i)(len * (1 - g->node_max_conflict))){
		return 0;
	} else {
		if(0){
			sort_array(pass->buffer, pass->size, u4i, num_cmpgt(regs->buffer[kbins[off + a].idx].rid, regs->buffer[kbins[off + b].idx].rid));
			for(i=1;i<pass->size;i++){
				r1 = ref_rbregv(regs, kbins[off + pass->buffer[i - 1]].idx);
				r2 = ref_rbregv(regs, kbins[off + pass->buffer[i]].idx);
				if(r1->rid == r2->rid){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
			}
		}
		return pass->size;
	}
}

thread_beg_def(mupd);
Graph *g;
rbregv *regs;
ptr_ref_t *kbins, *ksrts;
u8i off;
u4i len;
ptrrefv *rs;
thread_end_def(mupd);

thread_beg_func(mupd);
Graph *g;
rbregv *regs;
ptr_ref_t *kbins, *ksrts, *ks;
ptrrefv *rs;
u8i kboff;
u4i i, j;
g = mupd->g;
regs = mupd->regs;
kbins = mupd->kbins;
ksrts = mupd->ksrts;
rs = mupd->rs;
thread_beg_loop(mupd);
kboff = ksrts[mupd->off].idx;
clear_ptrrefv(rs);
for(i=0;i<mupd->len;i++){
	ks = ksrts + mupd->off + i;
	for(j=0;j<ks->cnt;j++){
		affact_array_rbregs_graph_core(g, regs, kbins + ks->idx + j, next_ref_ptrrefv(rs));
	}
}
thread_end_loop(mupd);
thread_end_func(mupd);

static inline u8i mul_update_rbregs_graph(Graph *g, rbregv *regs, ptr_ref_t *kbins, ptr_ref_t *ksrts, u8i kscnt, u4i ncpu, u4i batch_size){
	BitVec *regbits;
	ptrrefv *rs;
	u4v *pass;
	node_t *v;
	ptr_ref_t *ks, *rb;
	rb_reg_t *r;
	u8i ksidx, idx, nreg;
	u4i j, m, len, num;
	thread_prepare(mupd);
	regbits = init_bitvec(regs->size);
	thread_beg_init(mupd, ncpu);
	mupd->g = g;
	mupd->regs = regs;
	mupd->kbins = kbins;
	mupd->ksrts = ksrts;
	mupd->rs = init_ptrrefv(32);
	mupd->off = 0;
	mupd->len = 0;
	thread_end_init(mupd);
	nreg = 0;
	clear_nodev(g->nodes);
	pass = init_u4v(1024);
	for(ksidx=0;ksidx<kscnt+ncpu;){
		thread_wait_next(mupd);
		if(ksidx < kscnt){
			len = num_min(batch_size, kscnt - ksidx);
			ksidx += len;
		} else {
			len = 0;
			ksidx ++;
		}
		rs = mupd->rs;
		num = 0;
		for(m=0;m<mupd->len;m++){
			ks = ksrts + mupd->off + m;
			clear_u4v(pass);
			for(j=0;j<ks->cnt;j++){
				rb = ref_ptrrefv(rs, num ++);
				if(rb->cnt == 0) continue;
				if(kbins[ks->idx + j].idx < rb->idx || kbins[ks->idx + j].idx >= rb->idx + rb->cnt){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				if(reg_count_bitvec(regbits, rb->idx, rb->idx + rb->cnt) == 0){
					push_u4v(pass, j);
				}
			}
			if(pass->size < g->min_node_cov || pass->size < (u4i)(ks->cnt * (1 - g->node_max_conflict))){
			} else {
				for(j=0;j<pass->size;j++){
					idx = kbins[ks->idx + pass->buffer[j]].idx;
					one_bitvec(regbits, idx);
					r = ref_rbregv(regs, idx);
					r->bidx = g->nodes->size; // translate bidx into node_idx
					r->closed = 0;
				}
				v = next_ref_nodev(g->nodes);
				ZEROS(v);
				v->regs.idx = nreg;
				v->regs.cnt = 0;
				nreg += pass->size;
			}
		}
		mupd->off = ksidx - len;
		mupd->len = len;
		if(len){
			thread_wake(mupd);
		}
	}
	free_bitvec(regbits);
	free_u4v(pass);
	thread_beg_close(mupd);
	free_ptrrefv(mupd->rs);
	thread_end_close(mupd);
	return nreg;
}

static inline void build_nodes_graph(Graph *g, u8i maxbp, int ncpu, FileReader *pws, char *prefix, char *dump_kbm){
	kbm_bin_t *kb;
	read_t *rd;
	rb_reg_t *r;
	rd_reg_t *rg;
	node_t *v;
	reg_t *reg;
	ptr_ref_t *kbins, *ksrts;
	u4v *maps[3], *pass;
	u8i idx, bidx, bcnt, cnt, nreg;
	u4i rid, i, num;
	u2i *kbcovs, *kboris;
	rbregv *regs;
	regs = init_rbregv(1024);
	maps[0] = init_u4v(4);
	maps[1] = init_u4v(4);
	maps[2] = init_u4v(4);
	kbcovs = calloc(g->kbm->bins->size, 2);
	if(pws){
		cnt = load_alignments_core(g, pws, regs, maps);
	} else {
		UNUSED(maxbp);
		cnt = proc_alignments_core(g, ncpu, regs, maps, prefix, dump_kbm);
	}
	free_u4v(maps[0]);
	free_u4v(maps[1]);
	free_u4v(maps[2]);
	print_proc_stat_info(0);
	kbcovs_read_rdregs(g, regs, kbcovs, 1);
	for(idx=0;idx<g->kbm->bins->size;idx++) kbcovs[idx] ++;
	kboris = malloc(g->kbm->bins->size * 2);
	memcpy(kboris, kbcovs, g->kbm->bins->size * 2);
	print_avg_kbcovs(g, kbcovs, kboris);
	if(1){
		BitVec *kbits;
		kbits  = init_bitvec(g->kbm->bins->size);
		mask_grouped_minor_kbins_by_kbcovs(g, regs, kbcovs, kbits);
		for(idx=cnt=0;idx<regs->size;idx++){
			r = ref_rbregv(regs, idx);
			if(get_bitvec(kbits, r->bidx) || kbcovs[r->bidx] < g->min_node_cov){
				r->closed = 1;
				cnt ++;
			}
		}
		for(idx=0;idx<g->kbm->bins->size;idx++){
			if(get_bitvec(kbits, idx) || kbcovs[idx] < g->min_node_cov) kbcovs[idx] = 0;
		}
		psort_array(regs->buffer, regs->size, rb_reg_t, ncpu, num_cmpgt(a.closed, b.closed));
		regs->size -= cnt;
		recap_rbregv(regs, regs->size);
		free_bitvec(kbits);
		fprintf(KBM_LOGF, "[%s] masked %llu minor kbins from %llu\n", date(), cnt, regs->size + cnt); fflush(KBM_LOGF);
		print_avg_kbcovs(g, kbcovs, kboris);
	}
	cnt = 0;
	for(idx=0;idx<g->kbm->bins->size;idx++){
		if(kbcovs[idx] >= g->min_node_cov){
			kb = ref_kbmbinv(g->kbm->bins, idx);
			push_rbregv(regs, (rb_reg_t){idx, 0, kb->off, g->reglen, 0, kb->ridx});
			cnt ++;
		}
	}
	ksrts = malloc(cnt * sizeof(ptr_ref_t));
	fprintf(KBM_LOGF, "[%s] sorting %llu/%llu kbins ... ", date(), regs->size, cnt); fflush(KBM_LOGF);
	psort_array(regs->buffer, regs->size, rb_reg_t, ncpu, num_cmpgtx(a.rid, b.rid, a.off, b.off));
	kbins = malloc(regs->size * sizeof(ptr_ref_t));
	rid = MAX_U4;
	for(idx=0;idx<regs->size;idx++){
		r = ref_rbregv(regs, idx);
		rd = ref_readv(g->reads, r->rid);
		if(r->rid != rid){
			rid = r->rid;
			rd->regs.idx = idx;
			rd->regs.cnt = 1;
		} else {
			rd->regs.cnt ++;
		}
	}
	thread_fast_run(kdup, ncpu, THREAD_EXPR(
		u4i i;
		u4v *idxs = init_u4v(32);
		for(i=TIDX;i<g->reads->size;i+=NCPU){
			del_dup_bidx_read_rbregs(g, regs, i, idxs);
		}
		free_u4v(idxs);
	));
	rid = 0;
	bcnt = cnt = 0;
	for(idx=0;idx<regs->size;idx++){
		r = ref_rbregv(regs, idx);
		if(r->closed) continue;
		r->closed = 1;
		rd = ref_readv(g->reads, r->rid);
		if(r->rid != rid){
			rid = r->rid;
			cnt = 0;
		} else {
			cnt ++;
		}
		kbins[bcnt].idx = idx;
		kbins[bcnt].cnt = cnt;
		bcnt ++;
	}
	psort_array(kbins, bcnt, ptr_ref_t, ncpu, num_cmpgtx(regs->buffer[a.idx].bidx, regs->buffer[b.idx].bidx, a.idx, b.idx));
	bidx = MAX_U8;
	cnt = MAX_U8;
	for(idx=0;idx<bcnt;idx++){
		r = ref_rbregv(regs, kbins[idx].idx);
		if(bidx != r->bidx){
			if(!(bidx != MAX_U8 && kbcovs[bidx] < g->min_node_cov)){
				cnt ++;
			}
			bidx = r->bidx;
			ksrts[cnt].idx = idx;
			ksrts[cnt].cnt = 1;
		} else {
			ksrts[cnt].cnt ++;
		}
	}
	if(!(bidx != MAX_U8 && kbcovs[bidx] < g->min_node_cov)){
		cnt ++;
	}
	psort_array(ksrts, cnt, ptr_ref_t, ncpu, num_cmpgtx(b.cnt, a.cnt, a.idx, b.idx));
	free(kbcovs);
	free(kboris);
	fprintf(KBM_LOGF, " Done\n"); fflush(KBM_LOGF);
	fprintf(KBM_LOGF, "[%s] generating nodes ... ", date()); fflush(KBM_LOGF);
	if(1){
		nreg = mul_update_rbregs_graph(g, regs, kbins, ksrts, cnt, ncpu, 100);
	} else {
		nreg = 0;
		pass = init_u4v(32);
		clear_nodev(g->nodes);
		for(idx=0;idx<cnt;idx++){
			num = update_rbregs_graph(g, regs, kbins, ksrts[idx].idx, ksrts[idx].cnt, pass);
			if(num){
				for(i=0;i<pass->size;i++){
					r = ref_rbregv(regs, kbins[ksrts[idx].idx + pass->buffer[i]].idx);
					r->bidx = g->nodes->size; // translate bidx into node_idx
					r->closed = 0;
				}
				v = next_ref_nodev(g->nodes);
				ZEROS(v);
				v->regs.idx = nreg;
				v->regs.cnt = 0;
				nreg += num;
			}
		}
		free_u4v(pass);
	}
	free(ksrts);
	free(kbins);
	fprintf(KBM_LOGF, " Done, %llu/%llu nodes\n", nreg, g->nodes->size); fflush(KBM_LOGF);
	fprintf(KBM_LOGF, "[%s] generating read paths ...", date());
	clear_rdregv(g->rdregs);
	for(rid=0;rid<g->reads->size;rid++){
		rd = ref_readv(g->reads, rid);
		idx = rd->regs.idx;
		rd->regs.idx = g->rdregs->size;
		for(i=0;i<rd->regs.cnt;i++){
			r = ref_rbregv(regs, idx + i);
			if(r->closed) continue;
			rg = next_ref_rdregv(g->rdregs);
			rg->node = r->bidx;
			rg->dir  = r->dir;
			rg->off  = r->off;
			rg->len  = r->len;
			rg->closed = 0;
		}
		rd->regs.cnt = g->rdregs->size - rd->regs.idx;
#if DEBUG
		if(rd->regs.cnt > WT_MAX_RDREG_AUX){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
#endif
	}
	free_rbregv(regs);
	fprintf(KBM_LOGF, " Done, %llu\n", g->rdregs->size); fflush(KBM_LOGF);
	fprintf(KBM_LOGF, "[%s] generating nodes regs ...", date());
	renew_regv(g->regs, g->rdregs->size);
	g->regs->size = g->rdregs->size;
	for(rid=0;rid<g->reads->size;rid++){
		rd = ref_readv(g->reads, rid);
		for(i=0;i<rd->regs.cnt;i++){
			rg = ref_rdregv(g->rdregs, rd->regs.idx + i);
			v = ref_nodev(g->nodes, rg->node);
			reg = ref_regv(g->regs, v->regs.idx + v->regs.cnt);
			reg->rid = rid;
			reg->dir = rg->dir;
			reg->off = rg->off;
			reg->len = rg->len;
			reg->closed = 0;
			reg->rdreg_aux = i;
			v->regs.cnt ++;
		}
	}
	if(0){
		for(idx=0;idx<g->nodes->size;idx++){
			v = ref_nodev(g->nodes, idx);
			for(i=0;i<v->regs.cnt;i++){
				reg = ref_regv(g->regs, v->regs.idx + i);
				rd  = ref_readv(g->reads, reg->rid);
				rg  = ref_rdregv(g->rdregs, rd->regs.idx + reg->rdreg_aux);
				if(rg->node != idx){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
				if(rg->off != reg->off){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					abort();
				}
			}
		}
	}
	fprintf(KBM_LOGF, " Done, %llu\n", g->regs->size); fflush(KBM_LOGF);
	print_proc_stat_info(0);
}

static inline u4i estimate_genome_size(Graph *g, unsigned long long tot_bp, FILE *out){
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
	fprintf(out, "[%s] median node depth = %d, average node depth = %0.1f\n", date(), mid, avg);
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
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		n->erefs[0] = n->erefs[1] = EDGE_REF_NULL;
	}
}

typedef struct {
	uint64_t idx:46; int off:18;
} edge_off_t;
define_list(edgeoffv, edge_off_t);

static inline int estimate_edge_length(edge_off_t *ps, uint32_t size, uint32_t idxs[2]){
	int64_t tot;
	uint32_t i, b, e, mi;
	int max, len, avg;
	 idxs[0] = 0; idxs[1] = size;
	if(size == 0){ return 0; }
	if(size <= 2){ return ps[size/2].off; }
	b = 0; e = size;
	for(i=b,tot=0;i<e;i++) tot += ps[i].off;
	len = tot / (e - b);
	while(b + 2 < e){
		max = 0; mi = 0;
		for(i=b+1;i<e;i++){
			if(ps[i].off - ps[i-1].off > max){
				max = ps[i].off - ps[i-1].off;
				mi = i;
			}
		}
		if(max < len * 0.5) break;
		else if(max < 2) break;
		if(mi - b > e - mi) e = mi;
		else b = mi;
		for(i=b,tot=0;i<e;i++) tot += ps[i].off;
		avg = tot / (e - b);
		if(num_diff(avg, len) < num_max(avg * 0.2, 1)) break;
		len = avg;
	}
	return len;
}

// s = src, d = dst
#define get_edge_sidx(e, k) ((u8i)((e)->es[k].node))
#define set_edge_sidx(e, k, v) (e)->es[k].node = (v)
#define get_edge_didx(e, k) ((u8i)(e)->es[!(k)].node)
#define set_edge_didx(e, k, v) (e)->es[!(k)].node = (v)
#define get_edge_sdir(e, k) (e)->es[k].dir
#define set_edge_sdir(e, k, v) (e)->es[k].dir = (v)
#define get_edge_ddir(e, k) !(e)->es[!(k)].dir
#define set_edge_ddir(e, k, v) (e)->es[!(k)].dir = !(v)
#define is_edge_closed(e) (e)->es[0].closed
#define get_edge_closed(e) is_edge_closed(e)
#define set_edge_closed(e, v) (e)->es[0].closed = (v)
#define get_edge_hard(e) (e)->es[1].closed
#define set_edge_hard(e, v) (e)->es[1].closed = (v)
#define get_edge_cov(e) (e)->es[0].cov
#define set_edge_cov(e, v) (e)->es[0].cov = (v)
#define inc_edge_cov(e, v) ((e)->es[0].cov = ((e)->es[0].cov + (v) < WT_MAX_EDGE_COV)? (e)->es[0].cov + (v) : WT_MAX_EDGE_COV)
#define dec_edge_cov(e, v) ((e)->es[0].cov = ((e)->es[0].cov < (v))? 0 : (e)->es[0].cov - (v))
#define get_edge_off(e) (e)->es[0].off
#define set_edge_off(e, v) (e)->es[0].off = (v)
#define get_edge_status(e) (e)->es[0].status
#define set_edge_status(e) (e)->es[0].status = (v)
#define get_edge_flag(e) (e)->es[0].flag
#define set_edge_flag(e, v) (e)->es[0].flag = (v)
// set edge be keeped, suggest not closed
#define get_edge_keep(e) (e)->es[1].closed
#define set_edge_keep(e, v) (e)->es[1].closed = (v)
// edge port, 0: default, > 0: cannot go in other port than 0
#define get_edge_port(e) ((((u8i)(e)->es[0].port) << 20) | (e)->es[1].port)
#define set_edge_port(e, v) { (e)->es[0].port = (v) >> 20; (e)->es[1].port = (v) & 0xFFFFFU; }
#define get_reg_port(reg) ((((u8i)(reg)->port1) << 16) | (reg)->port2)
#define set_reg_port(reg, v) { (reg)->port1 = (v) >> 16; (reg)->port2 = (v) & 0xFFFFU; }
#define get_reg_closed(reg) (reg)->closed
#define set_reg_closed(g, reg, v) { (reg)->closed = (v); (g)->rdregs->buffer[(g)->reads->buffer[(reg)->rid].regs.idx + (reg)->rdreg_aux].closed = (v); }

static inline u4i cal_edge_cov_graph_core(Graph *g, u8i node1, u4i dir1, u8i node2, u4i dir2, int all){
	node_t *vs[2];
	reg_t *rs[2], *re[2];
	u4i cov, k;
	vs[0] = ref_nodev(g->nodes, node1);
	vs[1] = ref_nodev(g->nodes, node2);
	if(vs[0]->regs.cnt == 0 || vs[1]->regs.cnt == 0) return 0;
	rs[0] = ref_regv(g->regs, vs[0]->regs.idx);
	rs[1] = ref_regv(g->regs, vs[1]->regs.idx);
	re[0] = rs[0] + vs[0]->regs.cnt;
	re[1] = rs[1] + vs[1]->regs.cnt;
	cov = 0;
	while(1){
		if(rs[0]->rid < rs[1]->rid){
			rs[0] ++;
			if(rs[0] >= re[0]) break;
		} else if(rs[0]->rid > rs[1]->rid){
			rs[1] ++;
			if(rs[1] >= re[1]) break;
		} else {
			if(all || (rs[0]->closed == 0 && rs[1]->closed == 0)){
				k = rs[0]->off > rs[1]->off;
				if(dir1 == (rs[0]->dir ^ k) && dir2 == (rs[1]->dir ^ k)){
					cov ++;
				}
			}
			rs[0] ++;
			if(rs[0] >= re[0]) break;
			rs[1] ++;
			if(rs[1] >= re[1]) break;
		}
	}
	return cov;
}

#define cal_edge_cov_graph(g, e) cal_edge_cov_graph_core(g, get_edge_sidx(e, 0), get_edge_sdir(e, 0), get_edge_didx(e, 0), get_edge_ddir(e, 0), 0)

static inline void add_edge_graph_core(nodev *nodes, edgev *edges, edge_t *e){
	node_t *n;
	int k;
	for(k=0;k<2;k++){
		n = ref_nodev(nodes, e->es[k].node);
		e->es[k].next = n->erefs[e->es[k].dir].idx;
		e->es[k].flg  = n->erefs[e->es[k].dir].flg;
		n->erefs[e->es[k].dir].idx = offset_edgev(edges, e);
		n->erefs[e->es[k].dir].flg = k;
		if(!is_edge_closed(e)){
			n->erefs[e->es[k].dir].cnt ++;
		}
	}
}

#define add_edge_graph(g, e) add_edge_graph_core((g)->nodes, (g)->edges, e)

static inline void cut_edge_graph_core(nodev *nodes, edgev *edges, edge_t *e, int closed_val){
	UNUSED(edges);
	if(!e->es[0].closed){
		ref_nodev(nodes, e->es[0].node)->erefs[e->es[0].dir].cnt --;
		ref_nodev(nodes, e->es[1].node)->erefs[e->es[1].dir].cnt --;
	}
	e->es[0].closed = closed_val;
}

#define cut_edge_graph(g, e) cut_edge_graph_core((g)->nodes, (g)->edges, e, WT_EDGE_CLOSED_MASK)

static inline void revive_edge_graph(Graph *g, edge_t *e){
	if(e->es[0].closed == WT_EDGE_CLOSED_NULL) return;
	e->es[0].closed = WT_EDGE_CLOSED_NULL;
	ref_nodev(g->nodes, e->es[0].node)->erefs[e->es[0].dir].cnt ++;
	ref_nodev(g->nodes, e->es[1].node)->erefs[e->es[1].dir].cnt ++;
}

#define count_living_edges_graph(g, n, dir) (n)->erefs[dir].cnt

static inline u4i beg_iter_edges_graph(node_t *n, int dir, edge_ref_t *ref){
	ref->idx  = n->erefs[dir].idx;
	ref->flg2 = n->erefs[dir].flg;
	ref->flg  = ref->flg2;
	ref->cnt  = 0;
	return n->erefs[dir].cnt;
}

static inline edge_t* ref_iter_edges_graph_core(edgev *edges, edge_ref_t *ref, int all){
	edge_t *e;
	while(ref->idx){
		e = ref_edgev(edges, ref->idx);
		ref->flg  = ref->flg2;
		ref->idx  = e->es[ref->flg].next;
		ref->flg2 = e->es[ref->flg].flg;
		if(!all && is_edge_closed(e)){
			continue;
		}
		ref->cnt ++;
		return e;
	}
	return NULL;
}

#define ref_iter_edges_graph(g, ref) ref_iter_edges_graph_core((g)->edges, ref, 0)

// sort edges by off(ASC) and then cov(DSC)
static inline void srt_node_edges_graph(Graph *g, u8i nid, int dir){
	node_t *n;
	edge_ref_t F1, F2, F3;
	edge_t *e1, *e2, *e3;
	u4i cnt;
	int xchg;
	n = ref_nodev(g->nodes, nid);
	cnt = 0;
	do {
		cnt ++;
		xchg = 0;
		F1.idx = 0;
		F1.flg = 0;
		F1.cnt = 0;
		e1 = NULL;
		F2.idx = n->erefs[dir].idx;
		F2.flg = n->erefs[dir].flg;
		F2.cnt = 0;
		if(F2.idx == 0) break;
		e2 = ref_edgev(g->edges, F2.idx);
		while(e2->es[F2.flg].next){
			F3.idx = e2->es[F2.flg].next;
			F3.flg = e2->es[F2.flg].flg;
			e3 = ref_edgev(g->edges, F3.idx);
			if((get_edge_off(e2) > get_edge_off(e3)) || (get_edge_off(e2) == get_edge_off(e3) && get_edge_cov(e2) < get_edge_cov(e3))){
				xchg ++;
				if(F1.idx){
					e1 = ref_edgev(g->edges, F1.idx);
					e1->es[F1.flg].next = F3.idx;
					e1->es[F1.flg].flg  = F3.flg;
				} else {
					n->erefs[dir].idx = F3.idx;
					n->erefs[dir].flg = F3.flg;
				}
				e2->es[F2.flg].next = e3->es[F3.flg].next;
				e2->es[F2.flg].flg  = e3->es[F3.flg].flg;
				e3->es[F3.flg].next = F2.idx;
				e3->es[F3.flg].flg  = F2.flg;
				F1 = F3;
			} else {
				F1 = F2;
				F2 = F3;
				e2 = e3;
			}
		}
	} while(xchg);
}

static inline void del_node_edges_graph(Graph *g, node_t *n){
	edge_ref_t eref;
	edge_t *e;
	u4i k;
	for(k=0;k<2;k++){
		beg_iter_edges_graph(n, k, &eref);
		while((e = ref_iter_edges_graph(g, &eref))){
			cut_edge_graph_core(g->nodes, g->edges, e, WT_EDGE_CLOSED_HARD);
		}
	}
}

static inline void del_node_graph(Graph *g, node_t *n){
	del_node_edges_graph(g, n);
	n->closed = 1;
}

// dir = 2 means either strand
static inline edge_t* edge_node2node_graph(Graph *g, u8i node1, int dir1, u8i node2, int dir2, int *flg){
	node_t *n;
	edge_ref_t eref;
	edge_t *e;
	int dire;
	n = ref_nodev(g->nodes, node1);
	if(dir1 > 1){
		dir1 = 0; dire = 2;
	} else {
		dire = dir1 + 1;
	}
	while(dir1 < dire){
		beg_iter_edges_graph(n, dir1, &eref);
		while((e = ref_iter_edges_graph_core(g->edges, &eref, 1))){
			if(get_edge_didx(e, eref.flg) == node2){
				if(dir2 > 1) break;
				if(get_edge_ddir(e, eref.flg) == dir2) break;
			}
		}
		if(e){
			if(flg) *flg = eref.flg;
			return e;
		}
		dir1 ++;
	}
	return NULL;
}

static inline u8i build_edges_graph(Graph *g, int ncpu, FILE *log){
	read_t *rd;
	edge_t *E, *e;
	rd_reg_t *r1, *r2;
	edgeoffv *offs;
	u8i idx, lst, ret, *u;
	u4i rid, i, k, idxs[2];
	int exists, eoff;
	UNUSED(log);
	clear_edgev(g->edges);
	E = next_ref_edgev(g->edges);
	memset(E, 0, sizeof(edge_t));
	set_edge_cov(E, 1);
	set_edge_closed(E, WT_EDGE_CLOSED_MASK);
	//cut_edge_graph(g, e, WT_EDGE_CLOSED_MASK);
	//E->cov = 1;
	//E->status = 0;
	clear_edgehash(g->ehash);
	offs = init_edgeoffv(32);
	ret = 0;
	for(rid=0;rid<g->kbm->reads->size;rid++){
		rd = ref_readv(g->reads, rid);
		if(rd->regs.cnt < 2) continue;
		r1 = NULL;
		for(i=0;i<rd->regs.cnt;i++){
			r2 = ref_rdregv(g->rdregs, rd->regs.idx + i);
			if(r2->closed) continue;
			if(g->nodes->buffer[r2->node].closed) continue;
			if(r1){
				encap_edgev(g->edges, 1);
				E = ref_edgev(g->edges, 0);
				k = r1->node > r2->node;
				set_edge_sidx(E, k, r1->node);
				set_edge_didx(E, k, r2->node);
				set_edge_sdir(E, k, r1->dir);
				set_edge_ddir(E, k, r2->dir);
				eoff = Int(r2->off) - Int(r1->off + r1->len);
				set_edge_off(E, eoff);
				u = prepare_edgehash(g->ehash, 0, &exists);
				if(exists){
					e = ref_edgev(g->edges, *u);
					if(get_edge_cov(e) < WT_MAX_EDGE_COV) set_edge_cov(e, get_edge_cov(e) + 1);
				} else {
					*u = g->edges->size;
					e = next_ref_edgev(g->edges);
					*e = *E;
				}
				if(get_edge_cov(e) == g->min_edge_cov){
					set_edge_closed(e, 0);
					add_edge_graph(g, e);
					ret ++;
				}
				push_edgeoffv(offs, (edge_off_t){*u, eoff});
			}
			r1 = r2;
		}
	}
	if(g->edges->size == 0){
		free_edgeoffv(offs);
		return ret;
	}
	psort_array(offs->buffer, offs->size, edge_off_t, ncpu, num_cmpgtx(a.idx, b.idx, a.off, b.off));
	lst = 0;
	for(idx=1;idx<=offs->size;idx++){
		if(idx < offs->size && offs->buffer[idx].idx == offs->buffer[lst].idx) continue;
		if(1){
			set_edge_off(ref_edgev(g->edges, offs->buffer[lst].idx), offs->buffer[(lst+idx)/2].off);
			set_edge_cov(ref_edgev(g->edges, offs->buffer[lst].idx), idx - lst);
		} else {
			set_edge_off(ref_edgev(g->edges, offs->buffer[lst].idx), estimate_edge_length(offs->buffer + lst, idx - lst, idxs));
			set_edge_cov(ref_edgev(g->edges, offs->buffer[lst].idx), idxs[1] - idxs[0]);
		}
		lst = idx;
	}
	free_edgeoffv(offs);
	thread_fast_run(srt_edges, ncpu, EXPR(u8i ni; for(ni=TIDX;ni<g->nodes->size;ni+=NCPU){ srt_node_edges_graph(g, ni, 0); srt_node_edges_graph(g, ni, 1); }));
	return ret;
}

static inline void sort_read_rdregs_graph(Graph *g, u4i rid){
	read_t *rd;
	rd = ref_readv(g->reads, rid);
	sort_array(g->rdregs->buffer + rd->regs.idx, rd->regs.cnt, rd_reg_t, num_cmpgt(a.off, b.off));
}

static inline void load_nodes_graph(Graph *g, FileReader *nds, int ncpu){
	read_t *rd;
	node_t *n;
	rd_reg_t *rg;
	reg_t *reg;
	u8i nid, nreg;
	u4i i, rid;
	char *str, *tok;
	int ncol, closed, nwarn;
	clear_nodev(g->nodes);
	clear_regv(g->regs);
	nwarn = 0;
	while((ncol = readtable_filereader(nds)) != -1){
		if(nds->line->string[0] == '#') continue;
		if(ncol < 2) continue;
		nreg = atoi(get_col_str(nds, 1));
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
		n->regs.cnt = 0;
		for(i=0;i<nreg;i++){
			reg = next_ref_regv(g->regs);
			reg->closed = 0;
			str = get_col_str(nds, 2 + i);
			{
				if(str[get_col_len(nds, 2 + i)-1] == '*'){
					reg->closed = 1;
				}
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->len = atoi(tok);
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->off = atoi(tok);
				tok = rindex(str, '_'); *tok = '\0'; tok ++;
				reg->dir = (tok[0] == 'R');
				if((rid = getval_cuhash(g->kbm->tag2idx, str)) == MAX_U4){
					g->regs->size --;
					if(nwarn < 10){
						fprintf(stderr, " -- Cannot find read \"%s\" in LINE:%llu in %s -- %s:%d --\n", str, nds->n_line, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						nwarn ++;
					}
					continue;
				}
				reg->rid = rid;
			}
			rd = ref_readv(g->reads, reg->rid);
			if(rd->regs.cnt == WT_MAX_RDREG_AUX){
				g->regs->size --;
				continue;
			}
			reg->rdreg_aux = rd->regs.cnt ++;
			n->regs.cnt ++;
		}
	}
	nreg = 0;
	for(rid=0;rid<g->reads->size;rid++){
		rd = ref_readv(g->reads, rid);
		rd->regs.idx = nreg;
		nreg += rd->regs.cnt;
		rd->regs.cnt = 0;
	}
	renew_rdregv(g->rdregs, nreg);
	g->rdregs->size = nreg;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		for(i=0;i<n->regs.cnt;i++){
			reg = ref_regv(g->regs, n->regs.idx + i);
			rd = ref_readv(g->reads, reg->rid);
			rg = ref_rdregv(g->rdregs, rd->regs.idx + rd->regs.cnt);
			rg->node = nid;
			rg->dir  = reg->dir;
			rg->off  = reg->off;
			rg->len  = reg->len;
			rg->closed = 0;
#if DEBUG
			if(reg->rdreg_aux != rd->regs.cnt){
				fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				abort();
			}
#endif
			rd->regs.cnt ++;
		}
	}
	thread_fast_run(rdsrt, ncpu, THREAD_EXPR(u4i ri; for(ri=TIDX;ri<g->reads->size;ri+=NCPU){ sort_read_rdregs_graph(g, ri); }));
}

static inline u8i mask_nodes_by_edge_cov_graph(Graph *g, u4i min_node_cov, float min_edge_cov_ratio, FILE *out){
	node_t *n;
	edge_t *e;
	edge_ref_t eref;
	u8i nid, ret;
	u4i max, k;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		//if(n->regs.cnt < min_node_cov) continue;
		if(n->cov < min_node_cov) continue;
		max = 0;
		for(k=0;k<2;k++){
			beg_iter_edges_graph(n, k, &eref);
			while((e = ref_iter_edges_graph(g, &eref))){
				if(get_edge_cov(e) > max) max = get_edge_cov(e);
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

static inline void get_subgraph_nodes_graph(Graph *g, ptrrefhash *nodes, u8v *stack, uint16_t max_step, u4i closed_val){
	node_t *n;
	edge_ref_t eref;
	edge_t *e;
	ptr_ref_t *p, *pp;
	u8i nid;
	u4i k, cnt;
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
			beg_iter_edges_graph(n, k, &eref);
			while((e = ref_iter_edges_graph_core(g->edges, &eref, 1))){
				if(is_edge_closed(e) >= closed_val) continue;
				nid = get_edge_didx(e, eref.flg);
				pp = prepare_ptrrefhash(nodes, (ptr_ref_t){nid, 0}, &exists);
				if(exists) continue;
				pp->idx = nid; pp->cnt = cnt + 1;
				push_u8v(stack, nid);
			}
		}
	}
}

static inline u8i print_local_dot_graph(Graph *g, u8i local_dot_node, u4i local_dot_step, char *filename){
	FILE *out;
	ptrrefhash *hash;
	u8v *stack;
	ptr_ref_t *p;
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t eref;
	edge_t *e;
	u8i i;
	u4i j, k, max;
	out = open_file_for_write(filename, NULL, 1);
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
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"%s]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->off, r->len, p->idx == local_dot_node? " style=filled fillcolor=yellow" : "");
	}
	reset_iter_ptrrefhash(hash);
	while((p = ref_iter_ptrrefhash(hash))){
		i = p->idx;
		n = ref_nodev(g->nodes, i);
		for(k=0;k<2;k++){
			beg_iter_edges_graph(n, k, &eref);
			while((e = ref_iter_edges_graph(g, &eref))){
				fprintf(out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", (u8i)get_edge_sidx(e, eref.flg), (u8i)get_edge_didx(e, eref.flg), "+-"[get_edge_sdir(e, eref.flg)], "+-"[get_edge_ddir(e, eref.flg)], get_edge_cov(e), get_edge_off(e), colors[get_edge_sdir(e, eref.flg)][get_edge_ddir(e, eref.flg)]);
			}
		}
	}
	fprintf(out, "}\n");
	fclose(out);
	return 0;
}

static inline u8i mask_possible_tip_nodes_graph(Graph *g){
	node_t *n;
	reg_t *r;
	read_t *rd;
	u8i ret, i;
	u4i j, cnt;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		cnt = 0;
		for(j=0;j<n->regs.cnt;j++){
			r = ref_regv(g->regs, n->regs.idx + j);
			rd = ref_readv(g->reads, r->rid);
			if(r->rdreg_aux == 0 || r->rdreg_aux + 1 == rd->regs.cnt) continue;
			cnt ++;
		}
		if(cnt < g->min_edge_cov){
			n->closed = 1;
			ret ++;
		}
	}
	return ret;
}

static inline void print_node_edges_graph(Graph *g, u8i nid, int dir, FILE *out){
	node_t *n;
	edge_ref_t eref;
	edge_t *e;
	n = ref_nodev(g->nodes, nid);
	beg_iter_edges_graph(n, dir, &eref);
	while((e = ref_iter_edges_graph(g, &eref))){
		fprintf(out, "N%llu\t%c\tN%llu\t%c\t%d\t%d\n", (u8i)get_edge_sidx(e, eref.flg), "+-"[get_edge_sdir(e, eref.flg)], (u8i)get_edge_didx(e, eref.flg), "+-"[get_edge_ddir(e, eref.flg)], get_edge_cov(e), get_edge_off(e));
	}
}

static inline int inclusive_linear_trace_graph(Graph *g, trace_t *_t, int max_len, u4i max_step, u8i check_node_visit, tracev *path, int *msg){
	trace_t *t;
	node_t *n;
	edge_t *e;
	edge_ref_t eref;
	u4i step;
	int len;
	t = _t;
	len = g->reglen;
	step = 1;
	n = ref_nodev(g->nodes, t->node);
	if(msg) *msg = WT_TRACE_MSG_ONE;
	while(max_len == 0 || len <= max_len){
		if(n->erefs[t->dir].cnt == 0){
			if(msg) *msg = WT_TRACE_MSG_ZERO;
			break;
		} else if(n->erefs[t->dir].cnt > 1){
			if(msg) *msg = WT_TRACE_MSG_MORE;
			break;
		}
		if(max_step && step > max_step) break;
		beg_iter_edges_graph(n, t->dir, &eref);
		e = ref_iter_edges_graph(g, &eref);
		if(check_node_visit){
			if(n->bt_visit == check_node_visit){ if(msg) *msg = WT_TRACE_MSG_VISITED; break; }
			else n->bt_visit = check_node_visit;
		}
		t->edges[t->dir].idx = offset_edgev(g->edges, e);
		t->edges[t->dir].flg = eref.flg;
		if(path) t = next_ref_tracev(path);
		t->node = get_edge_didx(e, eref.flg);
		t->dir  = get_edge_ddir(e, eref.flg);
		t->edges[!t->dir] = (edge_ref_t){offset_edgev(g->edges, e), !eref.flg, 0, 0};
		t->edges[t->dir] = EDGE_REF_NULL;
		len += g->reglen + get_edge_off(e);
		step ++;
		n = ref_nodev(g->nodes, t->node);
		if(n->erefs[!t->dir].cnt > 1){
			if(msg) *msg = - 1 - WT_TRACE_MSG_MORE;
			break;
		}
	}
	return len;
}

static inline int exclusive_linear_trace_graph(Graph *g, trace_t *_t, int max_len, u4i max_step, u8i check_node_visit, tracev *path, int *msg){
	trace_t *t;
	node_t *n;
	edge_t *e;
	edge_ref_t eref;
	u4i step;
	int len;
	t = _t;
	len = g->reglen;
	step = 1;
	n = ref_nodev(g->nodes, t->node);
	if(msg) *msg = WT_TRACE_MSG_ONE;
	while(max_len == 0 || len <= max_len){
		if(n->erefs[t->dir].cnt == 0){
			if(msg) *msg = WT_TRACE_MSG_ZERO;
			break;
		} else if(n->erefs[t->dir].cnt > 1){
			if(msg) *msg = WT_TRACE_MSG_MORE;
			break;
		}
		if(max_step && step > max_step) break;
		beg_iter_edges_graph(n, t->dir, &eref);
		e = ref_iter_edges_graph(g, &eref);
		n = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
		if(n->erefs[!get_edge_ddir(e, eref.flg)].cnt > 1){
			if(msg) *msg = - 1 - WT_TRACE_MSG_MORE;
			break;
		}
		if(check_node_visit){
			if(n->bt_visit == check_node_visit){ if(msg) *msg = WT_TRACE_MSG_VISITED; break; }
			else { n->bt_visit = check_node_visit; }
		}
		t->edges[t->dir].idx = offset_edgev(g->edges, e);
		t->edges[t->dir].flg = eref.flg;
		if(path) t = next_ref_tracev(path);
		t->node = get_edge_didx(e, eref.flg);
		t->dir  = get_edge_ddir(e, eref.flg);
		t->edges[!t->dir] = (edge_ref_t){offset_edgev(g->edges, e), !eref.flg, 0, 0};
		t->edges[t->dir] = EDGE_REF_NULL;
		len += g->reglen + get_edge_off(e);
		step ++;
	}
	return len;
}

static inline u8i del_isolated_nodes_graph(Graph *g, FILE *log){
	node_t *n;
	u8i ret, i;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		if(n->erefs[0].cnt || n->erefs[1].cnt) continue;
		n->closed = 1;
		if(log) fprintf(log, "DEL_ISO\tN%llu\n", i);
		ret ++;
	}
	return ret;
}

static inline u8i cut_binary_edges_graph(Graph *g){
	UUhash *hash;
	node_t *n;
	edge_ref_t eref;
	edge_t *e, *p;
	u8i idx, nid, ret;
	ret = 0;
	hash = init_UUhash(13);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		clear_UUhash(hash);
		beg_iter_edges_graph(n, 0, &eref);
		while((e = ref_iter_edges_graph_core(g->edges, &eref, 1))){
			if(is_edge_closed(e) < WT_EDGE_CLOSED_LESS){
				put_UUhash(hash, (UUhash_t){get_edge_didx(e, eref.flg), offset_edgev(g->edges, e)});
			}
		}
		beg_iter_edges_graph(n, 1, &eref);
		while((e = ref_iter_edges_graph_core(g->edges, &eref, 1))){
			if(is_edge_closed(e) < WT_EDGE_CLOSED_LESS){
				idx = getval_UUhash(hash, get_edge_didx(e, eref.flg));
				if(idx != MAX_U8){
					p = ref_edgev(g->edges, idx);
					cut_edge_graph_core(g->nodes, g->edges, e, WT_EDGE_CLOSED_HARD);
					cut_edge_graph_core(g->nodes, g->edges, p, WT_EDGE_CLOSED_HARD);
					ret += 2;
				}
			}
		}
	}
	free_UUhash(hash);
	return ret;
}

#endif
