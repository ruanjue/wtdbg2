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

#ifndef __DBGCNS_CNS_RJ_H
#define __DBGCNS_CNS_RJ_H

#include "dna.h"
#include "list.h"
#include "hashset.h"
#include "thread.h"
#include "ksw.h"
#include "chararray.h"

/*
#ifndef DBGCNS_DEBUG
#define DBGCNS_DEBUG	0
#endif
*/

static int DBGCNS_DEBUG = 0;

#define DBGCNS_KMER_MAX_SIZE	21
#define DBGCNS_KMER_MAX_NODE_COV	0x3FF
#define DBGCNS_KMER_MAX_EDGE_COV	0xFF
#define DBGCNS_MAX_UID	0xFFF

typedef struct {
	uint64_t mer:42, cov:10, visit:12;
	uint8_t  edges[2][4];
} dbgcns_kmer_t;
define_list(dbgcnskmerv, dbgcns_kmer_t);
#define UD(E) ((dbgcnskmerv*)set->userdata)->buffer[E].mer
#define kmer_hashcode(E) u64hashcode(UD(E))
#define kmer_hashequals(E1, E2) UD(E1) == UD(E2)
define_hashset(dbgcnskmerhash, uint32_t, kmer_hashcode, kmer_hashequals);
#undef UD

typedef struct {
	u4i ksize;
	u8i kmask;
	int hz;
	dbgcnskmerv    *kmers;
	dbgcnskmerhash *khash;
	u1v *zseq;
} DBG;

#define BT_SCORE_MIN	-0x0FFFFFFF

typedef struct {
	u8i mer:42, off:21, closed:1;
	u2i n_in, n_visit;
	u4i edges, bt_node, bt_edge;
	int bt_score;
	u8i ptr; // refer to dbg->kmers
} fbg_kmer_t;
#define fbgkmer_hashcode(E) u64hashcode((E).mer)
#define fbgkmer_equals(E1, E2) ((E1).mer == (E2).mer)
define_hashset(fbgkmerh, fbg_kmer_t, fbgkmer_hashcode, fbgkmer_equals);

typedef struct {
	u4i node, cov:28, most:1, key:2, select:1;
	u4i dist;
	u4i link;
	u4i next;
} fbg_edge_t;
define_list(fbgedgev, fbg_edge_t);

typedef struct {
	u4i ridx, roff:15, rlen:15, key:1, select:1;
	u4i next;
} fbg_link_t;
define_list(fbglinkv, fbg_link_t);

typedef struct {
	u8i ridx:16, kidx:30, koff:17, closed:1;
} rd_kmer_t;
define_list(rdkmerv, rd_kmer_t);

typedef struct {
	u4i lidx, len, gidx;
} link_grp_t;
define_list(linkgrpv, link_grp_t);

typedef struct {
	fbgkmerh *kmers;
	fbgedgev *edges;
	fbglinkv *links;
	linkgrpv *grps;
	rdkmerv  *mats;
	u1v *starseq;
} FBG;

#define DBGCNS_DP_SCORE_MIN	-0x7FFFFFFF
#define DBGCNS_PATH_M	0
#define DBGCNS_PATH_X	1
#define DBGCNS_PATH_I	2
#define DBGCNS_PATH_D	3
#define DBGCNS_CNS_NON	0
#define DBGCNS_CNS_TIP	1
#define DBGCNS_CNS_CUT	2
#define DBGCNS_CNS_EXT	3
#define DBGCNS_CNS_HIT	4

typedef struct {
	union {
		struct { uint32_t kidx:30, path:2; uint32_t qpos; };
		uint64_t identifier;
	};
	int score;
	uint32_t bt_idx;
} dbgcns_dp_t;
define_list(dbgcnsdpv, dbgcns_dp_t);
#define DD(E) ((dbgcnsdpv*)set->userdata)->buffer[E]
#define dp_hashcode(E) u64hashcode(DD(E).identifier)
#define dp_hashequals(E1, E2) DD(E1).identifier == DD(E2).identifier
define_hashset(dbgcnsdphash, uint32_t, dp_hashcode, dp_hashequals);
#undef DD

typedef struct {uint64_t off:40, len:23, solid:1;} blk_t;
define_list(blkv, blk_t);

typedef struct {
	DBG      *g;
	FBG      *fbg;
	int      C, M, X, I, D, E, H, L;
	int      Z, W;

	u1v      *qseqs;
	blkv     *qblks;

	uint8_t  *qry;
	uint32_t qlen;
	u4i      qidx;
	int      avg_cov;

	dbgcnsdpv      *dps;
	u4v      *heap;
	dbgcnsdphash   *hash;

	b4v      *qmaxs;
	uint32_t qtop;
	int      max_score;
	uint32_t best_idx;

	String   *seq;
	u1v      *cns;
	u1v      *cigars;
	int      alns[4];
} CNS;

static inline DBG* init_dbg(uint32_t ksize){
	DBG *g;
	if(ksize > DBGCNS_KMER_MAX_SIZE){
		fprintf(stderr, " -- ksize MUST be no greater than %d, but is %d in %s -- %s:%d --\n", DBGCNS_KMER_MAX_SIZE, ksize, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		exit(1);
	}
	g = malloc(sizeof(DBG));
	g->hz = 0;
	g->ksize = ksize;
	g->kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - ksize) << 1);
	g->kmers = init_dbgcnskmerv(32);
	next_ref_dbgcnskmerv(g->kmers);
	memset(g->kmers->buffer, 0, sizeof(dbgcns_kmer_t));
	g->khash = init_dbgcnskmerhash(1023);
	set_userdata_dbgcnskmerhash(g->khash, g->kmers);
	g->zseq = init_u1v(64);
	return g;
}

static inline void free_dbg(DBG *g){
	free_dbgcnskmerv(g->kmers);
	free_dbgcnskmerhash(g->khash);
	free_u1v(g->zseq);
	free(g);
}

static inline void clear_dbg(DBG *g){
	clear_dbgcnskmerv(g->kmers);
	next_ref_dbgcnskmerv(g->kmers);
	memset(g->kmers->buffer, 0, sizeof(dbgcns_kmer_t));
	clear_dbgcnskmerhash(g->khash);
}

static inline int add_kmer_dbg(DBG *g, u8i kmer, u4i ridx, uint8_t fbase, uint8_t rbase, int inc){
	dbgcns_kmer_t *k;
	u8i krev;
	uint32_t *u;
	int dir, exists;
	if(0){
		krev = dna_rev_seq(kmer, g->ksize);
		if(kmer < krev){ dir = 0; }
		else { kmer = krev; dir = 1; }
	} else {
		dir = 0;
	}
	g->kmers->buffer[0].mer = kmer;
	u = prepare_dbgcnskmerhash(g->khash, 0, &exists);
	if(exists){
		k = ref_dbgcnskmerv(g->kmers, *u);
	} else {
		*u = g->kmers->size;
		k = next_ref_dbgcnskmerv(g->kmers);
		memset(k, 0, sizeof(dbgcns_kmer_t));
		k->mer = kmer;
		k->cov = 0;
		k->visit = 0;
	}
	if(k->visit != ridx + 1 && k->cov < DBGCNS_KMER_MAX_NODE_COV && (inc || k->cov == 0)){
		k->cov ++;
	}
	k->visit = ridx + 1;
	if(dir){
		if(fbase < 4){
			fbase = (~fbase) & 0x03;
			if(k->edges[1][fbase] < DBGCNS_KMER_MAX_EDGE_COV && (inc || k->edges[1][fbase] == 0)) k->edges[1][fbase] ++;
		}
		if(rbase < 4){
			if(k->edges[0][rbase] < DBGCNS_KMER_MAX_EDGE_COV && (inc || k->edges[0][rbase] == 0)) k->edges[0][rbase] ++;
		}
	} else {
		if(fbase < 4){
			if(k->edges[0][fbase] < DBGCNS_KMER_MAX_EDGE_COV && (inc || k->edges[0][fbase] == 0)) k->edges[0][fbase] ++;
		}
		if(rbase < 4){
			rbase = (~rbase) & 0x03;
			if(k->edges[1][rbase] < DBGCNS_KMER_MAX_EDGE_COV && (inc || k->edges[1][rbase] == 0)) k->edges[1][rbase] ++;
		}
	}
	return exists;
}

static inline void homopolymer_compress_dbg(DBG *g, u1i *seq, u4i len){
	u4i i;
	u1i b;
	clear_u1v(g->zseq);
	b = 4;
	for(i=0;i<len;i++){
		if(seq[i] == b) continue;
		b = seq[i];
		push_u1v(g->zseq, b);
	}
}

static inline void add_seq_dbg(DBG *g, u4i ridx, uint8_t *seq, uint32_t len, int cov_inc){
	u8i kmer;
	uint32_t i;
	uint8_t b, f, r;
	if(g->hz){
		homopolymer_compress_dbg(g, seq, len);
		seq = g->zseq->buffer;
		len = g->zseq->size;
	}
	kmer = 0;
	for(i=0;i<len;){
		b = seq[i];
		kmer = ((kmer << 2) | b) & g->kmask;
		i ++;
		if(i < g->ksize) continue;
		f = (i < len)? seq[i] : 4;
		r = (i > g->ksize)? seq[i - g->ksize] : 4;
		add_kmer_dbg(g, kmer, ridx, f, r, cov_inc);
	}
}

static inline int kmer_cov_seq_dbg(DBG *g, u1i *seq, u4i len, u4i uid){
	dbgcns_kmer_t *k;
	u8i kmer;
	uint32_t i, *u;
	uint8_t b;
	int cov;
	uid = uid % DBGCNS_MAX_UID;
	if(g->hz){
		homopolymer_compress_dbg(g, seq, len);
		seq = g->zseq->buffer;
		len = g->zseq->size;
	}
	cov = 0;
	kmer = 0;
	for(i=0;i<len;){
		b = seq[i];
		kmer = ((kmer << 2) | b) & g->kmask;
		i ++;
		if(i < g->ksize) continue;
		g->kmers->buffer[0].mer = kmer;
		u = get_dbgcnskmerhash(g->khash, 0);
		if(u){
			k = ref_dbgcnskmerv(g->kmers, *u);
			if(k->visit == uid){ cov --; continue; }
			k->visit = uid;
			cov += k->cov > 1? k->cov + 1 : - 1;
		}
	}
	return cov;
}

static inline void print_kmers_dbg(DBG *g, FILE *out){
	dbgcns_kmer_t *k;
	uint64_t i;
	char seq[DBGCNS_KMER_MAX_SIZE + 1];
	for(i=1;i<g->kmers->size;i++){
		k = ref_dbgcnskmerv(g->kmers, i);
		kmer2seq(seq, k->mer, g->ksize);
		fprintf(out, "%s\t%d\t%d,%d,%d,%d\t%d,%d,%d,%d\n", seq, k->cov,
			k->edges[0][0], k->edges[0][1], k->edges[0][2], k->edges[0][3],
			k->edges[1][0], k->edges[1][1], k->edges[1][2], k->edges[1][3]);
	}
}

static inline CNS* init_cns(uint32_t ksize, int Z, int W, int M, int X, int I, int D, int E, int H, int L){
	CNS *cns;
	cns = malloc(sizeof(CNS));
	cns->g = init_dbg(ksize);
	cns->fbg = malloc(sizeof(FBG));
	cns->fbg->kmers = init_fbgkmerh(1023);
	cns->fbg->edges = init_fbgedgev(32);
	memset(next_ref_fbgedgev(cns->fbg->edges), 0, sizeof(fbg_edge_t));
	cns->fbg->links = init_fbglinkv(32);
	memset(next_ref_fbglinkv(cns->fbg->links), 0, sizeof(fbg_link_t));
	cns->fbg->grps = init_linkgrpv(32);
	cns->fbg->mats = init_rdkmerv(32);
	cns->fbg->starseq = init_u1v(32);
	cns->qseqs = init_u1v(32);
	cns->qblks = init_blkv(32);
	cns->Z = Z;
	cns->W = W;
	cns->C = 1;
	cns->M = M;
	cns->X = X;
	cns->I = I;
	cns->D = D;
	cns->E = E;
	cns->H = H;
	cns->L = L;
	cns->qlen = 0;
	cns->qidx = 0;
	cns->avg_cov = 0;
	cns->dps  = init_dbgcnsdpv(32);
	next_ref_dbgcnsdpv(cns->dps);
	memset(cns->dps->buffer, 0, sizeof(dbgcns_dp_t)); // no need, it is always zero-filled
	cns->heap = init_u4v(32);
	cns->hash = init_dbgcnsdphash(1023);
	set_userdata_dbgcnsdphash(cns->hash, cns->dps);
	cns->qmaxs = init_b4v(32);
	cns->qtop = 0;
	cns->max_score = DBGCNS_DP_SCORE_MIN;
	cns->best_idx  = 0;
	cns->seq = init_string(32);
	cns->cns = init_u1v(32);
	cns->cigars = init_u1v(32);
	return cns;
}

static inline void free_cns(CNS *cns){
	free_dbg(cns->g);
	free_fbgkmerh(cns->fbg->kmers);
	free_fbgedgev(cns->fbg->edges);
	free_fbglinkv(cns->fbg->links);
	free_linkgrpv(cns->fbg->grps);
	free_rdkmerv(cns->fbg->mats);
	free_u1v(cns->fbg->starseq);
	free(cns->fbg);
	free_u1v(cns->qseqs);
	free_blkv(cns->qblks);
	free_dbgcnsdpv(cns->dps);
	free_u4v(cns->heap);
	free_dbgcnsdphash(cns->hash);
	free_b4v(cns->qmaxs);
	free_string(cns->seq);
	free_u1v(cns->cns);
	free_u1v(cns->cigars);
	free(cns);
}

static inline void reset_cns(CNS *cns){
	clear_dbg(cns->g);
	clear_u1v(cns->qseqs);
	clear_blkv(cns->qblks);
	clear_fbgkmerh(cns->fbg->kmers);
	clear_fbgedgev(cns->fbg->edges);
	memset(next_ref_fbgedgev(cns->fbg->edges), 0, sizeof(fbg_edge_t));
	clear_fbglinkv(cns->fbg->links);
	memset(next_ref_fbglinkv(cns->fbg->links), 0, sizeof(fbg_link_t));
	cns->qry = NULL;
	cns->qlen = 0;
}

static inline void add_seq_cns(CNS *cns, char *seq, int len, int solid){
	blk_t *b;
	int i;
	b = next_ref_blkv(cns->qblks);
	b->off = cns->qseqs->size;
	b->len = len;
	b->solid = solid;
	for(i=0;i<len;i++) push_u1v(cns->qseqs, base_bit_table[(int)seq[i]]);
}

static inline void ready_cns(CNS *cns){
	UNUSED(cns);
	//blk_t *b;
	//u4i i;
	//for(i=0;i<cns->qblks->size;i++){
		//b = ref_blkv(cns->qblks, i);
		//add_seq_dbg(cns->g, cns->qseqs->buffer + b->off, b->len);
	//}
}

static inline int dbg_cns_core(CNS *cns){
	dbgcns_dp_t *dp, *dp2;
	dbgcns_kmer_t *k;
	u8i kmer, knew;
	uint32_t i, dp_idx, kidx, *u;
	int sum, cov, cut;
	uint8_t b, q;
	int exists, score, nadd, mat_only;
	if(cns->heap->size == 0) return DBGCNS_CNS_NON;
	dp_idx = cns->heap->buffer[0]; //array_heap_pop(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	encap_dbgcnsdpv(cns->dps, 9);
	dp = ref_dbgcnsdpv(cns->dps, dp_idx);
	mat_only = 0;
	if(dp->qpos >= cns->qlen) return DBGCNS_CNS_HIT;
	if(dp->qpos + cns->W < cns->qtop){
		mat_only = 1;
		//array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		//return DBGCNS_CNS_CUT;
	} else if(dp->qpos > cns->qtop){
		cns->qtop = dp->qpos;
		for(i=cns->heap->size;i>0;i--){
			if(cns->dps->buffer[cns->heap->buffer[i-1]].qpos + cns->W < cns->qtop){
				array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, i - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			}
		}
	}
	if(dp->score < cns->qmaxs->buffer[dp->qpos]){
		mat_only = 1;
		//array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		//return DBGCNS_CNS_CUT;
	} else if(dp->score + cns->Z * cns->X > cns->qmaxs->buffer[dp->qpos]){
		cns->qmaxs->buffer[dp->qpos] = dp->score + cns->Z * cns->X;
	}
	u = prepare_dbgcnsdphash(cns->hash, dp_idx, &exists);
	if(exists){
		dp2 = ref_dbgcnsdpv(cns->dps, *u);
		if(dp->score > dp2->score){
			*u = dp_idx;
		} else {
			array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			return DBGCNS_CNS_CUT;
		}
	} else {
		*u = dp_idx;
	}
	array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	k = ref_dbgcnskmerv(cns->g->kmers, dp->kidx);
	kmer = k->mer;
	q = cns->qry[dp->qpos];
	sum = k->edges[0][0] + k->edges[0][1] + k->edges[0][2] + k->edges[0][3];
	if(sum == 0) return DBGCNS_CNS_TIP;
	cut = num_max((cns->avg_cov + cns->L - 1) / cns->L, 1);
	//cut = num_max((cns->avg_cov + 2) / 3, 1);
	nadd = 0;
	for(b=0;b<4;b++){
		if((cov = k->edges[0][b]) == 0) continue;
		if(b == q){
			score = dp->score + cns->M;
		} else if(mat_only){
			continue;
		} else {
			score = dp->score + cns->X;
		}
		//score += (cov > 1)? (cov >= cut? 1 : 0) : -1;
		score += (cov > cut)? (cns->H + (b == q? 1 : 0)) : cov - cut;
		//score += (cov > cut)? cns->H : cov - cut;
		knew = ((kmer << 2) | b) & cns->g->kmask;
		cns->g->kmers->buffer[0].mer = knew;
		u = get_dbgcnskmerhash(cns->g->khash, 0);
		kidx = *u;
		dp2 = next_ref_dbgcnsdpv(cns->dps);
		dp2->kidx = kidx;
		dp2->path = (b == q)? DBGCNS_PATH_M : DBGCNS_PATH_X;
		dp2->qpos = dp->qpos + 1;
		dp2->score = score;
		dp2->bt_idx = dp_idx;
		array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		nadd ++;
		// deletion
		if(dp->path != DBGCNS_PATH_I){
			dp2 = next_ref_dbgcnsdpv(cns->dps);
			score = dp->score + (cns->E + (dp->path == DBGCNS_PATH_D? 0 : cns->D));
			if(dp->path != DBGCNS_PATH_D) score += (cov > cut)? 0 : cov - cut;
			else score += (cov > cut)? cns->H : cov - cut;
			dp2->kidx = kidx;
			dp2->path = DBGCNS_PATH_D;
			dp2->qpos = dp->qpos;
			dp2->score = score;
			dp2->bt_idx = dp_idx;
			array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			nadd ++;
		}
	}
	// insertion
	if(mat_only == 0 && dp->path != DBGCNS_PATH_D){
		dp2 = next_ref_dbgcnsdpv(cns->dps);
		score = dp->score + (cns->E + (dp->path == DBGCNS_PATH_I? 0 : cns->I));
		//if(dp->qpos && cns->qry[dp->qpos] == cns->qry[dp->qpos - 1]) score += 1; // homopolymer merge
		score += 1 - cut;
		dp2->kidx = dp->kidx;
		dp2->path = DBGCNS_PATH_I;
		dp2->qpos = dp->qpos + 1;
		dp2->score = score;
		dp2->bt_idx = dp_idx;
		array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
		nadd ++;
	}
	return nadd? DBGCNS_CNS_EXT : DBGCNS_CNS_CUT;
}

static inline void ready_core_cns(CNS *cns, int candidate_mode, int reflen){
	blk_t *b;
	u4i i;
	int tot;
	if(candidate_mode == 3) i = 1;
	else i = 0;
	tot = 0;
	for(;i<cns->qblks->size;i++){
		b = ref_blkv(cns->qblks, i);
		add_seq_dbg(cns->g, i, cns->qseqs->buffer + b->off, b->len, 1);
		tot += b->len;
	}
	if(candidate_mode == 3){
		b = ref_blkv(cns->qblks, 0);
		add_seq_dbg(cns->g, i, cns->qseqs->buffer + b->off, b->len, 0);
		tot += b->len;
	}
	cns->avg_cov = (tot + reflen - 1) / reflen;
}

static inline int run_core_cns(CNS *cns, uint8_t *qry, uint32_t qlen){
	dbgcns_dp_t *dp;
	dbgcns_kmer_t *k;
	uint32_t i, kmer, kidx, dp_idx, *u;
	int status;
	{
		cns->qry  = qry;
		cns->qlen = qlen;
		// reset auxiliaries
		clear_dbgcnsdpv(cns->dps);
		next_ref_dbgcnsdpv(cns->dps);
		memset(cns->dps->buffer, 0, sizeof(dbgcns_dp_t));
		clear_u4v(cns->heap);
		clear_dbgcnsdphash(cns->hash);
		clear_b4v(cns->qmaxs);
		for(i=0;i<qlen;i++) push_b4v(cns->qmaxs, DBGCNS_DP_SCORE_MIN);
		cns->qtop = 0;
		cns->max_score= DBGCNS_DP_SCORE_MIN;
		cns->best_idx = 0;
		clear_u1v(cns->cns);
		clear_string(cns->seq);
		clear_u1v(cns->cigars);
		cns->alns[0] = cns->alns[1] = cns->alns[2] = cns->alns[3] = 0;
	}
	// set first kmer
	kmer = 0;
	for(cns->qtop=0;cns->qtop<cns->g->ksize;cns->qtop++){
		kmer = (kmer << 2) | qry[cns->qtop];
	}
	cns->g->kmers->buffer[0].mer = kmer;
	u = get_dbgcnskmerhash(cns->g->khash, 0);
	if(u == NULL) return 0;
	kidx = *u;
	dp = next_ref_dbgcnsdpv(cns->dps);
	dp->kidx = kidx;
	dp->path = DBGCNS_PATH_M;
	dp->qpos = cns->qtop;
	dp->score = 0;
	dp->bt_idx = 0;
	array_heap_push(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, cns->dps->size - 1, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
	// dbg traversing
	while(cns->heap->size){
		status = dbg_cns_core(cns);
		if(status == DBGCNS_CNS_HIT){
			dp_idx = cns->heap->buffer[0];
			array_heap_remove(cns->heap->buffer, cns->heap->size, cns->heap->cap, uint32_t, 0, num_cmp(cns->dps->buffer[b].score, cns->dps->buffer[a].score));
			dp = ref_dbgcnsdpv(cns->dps, dp_idx);
			if(dp->score > cns->max_score){
				cns->max_score = dp->score;
				cns->best_idx = dp_idx;
			}
		}
	}
	if(cns->best_idx == 0) return 0;
	// traceback to get cns seq
	dp_idx = cns->best_idx;
	while(dp_idx){
		dp = ref_dbgcnsdpv(cns->dps, dp_idx);
		push_u1v(cns->cigars, dp->path);
		cns->alns[dp->path] ++;
		if(dp->path != DBGCNS_PATH_I){
			k = ref_dbgcnskmerv(cns->g->kmers, dp->kidx);
			push_u1v(cns->cns, k->mer & 0x03);
			if(k->cov == 1){
				add_char_string(cns->seq, "acgt"[k->mer & 0x03]);
			} else {
				add_char_string(cns->seq, "ACGT"[k->mer & 0x03]);
			}
		}
		dp_idx = dp->bt_idx;
	}
	// first ksize - 1 bases may be not corrected, truncated
	reverse_string(cns->seq);
	reverse_u1v(cns->cns);
	reverse_u1v(cns->cigars);
	return cns->seq->size;
}

static inline int hierarchical_clustering_edge_links(FBG *fbg, fbg_edge_t *e, linkgrpv *grps, double max_var){
	fbg_link_t *lnk;
	u4i lidx, gidx;
	u4i i, b;
	double sum, avg;
	clear_linkgrpv(grps);
	lidx = e->link;
	while(lidx){
		lnk = ref_fbglinkv(fbg->links, lidx);
		push_linkgrpv(grps, (link_grp_t){lidx, lnk->rlen, 0});
		lidx = lnk->next;
	}
	if(grps->size == 0) return 0;
	sort_array(grps->buffer, grps->size, link_grp_t, num_cmpgt(a.len, b.len));
	gidx = 0;
	b = 0;
	sum = avg = grps->buffer[0].len;
	for(i=1;i<grps->size;i++){
		if(grps->buffer[i].len - avg > avg * max_var){
			gidx ++;
			sum = avg = grps->buffer[i].len;
			b = i;
		} else {
			sum += grps->buffer[i].len;
			avg = sum / (i - b + 1);
			if(avg - grps->buffer[b].len > avg * max_var || grps->buffer[i].len - avg > avg * max_var){
				gidx ++;
				sum = avg = grps->buffer[i].len;
				b = i;
			}
		}
		grps->buffer[i].gidx = gidx;
	}
	return gidx + 1;
}

static inline int revise_edge_fbg(FBG *fbg, u4i kidx, u4i eidx, linkgrpv *grps, double max_var){
	fbg_kmer_t *k;
	fbg_edge_t *e, *p;
	fbg_link_t *lnk;
	u4i eidx2, eidx3;
	u4i i, b, j, ng, avg, gidx;
	u4i max_cov, max_eidx, key;
	double sum;
	e = ref_fbgedgev(fbg->edges, eidx);
	ng = hierarchical_clustering_edge_links(fbg, e, grps, max_var);
	encap_fbgedgev(fbg->edges, ng - 1);
	k = ref_fbgkmerh(fbg->kmers, kidx);
	e = ref_fbgedgev(fbg->edges, eidx);
	e->key = key = e->key? 2 : 0;
	e->most = 0;
	if(DBGCNS_DEBUG){
		if(ng > 1){
			fprintf(stderr, "REVISE K%d -> K%d cov = %d into %d edges --\n", k->off, ref_fbgkmerh(fbg->kmers, e->node)->off, e->cov, ng); fflush(stderr);
		}
	}
	ref_fbgkmerh(fbg->kmers, e->node)->n_in += ng - 1;
	p = e;
	p->link = 0;
	eidx3 = e->next;
	gidx = 0;
	sum = 0;
	max_cov = 0;
	max_eidx = 0;
	for(i=b=0;i<=grps->size;i++){
		if(i == grps->size || grps->buffer[i].gidx != gidx){
			avg = sum / (i - b) + 0.5;
			if(p == NULL){
				eidx2 = fbg->edges->size;
				p = next_ref_fbgedgev(fbg->edges);
				p->node = e->node;
				e->next = eidx2;
				p->next = eidx3;
			}
			p->link = 0;
			p->cov = i - b;
			p->dist = avg;
			p->key = key;
			p->select = 0;
			p->most = 0;
			if(p->cov > max_cov){
				max_cov = p->cov;
				max_eidx = offset_fbgedgev(fbg->edges, p);
			} else if(p->cov == max_cov){
				max_eidx = 0;
			}
			if(DBGCNS_DEBUG){
				if(ng > 1){
					fprintf(stderr, "+ %d %d --\n", p->cov, p->dist); fflush(stderr);
				}
			}
			for(j=b;j<i;j++){
				lnk = ref_fbglinkv(fbg->links, grps->buffer[j].lidx);
				lnk->next = p->link;
				p->link = grps->buffer[j].lidx;
				if(lnk->key) p->key = 1;
			}
			e = p;
			p = NULL;
			b = i;
			gidx ++;
			sum = 0;
		}
		if(i < grps->size) sum += grps->buffer[i].len;
	}
	if(max_eidx){
		ref_fbgedgev(fbg->edges, max_eidx)->most = 1;
	}
	return gidx;
}

static inline void revise_edge_cov_fbg(FBG *fbg, fbg_edge_t *e){
	fbg_link_t *lnk;
	u4i lidx, cnt;
	double sum, var, max, avg, std;
	sum = 0;
	var = 0;
	max = 0;
	cnt = 0;
	lidx = e->link;
	while(lidx){
		lnk = ref_fbglinkv(fbg->links, lidx);
		if((double)lnk->rlen > max) max = lnk->rlen;
		sum += lnk->rlen;
		var += lnk->rlen * lnk->rlen;
		cnt ++;
		lidx = lnk->next;
	}
	if(cnt == 0) return;
	avg = sum / cnt;
	std = sqrt(var / cnt - avg * avg);
	if(std > 10 && std > 0.2 * avg) std = 0.2 * max;
	if(std < 1) std = 1;
	cnt = 0;
	sum = 0;
	lidx = e->link;
	while(lidx){
		lnk = ref_fbglinkv(fbg->links, lidx);
		if(num_diff((double)lnk->rlen, avg) <= std){
			sum += lnk->rlen;
			cnt ++;
		}
		lidx = lnk->next;
	}
	if(cnt == 0) cnt = 1;
	avg = sum / cnt;
	e->cov = cnt;
	e->dist = avg + 0.5;
}

static inline void build_DirectFuzzyBruijnGraph(CNS *cns, u4i ridx, u4i cov_cutoff, double max_dist_var){
	DBG *g;
	dbgcns_kmer_t *k;
	fbg_kmer_t A, *a;
	fbg_edge_t *e;
	fbg_link_t *l;
	rd_kmer_t *rk;
	u8i kmer;
	u4i r, rr, i, j, beg, end, c, len, idx1, off1, idx2, off2, *u;
	u4i eidx;
	u1i b, *seq;
	int exists;
	g = cns->g;
	// select high cov kmers
	if(cov_cutoff){
		seq = cns->qseqs->buffer + cns->qblks->buffer[ridx].off;
		len = cns->qblks->buffer[ridx].len;
		memset(&A, 0, sizeof(fbg_kmer_t));
		kmer = 0;
		j = 0x0000FFFFU;
		for(i=0;i<len;){
			b = seq[i];
			kmer = ((kmer << 2) | b) & g->kmask;
			i ++;
			if(i < g->ksize) continue;
			g->kmers->buffer[0].mer = kmer;
			u = get_dbgcnskmerhash(g->khash, 0);
			if(u == NULL) continue;
			k = ref_dbgcnskmerv(g->kmers, *u);
			if(k->cov < cov_cutoff) continue;
			//if(j + 1 == i){ j = i; continue; }
			A.mer = kmer;
			a = prepare_fbgkmerh(cns->fbg->kmers, A, &exists);
			if(exists){
				a->closed = 1;
			} else {
				a->mer = kmer;
				a->off = i - g->ksize;
				a->closed = 0;
				a->n_in = 0;
				a->n_visit = 0;
				a->edges = 0;
				a->ptr = *u;
				j = i;
			}
		}
	} else {
		// select best cov per 10 bp, but cov >= 4
		clear_rdkmerv(cns->fbg->mats);
		seq = cns->qseqs->buffer + cns->qblks->buffer[ridx].off;
		len = cns->qblks->buffer[ridx].len;
		memset(&A, 0, sizeof(fbg_kmer_t));
		kmer = 0;
		for(i=0;i<len;){
			b = seq[i];
			kmer = ((kmer << 2) | b) & g->kmask;
			i ++;
			if(i < g->ksize) continue;
			g->kmers->buffer[0].mer = kmer;
			u = get_dbgcnskmerhash(g->khash, 0);
			if(u == NULL) continue;
			if(g->kmers->buffer[*u].cov < 4) continue;
			push_rdkmerv(cns->fbg->mats, (rd_kmer_t){ridx, *u, i - g->ksize, 0});
		}
		sort_array(cns->fbg->mats->buffer, cns->fbg->mats->size, rd_kmer_t, num_cmpgtx(g->kmers->buffer[b.kidx].cov, g->kmers->buffer[a.kidx].cov, a.kidx, b.kidx));
		for(i=1;i<cns->fbg->mats->size;i++){
			rk = ref_rdkmerv(cns->fbg->mats, i);
			if(rk->kidx == ref_rdkmerv(cns->fbg->mats, i - 1)->kidx){
				rk->closed = 1;
				ref_rdkmerv(cns->fbg->mats, i - 1)->closed = 1;
			}
		}
		sort_array(cns->fbg->mats->buffer, cns->fbg->mats->size, rd_kmer_t, num_cmpgt(a.koff, b.koff));
		beg = 0;
		j = 0xFFFFFFFFU;
		for(i=0;i<cns->fbg->mats->size;i++){
			rk = ref_rdkmerv(cns->fbg->mats, i);
			if(rk->closed) continue;
			if(j == 0xFFFFFFFFU) j = i;
			if(rk->koff - beg < 10){
				if(g->kmers->buffer[rk->kidx].cov > g->kmers->buffer[cns->fbg->mats->buffer[j].kidx].cov){
					j = i;
				}
				continue;
			}
			rk = ref_rdkmerv(cns->fbg->mats, j);
			beg = rk->koff;
			j = i;
			i --;
			k = ref_dbgcnskmerv(g->kmers, rk->kidx);
			kmer = k->mer;
			A.mer = kmer;
			a = prepare_fbgkmerh(cns->fbg->kmers, A, &exists);
			if(exists){
				a->closed = 1;
			} else {
				a->mer = kmer;
				a->off = rk->koff;
				a->closed = 0;
				a->n_in = 0;
				a->n_visit = 0;
				a->edges = 0;
				a->ptr = rk->kidx;
			}
		}
	}
	clear_rdkmerv(cns->fbg->mats);
	for(r=0;r<cns->qblks->size;r++){
		seq = cns->qseqs->buffer + cns->qblks->buffer[r].off;
		len = cns->qblks->buffer[r].len;
		memset(&A, 0, sizeof(fbg_kmer_t));
		kmer = 0;
		for(i=0;i<len;){
			b = seq[i];
			kmer = ((kmer << 2) | b) & g->kmask;
			i ++;
			if(i < g->ksize) continue;
			A.mer = kmer;
			a = get_fbgkmerh(cns->fbg->kmers, A);
			if(a == NULL) continue;
			if(a->closed) continue;
			push_rdkmerv(cns->fbg->mats, (rd_kmer_t){r, offset_fbgkmerh(cns->fbg->kmers, a), i - g->ksize, 0});
		}
	}
	sort_array(cns->fbg->mats->buffer, cns->fbg->mats->size, rd_kmer_t, num_cmpgtx(a.ridx, b.ridx, a.kidx, b.kidx));
	for(i=1;i<cns->fbg->mats->size;i++){
		if(cns->fbg->mats->buffer[i-1].ridx == cns->fbg->mats->buffer[i].ridx && cns->fbg->mats->buffer[i-1].kidx == cns->fbg->mats->buffer[i].kidx){
			ref_fbgkmerh(cns->fbg->kmers, cns->fbg->mats->buffer[i].kidx)->closed = 1;
		}
	}
	for(i=0;i<cns->fbg->mats->size;i++){
		if(ref_fbgkmerh(cns->fbg->kmers, cns->fbg->mats->buffer[i].kidx)->closed){
			cns->fbg->mats->buffer[i].closed = 1;
		}
	}
	sort_array(cns->fbg->mats->buffer, cns->fbg->mats->size, rd_kmer_t, num_cmpgtx(a.ridx, b.ridx, a.koff, b.koff));
	// add edges
	for(i=0;i+1<cns->fbg->mats->size;i++){
		if(cns->fbg->mats->buffer[i].closed) continue;
		idx1 = cns->fbg->mats->buffer[i].kidx;
		off1 = cns->fbg->mats->buffer[i].koff;
		a = ref_fbgkmerh(cns->fbg->kmers, idx1);
		for(j=i+1,c=0;j<cns->fbg->mats->size&&c<1;j++){
			if(cns->fbg->mats->buffer[j].ridx != cns->fbg->mats->buffer[i].ridx) break;
			if(cns->fbg->mats->buffer[j].closed) continue;
			idx2 = cns->fbg->mats->buffer[j].kidx;
			off2 = cns->fbg->mats->buffer[j].koff;
			if(a->off >= ref_fbgkmerh(cns->fbg->kmers, idx2)->off) continue;
			c ++;
			eidx = a->edges;
			while(eidx){
				e = ref_fbgedgev(cns->fbg->edges, eidx);
				if(e->node == idx2) break;
				eidx = e->next;
			}
			if(eidx == 0){
				eidx = cns->fbg->edges->size;
				e = next_ref_fbgedgev(cns->fbg->edges);
				e->node = idx2;
				e->link = 0;
				e->next = a->edges;
				e->cov  = 0;
				e->dist = 0;
				e->key = 0;
				e->select = 0;
				a->edges = eidx;
				ref_fbgkmerh(cns->fbg->kmers, idx2)->n_in ++;
			}
			e = ref_fbgedgev(cns->fbg->edges, eidx);
			if(cns->fbg->mats->buffer[i].ridx == ridx) e->key = 1;
			e->cov ++;
			l = next_ref_fbglinkv(cns->fbg->links);
			l->ridx = cns->fbg->mats->buffer[i].ridx;
			l->roff = off1;
			l->rlen = off2 - off1;
			l->key = (cns->fbg->mats->buffer[i].ridx == ridx);
			l->select = 0;
			l->next = e->link;
			e->link = cns->fbg->links->size - 1;
		}
	}
	rr = 0xFFFFFFFFU;
	for(beg=end=0;beg<cns->fbg->mats->size;beg=end){
		r = rr;
		for(;end<cns->fbg->mats->size;end++){
			if(cns->fbg->mats->buffer[end].closed) continue;
			if(r == 0xFFFFFFFFU){
				r = cns->fbg->mats->buffer[end].ridx;
			} else if(cns->fbg->mats->buffer[end].ridx == r) continue;
			else { rr = cns->fbg->mats->buffer[end].ridx; break; }
		}
		if(r == ridx) continue; // reference seq
		for(i=beg;i+2<end;i++){
			if(cns->fbg->mats->buffer[i].closed) continue;
			idx1 = cns->fbg->mats->buffer[i].kidx;
			off1 = cns->fbg->mats->buffer[i].koff;
			a = ref_fbgkmerh(cns->fbg->kmers, idx1);
			c = 0;
			for(j=i+1;j<end;j++){
				if(cns->fbg->mats->buffer[j].closed) continue;
				idx2 = cns->fbg->mats->buffer[j].kidx;
				off2 = cns->fbg->mats->buffer[j].koff;
				if(a->off >= ref_fbgkmerh(cns->fbg->kmers, idx2)->off) continue;
				c ++;
				if(c < 2) continue;
				if(c > 5) break;
				eidx = a->edges;
				while(eidx){
					e = ref_fbgedgev(cns->fbg->edges, eidx);
					if(e->node == idx2) break;
					eidx = e->next;
				}
				if(eidx == 0) continue;
				e = ref_fbgedgev(cns->fbg->edges, eidx);
				e->cov ++;
				l = next_ref_fbglinkv(cns->fbg->links);
				l->ridx = r;
				l->roff = off1;
				l->rlen = off2 - off1;
				l->key  = 0;
				l->select = 0;
				l->next = e->link;
				e->link = cns->fbg->links->size - 1;
			}
		}
	}
	if(0){
		for(i=1;i<cns->fbg->edges->size;i++){
			revise_edge_cov_fbg(cns->fbg, ref_fbgedgev(cns->fbg->edges, i));
		}
	} else {
		reset_iter_fbgkmerh(cns->fbg->kmers);
		while((a = ref_iter_fbgkmerh(cns->fbg->kmers))){
			eidx = a->edges;
			while(eidx){
				e = ref_fbgedgev(cns->fbg->edges, eidx);
				eidx = e->next;
				revise_edge_fbg(cns->fbg, offset_fbgkmerh(cns->fbg->kmers, a), offset_fbgedgev(cns->fbg->edges, e), cns->fbg->grps, max_dist_var);
			}
		}
	}
}

static inline void print_dot_DirectFuzzyBruijnGraph(CNS *cns, FILE *out){
	fbg_kmer_t *k, *n;
	fbg_edge_t *e;
	fbg_link_t *l;
	u4i eidx, lidx;
	if(out == NULL) return;
	fprintf(out, "digraph {\n\trankdir=LR\n");
	reset_iter_fbgkmerh(cns->fbg->kmers);
	while((k = ref_iter_fbgkmerh(cns->fbg->kmers))){
		if(k->closed) continue;
		eidx = k->edges;
		while(eidx){
			e = ref_fbgedgev(cns->fbg->edges, eidx);
			if(e->cov >= 3 || e->key || e->select){
				n = ref_fbgkmerh(cns->fbg->kmers, e->node);
				lidx = e->link;
				while(lidx){
					l = ref_fbglinkv(cns->fbg->links, lidx);
					fprintf(out, "\tK%d -> K%d [label=\"R%04d_%d_%d(%d:%d)\" color=%s%s]\n", k->off, n->off, l->ridx, l->roff, l->rlen, e->cov, e->dist, l->key? "blue" : (e->most? "green" : "black"), l->select? " style=dashed" : "");
					lidx = l->next;
				}
			}
			eidx = e->next;
		}
	}
	fprintf(out, "}\n");
}

static inline void DP_best_path_DirectFuzzyBruijnGraph(CNS *cns){
	FBG *fbg;
	fbg_kmer_t *k, *n;
	fbg_edge_t *e;
	u4v *heap;
	u4i kidx, nb, ne, nboff, neoff;
	u4i eidx, etmp;
	int ref_score, ref_one, alt_one, most_score, max_dist_var, var, score;
	float dist_var;
	ref_score = 1;
	ref_one = -50;
	alt_one = -1000;
	most_score = 2;
	max_dist_var = 100;
	dist_var = -0.5;
	fbg = cns->fbg;
	reset_iter_fbgkmerh(fbg->kmers);
	nb = ne = 0xFFFFFFFFU;
	nboff = 0xFFFFFFFFU;
	neoff = 0;
	while((k = ref_iter_fbgkmerh(fbg->kmers))){
		k->n_visit = 0;
		k->bt_node = 0xFFFFFFFFU;
		k->bt_edge = 0;
		k->bt_score = BT_SCORE_MIN;
		if(k->closed) continue;
		if(k->off < nboff){ nb = offset_fbgkmerh(fbg->kmers, k); nboff = k->off; }
		if(k->off >= neoff){ ne = offset_fbgkmerh(fbg->kmers, k); neoff = k->off; }
	}
	if(nb == 0xFFFFFFFFU) return;
	heap = init_u4v(32);
	k = ref_fbgkmerh(fbg->kmers, nb);
	k->bt_score = 0;
	push_u4v(heap, nb);
	while(heap->size){
		kidx = heap->buffer[--heap->size];
		k = ref_fbgkmerh(fbg->kmers, kidx);
		eidx = k->edges;
		while(eidx){
			e = ref_fbgedgev(fbg->edges, eidx);
			n = ref_fbgkmerh(fbg->kmers, e->node);
			score = k->bt_score + e->most * most_score;
			var = n->off - k->off;
			var = num_diff(var, (int)e->dist);
			if(e->key == 1){
				score += ref_score + (e->cov <= 1? ref_one : 0);
			} else if(e->key == 2){
				score += (e->cov <= 1? alt_one : 0);
			} else {
				if(var > max_dist_var) score = -1000000;
				else score += var * dist_var + (e->cov <= 1? alt_one : 0);
			}
			if(score > n->bt_score){
				n->bt_score = score;
				n->bt_node = kidx;
				n->bt_edge = eidx;
			}
			n->n_visit ++;
			if(n->n_visit >= n->n_in){
				push_u4v(heap, e->node);
			}
			eidx = e->next;
		}
	}
	free_u4v(heap);
	k = ref_fbgkmerh(fbg->kmers, ne);
	if(ne != nb && k->bt_edge == 0){
		FILE *out = open_file_for_write("debug.dot", NULL, 1);
		print_dot_DirectFuzzyBruijnGraph(cns, out);
		fclose(out);
		fprintf(stderr, " -- something wrong, nb = %d(K%04d), ne = %d(K%04d) in %s -- %s:%d --\n", nb, fbg->kmers->array[nb].off, ne, fbg->kmers->array[ne].off, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		exit(1);
	}
	kidx = 0xFFFFFFFFU;
	eidx = 0;
	while(kidx != nb){
		k = ref_fbgkmerh(fbg->kmers, ne);
		ne = k->bt_node;
		swap_tmp(eidx, k->bt_edge, etmp);
		k->bt_node = kidx;
		kidx = offset_fbgkmerh(fbg->kmers, k);
	}
}

static inline int correct_struct_DirectFuzzyBruijnGraph(CNS *cns, u4i ridx){
	FBG *fbg;
	fbg_kmer_t *k, *n;
	fbg_edge_t *e;
	fbg_link_t *lnk;
	int chg;
	u4i off, kidx, koff, eidx, lidx, key, upd;
	fbg = cns->fbg;
	reset_iter_fbgkmerh(fbg->kmers);
	kidx = 0xFFFFFFFFU;
	koff = 0xFFFFFFFFU;
	while((k = ref_iter_fbgkmerh(fbg->kmers))){
		if(k->closed) continue;
		if(k->off < koff){ kidx = offset_fbgkmerh(fbg->kmers, k); koff = k->off; }
	}
	clear_u1v(fbg->starseq);
	if(DBGCNS_DEBUG){
		fprintf(stderr, " -- select seq[%d] len=%d in %s -- %s:%d --\n", ridx, cns->qblks->buffer[ridx].len, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	}
	if(kidx == 0xFFFFFFFFU){
		append_array_u1v(fbg->starseq, cns->qseqs->buffer + cns->qblks->buffer[ridx].off, cns->qblks->buffer[ridx].len);
		return 0;
	}
	chg = 0;
	off = 0;
	while(1){
		k = ref_fbgkmerh(fbg->kmers, kidx);
		if(off < k->off){
			if(DBGCNS_DEBUG){
				fprintf(stderr, "-- %d + %d bases from seq[%d], offset %d -> %d\n", (int)fbg->starseq->size, k->off - off, ridx, off, k->off);
			}
			append_array_u1v(fbg->starseq, cns->qseqs->buffer + cns->qblks->buffer[ridx].off + off, k->off - off);
			off = k->off;
		}
		if(k->bt_edge == 0) break;
		eidx = k->bt_edge;
		e = ref_fbgedgev(fbg->edges, eidx);
		e->select = 1;
		key = 0;
		double sum, var, cnt, avg, std, min;
		sum = 0;
		var = 0;
		cnt = 0;
		lidx = e->link;
		while(lidx){
			lnk = ref_fbglinkv(fbg->links, lidx);
			sum += lnk->rlen;
			var += lnk->rlen * lnk->rlen;
			cnt ++;
			if(lnk->key) key = lidx;
			lidx = lnk->next;
		}
		if(key){
			upd = key;
		} else {
			avg = sum / cnt;
			std = sqrt(var / cnt - avg * avg);
			if(std < 1) std = 1;
			lidx = e->link;
			min = 100000;
			upd = key;
			while(lidx){
				lnk = ref_fbglinkv(fbg->links, lidx);
				if(num_diff((double)lnk->rlen, avg) < min){
					upd = lidx;
					min = num_diff((double)lnk->rlen, avg);
				}
				lidx = lnk->next;
			}
		}
		lnk = ref_fbglinkv(fbg->links, upd);
		lnk->select = 1;
		kidx = k->bt_node;
		n = ref_fbgkmerh(fbg->kmers, kidx);
		if(DBGCNS_DEBUG){
			fprintf(stderr, "-- %d + %d bases from seq[%d], offset %d + %d -> %d\n", (int)fbg->starseq->size, lnk->rlen, lnk->ridx, off, n->off - off, n->off);
		}
		append_array_u1v(fbg->starseq, cns->qseqs->buffer + cns->qblks->buffer[lnk->ridx].off + lnk->roff, lnk->rlen);
		off = n->off;
	}
	k = NULL;
	if(off < cns->qblks->buffer[ridx].len){
		if(DBGCNS_DEBUG){
			fprintf(stderr, "-- %d + %d bases from seq[%d], offset %d -> %d\n", (int)fbg->starseq->size, cns->qblks->buffer[ridx].len - off, ridx, off, cns->qblks->buffer[ridx].len);
		}
		append_array_u1v(fbg->starseq, cns->qseqs->buffer + cns->qblks->buffer[ridx].off + off, cns->qblks->buffer[ridx].len - off);
	}
	return chg;
}

static inline int homopolymer_analysis_cns(CNS *cns){
	u8i kmer, xmask[2];
	u4i chg, i, j, l, r, c, brun, *u, kcnts[3];
	u1i b;
	char kstr[64];
	UNUSED(kstr); // only for compile warning
	xmask[1] = 0xFFFFFFFFFFFFFFFFLLU >> (64 - (cns->g->ksize / 2 * 2));
	xmask[0] = (~xmask[1]) & cns->g->kmask;
	b = 4; brun = 0;
	clear_u1v(cns->g->zseq);
	chg = 0;
	for(i=0;i<cns->g->ksize;i++) push_u1v(cns->g->zseq, cns->cns->buffer[i]);
	for(;i+cns->g->ksize<cns->cns->size;i++){
		if(cns->cns->buffer[i] == b){
			brun ++;
		} else {
			if(brun >= 3 && brun + 2 <= cns->g->ksize){
				if(DBGCNS_DEBUG){
					fprintf(stderr, "POLY %c(%d) at pos %d\n", "ACGT"[b], brun, i);
				}
				r = i - brun;
				c = (cns->g->ksize - brun) / 2;
				l = r - c;
				kmer = 0;
				for(j=l;j<r;j++) kmer = (kmer << 2) | cns->cns->buffer[j];
				for(j=1;j<brun;j++) kmer = (kmer << 2) | b;
				l = i;
				r = l + (cns->g->ksize - brun - c) + 1;
				for(j=l;j<r;j++) kmer = (kmer << 2) | cns->cns->buffer[j];
				cns->g->kmers->buffer[0].mer = kmer;
				u = get_dbgcnskmerhash(cns->g->khash, 0);
				kcnts[0] = u? ref_dbgcnskmerv(cns->g->kmers, *u)->cov : 0;
				if(DBGCNS_DEBUG){
					kmer2seq(kstr, kmer, cns->g->ksize);
					fprintf(stderr, "%s\t%d\n", kstr, kcnts[0]);
				}
				kmer = (kmer & xmask[0]) | (((u8i)b) << (cns->g->ksize / 2 * 2 - 2)) | ((kmer & xmask[1]) >> 2);
				cns->g->kmers->buffer[0].mer = kmer;
				u = get_dbgcnskmerhash(cns->g->khash, 0);
				kcnts[1] = u? ref_dbgcnskmerv(cns->g->kmers, *u)->cov : 0;
				if(DBGCNS_DEBUG){
					kmer2seq(kstr, kmer, cns->g->ksize);
					fprintf(stderr, "%s\t%d\n", kstr, kcnts[1]);
				}
				kmer = (kmer & xmask[0]) | (((u8i)b) << (cns->g->ksize / 2 * 2 - 2)) | ((kmer & xmask[1]) >> 2);
				cns->g->kmers->buffer[0].mer = kmer;
				u = get_dbgcnskmerhash(cns->g->khash, 0);
				kcnts[2] = u? ref_dbgcnskmerv(cns->g->kmers, *u)->cov : 0;
				if(DBGCNS_DEBUG){
					kmer2seq(kstr, kmer, cns->g->ksize);
					fprintf(stderr, "%s\t%d\n", kstr, kcnts[2]);
				}
				if(kcnts[0] > kcnts[1] && kcnts[1] >= kcnts[2]){ // there is a insertion base
					brun --;
					chg ++;
					if(DBGCNS_DEBUG){
						fprintf(stderr, "#HOMO ins\n");
					}
				} else if(kcnts[2] > kcnts[1] && kcnts[1] >= kcnts[0]){ // deletion
					brun ++;
					chg ++;
					if(DBGCNS_DEBUG){
						fprintf(stderr, "#HOMO del\n");
					}
				}
			}
			for(j=0;j<brun;j++) push_u1v(cns->g->zseq, b);
			b = cns->cns->buffer[i];
			brun = 1;
		}
	}
	for(j=0;j<brun;j++) push_u1v(cns->g->zseq, b);
	for(;i<cns->cns->size;i++) push_u1v(cns->g->zseq, cns->cns->buffer[i]);
	if(chg){
		clear_u1v(cns->cns);
		append_u1v(cns->cns, cns->g->zseq);
		clear_string(cns->seq);
		encap_string(cns->seq, cns->g->zseq->size);
		for(i=0;i<cns->g->zseq->size;i++){
			cns->seq->string[i] = bit_base_table[cns->g->zseq->buffer[i]];
		}
		cns->seq->size = i;
		cns->seq->string[i] = '\0';
	}
	return chg;
}

static inline int run_cns(CNS *cns, int candidate_mode, int corr_struct){
	u4i i;
	if(cns->qblks->size == 0) return 0;
	cns->qidx = 0;
	if(candidate_mode == 4){ // longest
		for(i=0;i<cns->qblks->size;i++){
			if(cns->qblks->buffer[i].solid == 0) continue;
			if(cns->qblks->buffer[i].len > cns->qblks->buffer[cns->qidx].len) cns->qidx = i;
		}
	} else if(candidate_mode == 5){ // shortest
		for(i=0;i<cns->qblks->size;i++){
			if(cns->qblks->buffer[i].solid == 0) continue;
			if(cns->qblks->buffer[i].len < cns->qblks->buffer[cns->qidx].len) cns->qidx = i;
		}
	} else if(candidate_mode == 3){ // first and but not increase coverage
		cns->qidx = 0;
	} else if(candidate_mode == 0){ // kmer coverage
		cns->qidx = 0;
	} else if(candidate_mode == 2){ // first and include
		cns->qidx = 0;
	} else if(candidate_mode < 0){
		cns->qidx = num_min(-1 - candidate_mode, (int)(cns->qblks->size - 1));
	} else if(candidate_mode == 1){ // median
		u2v *idxs;
		idxs = init_u2v(cns->qblks->size);
		for(i=0;i<cns->qblks->size;i++){
			if(cns->qblks->buffer[i].solid == 0) continue;
			push_u2v(idxs, i);
		}
		cns->qidx = quick_median_array(idxs->buffer, idxs->size, u2i, num_cmpgt(cns->qblks->buffer[a].len, cns->qblks->buffer[b].len));
	} else {
		fprintf(stderr, " -- Unknown candidate mode %d in %s -- %s:%d --\n", candidate_mode, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	}
	ready_core_cns(cns, candidate_mode, cns->qblks->buffer[cns->qidx].len);
	if(candidate_mode == 0){
		double max, cov;
		max = - 1000000;
		for(i=0;i<cns->g->kmers->size;i++) cns->g->kmers->buffer[i].visit = 0;
		for(i=cns->qblks->size;i>0;i--){
			if(cns->qblks->buffer[i - 1].solid == 0) continue;
			cov = kmer_cov_seq_dbg(cns->g, cns->qseqs->buffer + cns->qblks->buffer[i - 1].off, cns->qblks->buffer[i - 1].len, (i - 1 + 1) % 255) * 1.0 / cns->qblks->buffer[i - 1].len;
			if(cov > max){
				max = cov;
				cns->qidx = i - 1;
			}
		}
		// revise cns->avg_cov
		cns->avg_cov = cns->avg_cov * cns->qblks->buffer[0].len / cns->qblks->buffer[cns->qidx].len;
	}
	if(corr_struct){
		build_DirectFuzzyBruijnGraph(cns, cns->qidx, 5, 0.2);
		//build_DirectFuzzyBruijnGraph(cns, cns->qidx, 0, 0.2);
		DP_best_path_DirectFuzzyBruijnGraph(cns);
		correct_struct_DirectFuzzyBruijnGraph(cns, cns->qidx);
		if(DBGCNS_DEBUG){
			FILE *dotf = open_file_for_write("debug.dot", NULL, 1);
			print_dot_DirectFuzzyBruijnGraph(cns, dotf);
			fclose(dotf);
		}
		run_core_cns(cns, cns->fbg->starseq->buffer, cns->fbg->starseq->size);
	} else {
		run_core_cns(cns, cns->qseqs->buffer + cns->qblks->buffer[cns->qidx].off, cns->qblks->buffer[cns->qidx].len);
	}
	homopolymer_analysis_cns(cns);
	return cns->seq->size;
}

#endif
