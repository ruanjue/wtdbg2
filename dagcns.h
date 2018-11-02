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

#include "dna.h"
#include "chararray.h"
#include "list.h"
#include "queue.h"
#include "hashset.h"
#include "general_graph.h"
#include <math.h>
#include <float.h>

#ifndef PACBIO_PROBS_DAGCNS_RJ_H
#define PACBIO_PROBS_DAGCNS_RJ_H

#define DAGCNS_MAX_LEN	0x3FFF // 16k

static int dagcns_debug = 0;

typedef struct {
	int x, y;
	u4i dnidx;
	f8i probs[2];
} bdp_node_t;
define_list(bdpnodev, bdp_node_t);

typedef struct {
	f8i prob;
	u1i cigar, base;
} bdp_edge_t;
define_list(bdpedgev, bdp_edge_t);

define_simple_geg_callback(bdp, bdpnodev, bdp_node_t, bdpedgev, bdp_edge_t);

typedef struct {
	u4i gnidx, dnidx;
} bdp_link_t;
define_list(bdplinkv, bdp_link_t);

typedef struct {
	uint32_t nodes[2];
	uint32_t links[2];
	uint32_t cov:28, visit:1, closed:1, cns:1, alt:1;
	double   score;
} dagedge_t;
define_list(dagedgev, dagedge_t);

#define NODE_MAX_FW_EDGE	0xFFFFFFFFU
typedef struct dagnode_t {
	uint32_t pos:28, base:2, cns:1, visit:1;
	uint32_t fw_edge;
	uint32_t edges[2];
	f8i      aux;
} dagnode_t;
define_list(dagnodev, dagnode_t);

typedef struct {
	uint32_t pos;
	uint32_t bases[4];
} dagsnp_t;
define_list(dagsnpv, dagsnp_t);

typedef struct {
	u8list   *cns;
	u32list  *deps;
	dagnodev *nodes;
	dagedgev *edges;
	u32list  *trash;
	String   *alns[2];
	int      W, M, X, I, D, E;
	f4i      pM, pX, pI, pD; // log(prob_M)
	double   ref_penalty, alt_penalty; // 0.5 and 0.2
	double   cns_score;
	uint32_t cns_head, backbone_size;
} DAGCNS;

static inline DAGCNS* init_dagcns(int W, int M, int X, int I, int D, int E, f4i pM, f4i pX, f4i pI, f4i pD){
	DAGCNS *g;
	g = malloc(sizeof(DAGCNS));
	g->cns = init_u8list(1024);
	g->deps = init_u32list(1024);
	g->nodes = init_dagnodev(1024);
	g->edges = init_dagedgev(1024);
	g->trash = init_u32list(1024);
	g->alns[0] = init_string(1024);
	g->alns[1] = init_string(1024);
	g->cns_score = 0;
	g->cns_head = 0xFFFFFFFFU;
	g->backbone_size = 0;
	g->W = W;
	g->M = M;
	g->X = X;
	g->I = I;
	g->D = D;
	g->E = E;
	g->pM = pM;
	g->pX = pX;
	g->pI = pI;
	g->pD = pD;
	g->ref_penalty = 0.5;
	g->alt_penalty = 0.2;
	return g;
}

static inline void free_dagcns(DAGCNS *g){
	free_dagnodev(g->nodes);
	free_dagedgev(g->edges);
	free_u32list(g->trash);
	free_u8list(g->cns);
	free_u32list(g->deps);
	free_string(g->alns[0]);
	free_string(g->alns[1]);
	free(g);
}

static inline void reset_dagcns(DAGCNS *g){
	clear_dagnodev(g->nodes);
	clear_dagedgev(g->edges);
	clear_u32list(g->trash);
	clear_u8list(g->cns);
	clear_u32list(g->deps);
	g->cns_score = 0;
	g->cns_head = 0xFFFFFFFFU;
	g->backbone_size = 0;
}

static uint32_t prepare_node_dagcns(DAGCNS *g, uint32_t pos, uint8_t base){
	dagnode_t *n;
	n = next_ref_dagnodev(g->nodes);
	n->pos = pos;
	n->base = base;
	n->cns = 0;
	n->aux  = 0;
	n->visit = 0;
	n->edges[0] = 0xFFFFFFFFU;
	n->edges[1] = 0xFFFFFFFFU;
	n->fw_edge = NODE_MAX_FW_EDGE;
	return g->nodes->size - 1;
}

static inline dagedge_t* find_edge_by_node_dagcns(DAGCNS *g, uint32_t nid1, uint32_t nid2, int dir){
	dagnode_t *n;
	dagedge_t *e;
	n = ref_dagnodev(g->nodes, nid1);
	if(n->edges[dir] != 0xFFFFFFFFU){
		e = ref_dagedgev(g->edges, n->edges[dir]);
		while(1){
			if(e->nodes[dir] == nid2) return e;
			if(e->links[dir] == 0xFFFFFFFFU) break;
			e = ref_dagedgev(g->edges, e->links[dir]);
		}
	}
	return NULL;
}

static inline dagedge_t* add_edge_dagcns(DAGCNS *g, uint32_t nid1, uint32_t nid2, int dir){
	dagnode_t *n;
	dagedge_t *e;
	uint32_t eid;
	n = ref_dagnodev(g->nodes, nid1);
	if(pop_u32list(g->trash, &eid)){
		e = ref_dagedgev(g->edges, eid);
	} else {
		eid = g->edges->size;
		e = next_ref_dagedgev(g->edges);
	}
	e->nodes[!dir] = nid1;
	e->nodes[dir] = nid2;
	e->cov = 1;
	e->score = 0;
	e->visit = 0;
	e->closed = 0;
	e->cns = 0;
	e->alt = 0;
	e->links[dir] = n->edges[dir];
	n->edges[dir] = eid;
	n = ref_dagnodev(g->nodes, nid2);
	e->links[!dir] = n->edges[!dir];
	n->edges[!dir] = eid;
	return e;
}

static inline dagedge_t* prepare_edge_dagcns(DAGCNS *g, uint32_t nid1, uint32_t nid2, int dir){
	dagedge_t *e;
	e = find_edge_by_node_dagcns(g, nid1, nid2, dir);
	if(e){ e->cov ++; return e; }
	return add_edge_dagcns(g, nid1, nid2, dir);
}

static inline void gen_pregraph_dagcns(DAGCNS *g){
	dagedge_t *e;
	uint32_t i;
	clear_dagnodev(g->nodes);
	clear_dagedgev(g->edges);
	clear_u32list(g->trash);
	clear_u32list(g->deps);
	g->backbone_size = g->cns->size;
	for(i=0;i<g->cns->size;i++){
		push_u32list(g->deps, 0);
		prepare_node_dagcns(g, i, g->cns->buffer[i]);
		if(i){ // make sure the graph is conntective even the alignment is partial
			e = add_edge_dagcns(g, i - 1, i, 0);
			e->cov = 0;
		}
	}
}

static inline int remove_edge_dagcns(DAGCNS *g, uint32_t eid){
	dagnode_t *n;
	dagedge_t *e;
	uint32_t i, lst;
	for(i=0;i<2;i++){
		e = ref_dagedgev(g->edges, eid);
		lst = e->links[i];
		n = ref_dagnodev(g->nodes, e->nodes[!i]);
		if(n->edges[i] == eid){
			n->edges[i] = lst;
		} else if(n->edges[i] == 0xFFFFFFFFU){
			//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
			return 0;
		} else {
			e = ref_dagedgev(g->edges, n->edges[i]);
			while(1){
				if(e->links[i] == eid){
					e->links[i] = lst; break;
				} else if(e->links[i] == 0xFFFFFFFFU){
					//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
					return 0;
				} else {
					e = ref_dagedgev(g->edges, e->links[i]);
				}
			}
		}
	}
	push_u32list(g->trash, eid);
	return 1;
}

#define MIN_SCORE	-0x0FFFFFFF

static inline f8i log_sum(f8i a, f8i b){
	f8i c;
	c = num_max(a, b);
	if(c - a >= 10 || c - b >= 10) return c;
	return logl(expl(a - c) + expl(b - c)) + c;
}

typedef struct {
	int x, y, d;
} bdp_beg_t;
define_list(bdpbegv, bdp_beg_t);

static inline void fprint_dot_bdpgraph(GEGraph *g, bdpnodev *bnodes, bdpedgev *bedges, char *prefix, char *suffix){
	static const char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};
	FILE *out;
	ge_node_t *n;
	bdp_node_t *bn;
	ge_edge_t *e;
	bdp_edge_t *be;
	u8i i;
	out = open_file_for_write(prefix, suffix, 1);
	fprintf(out, "digraph {\nnode [shape=record]\nrankdir=LR\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_genodev(g->nodes, i);
		if(n->closed) continue;
		bn = ref_bdpnodev(bnodes, offset_genodev(g->nodes, n));
		fprintf(out, " N%llu [label=\"{N%llu|%d|%d}|{%.4Lf|%.4Lf}\"]\n", i, i, bn->x, bn->y, bn->probs[0], bn->probs[1]);
	}
	for(i=1;i<g->edges->size;i++){
		e = ref_geedgev(g->edges, i);
		if(e->closed) continue;
		be = ref_bdpedgev(bedges, offset_geedgev(g->edges, e));
		fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d:%c:%.4Lf\" color=%s]\n", (u8i)e->node1, (u8i)e->node2, "+-"[e->dir1], "+-"[e->dir2], e->cov, "MIDX"[be->cigar], be->prob, colors[e->dir1][e->dir2]);
		fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d:%c:%.4Lf\" color=%s]\n", (u8i)e->node2, (u8i)e->node1, "-+"[e->dir2], "-+"[e->dir1], e->cov, "MIDX"[be->cigar], be->prob, colors[!e->dir2][!e->dir1]);
	}
	fprintf(out, "}\n");
	fclose(out);
}

static inline u4i dp_matrix2alignment_graph(DAGCNS *dag, u1i *query, u1i *z, int x, int y, GEGraph *g, bdpnodev *bnodes, bdpedgev *bedges){
	ge_node_t *gn, *gn2;
	ge_edge_ref_t *gf;
	ge_edge_t *ge;
	bdp_node_t *bn, *bn2;
	bdp_edge_t *be;
	bdp_beg_t BEG;
	UUhash *bnhash;
	UUhash_t *u;
	bdpbegv *stack;
	f8i cigar_probs[4], curator;
	u8v *idxs;
	u1i *target;
	u4i i, k, nidx, nlst, nbegs[2], visit, ret;
	int d, f, n_col, ops[2], base, exists;
	cigar_probs[0] = dag->pM;
	cigar_probs[1] = dag->pI;
	cigar_probs[2] = dag->pD;
	cigar_probs[3] = dag->pX;
	n_col = dag->cns->size;
	target = dag->cns->buffer;
	bdp_set_callbacks_gegraph(g, bnodes, bedges);
	reset_gegraph(g);
	stack = init_bdpbegv(4);
	push_bdpbegv(stack, (bdp_beg_t){x, y, 0});
	idxs = init_u8v(4);
	nbegs[0] = nbegs[1] = 0;
	bnhash = init_UUhash(1023);
	ret = 0;
	while(stack->size){
		BEG = stack->buffer[--stack->size];
		u = prepare_UUhash(bnhash, (((b8i)BEG.x) << 32) | BEG.y, &exists);
		if(exists){
			nidx = u->val;
		} else {
			gn = add_node_gegraph(g);
			nidx = offset_genodev(g->nodes, gn);
			u->key = (((b8i)BEG.x) << 32) | BEG.y;
			u->val = nidx;
			bn = ref_bdpnodev(bnodes, nidx);
			bn->x = BEG.x; bn->y = BEG.y;
			bn->dnidx = MAX_VALUE_U4;
		}
		nlst = nidx;
		x = BEG.x; y = BEG.y;
		d = BEG.d;
		f = 0;
		clear_u8v(idxs);
		while(x >= 0 && y >= 0){
			d = (z[x * n_col + y] >> (d << 1)) & 0x03;
			if(d == 0){
				if(query[x] == target[y]){
					if(z[x * n_col + y] & (1 << 6)){
						f = 1;
						break;
					} else if(BEG.d == 0){
						z[x * n_col + y] |= 1 << 6;
						push_bdpbegv(stack, (bdp_beg_t){x, y, 1});
						push_bdpbegv(stack, (bdp_beg_t){x, y, 2});
					}
					ops[0] = 0; ops[1] = 3;
				} else {
					ops[0] = 1; ops[1] = 2;
				}
				//x --; y --; 
			} else if(d == 1){
				ops[0] = 1; ops[1] = 3;
				//x --;
			} else {
				ops[0] = 2; ops[1] = 3;
				//y --;
			}
			for(i=0;i<2;i++){
				if(ops[i] == 3) break;
				base = query[x];
				if(ops[i] == 0){ x --; y --; }
				else if(ops[i] == 1){ x --; }
				else { y --; }
				u = prepare_UUhash(bnhash, (((b8i)x) << 32) | y, &exists);
				if(exists){
					nidx = u->val;
				} else {
					gn = add_node_gegraph(g);
					nidx = offset_genodev(g->nodes, gn);
					u->key = (((b8i)x) << 32) | y;
					u->val = nidx;
					bn = ref_bdpnodev(bnodes, nidx);
					bn->x = x; bn->y = y;
					bn->dnidx = MAX_VALUE_U4;
				}
				ge = prepare_edge_gegraph(g, nlst, 1, nidx, 1, &exists);
				if(exists){
					ge->cov ++;
				} else {
					ge->cov = 1;
					be = ref_bdpedgev(bedges, offset_geedgev(g->edges, ge));
					be->cigar = ops[i];
					be->base  = base;
				}
				nlst = nidx;
				if(BEG.d) push_u8v(idxs, offset_geedgev(g->edges, ge));
			}
		}
		if(BEG.d){
			if(f == 0){
				// not a bubble
				for(i=0;i<idxs->size;i++){
					ge = ref_geedgev(g->edges, idxs->buffer[i]);
					ge->cov --;
					if(ge->cov == 0) cut_edge_gegraph(g, ge);
				}
			} else ret ++;
		} else {
			ret ++;
			nbegs[0] = g->nodes->size > 1? g->nodes->size - 2 : 0;
			nbegs[1] = 0;
		}
	}
	free_bdpbegv(stack);
	free_UUhash(bnhash);
	u4i cnt = 0;
	for(i=0;i<g->nodes->size;i++){
		gn = ref_genodev(g->nodes, i);
		if(gn->edges[0].cnt || gn->edges[1].cnt){ cnt++; continue; }
		gn->closed = 1;
	}
	// calculate probabilities (forward + backward)
	if(nbegs[0] == nbegs[1]){
		free_u8v(idxs);
		return ret;
	}
	for(k=0;k<2;k++){
		bn = ref_bdpnodev(bnodes, nbegs[k]);
		bn->probs[k] = 0;
		clear_u8v(idxs);
		push_u8v(idxs, nbegs[k]);
		visit = k + 1;
		while(idxs->size){
			gn = ref_genodev(g->nodes, idxs->buffer[idxs->size - 1]);
			bn = ref_bdpnodev(bnodes, idxs->buffer[idxs->size - 1]);
			idxs->size --;
			geg_beg_iter_edges(g, gn, k, gf, ge);
			if(ge->closed) continue;
			be = ref_bdpedgev(bedges, offset_geedgev(g->edges, ge));
			gn2 = ref_genodev(g->nodes, gf->flg? ge->node1 : ge->node2);
			bn2 = ref_bdpnodev(bnodes, gf->flg? ge->node1 : ge->node2);
			if(gn2->bt_visit == visit){
				bn2->probs[k] = log_sum(bn2->probs[k], bn->probs[k] + cigar_probs[be->cigar]); // p1 + p2
			} else {
				gn2->bt_visit = visit;
				gn2->unvisit = gn2->edges[!k].cnt;
				bn2->probs[k] = bn->probs[k] + cigar_probs[be->cigar]; // p1 * p2
			}
			if(gn2->unvisit) gn2->unvisit --;
			if(gn2->unvisit == 0){
				push_u8v(idxs, offset_genodev(g->nodes, gn2));
			}
			geg_end_iter_edges();
		}
	}
	free_u8v(idxs);
	// calculate edge prob
	curator = -FLT_MAX;
	for(i=1;i<g->edges->size;i++){
		ge = ref_geedgev(g->edges, i);
		if(ge->closed) continue;
		be = ref_bdpedgev(bedges, i);
		be->prob = cigar_probs[be->cigar] + ref_bdpnodev(bnodes, ge->node1)->probs[ge->dir1] + ref_bdpnodev(bnodes, ge->node2)->probs[!ge->dir2];
		if(be->prob > curator) curator = be->prob;
	}
	for(i=1;i<g->edges->size;i++){
		ge = ref_geedgev(g->edges, i);
		if(ge->closed) continue;
		be = ref_bdpedgev(bedges, i);
		be->prob = expl(be->prob - curator);
	}
	//fprint_dot_bdpgraph(g, bnodes, bedges, "geg.dot", NULL);
	return nbegs[0];
}

static inline u4i branched_dynamic_programming_alignment(DAGCNS *g, u1i *query, int ql, GEGraph *geg, bdpnodev *bnodes, bdpedgev *bedges, u1v *mem_buffer){
	u4i nbeg;
	int *rh, *re;
	u1i *z, *zi, d, *target;
	int tl, n_col;
	int i, j, jc, jb, je, h1, h, m, e, f, t;
	int mi, mj, max;
	int W, M, X, I, D, E;
	tl = g->cns->size;
	target = g->cns->buffer;
	if(ql <= 0 || tl <= 0){ return 0; }
	n_col = tl;
	encap_u1v(mem_buffer, kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x(((long long)ql) * n_col));
	rh = (int*)(mem_buffer->buffer + mem_buffer->size);
	re = (int*)(mem_buffer->buffer + mem_buffer->size + kswx_roundup8x((tl + 2) * sizeof(int)));
	z  =       (mem_buffer->buffer + mem_buffer->size + kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x((tl + 2) * sizeof(int)));
	W = g->W;
	M = g->M;
	X = g->X;
	I = g->I;
	D = g->D;
	E = g->E;
	// banded DP, global alignment
	rh[0] = 0;
	re[1] = 0 + D + E;
	for(j=2;j<=tl&&j<=W;j++) rh[j] = rh[j-1] + E;
	for(;j<tl;j++) rh[j] = MIN_SCORE;
	for(j=0;j<=tl;j++) re[j] = MIN_SCORE;
	mi = mj = 0;
	max = MIN_SCORE;
	for(i=0;i<ql;i++){
		f = MIN_SCORE;
		jc = (i * tl / ql); // in case of biased gaps
		jb = jc > W? jc - W : 0;
		je = jc + W + 1 < tl? jc + W + 1 : tl;
		h1 = jb == 0? (I + E * (i + 1)) : MIN_SCORE;
		zi = &z[i * n_col];
		for(j=jb;j<je;j++){
			m = rh[j] + ((query[i] == target[j])? M : X);
			rh[j] = h1;
			e = re[j];
			d = m >= e? 0 : 1;
			h = m >= e? m : e;
			d = h >= f? d : 2;
			h = h >= f? h : f;
			h1 = h;
			t = m + I + E;
			e = e + E;
			d |= e > t? 1<<2 : 0;
			e = e > t? e : t;
			re[j] = e;
			t = m + D + E;
			f = f + E;
			d |= f > t? 2<<4 : 0;
			f = f > t? f : t;
			zi[j] = d;
		}
		rh[j] = h1; re[j] = MIN_SCORE;
		if(i + 1 == ql){
			for(j=jb;j<je;j++){
				if(rh[j + 1] > max){
					max = rh[j + 1];
					mi = i; mj = j;
				}
			}
		} else if(je == tl){
			if(h1 > max){
				max = h1;
				mi = i; mj = tl - 1;
			}
		}
	}
	if(max == MIN_SCORE) return 0;
	nbeg = dp_matrix2alignment_graph(g, query, z, mi, mj, geg, bnodes, bedges);
	return nbeg;
}

static inline void bdpgraph2dagcns(DAGCNS *dg, GEGraph *gg, bdpnodev *bnodes, bdpedgev *bedges, u4i nbeg, bdplinkv *stack){
	dagnode_t *dn, *dn2;
	dagedge_t *de;
	ge_node_t *gn, *gn2;
	ge_edge_ref_t *gf;
	ge_edge_t *ge;
	bdp_node_t *bn, *bn2;
	bdp_edge_t *be;
	bdp_link_t T;
	u4i i, j, beg, end;
	int open;
	for(i=0;i<gg->nodes->size;i++) gg->nodes->buffer[i].bt_visit = 0;
	clear_bdplinkv(stack);
	push_bdplinkv(stack, (bdp_link_t){nbeg, 0xFFFFFFFFU});
	beg = bnodes->buffer[nbeg].y;
	end = beg;
	open = 0;
	while(stack->size){
		T = stack->buffer[--stack->size];
		gn = ref_genodev(gg->nodes, T.gnidx);
		bn = ref_bdpnodev(bnodes, T.gnidx);
		if(bn->y > (int)end) end = bn->y;
		dn = (T.dnidx == 0xFFFFFFFFU)? NULL : ref_dagnodev(dg->nodes, T.dnidx);
		geg_beg_iter_edges(gg, gn, 0, gf, ge);
		if(ge->closed) continue;
		be = ref_bdpedgev(bedges, offset_geedgev(gg->edges, ge));
		gn2 = ref_genodev(gg->nodes, gf->flg? ge->node1 : ge->node2);
		bn2 = ref_bdpnodev(bnodes, gf->flg? ge->node1 : ge->node2);
		if(gn2->bt_visit == 0){
			gn2->bt_visit = 1;
			open ++;
			gn2->unvisit = gn2->edges[1].cnt;
		}
		if(gn2->unvisit) gn2->unvisit --;
		if(gn2->unvisit == 0){
			open --;
		}
		if(open < 0){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
		if(bn2->dnidx == MAX_VALUE_U4){
			if(be->cigar == 2){
				dn2 = dn;
			} else if(dn && (open || be->cigar == 1)){
				dn2 = ref_dagnodev(dg->nodes, prepare_node_dagcns(dg, bn2->y, be->base));
				dn  = dn? ref_dagnodev(dg->nodes, T.dnidx) : NULL;
			} else if((dn || open == 0) && be->cigar == 0){ // cigar == 0 && open == 0
				dn2 = ref_dagnodev(dg->nodes, bn2->y);
			} else dn2 = NULL;
			bn2->dnidx = dn2? offset_dagnodev(dg->nodes, dn2) : MAX_VALUE_U4;
		} else {
			dn2 = ref_dagnodev(dg->nodes, bn2->dnidx);
		}
		if(dn && dn2 && dn != dn2){
			de = prepare_edge_dagcns(dg, offset_dagnodev(dg->nodes, dn), offset_dagnodev(dg->nodes, dn2), 0);
			de->score += be->prob;
		}
		if(gn2->unvisit == 0){
			push_bdplinkv(stack, (bdp_link_t){offset_genodev(gg->nodes, gn2), dn2? offset_dagnodev(dg->nodes, dn2) : 0xFFFFFFFFU});
		}
		geg_end_iter_edges();
	}
	for(j=beg;j<end;j++) dg->deps->buffer[j] ++;
	return;
}

static inline void merge_nodes_core_dagcns(DAGCNS *g, uint32_t nid, u32list *stack, u32list *cache[4], int dir){
	dagnode_t *n0, *n2, *n;
	dagedge_t *e, *e2, *e1;
	uint32_t base, eid, nid1, i, ret;
	clear_u32list(stack);
	push_u32list(stack, nid);
	ret = 0;
	while(pop_u32list(stack, &nid)){
		n0 = ref_dagnodev(g->nodes, nid);
		if((eid = n0->edges[dir]) == 0xFFFFFFFFU) continue;
		clear_u32list(cache[0]);
		clear_u32list(cache[1]);
		clear_u32list(cache[2]);
		clear_u32list(cache[3]);
		while(1){
			e = ref_dagedgev(g->edges, eid);
			n = ref_dagnodev(g->nodes, e->nodes[dir]);
			e2 = ref_dagedgev(g->edges, n->edges[!dir]);
			if(e2->links[!dir] == 0xFFFFFFFFU){ // check whether there is only one edge from n -(!dir)-> n0
				push_u32list(cache[n->base], eid);
			}
			if((eid = e->links[dir]) == 0xFFFFFFFFU) break;
		}
		for(base=0;base<4;base++){
			if(cache[base]->size < 2) continue;
			for(i=0;i<cache[base]->size;i++) ref_dagedgev(g->edges, cache[base]->buffer[i])->visit = 1;
			e1 = ref_dagedgev(g->edges, cache[base]->buffer[0]);
			n = ref_dagnodev(g->nodes, e1->nodes[dir]);
			eid = n->edges[dir];
			nid1 = e1->nodes[dir];
			while(eid != 0xFFFFFFFFU){
				e = ref_dagedgev(g->edges, eid);
				e->visit = 1;
				eid = e->links[dir];
			}
			for(i=1;i<cache[base]->size;i++){
				e2 = ref_dagedgev(g->edges, cache[base]->buffer[i]);
				n2 = ref_dagnodev(g->nodes, e2->nodes[dir]);
				e1->cov += e2->cov;
				e1->score += e2->score;
				remove_edge_dagcns(g, cache[base]->buffer[i]);
				eid = n2->edges[dir];
				while(eid != 0xFFFFFFFFU){
					e2 = ref_dagedgev(g->edges, eid);
					e  = prepare_edge_dagcns(g, nid1, e2->nodes[dir], dir);
					{
						e1 = ref_dagedgev(g->edges, cache[base]->buffer[0]); // memory referred by e1 may be freed in prepare_edge_dagcns
					}
					e->cov = e->cov - 1 + e2->cov;
					e->score = e->score + e2->score;
					e->visit = 1;
					eid = e2->links[dir];
				}
				eid = n2->edges[dir];
				while(eid != 0xFFFFFFFFU){
					e2 = ref_dagedgev(g->edges, eid);
					remove_edge_dagcns(g, eid); // e2->links retain the same values after removing
					eid = e2->links[dir];
				}
			}
			//n = ref_dagnodev(g->nodes, e1->nodes[dir]);
			//eid = n->edges[dir];
			//if(eid != 0xFFFFFFFFU && g->edges->buffer[eid].links[dir] == 0xFFFFFFFFU) continue; // we had merged a bubble branch1:A->C->T, branch2:A->C->T
			push_u32list(stack, nid1);
		}
	}
}

static inline int has_non_visited_edge_dagcns(DAGCNS *g, uint32_t nid, int dir){
	dagnode_t *n;
	dagedge_t *e;
	uint32_t eid;
	n = ref_dagnodev(g->nodes, nid);
	eid = n->edges[dir];
	while(eid != 0xFFFFFFFFU){
		e = ref_dagedgev(g->edges, eid);
		if(e->visit == 0) return 1;
		eid = e->links[dir];
	}
	return 0;
}

static inline void print_local_dot_dagcns(DAGCNS *g, uint32_t nid, int distance, FILE *out){
	u32list *stack;
	u32hash *hash;
	dagnode_t *n, *n1, *n2;
	dagedge_t *e;
	uint32_t id1, id2, eid, *u;
	int lo, hi, dir, exists;
	n = ref_dagnodev(g->nodes, nid);
	stack = init_u32list(32);
	hash = init_u32hash(1023);
	lo = n->pos - distance;
	hi = n->pos + distance;
	push_u32list(stack, nid);
	put_u32hash(hash, nid);
	fprintf(out, "digraph {\nrankdir=LR\n");
	while(stack->size){
		id1 = stack->buffer[--stack->size];
		n1 = ref_dagnodev(g->nodes, id1);
		for(dir=0;dir<1;dir++){
			eid = n1->edges[dir];
			while(eid != 0xFFFFFFFFU){
				e = ref_dagedgev(g->edges, eid);
				id2 = e->nodes[dir];
				n2 = ref_dagnodev(g->nodes, id2);
				fprintf(out, "N%d_%d_%c -> N%d_%d_%c [label=\"%d:%0.6f\" color=%s]\n", id1, n1->pos, "ACGT"[n1->base], id2, n2->pos, "ACGT"[n2->base], e->cov, e->score, e->visit? "blue" : "black");
				if(n2->pos >= lo && n2->pos <= hi){
					u = prepare_u32hash(hash, id2, &exists);
					if(exists){
					} else {
						*u = id2;
						push_u32list(stack, id2);
					}
				}
				eid = e->links[dir];
			}
		}
	}
	fprintf(out, "}\n");
	fflush(out);
	free_u32list(stack);
	free_u32hash(hash);
}

static inline void fprint_local_dot_dagcns(DAGCNS *g, u4i nid, int level, char *prefix, char *suffix){
	FILE *out;
	out = open_file_for_write(prefix, suffix, 1);
	print_local_dot_dagcns(g, nid, level, out);
	fclose(out);
}

static inline void merge_nodes_dagcns(DAGCNS *g){
	dagnode_t *n;
	dagedge_t *e;
	u32list *stack, *cache[4];
	u32fifo *queue;
	uint32_t i, nid, eid;
	stack = init_u32list(1024);
	cache[0] = init_u32list(4);
	cache[1] = init_u32list(4);
	cache[2] = init_u32list(4);
	cache[3] = init_u32list(4);
	for(i=0;i<g->edges->size;i++) g->edges->buffer[i].visit = 0;
	for(i=0;i<g->nodes->size;i++) g->nodes->buffer[i].visit = 0;
	queue = init_u32fifo();
	for(i=0;i<g->nodes->size;i++){
		n = ref_dagnodev(g->nodes, i);
		if(n->edges[1] != 0xFFFFFFFFU) continue;
		push_u32fifo(queue, i);
	}
	//dagcns_debug = 2;
	while(pop_u32fifo(queue, &nid)){
		if(dagcns_debug > 1) fprintf(stdout, "\npop(%u) %u\n", (u4i)queue->size, nid);
		n = ref_dagnodev(g->nodes, nid);
		if(n->visit) continue;
		n->visit = 1;
		merge_nodes_core_dagcns(g, nid, stack, cache, 1);
		merge_nodes_core_dagcns(g, nid, stack, cache, 0);
		n = ref_dagnodev(g->nodes, nid);
		eid = n->edges[0];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			e->visit = 1;
			eid = e->links[0];
		}
		eid = n->edges[0];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			if(!has_non_visited_edge_dagcns(g, e->nodes[0], 1)){
				if(dagcns_debug > 1) fprintf(stdout, "push %u\n", e->nodes[0]);
				push_u32fifo(queue, e->nodes[0]);
			}
			eid = e->links[0];
		}
	}
	free_u32fifo(queue);
	free_u32list(stack);
	free_u32list(cache[0]);
	free_u32list(cache[1]);
	free_u32list(cache[2]);
	free_u32list(cache[3]);
}

static inline void print_seq_dagcns(DAGCNS *g, FILE *out){
	char buffer[100];
	uint32_t i, j;
	for(i=j=0;i<g->cns->size;i++){
		buffer[j++] = bit_base_table[g->cns->buffer[i]];
		if(j == 99){
			buffer[j] = '\0';
			fprintf(out, "%s", buffer);
			j = 0;
		}
	}
	buffer[j] = '\0';
	fprintf(out, "%s", buffer);
}

static inline void gen_consensus_dagcns(DAGCNS *g, u32list *map){
	dagnode_t *n1, *n2;
	dagedge_t *e;
	u32fifo *queue;
	uint32_t i, lst, nid, eid, best_e;
	f8i best_s, score;
	queue = init_u32fifo();
	if(queue == NULL){ // un-reachable, but is used to call fprint_local_dot_dagcns in gdb Debug
		fprint_local_dot_dagcns(g, 0, 10, "test.dot", NULL);
	}
	for(i=0;i<g->nodes->size;i++){
		n1 = ref_dagnodev(g->nodes, i);
		if(n1->edges[0] == 0xFFFFFFFFU && n1->edges[1] != 0xFFFFFFFFU){
			push_u32fifo(queue, i);
			n1->fw_edge = NODE_MAX_FW_EDGE;
			n1->aux = 0;
		}
	}
	for(i=0;i<g->edges->size;i++) g->edges->buffer[i].visit = 0;
	while(pop_u32fifo(queue, &nid)){
		best_s = - FLT_MAX;
		best_e = 0xFFFFFFFFU;
		n1 = ref_dagnodev(g->nodes, nid);
		eid = n1->edges[0];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			n2 = ref_dagnodev(g->nodes, e->nodes[0]);
			if(e->nodes[0] < g->backbone_size){
				//score = n2->aux + e->cov - g->ref_penalty * g->deps->buffer[n1->pos];
				score = n2->aux + e->score - g->ref_penalty * g->deps->buffer[n1->pos];
			} else {
				//score = n2->aux + e->cov - g->alt_penalty * g->deps->buffer[n1->pos];
				score = n2->aux + e->score - g->alt_penalty * g->deps->buffer[n1->pos];
			}
			if(score > best_s){
				best_s = score;
				best_e = eid;
			}
			eid = e->links[0];
		}
		if(best_s > - FLT_MAX) n1->aux = best_s;
		n1->fw_edge = best_e;
		eid = n1->edges[1];
		while(eid != 0xFFFFFFFFU){
			e = ref_dagedgev(g->edges, eid);
			e->visit = 1;
			if(!has_non_visited_edge_dagcns(g, e->nodes[1], 0)){
				push_u32fifo(queue, e->nodes[1]);
			}
			eid = e->links[1];
		}
	}
	free_u32fifo(queue);
	clear_u8list(g->cns);
	clear_u32list(g->deps);
	if(map) clear_u32list(map);
	g->cns_head = 0;
	if(g->nodes->size == 0) return;
	n1 = ref_dagnodev(g->nodes, g->cns_head);
	g->cns_score = n1->aux;
	n1->cns = 1;
	lst = 0;
	if(map && g->cns_head < g->backbone_size){
		while(lst < g->cns_head){ push_u32list(map, g->cns->size); lst ++; }
	}
	push_u8list(g->cns, n1->base);
	push_u32list(g->deps, 0);
	while(n1->fw_edge != NODE_MAX_FW_EDGE){
		e = ref_dagedgev(g->edges, n1->fw_edge);
		e->cns = 1;
		if(map && e->nodes[0] < g->backbone_size){
			while(lst < e->nodes[0]){ push_u32list(map, g->cns->size); lst ++; }
		}
		n1 = ref_dagnodev(g->nodes, e->nodes[0]);
		n1->cns = 1;
		push_u8list(g->cns, n1->base);
		push_u32list(g->deps, 0);
	}
	if(map) while(lst <= g->backbone_size){ push_u32list(map, g->cns->size); lst ++; }
}

#endif
