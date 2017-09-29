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
#include "string.h"
#include "list.h"
#include "hashset.h"
#include "kswx.h"
#include "hzm_aln.h"
#include <float.h>

#ifndef PO_MSA_RJ_H
#define PO_MSA_RJ_H

#define BLOCK_SIZE	32
#define NODE_NULL	0
#define EDGE_NULL	0
#define LINK_NULL	0
#define POS_NULL	0x0FFFFFFFU
#define DP_NULL	0
#define NODE_MAX_EDGES	32
#define SCORE_MIN	-0x3FFFFFFF
#define ROW_MAX_OFF	2
#define SC_IDX_NULL	0xFFFFFFFFU

typedef struct {
	uint32_t pos:28, base:2, cns:1, ref:1;
	uint32_t n_in:16, n_ot:16;
	uint32_t edges;
	uint32_t aligned[2]; // nodes are aligned, but with different base
	uint32_t cnt:16, cov:16;
	uint32_t bt;
	uint32_t aux;
	int score;
} pognode_t;
define_list(pognodev, pognode_t);

typedef struct {
	uint32_t node;
	uint32_t cov:28, visit:1, closed:1, cns:1, ref:1;
	uint32_t next;
	uint32_t link;
} pogedge_t;
define_list(pogedgev, pogedge_t);

typedef struct {
	uint32_t read;
	uint32_t next;
} poglink_t;
define_list(poglinkv, poglink_t);

typedef struct {
	int h, e;
	uint64_t hbt:26, ebt:26, d:8, mg:4;
} pocell_t;
define_list(pocellv, pocell_t);

typedef struct {
	uint32_t node;
	uint32_t qx, qy;
	int score;
	uint64_t coff:43, woff:2, winc:2, n_in:16, ready:1;
} podp_t;
define_list(podpv, podp_t);

typedef struct {
	uint32_t node;
	uint32_t pos;
	uint32_t btn, bte;
} pomat_t;
define_list(pomatv, pomat_t);

typedef struct {
	u1v *cns;
	int cns_score;
	u1v *ref;
	u4v *map;
	pognodev *nodes;
	pogedgev *edges;
	poglinkv *links;
	HZMAux   *aux;
	uint32_t beginning;
	uint32_t W, MW;
	int      has_merge; // whether consider homopolymer merge error
	pocellv  *cells;
	podpv    *dps;
	pomatv   *mats;
	u64hash  *hash;
	u4v      *heap;
	uint32_t aln_score;
	uint32_t counter;
} POMSA;

static inline POMSA* init_pomsa(){
	POMSA *g;
	g = malloc(sizeof(POMSA));
	g->cns = init_u1v(1024);
	g->ref = init_u1v(1024);
	g->cns_score = SCORE_MIN;
	g->map = init_u4v(1024);
	g->nodes = init_pognodev(1024);
	g->edges = init_pogedgev(1024);
	g->links = init_poglinkv(1024);
	g->aux = init_hzmaux();
	g->aux->has_alignment = 0;
	g->beginning = 0;
	g->W  = 100;
	g->MW = 800;
	g->has_merge = 0;
	g->cells = init_pocellv(1024);
	g->dps   = init_podpv(1024);
	g->mats  = init_pomatv(1024);
	g->hash  = init_u64hash(1023);
	g->heap  = init_u4v(1024);
	g->aln_score = 0;
	g->counter = 0;
	return g;
}

// To change alignment parameters, modify g->aux->***

static inline void free_pomsa(POMSA *g){
	free_u1v(g->cns);
	free_u1v(g->ref);
	free_u4v(g->map);
	free_pognodev(g->nodes);
	free_pogedgev(g->edges);
	free_poglinkv(g->links);
	free_hzmaux(g->aux);
	free_pocellv(g->cells);
	free_podpv(g->dps);
	free_pomatv(g->mats);
	free_u64hash(g->hash);
	free_u4v(g->heap);
	free(g);
}

static inline uint32_t add_node_pomsa(POMSA *g, uint32_t pos, uint8_t base, int is_ref){
	pognode_t *n;
	pogedge_t *e;
	uint32_t nid;
	nid = g->nodes->size;
	n = next_ref_pognodev(g->nodes);
	n->pos = pos;
	n->base = base;
	n->cns = 0;
	n->ref = is_ref;
	n->aux = 0;
	n->cnt = 0;
	n->cov = 0;
	n->bt  = NODE_NULL;
	n->n_in = 0;
	n->n_ot = 0;
	n->edges = g->edges->size;
	e = next_ref_pogedgev(g->edges);
	memset(e, 0, sizeof(pogedge_t));
	n->aligned[0] = NODE_NULL;
	n->aligned[1] = NODE_NULL;
	return nid;
}

static inline uint32_t add_edge_pomsa(POMSA *g, uint32_t node1, uint32_t node2, int is_ref){
	pognode_t *n;
	pogedge_t *e;
	if(node2 == NODE_NULL){
		fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	if(node1 == node2) return EDGE_NULL;
	encap_pogedgev(g->edges, 1); // make sure ptr-"e" is not broken down by realloc
	n = ref_pognodev(g->nodes, node1);
	//ref_pognodev(g->nodes, node1)->cov ++;
	//ref_pognodev(g->nodes, node2)->cov ++;
	e = ref_pogedgev(g->edges, n->edges);
	while(e->node != NODE_NULL){
		if(e->node == node2){
			e->cov ++;
			e->ref |= is_ref;
			return e - g->edges->buffer;
		}
		e = ref_pogedgev(g->edges, e->next);
	}
	if(ref_pognodev(g->nodes, node1)->n_ot == NODE_MAX_EDGES) return EDGE_NULL;
	if(ref_pognodev(g->nodes, node2)->n_in == NODE_MAX_EDGES) return EDGE_NULL;
	ref_pognodev(g->nodes, node1)->n_ot ++;
	ref_pognodev(g->nodes, node2)->n_in ++;
	e->node = node2;
	e->cov = 1;
	e->visit = 0;
	e->closed = 0;
	e->ref = is_ref;
	e->cns = 0;
	e->link = LINK_NULL;
	e->next = g->edges->size;
	memset(next_ref_pogedgev(g->edges), 0, sizeof(pogedge_t));
	return e - g->edges->buffer;
}

static inline pogedge_t* find_edge_pomsa(POMSA *g, uint32_t node1, uint32_t node2){
	pognode_t *n;
	pogedge_t *e;
	if(node2 == NODE_NULL){
		fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	if(node1 == node2) return NULL;
	n = ref_pognodev(g->nodes, node1);
	e = ref_pogedgev(g->edges, n->edges);
	while(e->node != NODE_NULL){
		if(e->node == node2) return e;
		e = ref_pogedgev(g->edges, e->next);
	}
	return NULL;
}

static inline int add_edge_read_pomsa(POMSA *g, pogedge_t *e, uint32_t read){
	poglink_t *l, *p;
	uint32_t link;
	l = next_ref_poglinkv(g->links);
	l->read = read;
	link = e->link;
	p = NULL;
	while(link != LINK_NULL){
		if(ref_poglinkv(g->links, link)->read == read) return 0;
		if(ref_poglinkv(g->links, link)->read > read) break;
		p = ref_poglinkv(g->links, link);
		link = p->next;
	}
	if(p){
		l->next = p->next;
		p->next = g->links->size - 1;
	} else {
		l->next = e->link;
		e->link = g->links->size - 1;
	}
	return 1;
}

static inline void print_local_dot_pomsa(POMSA *g, char *dot_name, uint32_t beg, uint32_t end, FILE *out){
	pognode_t *n1, *n2;
	pogedge_t *e;
	poglink_t *l;
	u4v *stack;
	BitVec *flags;
	char *colors[2][2] = {{"black", "yellow"}, {"blue", "yellow"}};
	uint32_t idx, link;
	if(beg > g->ref->size) return;
	if(end < beg) end = beg;
	fprintf(out, "digraph %s {\n", dot_name? dot_name : "");
	fprintf(out, "\trankdir=LR\n");
	stack = init_u4v(64);
	flags = init_bitvec(g->nodes->size);
	push_u4v(stack, beg + 1);
	one_bitvec(flags, beg + 1);
	while(stack->size){
		idx = stack->buffer[--stack->size];
		n1 = ref_pognodev(g->nodes, idx);
		fprintf(out, "\t%c%d_%d [label=\"%u:%d:%c:%d:%d\" color=\"%s\"]\n", "ACGT"[n1->base], idx, n1->pos, idx, n1->pos, "ACGT"[n1->base], n1->n_in, n1->n_ot, colors[n1->ref][n1->cns]);
		if(n1->pos > end) continue;
		e  = ref_pogedgev(g->edges, n1->edges);
		while(e->node != NODE_NULL){
			n2 = ref_pognodev(g->nodes, e->node);
			if(0){
				link = e->link;
				while(link != LINK_NULL){
					l = ref_poglinkv(g->links, link);
					link = l->next;
					fprintf(out, "\t%c%d_%d -> %c%d_%d [label=\"%03u:%u\" color=%s]\n", "ACGT"[n1->base], idx, n1->pos, "ACGT"[n2->base], e->node, n2->pos, l->read, e->cov, colors[e->ref][e->cns]);
				}
			} else {
				fprintf(out, "\t%c%d_%d -> %c%d_%d [label=\"%d\" color=%s]\n", "ACGT"[n1->base], idx, n1->pos, "ACGT"[n2->base], e->node, n2->pos, e->cov, colors[e->ref][e->cns]);
			}
			if(get_bitvec(flags, e->node) == 0){
				one_bitvec(flags, e->node);
				push_u4v(stack, e->node);
			}
			e = ref_pogedgev(g->edges, e->next);
		}
	}
	fprintf(out, "}\n");
	free_u4v(stack);
	free_bitvec(flags);
}

static inline void printf_local_dot_pomsa(POMSA *g, uint32_t beg, uint32_t end, char *filename){
	FILE *out;
	if((out = fopen(filename, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s for write in %s -- %s:%d --\n", filename, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	}
	print_local_dot_pomsa(g, NULL, beg, end, out);
	fclose(out);
}

static inline void beg_update_pomsa(POMSA *g, uint8_t *ref, uint32_t len){
	pognode_t *n;
	pogedge_t *e;
	poglink_t *l;
	uint32_t i;
	clear_u1v(g->ref);
	clear_pognodev(g->nodes);
	clear_pogedgev(g->edges);
	clear_poglinkv(g->links);
	append_array_u1v(g->ref, ref, len);
	n = next_ref_pognodev(g->nodes);
	memset(n, 0, sizeof(pognode_t));
	e = next_ref_pogedgev(g->edges);
	memset(e, 0, sizeof(pogedge_t));
	l = next_ref_poglinkv(g->links);
	memset(l, 0, sizeof(poglink_t));
	reset_hzmaux(g->aux);
	app_tseq_hzmaux(g->aux, ref, len);
	ready_hzmaux(g->aux);
	for(i=0;i<len;i++) add_node_pomsa(g, i, ref[i], 1);
	for(i=0;i+1<len;i++) add_edge_pomsa(g, i + 1, i + 2, 1);
	g->aln_score = 0;
	g->counter = 0;
}

// If seq already aligned, set hit.aln to positive value, hit.qb/qe/tb/te will be used directly
// Else, set hit.aln to zero or negative value, hit.tb/te will be used to restrict the new alignment
// single thread only
static inline int update_pomsa(POMSA *g, uint8_t *seq, uint32_t len, float min_sm, kswx_t *hit){
	pognode_t *n, *n2;
	pogedge_t *e, *ee;
	podp_t *dp, *dp1, *dp2;
	u4v *stack;
	pocell_t *cs1, *cs2, *c1, *c2;
	uint32_t qb, qe, tb, te, i, j, beg, end, base, mg_base;
	uint32_t idx, max_idx, local_idx, y, z, pos, cur, nxt;
	int max_score, max_local, ms, mi, mz, x;
	int8_t *qps, *qp;
	g->counter ++;
	if(hit->aln > 0){
		g->aux->hit = *hit;
	} else {
		if(hit->tb <= 0) hit->tb = 0;
		if(hit->te <= 0) hit->te = len;
		if(align_hzmaux(g->aux, 0, seq, NULL, len, hit->tb, hit->te, 0, min_sm) == 0) return 0;
		*hit = g->aux->hit;
	}
	qb = g->aux->hit.qb;
	qe = g->aux->hit.qe;
	tb = g->aux->hit.tb;
	te = g->aux->hit.te;
	stack = init_u4v(1024);
	// set node->aux to adjusted n_in (may less than n_in, for tb greater than 0)
	{
		for(i=0;i<g->nodes->size;i++){ g->nodes->buffer[i].cnt = 0; g->nodes->buffer[i].aux = 0; }
		push_u4v(stack, tb + 1);
		while(stack->size){
			idx = stack->buffer[--stack->size];
			n = ref_pognodev(g->nodes, idx);
			n->aux ++;
			if(n->pos > te) continue;
			if(n->aux > 1) continue;
			e = ref_pogedgev(g->edges, n->edges);
			while(e->node != NODE_NULL){
				push_u4v(stack, e->node);
				e = ref_pogedgev(g->edges, e->next);
			}
		}
	}
	clear_podpv(g->dps);
	dp = next_ref_podpv(g->dps);
	memset(dp, 0, sizeof(podp_t));
	clear_pocellv(g->cells);
	// init qprof
	{
		qps = malloc(len * 4);
		for(i=0;i<4;i++){
			qp = qps + i * len;
			for(j=0;j<len;j++){
				qp[j] = (i == seq[j])? g->aux->M : g->aux->X;
			}
		}
	}
	dp = next_ref_podpv(g->dps);
	dp->node = tb + 1; // ref node id
	dp->qx   = qb;
	dp->qy   = dp->qx + g->W;
	if(dp->qy > qe) dp->qy = qe;
	if(dp->qx >= dp->qy) return 0;
	dp->coff = g->cells->size;
	encap_and_inc_pocellv(g->cells, dp->qy - dp->qx);
	dp->woff = 0;
	dp->winc = 0;
	ref_pognodev(g->nodes, dp->node)->cnt = 1;
	dp->ready = 1;
	base = ref_pognodev(g->nodes, dp->node)->base;
	// setup first row
	{
		cs1 = g->cells->buffer + dp->coff;
		qp = qps + base * len + dp->qx;
		for(i=0;i<g->W;i++){
			cs1[i].h = qp[i];
			cs1[i].e = SCORE_MIN;
			cs1[i].d = 0;
			cs1[i].mg = 4; // No homopolymer merge possiblity
			cs1[i].hbt = DP_NULL;
			cs1[i].ebt = DP_NULL;
		}
	}
	dp->score = g->aux->M;
	push_u4v(stack, g->dps->size - 1);
	max_score = SCORE_MIN;
	max_idx = DP_NULL;
	max_local = SCORE_MIN;
	local_idx = DP_NULL;
	while(stack->size){
		idx = stack->buffer[--stack->size];
		encap_podpv(g->dps, NODE_MAX_EDGES);
		dp = ref_podpv(g->dps, idx);
		if(dp->qx >= qe){
			continue;
		}
		if(dp->qy >= qe && dp->qx < qe){
			c1 = ref_pocellv(g->cells, dp->coff + qe - 1 - dp->qx);
			if(c1->h > max_score){
				max_score = c1->h;
				max_idx   = idx;
			}
		}
		if(dp->qx + dp->woff >= qe){
			continue;
		}
		n  = ref_pognodev(g->nodes, dp->node);
		if(n->pos >= te){
			continue;
		}
		ee = ref_pogedgev(g->edges, n->edges);
		while(ee->node != NODE_NULL){
			e = ee;
			ee = ref_pogedgev(g->edges, e->next);
			dp2 = next_ref_podpv(g->dps);
			dp2->node = e->node;
			n2 = ref_pognodev(g->nodes, dp2->node);
			dp2->qx   = dp->qx + dp->woff;
			dp2->qy   = dp2->qx + g->W;
			if(dp2->qy < dp->qy) dp2->qy = dp->qy;
			if(dp->winc && dp2->qy - dp2->qx < dp->qy - dp->qx + dp->winc) dp2->qy = dp2->qx + dp->qy - dp->qx + dp->winc;
			if(dp2->qy > dp2->qx + g->MW) dp2->qy = dp2->qx + g->MW;
			//dp2->qy += dp->winc;
			if(dp2->qy > qe) dp2->qy = qe;
			dp2->coff = g->cells->size;
			if(dp2->qx >= dp2->qy) dp2->qy = dp2->qx;
			encap_and_inc_pocellv(g->cells, dp2->qy - dp2->qx);
			dp2->ready = 0;
			dp2->score = SCORE_MIN;
			base = n2->base;
			qp = qps + base * len + dp2->qx;
			// call h, e
			{
				cs1 = g->cells->buffer + dp->coff + dp->woff;
				cs2 = g->cells->buffer + dp2->coff;
				beg = 0;
				if(dp->qy < dp2->qx) end = 0; else end = dp->qy - dp2->qx;
				if(dp->woff == 0){
					cs2[0].h = ((dp2->qx == qb)? 0 : SCORE_MIN);
					cs2[0].e = cs1[0].e + g->aux->E;
					cs2[0].mg = cs1[0].mg;
					beg = 1;
				} else beg = 0;
				for(i=beg;i<end;i++){
					cs2[i].h = cs1[(int)i-1].h;
					cs2[i].e = cs1[i].e;
					cs2[i].mg = cs1[i].mg;
				}
				for(;i+dp2->qx<dp2->qy;i++){
					cs2[i].h = cs2[i].e = SCORE_MIN;
					cs2[i].mg = 4;
				}
			}
			// merge multiple in-edges
			if(n2->aux > 1){
				if(n2->cnt){
					dp1 = ref_podpv(g->dps, n2->bt);
					//merge dp2 into dp1
					while(1){
						uint64_t coff;
						beg = num_min(dp1->qx, dp2->qx);
						end = num_max(dp1->qy, dp2->qy);
						if(beg >= end) break;
						// alloc a larger alignment row
						if(beg < dp1->qx || end > dp1->qy){
							coff = g->cells->size;
							encap_and_inc_pocellv(g->cells, end - beg);
							// copy cells
							memcpy(g->cells->buffer + coff + dp1->qx - beg, g->cells->buffer + dp1->coff, sizeof(pocell_t) * (dp1->qy - dp1->qx));
							dp1->coff = coff;
							cs1 = g->cells->buffer + dp1->coff;
							// init new cells
							for(i=beg;i<dp1->qx;i++){
								c1 = cs1 + i - beg;
								c1->h = c1->e = SCORE_MIN;
								c1->hbt = c1->ebt = DP_NULL;
							}
							for(i=dp1->qy;i<end;i++){
								c1 = cs1 + i - beg;
								c1->h = c1->e = SCORE_MIN;
								c1->hbt = c1->ebt = DP_NULL;
							}
							dp1->qx = beg;
							dp1->qy = end;
						}
						cs1 = g->cells->buffer + dp1->coff;
						cs2 = g->cells->buffer + dp2->coff;
						for(i=dp2->qx;i<dp2->qy;i++){
							c1 = cs1 + i - dp1->qx;
							c2 = cs2 + i - dp2->qx;
							if(c2->h > c1->h){
								c1->h = c2->h;
								c1->hbt = idx;
							}
							if(c2->e > c1->e){
								c1->e   = c2->e;
								c1->mg  = c2->mg;
								c1->ebt = idx;
							}
						}
						break;
					}
					n2->cnt ++;
					if(n2->cnt < n2->aux) continue; // still need further merging
					// ready for fully calculation
					dp2 = dp1;
				} else {
					n2->cnt ++;
					n2->bt = dp2 - g->dps->buffer;
					end = dp2->qy - dp2->qx;
					cs2 = g->cells->buffer + dp2->coff;
					for(i=0;i<end;i++) cs2[i].hbt = cs2[i].ebt = idx;
					continue;
				}
			} else {
				n2->cnt ++;
				n2->bt = dp2 - g->dps->buffer;
				// set xbt
				end = dp2->qy - dp2->qx;
				cs2 = g->cells->buffer + dp2->coff;
				for(i=0;i<end;i++) cs2[i].hbt = cs2[i].ebt = idx;
			}
			// fully calculation
			ms = SCORE_MIN;
			mi = 0;
			mz = 0;
			base = ref_pognodev(g->nodes, dp2->node)->base;
			mg_base = g->has_merge? base : 4;
			qp = qps + base * len + dp2->qx;
			{
				register int rh, re, rm, rf, t, d;
				beg = 0;
				end = dp2->qy - dp2->qx;
				rf  = SCORE_MIN;
				cs2 = g->cells->buffer + dp2->coff;
				for(i=beg;i<end;i++){
					c2 = cs2 + i;
					rm = c2->h + qp[i];
					re = c2->e;
					if(rm >= re){ d = 0; rh = rm; }
					else        { d = 1; rh = re; }
					if(rh < rf){ d = 2; rh = rf; }
					if(rh > ms){ ms = rh; mi = i; mz = d; }
					c2->h = rh;

					t  = rm + g->aux->D + g->aux->E;
					if(c2->mg == base){
						// Homopolymer merge
					} else {
						re = re + g->aux->E;
					}
					if(re > t){ d |= 1 << 2; c2->mg = 4; }
					else      { re = t; c2->mg = mg_base;  }
					c2->e  = re;

					t = rm + g->aux->I + g->aux->E;
					rf = rf + g->aux->E;
					if(rf > t){ d |= 2 << 4; }
					else      { rf = t; }

					c2->d = d;
				}
				//{ fprintf(stdout, "q=%d[%d~%d=%d] t=%d z=%d score=%d->%d\n", dp2->qx + mi, dp2->qx, dp2->qy, dp2->qy - dp2->qx, dp2->node, mz, dp->score, ms); fflush(stdout); }
				/*
				if(dp2->node >= 32 && dp2->node <= 34 && dp2->qx == 0){
					for(i=30;i<36;i++){
						fprintf(stdout, "%d H=[%d,%u] E=[%d,%u] Z=%d\n", i, cs2[i].h, cs2[i].hbt, cs2[i].e, cs2[i].ebt, cs2[i].d);
					}
				}
				*/
			}
			{
				if(ms > max_local){
					max_local = ms;
					local_idx = dp2 - g->dps->buffer;
				}
			}
			{
				dp2->score = ms;
				if(mi >= (int)(3 * (dp2->qy - dp2->qx) / 5)) dp2->woff = 2;
				else if(mi <= (int)(2 * (dp2->qy - dp2->qx) / 5)) dp2->woff = 0;
				else dp2->woff = 1;
				dp2->winc = (mz || dp2->score <= dp->score)? 2 : 0;
				// when band is larger than g->W, reduce it
				//if(mz && dp2->qx + g->W < dp2->qy && dp2->qy < qe){
					//if(dp2->woff == 0) dp2->qy --;
					//else if(dp2->woff == 1) dp2->qx ++;
					//else { dp2->qx ++; dp2->qy --; }
				//}
			}
			dp2->ready = 1;
			push_u4v(stack, dp2 - g->dps->buffer);
		}
	}
	free(qps);
	free_u4v(stack);
	if(max_idx == DP_NULL){ return 0; }
	g->aln_score += max_score;
	g->aux->hit.score = max_score;
	// update graph
	nxt = NODE_NULL;
	x = qe - 1;
	y = max_idx;
	pos = g->nodes->buffer[g->dps->buffer[y].node].pos;
	g->aux->hit.qe = x + 1;
	g->aux->hit.te = pos + 1;
	g->aux->hit.aln = 0;
	g->aux->hit.mat = 0;
	g->aux->hit.mis = 0;
	g->aux->hit.ins = 0;
	g->aux->hit.del = 0;
	z = 0;
	//fprintf(stdout, " -- BEG in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stdout);
	while(y != DP_NULL && x >= (int)qb){
		dp = ref_podpv(g->dps, y);
		c1 = ref_pocellv(g->cells, dp->coff + x - dp->qx);
		z = (c1->d >> (z << 1)) & 0x03;
		encap_pognodev(g->nodes, 1);
		n = ref_pognodev(g->nodes, dp->node);
		pos = n->pos;
		g->aux->hit.aln ++;
		//fprintf(stdout, "X=%d Y=%d Z=%d H=%d E=%d\n", x, y, z, c1->h, c1->e);
		if(z == 0){
			if(seq[x] == n->base){
				cur = dp->node;
			} else {
				// try find the aligned node with idential base
				uint32_t dir, f, nid;
				pognode_t *nn;
				for(f=dir=0;dir<2;dir++){
					nid = n->aligned[dir];
					while(nid != NODE_NULL){
						nn = ref_pognodev(g->nodes, nid);
						if(nn->base == seq[x]){ f = nid; break; }
						nid = nn->aligned[dir];
					}
					if(f) break;
				}
				if(f){
					cur = f;
				} else {
					cur = add_node_pomsa(g, pos, seq[x], 0);
					nn = ref_pognodev(g->nodes, cur);
					nn->aligned[1] = dp->node;
					nn->aligned[0] = n->aligned[0];
					n->aligned[0]  = cur;
				}
			}
			if(nxt != DP_NULL) add_edge_read_pomsa(g, ref_pogedgev(g->edges, add_edge_pomsa(g, cur, nxt, 0)), g->counter);
			nxt = cur;
			g->aux->hit.mat ++;
			x --;
			y = c1->hbt;
		} else if(z == 1){
			y = c1->ebt;
			g->aux->hit.del ++;
		} else {
			cur = add_node_pomsa(g, pos, seq[x], 0);
			if(nxt != DP_NULL) add_edge_read_pomsa(g, ref_pogedgev(g->edges, add_edge_pomsa(g, cur, nxt, 0)), g->counter);
			nxt = cur;
			x --;
			g->aux->hit.ins ++;
		}
	}
	g->beginning = nxt;
	g->aux->hit.qb = g->aux->hit.qe - (g->aux->hit.mat + g->aux->hit.mis + g->aux->hit.ins);
	g->aux->hit.tb = pos;
	//fprintf(stdout, " -- END in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stdout);
	// node coverage
	{
		for(i=g->aux->hit.tb;(int)i<g->aux->hit.te;i++) g->nodes->buffer[i + 1].cov ++;
	}
	return 1;
}

// reset the coverage of all edges to zero
static inline void beg_match_pomsa(POMSA *g){
	pogedge_t *e;
	uint32_t i;
	g->links->size = 1;
	for(i=0;i<g->edges->size;i++){
		e = ref_pogedgev(g->edges, i);
		e->cov = 0;
		e->link = LINK_NULL;
	}
}

// single thread only
static inline int match_pomsa(POMSA *g, uint32_t id, uint8_t *seq, uint32_t len, uint32_t beginning, u4v *edges){
	pognode_t *n;
	pogedge_t *e;
	pomat_t *m, *v;
	uint64_t *u, U;
	uint32_t idx, end;
	int exists;
	clear_pomatv(g->mats);
	clear_u64hash(g->hash);
	clear_u4v(g->heap);
	next_ref_pomatv(g->mats); // NULL
	m = next_ref_pomatv(g->mats);
	m->node = beginning;
	m->pos  = 0;
	m->btn  = NODE_NULL;
	m->bte  = EDGE_NULL;
	put_u64hash(g->hash, (((uint64_t)m->node) << 32) | m->pos);
	end = 0;
	push_u4v(g->heap, 1);
	if(0){
		uint32_t i;
		fprintf(stderr, "BEG: node=%u\n", beginning);
		for(i=0;i<50;i++){
			fprintf(stderr, "%u:%c,", i, bit_base_table[seq[i]]);
		}
		fprintf(stderr, "\n");
	}
	while(g->heap->size){
		idx = g->heap->buffer[--g->heap->size]; // stack not heap
		//idx = array_heap_pop(g->heap->buffer, g->heap->size, g->heap->cap, uint32_t, num_cmp(g->mats->buffer[a].pos, g->mats->buffer[b].pos));
		encap_pomatv(g->mats, NODE_MAX_EDGES);
		m = ref_pomatv(g->mats, idx);
		//fprintf(stderr, "POP: %u\t%u\t%u\n", m->node, m->pos, (uint32_t)g->heap->size);
		if(m->pos + 1 == len){ end = idx; break; }
		n = ref_pognodev(g->nodes, m->node);
		e = ref_pogedgev(g->edges, n->edges);
		while(e->node != NODE_NULL){
			//fprintf(stderr, "MEET: %u\t%c\t%c\n", e->node, bit_base_table[g->nodes->buffer[e->node].base], bit_base_table[seq[m->pos + 1]]);
			if(g->nodes->buffer[e->node].base == seq[m->pos + 1]){
				U = (((uint64_t)e->node) << 32) | (m->pos + 1);
				u = prepare_u64hash(g->hash, U, &exists);
				if(!exists){
					*u = U;
					v = next_ref_pomatv(g->mats);
					v->node = e->node;
					v->pos  = m->pos + 1;
					v->btn  = idx;
					v->bte  = e - g->edges->buffer;
					//fprintf(stderr, "PUSH: %u\t%u\t%c\n", v->node, v->pos, bit_base_table[seq[m->pos + 1]]);
					//array_heap_push(g->heap->buffer, g->heap->size, g->heap->cap, uint32_t, g->mats->size - 1, num_cmp(g->mats->buffer[a].pos, g->mats->buffer[b].pos));
					push_u4v(g->heap, g->mats->size - 1);
				} else {
					//fprintf(stderr, "SKIP: %u\t%u\t%c\n", e->node, m->pos + 1, bit_base_table[seq[m->pos + 1]]);
				}
			}
			e = ref_pogedgev(g->edges, e->next);
		}
	}
	if(end == 0) return 0;
	m = ref_pomatv(g->mats, end);
	if(edges) clear_u4v(edges);
	while(m->btn != NODE_NULL){
		e = ref_pogedgev(g->edges, m->bte);
		e->cov ++;
		add_edge_read_pomsa(g, e, id);
		if(edges) push_u4v(edges, m->bte);
		m = ref_pomatv(g->mats, m->btn);
	}
	if(edges) reverse_u4v(edges);
	return 1;
}

static inline uint32_t intersection_edge_cov_pomsa(POMSA *g, pogedge_t *p, pogedge_t *e){
	uint32_t link1, link2, cov;
	if(p == NULL || e == NULL) return 0;
	link1 = p->link;
	link2 = e->link;
	cov = 0;
	while(link1 != LINK_NULL && link2 != LINK_NULL){
		if(g->links->buffer[link1].read == g->links->buffer[link2].read){
			cov ++;
			link1 = g->links->buffer[link1].next;
			link2 = g->links->buffer[link2].next;
		} else if(g->links->buffer[link1].read < g->links->buffer[link2].read){
			link1 = g->links->buffer[link1].next;
		} else {
			link2 = g->links->buffer[link2].next;
		}
	}
	return cov;
}

static inline void call_consensus_pomsa(POMSA *g){
	u4v *stack;
	pognode_t *n, *n2;
	pogedge_t *e;
	uint32_t i, idx, max_idx, cnt, pos, cov;
	int max_score, score, ncov;
	stack = init_u4v(64);
	for(i=0;i<g->nodes->size;i++){
		n = ref_pognodev(g->nodes, i);
		n->cnt = 0;
		n->score = SCORE_MIN;
		n->bt  = NODE_NULL;
		n->aux = EDGE_NULL; // last edge
		if(n->n_in == 0){
			n->score = 0;
			push_u4v(stack, i);
		}
	}
	max_score = SCORE_MIN;
	max_idx   = NODE_NULL;
	while(stack->size){
		idx = stack->buffer[--stack->size];
		n = ref_pognodev(g->nodes, idx);
		if((n->score > max_score) || (n->score == max_score && n->pos > pos)){
			max_score = n->score;
			max_idx = idx;
			pos = n->pos;
		}
		ncov = ref_pognodev(g->nodes, n->pos + 1)->cov / 3; // e->cov less than 1/3 will be punished
		e = ref_pogedgev(g->edges, n->edges);
		while(e->node != NODE_NULL){
			n2 = ref_pognodev(g->nodes, e->node);
			n2->cnt ++;
			cov = intersection_edge_cov_pomsa(g, g->edges->buffer + n->aux, e);
			if(1){
				score = n->score - ncov + cov;
			} else {
				if(n2->ref){
					score = n->score - ncov + (cov - 1);
				} else {
					score = n->score - ncov + cov;
				}
			}
			if(score > n2->score){
				n2->score = score;
				n2->bt  = idx;
				n2->aux = e - g->edges->buffer;
			}
			if(n2->cnt == n2->n_in){
				push_u4v(stack, e->node);
			}
			e = ref_pogedgev(g->edges, e->next);
		}
	}
	free_u4v(stack);
	clear_u1v(g->cns);
	g->cns_score = max_score;
	if(max_idx == NODE_NULL) return;
	clear_and_encap_u4v(g->map, g->ref->size + 1);
	g->map->size = g->ref->size + 1;
	for(i=0;i<g->map->size;i++) g->map->buffer[i] = 0;
	cnt = 0;
	idx = max_idx;
	while(idx != NODE_NULL){
		n = ref_pognodev(g->nodes, idx);
		n->cns = 1;
		g->edges->buffer[n->aux].cns = 1;
		push_u1v(g->cns, n->base);
		cnt ++;
		if(n->ref){
			g->map->buffer[n->pos] = cnt;
			cnt = 0;
		}
		idx = n->bt;
	}
	g->map->buffer[0] += cnt;
	reverse_u1v(g->cns);
	cnt = 0;
	for(i=0;i<g->map->size;i++){
		g->map->buffer[i] += cnt;
		cnt = g->map->buffer[i];
	}
}

static inline void print_cnsseq_pomsa(POMSA *g, FILE *out){
	uint32_t i;
	for(i=0;i<g->cns->size;i++){
		fputc(bit_base_table[g->cns->buffer[i]], out);
		if((i % 100) == 99) fputc('\n', out);
	}
	if((i % 100)) fputc('\n', out);
}

#endif
