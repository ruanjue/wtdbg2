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

#ifndef __GENERAL_GRAPH_RJ_H
#define __GENERAL_GRAPH_RJ_H

#include "list.h"
#include "hashset.h"

#define GEG_MAX_NODE	0xFFFFFFFFFFLLU
#define GEG_MAX_EDGE_CNT	0x3FFFF
#define GEG_MAX_EDGE_COV	0x7FFFF
#define GEG_MAX_EDGE_OFF	0x7FFFFF
#define GEG_MIN_EDGE_OFF	-0x7FFFFF

typedef struct {
	u8i node1:40, dir1:1, dir2:1, closed:2, cov:19, visit:1;
	u8i node2:40; b8i off:24;
} ge_edge_t;
define_list(geedgev, ge_edge_t);

static inline uint64_t _ge_edge_hashcode(ge_edge_t e){
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
#define GEEDGEHASH(idx) ((geedgev *)set->userdata)->buffer[idx]
#define ge_edge_hashcode(E) _ge_edge_hashcode(GEEDGEHASH(E))
#define ge_edge_hashequals(E1, E2) (GEEDGEHASH(E1).node1 == GEEDGEHASH(E2).node1 && GEEDGEHASH(E1).node2 == GEEDGEHASH(E2).node2 \
			&& GEEDGEHASH(E1).dir1 == GEEDGEHASH(E2).dir1 && GEEDGEHASH(E1).dir2 == GEEDGEHASH(E2).dir2)
define_hashset(geedgehash, u8i, ge_edge_hashcode, ge_edge_hashequals);

typedef struct {
	u8i idx:63, flg:1;
	u8i next;
} ge_edge_ref_t;
define_list(geedgerefv, ge_edge_ref_t);

typedef struct { uint64_t idx:46, cnt:18; } ge_ptr_ref_t;
static const ge_ptr_ref_t GE_PTR_REF_NULL = (ge_ptr_ref_t){0, 0};
define_list(geptrrefv, ge_ptr_ref_t);
typedef struct { uint64_t idx:46, cnt:18; } ge_vec_ref_t;
static const ge_vec_ref_t GE_VEC_REF_NULL = (ge_vec_ref_t){0, 0};
define_list(gevecrefv, ge_vec_ref_t);

typedef struct {
	u8i closed:1, bt_visit:40, bt_dir:1, bt_idx:18, status:4;
	u8i unvisit:18, aux:46;
	ge_ptr_ref_t edges[2];
} ge_node_t;
define_list(genodev, ge_node_t);

#define GEG_TRACE_MSG_ZERO	0
#define GEG_TRACE_MSG_ONE	1
#define GEG_TRACE_MSG_MORE	2
#define GEG_TRACE_MSG_VISITED	3
#define GEG_TRACE_MSG_UNDEF	4

typedef void (*geg_clr_node_callback)(void *userdata);
typedef void (*geg_add_node_callback)(void *userdata, u8i nidx);
typedef void (*geg_del_node_callback)(void *userdata, u8i nidx);
typedef void (*geg_clr_edge_callback)(void *userdata);
typedef void (*geg_add_edge_callback)(void *userdata, u8i eidx);
typedef void (*geg_del_edge_callback)(void *userdata, u8i eidx);

#define define_simple_geg_callback(tag, node_tag, node_tag_t, edge_tag, edge_tag_t)	\
void tag##nodeclr(void *aux){ clear_##node_tag((node_tag*)aux); }	\
void tag##nodeadd(void *aux, u8i idx){	\
	node_tag *nodes = (node_tag*)aux;	\
	if(idx < nodes->size){	\
		memset(ref_##node_tag(nodes, idx), 0, sizeof(node_tag_t));	\
	} else {	\
		memset(next_ref_##node_tag(nodes), 0, sizeof(node_tag_t));	\
	}	\
}	\
void tag##nodedel(void *aux, u8i idx){ UNUSED(aux); UNUSED(idx); }	\
void tag##edgeclr(void *aux){ clear_##edge_tag((edge_tag*)aux); }	\
void tag##edgeadd(void *aux, u8i idx){	\
	edge_tag *edges = (edge_tag*)aux;	\
	if(idx < edges->size){	\
		memset(ref_##edge_tag(edges, idx), 0, sizeof(edge_tag_t));	\
	} else {	\
		memset(next_ref_##edge_tag(edges), 0, sizeof(edge_tag_t));	\
	}	\
}	\
void tag##edgedel(void *aux, u8i idx){ UNUSED(aux); UNUSED(idx); }	\
static inline void tag##_set_callbacks_gegraph(GEGraph *g, node_tag *nodeaux, edge_tag *edgeaux){	\
	set_callbacks_gegraph(g, (void*)nodeaux, (void*)edgeaux, tag##nodeclr, tag##nodeadd, tag##nodedel, tag##edgeclr, tag##edgeadd, tag##edgedel);	\
}

typedef struct {
	genodev    *nodes;
	geedgev    *edges;
	geedgehash *ehash;
	geedgerefv *erefs;
	void       *nodeaux;
	void       *edgeaux;
	geg_clr_node_callback nodeclr;
	geg_add_node_callback nodeadd;
	geg_del_node_callback nodedel;
	geg_clr_edge_callback edgeclr;
	geg_add_edge_callback edgeadd;
	geg_del_edge_callback edgedel;
} GEGraph;

static inline GEGraph* init_gegraph(){
	GEGraph *g;
	g = malloc(sizeof(GEGraph));
	g->nodes = init_genodev(32);
	g->edges = init_geedgev(32);
	g->ehash = init_geedgehash(1023);
	set_userdata_geedgehash(g->ehash, g->edges);
	g->erefs = init_geedgerefv(32);
	g->nodeaux  = NULL;
	g->edgeaux  = NULL;
	g->nodeclr  = NULL;
	g->nodeadd  = NULL;
	g->nodedel  = NULL;
	g->edgeclr  = NULL;
	g->edgeadd  = NULL;
	g->edgedel  = NULL;
	return g;
}

static inline void set_callbacks_gegraph(GEGraph *g, void *nodeaux, void *edgeaux, geg_clr_node_callback nodeclr, geg_add_node_callback nodeadd, geg_del_node_callback nodedel, geg_clr_edge_callback edgeclr, geg_add_edge_callback edgeadd, geg_del_edge_callback edgedel){
	g->nodeaux  = nodeaux;
	g->edgeaux  = edgeaux;
	g->nodeclr  = nodeclr;
	g->nodeadd  = nodeadd;
	g->nodedel  = nodedel;
	g->edgeclr  = edgeclr;
	g->edgeadd  = edgeadd;
	g->edgedel  = edgedel;
}

static inline void free_gegraph(GEGraph *g){
	free_genodev(g->nodes);
	free_geedgev(g->edges);
	free_geedgehash(g->ehash);
	free_geedgerefv(g->erefs);
	free(g);
}

static inline void reset_gegraph(GEGraph *g){
	clear_genodev(g->nodes);
	if(g->nodeclr) g->nodeclr(g->nodeaux);
	clear_geedgev(g->edges);
	memset(next_ref_geedgev(g->edges), 0, sizeof(ge_edge_t));
	if(g->edgeclr){
		g->edgeclr(g->edgeaux);
		g->edgeadd(g->edgeaux, 0);
	}
	clear_geedgehash(g->ehash);
	clear_geedgerefv(g->erefs);
	memset(next_ref_geedgerefv(g->erefs), 0, sizeof(ge_edge_ref_t));
}

static inline ge_node_t* add_node_gegraph(GEGraph *g){
	ge_node_t *n;
	n = next_ref_genodev(g->nodes);
	memset(n, 0, sizeof(ge_node_t));
	if(g->nodeadd) g->nodeadd(g->nodeaux, offset_genodev(g->nodes, n));
	return n;
}

static inline ge_edge_t* prepare_edge_gegraph(GEGraph *g, u8i node1, int dir1, u8i node2, int dir2, int *exists){
	ge_node_t *n;
	ge_edge_t *e;
	ge_edge_ref_t *f;
	u8i *u;
	e = ref_geedgev(g->edges, 0);
	if(node1 <= node2){
		e->node1  = node1;
		e->dir1   = dir1;
		e->node2  = node2;
		e->dir2   = dir2;
	} else {
		e->node1  = node2;
		e->dir1   = !dir2;
		e->node2  = node1;
		e->dir2   = !dir1;
	}
	e->cov    = 0;
	e->off    = 0;
	e->visit = 0;
	e->closed = 0;
	u = prepare_geedgehash(g->ehash, 0, exists);
	if(*exists){
		return g->edges->buffer + *u;
	} else {
		*u = g->edges->size;
		e = next_ref_geedgev(g->edges);
		*e = g->edges->buffer[0];
		n = g->nodes->buffer + e->node1;
		f = next_ref_geedgerefv(g->erefs);
		f->idx = *u;
		f->flg = 0;
		f->next = n->edges[e->dir1].idx;
		n->edges[e->dir1].idx = g->erefs->size - 1;
		n->edges[e->dir1].cnt ++;
		n = g->nodes->buffer + e->node2;
		f = next_ref_geedgerefv(g->erefs);
		f->idx = *u;
		f->flg = 1;
		f->next = n->edges[!e->dir2].idx;
		n->edges[!e->dir2].idx = g->erefs->size - 1;
		n->edges[!e->dir2].cnt ++;
		if(g->edgeadd) g->edgeadd(g->edgeaux, offset_geedgev(g->edges, e));
		return e;
	}
}

static inline void cut_edge_core_gegraph(GEGraph *g, ge_edge_t *e, int closed_val){
	if(e->closed) return;
	e->closed = closed_val;
	ref_genodev(g->nodes, e->node1)->edges[e->dir1].cnt --;
	ref_genodev(g->nodes, e->node2)->edges[!e->dir2].cnt --;
	if(g->edgedel) g->edgedel(g->edgeaux, offset_geedgev(g->edges, e));
}

#define cut_edge_gegraph(g, e) cut_edge_core_gegraph(g, e, 1)

static inline void revive_edge_gegraph(GEGraph *g, ge_edge_t *e){
	if(e->closed == 0) return;
	e->closed = 0;
	ref_genodev(g->nodes, e->node1)->edges[e->dir1].cnt ++;
	ref_genodev(g->nodes, e->node2)->edges[!e->dir2].cnt ++;
	if(g->edgedel) g->edgeadd(g->edgeaux, offset_geedgev(g->edges, e));
}

static inline ge_edge_ref_t* single_edge_gegraph(GEGraph *g, ge_node_t *n, int dir, int *info){
	ge_edge_ref_t *f, *ret;
	uint64_t idx;
	ret = NULL;
	if(info){
		*info = GEG_TRACE_MSG_ZERO;
		if(n->edges[dir].cnt == 0) return NULL;
		idx = n->edges[dir].idx;
		while(idx){
			f = ref_geedgerefv(g->erefs, idx);
			idx = f->next;
			if(g->edges->buffer[f->idx].closed) continue;
			if(ret){ *info = GEG_TRACE_MSG_MORE; return NULL; }
			else { *info = GEG_TRACE_MSG_ONE; ret = f; }
		}
	} else {
		if(n->edges[dir].cnt == 0) return NULL;
		idx = n->edges[dir].idx;
		while(idx){
			f = ref_geedgerefv(g->erefs, idx);
			idx = f->next;
			if(g->edges->buffer[f->idx].closed) continue;
			if(ret){ return NULL; }
			else { ret = f; }
		}
	}
	return ret;
}

#define count_edges_gegraph(g, n, dir) (n)->edges[dir].cnt

// dir = 2 means either strand
static inline ge_edge_ref_t* edge_node2node_gegraph(GEGraph *g, u8i node1, int dir1, u8i node2, int dir2){
	ge_node_t *n;
	ge_edge_ref_t *f;
	ge_edge_t *e;
	uint64_t idx;
	int dire;
	n = ref_genodev(g->nodes, node1);
	if(dir1 > 1){
		dir1 = 0; dire = 2;
	} else {
		dire = dir1 + 1;
	}
	while(dir1 < dire){
		idx = n->edges[dir1].idx;
		while(idx){
			f = ref_geedgerefv(g->erefs, idx);
			idx = f->next;
			e = ref_geedgev(g->edges, f->idx);
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

static inline void del_node_edges_gegraph(GEGraph *g, ge_node_t *n, int closed_val){
	ge_edge_ref_t *f;
	ge_edge_t *e;
	uint64_t idx;
	uint32_t k;
	for(k=0;k<2;k++){
		idx = n->edges[k].idx;
		while(idx){
			f = ref_geedgerefv(g->erefs, idx);
			idx = f->next;
			e = g->edges->buffer + f->idx;
			cut_edge_core_gegraph(g, e, closed_val);
		}
	}
}

static inline void del_node_gegraph(GEGraph *g, ge_node_t *n){
	del_node_edges_gegraph(g, n, 1);
	n->closed = 1;
	if(g->nodedel) g->nodeadd(g->nodeaux, offset_genodev(g->nodes, n));
}

#define geg_beg_iter_edges(g, n, dir, f, e)	\
{	\
	u8i _geg_iter_idx;	\
	_geg_iter_idx = (n)->edges[dir].idx;	\
	while(_geg_iter_idx){	\
		(f) = ref_geedgerefv((g)->erefs, _geg_iter_idx);	\
		_geg_iter_idx = (f)->next;	\
		(e) = (g)->edges->buffer + (f)->idx

#define geg_end_iter_edges()	\
	}	\
}

static inline void print_dot_gegraph(GEGraph *g, FILE *out){
	static const char *colors[2][2] = {{"blue", "green"}, {"red", "gray"}};
	ge_node_t *n;
	ge_edge_t *e;
	u8i i;
	fprintf(out, "digraph {\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_genodev(g->nodes, i);
		if(n->closed) continue;
		fprintf(out, " N%llu\n", i);
	}
	for(i=1;i<g->edges->size;i++){
		e = ref_geedgev(g->edges, i);
		if(e->closed) continue;
		fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", (u8i)e->node1, (u8i)e->node2, "+-"[e->dir1], "+-"[e->dir2], e->cov, e->off, colors[e->dir1][e->dir2]);
		fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", (u8i)e->node2, (u8i)e->node1, "-+"[e->dir2], "-+"[e->dir1], e->cov, e->off, colors[!e->dir2][!e->dir1]);
	}
	fprintf(out, "}\n");
	fflush(out);
}

static inline void fprint_dot_gegraph(GEGraph *g, char *prefix, char *suffix){
	FILE *out;
	out = open_file_for_write(prefix, suffix, 1);
	print_dot_gegraph(g, out);
	fclose(out);
}

#endif
