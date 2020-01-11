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

#ifndef __WTDBG_GRAPH_RJ_H
#define __WTDBG_GRAPH_RJ_H

#include "wtdbg.h"
#include "filewriter.h"
#include "pgzf.h"

static inline u8i print_local_dot_graph(Graph *g, char *prefix, char *suffix);

static inline void print_node_edges_cov_graph(Graph *g, FILE *out){
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

static inline void del_node_edges_graph(Graph *g, node_t *n){
	edge_ref_t *f;
	edge_t *e;
	u8i idx;
	u4i k;
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

static inline void del_frg_lnks_graph(Graph *g, frg_t *n){
	edge_ref_t *f;
	lnk_t *e;
	u8i idx;
	u4i k;
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

static inline void del_node_graph(Graph *g, node_t *n){
	del_node_edges_graph(g, n);
	n->closed = 1;
}

static inline u8i mask_nodes_by_edge_cov_graph(Graph *g, u4i min_node_cov, float min_edge_cov_ratio, FILE *out){
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
define_list(subnodev, subnode_t);
#define subnode_hashcode(E) u64hashcode((E).node)
#define subnode_hashequals(E1, E2) (E1).node == (E2).node
define_hashset(subnodehash, subnode_t, subnode_hashcode, subnode_hashequals);

typedef struct {
	subnode_t *node;
	u4i cov:28, visit:1, fwd:1, dir:1, closed:1;
	int off;
	u4i next;
} subedge_t;
define_list(subedgev, subedge_t);

static inline int evaluate_node_connectivity_graph(Graph *g, u8i nid, u4v *rds, subnodehash *nodes, subedgev *edges, ptrrefv *stack){
	node_t *nd;
	read_t *rd;
	reg_t  *rg;
	subnode_t N, *n, *n1, *n2;
	subedge_t *e;
	ptr_ref_t *p;
	u8i idx, edx, aim;
	u4i i, k, k1, k2, cnt;
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

static inline void print_subgraph_dot(Graph *g, u8i id, subnodehash *nodes, subedgev *edges, FILE *out){
	subnode_t *n;
	subedge_t *e;
	u8i idx;
	int k;
	UNUSED(g);
	fprintf(out, "digraph N%llu {\n", id);
	fprintf(out, " N%llu [style=filled fillcolor=yellow]\n", id);
	reset_iter_subnodehash(nodes);
	while((n = ref_iter_subnodehash(nodes))){
		if(n->closed) continue;
		fprintf(out, "N%llu [label=\"N%llu%s(%u)\"]\n", (u8i)n->node, (u8i)n->node, n->closed? "*" : "", n->cov);
		for(k=0;k<2;k++){
			idx = n->edges[k].idx;
			while(idx){
				e = ref_subedgev(edges, idx);
				idx = e->next;
				fprintf(out, " N%llu -> N%llu [label=\"%c%c:%d:%d\"]\n", (u8i)n->node, (u8i)e->node->node, "+-"[k], "+-"[e->dir], e->cov, e->off);
			}
		}
	}
	fprintf(out, "}\n");
}

static inline void fprint_subgraph_dot(Graph *g, u8i id, subnodehash *nodes, subedgev *edges, char *filename){
	FILE *out;
	out = open_file_for_write(filename, NULL, 1);
	print_subgraph_dot(g, id, nodes, edges, out);
	fclose(out);
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

static inline u8i mask_nodes_by_connectivity_graph(Graph *g, int ncpu, FILE *out){
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

static inline u8i mask_read_weak_regs_graph(Graph *g, int ncpu){
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

static inline u8i mask_possible_tip_nodes_graph(Graph *g){
	node_t *n;
	reg_t *r;
	u8i ret, i;
	u4i j, cnt;
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

static inline void print_node_edges_graph(Graph *g, u8i nid, int dir, FILE *out){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	u8i idx;
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

static inline edge_ref_t* first_living_edge_graph(Graph *g, node_t *n, int dir, int *info){
	edge_ref_t *f, *ret;
	u8i idx;
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

static inline edge_ref_t* first_living_lnk_graph(Graph *g, frg_t *n, int dir, int *info){
	edge_ref_t *f, *ret;
	u8i idx;
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
static inline edge_ref_t* edge_node2node_graph(Graph *g, u8i node1, int dir1, u8i node2, int dir2){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	u8i idx;
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

static inline u8i linear_trace_graph(Graph *g, tracev *path, u8i max_step, int *msg){
	trace_t *t;
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	u8i step;
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

static inline u8i linear_path_graph(Graph *g, pathv *path, int max_len, int *msg){
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

static inline int wander_path_graph(Graph *g, u4i frg_idx, int frg_dir, pathv *heap, u8i vst, int max_len){
	path_t T, S;
	frg_t *n, *w;
	lnk_t *e;
	edge_ref_t *f;
	u8i idx;
	int len;
	clear_pathv(heap);
	T.off = 0;
	T.frg = frg_idx;
	T.dir = frg_dir;
	ref_frgv(g->frgs, frg_idx)->bt_visit = vst;
	push_pathv(heap, T);
	while(heap->size){
		S = heap->buffer[0];
		array_heap_remove(heap->buffer, heap->size, heap->cap, path_t, 0, num_cmp(a.off, b.off));
		n = ref_frgv(g->frgs, S.frg);
		idx = n->lnks[S.dir].idx;
		while(idx){
			f = ref_edgerefv(g->lrefs, idx);
			idx = f->next;
			e = ref_lnkv(g->lnks, f->idx);
			if(e->closed) continue;
			if(e->weak) continue;
			w = ref_frgv(g->frgs, f->flg? e->frg1: e->frg2);
			if(w->bt_visit == vst) continue;
			w->bt_visit = vst;
			len = S.off + e->off + w->len;
			if(len >= max_len){
				return len;
			}
			T.frg = f->flg? e->frg1: e->frg2;
			T.dir = f->flg? !e->dir1 : e->dir2;
			T.off = len;
			array_heap_push(heap->buffer, heap->size, heap->cap, path_t, T, num_cmp(a.off, b.off));
		}
	}
	return 0;
}

static inline int cal_offset_traces_graph(Graph *g, tracev *path, u8i beg, u8i end, int offset){
	trace_t *t;
	node_t *n;
	reg_t *r;
	edge_t *e;
	u8i i;
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

static inline int cal_offset_paths_graph(Graph *g, pathv *path, u8i beg, u8i end){
	path_t *t;
	frg_t *n;
	lnk_t *e;
	u8i i;
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

static inline u8i true_linear_unique_trace_graph(Graph *g, tracev *path, u8i max_step, u8i visit, int *msg){
	trace_t *t;
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	u8i step;
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

static inline u8i true_linear_unique_path_graph(Graph *g, pathv *path, u8i max_step, u8i visit, int *msg){
	path_t *t;
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	u8i step;
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

static inline u8i count_linear_trace_graph(Graph *g, trace_t *t, u8i max_step, int *msg){
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	u8i step;
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

static inline int count_linear_path_graph(Graph *g, path_t *t, int max_len, int *msg){
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

static inline u4i del_isolated_nodes_graph(Graph *g, FILE *log){
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

static inline u8i cut_binary_edges_graph(Graph *g){
	UUhash *hash;
	UUhash_t *u;
	node_t *n;
	edge_ref_t *f;
	edge_t *e, *p;
	u8i idx, nid, ret;
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

static inline u8i cut_binary_lnks_graph(Graph *g, FILE *info){
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

static inline u8i cut_low_cov_lnks_graph(Graph *g, int low_cov){
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

static inline u4i rescue_low_cov_transitive_edges_graph(Graph *g, u8i nid, u8v *edges, UUhash *hash){
	node_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	edge_t *e, *e1, *e2;
	reg_t *r;
	UUhash_t *u;
	u8i idx, nid2, nid3;
	u4i i, k, k2, k3, k4, ret;
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

static inline  u8i rescue_low_cov_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	u8i i, nid, ret;
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

static inline u8i rescue_high_cov_edges_graph(Graph *g, u4i max_step, u4i cov_cutoff){
	tracev *path;
	trace_t *t;
	node_t *v, *w;
	edge_t *e;
	edge_ref_t *f;
	u8i node, vst, idx, fidx[3], ret;
	u4i k, d, i, cov[2];
	for(node=0;node<g->nodes->size;node++){
		v = ref_nodev(g->nodes, node);
		v->bt_visit = 0;
	}
	ret = 0;
	vst = 0;
	path = init_tracev(8);
	for(node=0;node<g->nodes->size;node++){
		v = ref_nodev(g->nodes, node);
		if(v->closed) continue;
		for(k=0;k<2;k++){
			if(v->edges[k].cnt != 1) continue;
			idx = v->edges[k].idx;
			cov[0] = cov[1] = 0;
			fidx[0] = fidx[1] = fidx[2] = 0;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = ref_edgev(g->edges, f->idx);
				if(e->closed){
					if(e->cov > cov[1]){
						cov[1] = e->cov;
						fidx[1] = offset_edgerefv(g->erefs, f);
					}
				} else {
					cov[0] = e->cov;
					fidx[0] = offset_edgerefv(g->erefs, f);
					if(cov[0] >= cov_cutoff) break;
				}
			}
			if(cov[0] >= cov_cutoff){
				continue;
			}
			if(cov[0] >= cov[1]){
				continue;
			}
			clear_tracev(path);
			t = next_ref_tracev(path);
			t->node = node;
			t->dir  = k;
			linear_trace_graph(g, path, max_step, NULL);
			if(path->size <= 2){
				continue;
			}
			vst ++;
			for(i=1;i<path->size;i++){
				t = ref_tracev(path, i);
				w = ref_nodev(g->nodes, t->node);
				w->bt_visit = vst;
			}
			cov[1] = 0;
			fidx[1] = 0;
			idx = v->edges[k].idx;
			while(idx){
				f = ref_edgerefv(g->erefs, idx);
				idx = f->next;
				e = ref_edgev(g->edges, f->idx);
				if(e->closed == 0) continue;
				w = ref_nodev(g->nodes, f->flg? e->node1 : e->node2);
				if(w->bt_visit == vst){
					if(e->cov > cov[0]){
						if(e->cov > cov[1]){
							cov[1] = e->cov;
							fidx[1] = offset_edgerefv(g->erefs, f);
						}
					}
				}
			}
			if(fidx[1] == 0){
				continue;
			}
			f = ref_edgerefv(g->erefs, fidx[1]);
			e = ref_edgev(g->edges, f->idx);
			if(f->flg){
				w = ref_nodev(g->nodes, e->node1);
				d = e->dir1;
			} else {
				w = ref_nodev(g->nodes, e->node2);
				d = !e->dir2;
			}
			f = first_living_edge_graph(g, w, d, NULL); // assert f != NULL
			if(f == NULL){
				continue;
			}
			fidx[2] = offset_edgerefv(g->erefs, f);
			{
				e = ref_edgev(g->edges, ref_edgerefv(g->erefs, fidx[0])->idx);
				cut_edge_graph(g, e);
				e = ref_edgev(g->edges, ref_edgerefv(g->erefs, fidx[2])->idx);
				cut_edge_graph(g, e);
				e = ref_edgev(g->edges, ref_edgerefv(g->erefs, fidx[1])->idx);
				revive_edge_graph(g, e);
				ret ++;
			}
		}
	}
	free_tracev(path);
	return ret;
}

static inline u8i cut_relative_low_cov_edges_graph(Graph *g, float lowf){
	node_t *n;
	edge_t *e;
	edge_ref_t *f;
	u4v *ecovs, *emaxs;
	u8i idx, x, ret;
	u4i k, max, min, val, i;
	for(idx=0;idx<g->edges->size;idx++){
		e = ref_edgev(g->edges, idx);
		e->flag = 0;
	}
	ecovs = init_u4v(64);
	emaxs = init_u4v(64);
	ret = 0;
	for(idx=0;idx<g->nodes->size;idx++){
		n = ref_nodev(g->nodes, idx);
		if(n->closed) continue;
		clear_u4v(ecovs);
		max = 0;
		min = WT_MAX_EDGE_COV;
		for(k=0;k<2;k++){
			x = n->edges[k].idx;
			while(x){
				f = ref_edgerefv(g->erefs, x);
				x = f->next;
				e = ref_edgev(g->edges, f->idx);
				if(e->closed) continue;
				if(e->cov > max) max = e->cov;
				if(e->cov < min) min = e->cov;
				push_u4v(ecovs, (Int(e->off) << 12) | (e->cov));
			}
		}
		if(min >= max * lowf){
			continue;
		}
		min = max * lowf;
		sort_array(ecovs->buffer, ecovs->size, u4i, num_cmpgt(Int(a) >> 12, Int(b) >> 12));
		clear_u4v(emaxs);
		max = 0;
		for(i=ecovs->size;i>0;i--){
			val = ecovs->buffer[i - 1] & 0xFFF;
			if(val > max) max = val;
			push_u4v(emaxs, max);
		}
		reverse_u4v(emaxs);
		for(k=0;k<2;k++){
			x = n->edges[k].idx;
			while(x){
				f = ref_edgerefv(g->erefs, x);
				x = f->next;
				e = ref_edgev(g->edges, f->idx);
				if(e->closed) continue;
				if(e->cov >= min) continue;
				for(i=0;i<ecovs->size;i++){
					if(e->off >= (Int(ecovs->buffer[i]) >> 12)){
						break;
					}
				}
				if(i == ecovs->size) continue; // Never happen
				if(e->cov < lowf * emaxs->buffer[i]){
					cut_edge_graph(g, e);
					ret ++;
				}
			}
		}
	}
	free_u4v(ecovs);
	free_u4v(emaxs);
	return ret;
}

static inline u4i rescue_low_cov_tip_edges_core(Graph *g, u8i nid){
	node_t *n, *w, *ww;
	edge_t *e, *ee;
	edge_ref_t *f;
	u8i idx, wid;
	u4i k, dir, ret;
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

static inline  u8i rescue_low_cov_tip_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	u8i nid, ret;
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

static inline int rescue_mercy_edge_core_graph(Graph *g, u4i rid, BitVec *tips[2]){
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

static inline u8i rescue_mercy_edges_graph(Graph *g){
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

static inline u4i rescue_weak_tip_lnks_core(Graph *g, u8i nid){
	frg_t *n, *w, *ww;
	lnk_t *e, *ee;
	edge_ref_t *f;
	u8i idx, wid;
	u4i k, dir, ret;
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
			if(e->closed) continue;
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

static inline u8i rescue_weak_tip_lnks2_graph(Graph *g){
	u8i nid, ret;
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		ret += rescue_weak_tip_lnks_core(g, nid);
	}
	return ret;
}

static inline u8i rescue_weak_tip_lnks_graph(Graph *g){
	pathv *heap;
	u8v *weaks[2];
	frg_t *n;
	lnk_t *e;
	edge_ref_t *f;
	u8i nid, i, vst, idx, eidx, ret;
	u4i k;
	int max_len;
	weaks[0] = init_u8v(g->frgs->size);
	weaks[1] = init_u8v(g->frgs->size);
	ret = 0;
	for(nid=0;nid<g->frgs->size;nid++){
		n = ref_frgv(g->frgs, nid);
		n->bt_visit = 0;
	}
	heap = init_pathv(32);
	vst = 1;
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
					if(e->closed && e->weak){
						if(eidx == 0) eidx = f->idx;
						else eidx = MAX_VALUE_U8;
					}
					idx = f->next;
				}
			}
			if(eidx != 0 && eidx != MAX_VALUE_U8){
				e = ref_lnkv(g->lnks, eidx);
				max_len = wander_path_graph(g, nid, k, heap, vst, e->off);
				if(max_len >= e->off){
					eidx |= 1LLU << 63;
				}
			}
			push_u8v(weaks[k], eidx);
		}
	}
	free_pathv(heap);
	for(k=0;k<2;k++){
		for(i=0;i<weaks[k]->size;i++){
			if(weaks[k]->buffer[i] == 0 || weaks[k]->buffer[i] == MAX_VALUE_U8) continue;
			e = ref_lnkv(g->lnks, weaks[k]->buffer[i] & (MAX_VALUE_U8 >> 1));
			if(i != e->frg1) continue;
			if((weaks[0]->buffer[e->frg2] & (MAX_VALUE_U8 >> 1)) == (weaks[k]->buffer[i] & (MAX_VALUE_U8 >> 1))){
				if(!(weaks[0]->buffer[e->frg2] & (1LLU < 63)) || !(weaks[k]->buffer[i] & (1LLU << 63))){
					ret ++;
					revive_lnk_graph(g, e);
				}
			}
			if((weaks[1]->buffer[e->frg2] & (MAX_VALUE_U8 >> 1)) == (weaks[k]->buffer[i] & (MAX_VALUE_U8 >> 1))){
				if(!(weaks[1]->buffer[e->frg2] & (1LLU < 63)) || !(weaks[k]->buffer[i] & (1LLU << 63))){
					ret ++;
					revive_lnk_graph(g, e);
				}
			}
		}
	}
	free_u8v(weaks[0]);
	free_u8v(weaks[1]);
	return ret;
}

static inline int _scoring_edge_orders(Graph *g, u8i fidx){
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


static inline u4i reduce_transitive_edges_core_graph(Graph *g, u8i nid, u8v *edges, UUhash *hash, u4i closed_val){
	node_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	edge_t *e, *e1, *e2;
	UUhash_t *u;
	u8i idx, nid2, nid3;
	u4i i, k, k2, k3, k4, ret;
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
		//sort_array(edges->buffer, edges->size, u8i, num_cmpgt(_scoring_edge_orders(g, a), _scoring_edge_orders(g, b)));
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

static inline u8i reduce_transitive_edges_graph(Graph *g){
	u8v *edges;
	UUhash *hash;
	u8i nid, ret;
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

static inline void set_init_ends_graph(Graph *g){
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

static inline u4i myers_transitive_reduction_core_graph(Graph *g, u8i nid, float _fuzz, int closed){
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

static inline u8i myers_transitive_reduction_graph(Graph *g, float fuzz){
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

static inline u4i myers_transitive_reduction_core_frg_graph(Graph *g, u8i nid, float _fuzz){
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

static inline u8i myers_transitive_reduction_frg_graph(Graph *g, float fuzz){
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

static inline u4i detach_repetitive_frg_core_graph(Graph *g, u8i nid, u4i max_dist, u8i visit, u8v *nds, u8v *heap){
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

static inline u4i detach_repetitive_frg_graph(Graph *g, u4i max_dist){
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

static inline u4i reduce_transitive_lnks_core_graph(Graph *g, u8i nid, u8v *lnks, UUhash *hash, u4i closed_val){
	frg_t *n, *w, *v;
	edge_ref_t *f, *f2, *f3;
	lnk_t *e, *e1, *e2;
	UUhash_t *u;
	u8i idx, nid2, nid3;
	u4i i, k, k2, k3, k4, ret;
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
		//sort_array(lnks->buffer, lnks->size, u8i, num_cmpgt(_scoring_lnk_orders(g, a), _scoring_lnk_orders(g, b)));
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

static inline u8i reduce_transitive_lnks_graph(Graph *g){
	u8v *lnks;
	UUhash *hash;
	u8i nid, ret;
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

static inline u8i trim_tip_core_graph(Graph *g, uint16_t max_step, tracev *path, u8i nid, int hard_trim){
	trace_t *t, T;
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	u8i ret, idx;
	u4i i, dir, step, step2, found, n_in;
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
	msg1 = WT_TRACE_MSG_ZERO;
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

static inline u8i trim_tips_graph(Graph *g, uint16_t max_step, int hard_trim){
	tracev *path;
	node_t *n;
	u8i ret, nid;
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
static inline u4i trim_blunt_tip_core_graph(Graph *g, u8i nid){
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

static inline u8i trim_blunt_tips_graph(Graph *g){
	node_t *n;
	u8i ret, nid;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(trim_blunt_tip_core_graph(g, nid)) ret ++;
	}
	return ret;
}

static inline u8i trim_frgtip_core_graph(Graph *g, int max_len, pathv *path, u8i nid){
	path_t *t, T;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e;
	u8i ret, idx;
	u4i i, dir, found, n_in;
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
		return 0;
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

static inline u8i trim_frgtips_graph(Graph *g, int max_len){
	pathv *path;
	frg_t *n;
	u8i ret, nid;
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
	u8i dir:1, ind:1, step:8, bt:16, ending:16, keep:2;
	int score:20;
} bt_t;
define_list(btv, bt_t);
#define WT_MAX_BTIDX	0xFFFF

static inline u4i pop_bubble_backtrace_graph(Graph *g, btv *bts, u4i idx){
	bt_t *bt;
	u4i ret;
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
	u8i dir:1, ind:1, step:8, bt:16, ending:16, score:22;
} frg_bt_t;
define_list(frgbtv, frg_bt_t);
u4i pop_frg_bubble_backtrace_graph(Graph *g, frgbtv *bts, u4i idx){
	frg_bt_t *bt;
	u4i ret;
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

static inline u4i pop_bubble2_backtrace_graph(Graph *g, btv *bts, u4i _idx){
	bt_t *bt;
	u4i ret, i, idx;
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

static inline u4i pop_frg_bubble2_backtrace_graph(Graph *g, frgbtv *bts, u4i _idx){
	frg_bt_t *bt;
	u4i ret, i, idx;
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

static inline u4i safe_cut_redundant_edges_graph(Graph *g, btv *bts, bt_t *b1, bt_t *b2){
	u4i ret;
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

static inline u4i safe_cut_redundant_lnks_graph(Graph *g, frgbtv *bts, frg_bt_t *b1, frg_bt_t *b2){
	u4i ret;
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

static inline u4i pop_bubble_core_graph(Graph *g, uint16_t max_step, btv *bts, u4v *heap, u8i nid, u4i dir, u8i visit, int safe){
	bt_t *bt, *tb;
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	u8i ret, idx;
	u4i bidx, i, lst, unclosed;
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
	array_heap_push(heap->buffer, heap->size, heap->cap, u4i, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	n->single_in = 1;
	unclosed = 0;
	while(heap->size && heap->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, u4i, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
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
			//tb->score = bt->score + num_min(0, e->cov - 20); // set normal e->cov = 20
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
				//if(num_cmpgt(tb->score, bts->buffer[tb->n->bt_idx].score)){
				if(tb->step > bts->buffer[tb->n->bt_idx].step){
					bts->buffer[tb->n->bt_idx].ending = i;
					tb->n->bt_idx = i;
					tb->ending = 0;
				}
				if(count_living_edges_graph(g, tb->n, tb->dir)){
					array_heap_push(heap->buffer, heap->size, heap->cap, u4i, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
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

static inline u8i pop_bubbles_graph(Graph *g, uint16_t max_step, int safe){
	btv *bts;
	u4v *heap;
	node_t *n;
	u8i nid, visit, ret, _ret;
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

static inline u4i pop_frg_bubble_core_graph(Graph *g, uint16_t max_step, frgbtv *bts, u4v *heap, u8i nid, u4i dir, u8i visit){
	frg_bt_t *bt, *tb;
	frg_t *n;
	edge_ref_t *f;
	lnk_t *e;
	u8i ret, idx;
	u4i bidx, i, lst, unclosed;
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
	array_heap_push(heap->buffer, heap->size, heap->cap, u4i, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	n->single_in = 1;
	unclosed = 0;
	while(heap->size && heap->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, u4i, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
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
					array_heap_push(heap->buffer, heap->size, heap->cap, u4i, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
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

static inline u8i pop_frg_bubbles_graph(Graph *g, uint16_t max_step){
	frgbtv *bts;
	u4v *heap;
	frg_t *n;
	u8i nid, visit, ret, _ret;
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

static inline u4i resolve_yarn_core_graph(Graph *g, u4i max_step, btv *bts, u4v *heap, u8i nid, u4i dir, u8i visit){
	bt_t *bt, *tb;
	node_t *n, *m;
	edge_ref_t *f;
	edge_t *e;
	u8i ret, idx, tip_idx;
	u4i bidx, i, lst, tip;
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
	array_heap_push(heap->buffer, heap->size, heap->cap, u4i, bts->size - 1, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
	n->bt_visit = visit;
	n->bt_idx = bts->size - 1;
	tip = 0; tip_idx = WT_MAX_BTIDX;
	n_in = 1;
	while(heap->size && bts->size < WT_MAX_BTIDX){
		bidx = array_heap_pop(heap->buffer, heap->size, heap->cap, u4i, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
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
			//if(tb->n == bts->buffer[1].n){
				//return 0;
			//}
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
					array_heap_push(heap->buffer, heap->size, heap->cap, u4i, tb->n->bt_idx, num_cmpx(bts->buffer[a].step, bts->buffer[b].step, bts->buffer[b].score, bts->buffer[a].score));
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
static inline u8i resolve_yarns_graph(Graph *g, u4i max_step){
	btv *bts;
	u4v *heap;
	tracev *path;
	trace_t *t;
	node_t *n;
	u8i nid, visit, ret, _ret;
	int dir;
	ret = 0;
	for(nid=0;nid<g->nodes->size;nid++) g->nodes->buffer[nid].bt_visit = 0;
	bts = init_btv(32);
	heap = init_u4v(32);
	path = init_tracev(4);
	visit = 0;
#if DEBUG
	if(max_step == 1000000){ // never happen
		print_local_dot_graph(g, "1.dot", NULL);
	}
#endif
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(n->edges[0].cnt <= 1 && n->edges[1].cnt > 1){
			if(n->edges[0].cnt == 1){
				clear_tracev(path);
				t = next_ref_tracev(path);
				t->node = nid;
				t->dir = 0;
				if(linear_trace_graph(g, path, 4, NULL) < 4){
					// solve yarn from linear nodes
					continue;
				}
			}
			dir = 1;
		} else if(n->edges[1].cnt <= 1 && n->edges[0].cnt > 1){
			if(n->edges[0].cnt == 1){
				clear_tracev(path);
				t = next_ref_tracev(path);
				t->node = nid;
				t->dir = 1;
				if(linear_trace_graph(g, path, 4, NULL) < 4){
					// solve yarn from linear nodes
					continue;
				}
			}
			dir = 0;
		} else continue;
		_ret = resolve_yarn_core_graph(g, max_step, bts, heap, nid, dir, ++visit);
		if(_ret) ret ++;
	}
	free_btv(bts);
	free_u4v(heap);
	free_tracev(path);
	return ret;
}

static inline u8i remove_boomerangs_frg_graph(Graph *g, u4i max_frg_len){
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

static inline u8i cut_weak_branches_frg_graph(Graph *g){
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

static inline u8i mask_all_branching_nodes_graph(Graph *g){
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

static inline u8i gen_unitigs_graph(Graph *g){
	tracev *path;
	u4v *lens;
	trace_t *t;
	node_t *n;
	u8i nid, nutg, i;
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
		push_u4v(lens, cal_offset_traces_graph(g, path, 0, path->size, 0) * KBM_BIN_SIZE);
		for(i=0;i<path->size;i++){
			ref_nodev(g->nodes, path->buffer[i].node)->rep_idx = g->utgs->size;
		}
		push_vplist(g->utgs, path);
	}
	fprintf(KBM_LOGF, "[%s] ", date()); num_n50(lens, KBM_LOGF); fprintf(KBM_LOGF, "\n");
	free_u4v(lens);
	return nutg;
}

static inline seqletv* path2seqlets_graph(Graph *g, pathv *path){
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
			//for(j=p2->tx+1;j<=p2->ty;j++){
				//t1 = ref_tracev(g->traces, frg2->toff + frg2->tcnt - 1 - j);
				//t2 = ref_tracev(g->traces, frg2->toff + frg2->tcnt - 0 - j);
			for(j=p2->ty;j>p2->tx;j--){
				t1 = ref_tracev(g->traces, frg2->toff + j - 1);
				t2 = ref_tracev(g->traces, frg2->toff + j);
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

static inline u8i gen_contigs_graph(Graph *g, FILE *out){
	pathv *path;
	seqletv *qs;
	seqlet_t *q;
	path_t *t;
	frg_t *n;
	u8i nid, nctg, i, off;
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
		if(((int)q->off + (int)q->len) * KBM_BIN_SIZE < (int)g->min_ctg_len){
			free_seqletv(qs);
			continue;
		}
		off = 0;
		for(i=0;i<qs->size;i++){
			q = ref_seqletv(qs, i);
			q->off = off;
			off += q->len - g->reglen;
		}
		if(out){
			for(i=0;i<path->size;i++){
				t = ref_pathv(path, i);
				fprintf(out, "ctg%d\tF%d\t%c\t%d\n", (int)g->ctgs->size, t->frg, "+-*@"[t->dir], t->off * KBM_BIN_SIZE);
			}
		}
		push_vplist(g->ctgs, qs);
	}
	free_pathv(path);
	g->major_nctg = g->ctgs->size;
	return nctg;
}

static inline u8i gen_complex_contigs_graph(Graph *g){
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

static inline void n50_stat_contigs_graph(Graph *g){
	seqletv *qs;
	seqlet_t *q;
	u4v *lens;
	int len;
	u8i i;
	lens = init_u4v(g->major_nctg + 1);
	for(i=0;i<g->major_nctg;i++){
		qs = (seqletv*)get_vplist(g->ctgs, i);
		q = ref_seqletv(qs, qs->size - 1);
		len = ((int)q->off + (int)q->len) * KBM_BIN_SIZE;
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
static inline u4i count_isolated_reads_graph(Graph *g){
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

static inline void print_dot_subgraph(Graph *g, subnodehash *nodes, subedgev *edges, FILE *out){
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

static inline void fprintf_dot_subgraph(Graph *g, subnodehash *nodes, subedgev *edges, char *name_prefix, char *name_suffix){
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

static inline subedge_t* find_edge_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
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

static inline int cut_edge_core_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
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

static inline int cut_edge_subgraph(subnodehash *nodes, subedgev *edges, u4i node1, int dir1, u4i node2, int dir2){
	return cut_edge_core_subgraph(nodes, edges, node1, dir1, node2, dir2)
		+ cut_edge_core_subgraph(nodes, edges, node2, !dir2, node1, !dir1);
}

static inline u4i unitigs2frgs_graph(Graph *g, int ncpu){
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

static inline int scan_rd_lnk_core(Graph *g, u4i rid, lnk_t *lnk, u8v *regids){
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

static inline int scan_nd_lnk_core(Graph *g, u8i nid, lnk_t *lnk){
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

static inline u4i gen_lnks_graph(Graph *g, int ncpu, FILE *log){
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
			if(l->cov < g->max_node_cov_sg){
				l->weak = 1;
				//l->closed = WT_EDGE_CLOSED_LESS;
				l->closed = 0;
			} else {
				l->weak = 0;
				l->closed = 0;
			}
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

static inline int gen_seq_traces_graph(Graph *g, tracev *path, String *seq){
	trace_t *t1, *t2;
	reg_t *reg, *r1, *r2;
	edge_t *e;
	u4i i;
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
							inc = (r2->beg - r1->end) * KBM_BIN_SIZE;
							if(inc <= 0) break;
							encap_string(seq, inc);
							seq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[r1->rid].binoff * KBM_BIN_SIZE+ (r1->end * KBM_BIN_SIZE), inc, seq->string + seq->size);
							seq->size += inc;
						} else {
							if(!(t1->dir ^ r1->dir)){ r1++; r2++; continue; }
							inc = (r1->beg - r2->end) * KBM_BIN_SIZE;
							if(inc <= 0) break;
							encap_string(seq, inc);
							revseq_basebank(g->kbm->rdseqs, (g->kbm->reads->buffer[r1->rid].seqoff + r2->end) * KBM_BIN_SIZE, inc, seq->string + seq->size);
							seq->size += inc;
						}
						inc = 0;
						found = 1; break;
					}
				}
				if(found == 0){ inc = e->off; break; }
			} while(0);
			if(inc > 0){ inc = 0; while(inc++ < e->off * KBM_BIN_SIZE) add_char_string(seq, 'N'); }
			else if(inc < 0){
				if(seq->size + inc < 0) seq->size = 0;
				else seq->size += inc;
				seq->string[seq->size] = '\0';
			}
		}
		t1 = t2;
		reg = ref_regv(g->regs, ref_nodev(g->nodes, t1->node)->regs.idx);
		inc = (reg->end - reg->beg) * KBM_BIN_SIZE;
		encap_string(seq, inc);
		if(t1->dir ^ reg->dir) revseq_basebank(g->kbm->rdseqs, (g->kbm->reads->buffer[reg->rid].seqoff + reg->beg) * KBM_BIN_SIZE, inc, seq->string + seq->size);
		else                   seq_basebank(g->kbm->rdseqs, (g->kbm->reads->buffer[reg->rid].seqoff + reg->beg) * KBM_BIN_SIZE, inc, seq->string + seq->size);
		seq->size += inc;
	}
	return seq->size;
}

typedef struct {
	u8i rid:26, dir:1, beg:18, end:18, view:1;
} lay_reg_t;
define_list(layregv, lay_reg_t);

typedef struct {
	seqlet_t edge;
	u8i roff:48, rcnt:16;
} lay_t;
define_list(layv, lay_t);

static inline void gen_lay_regs_core_graph(Graph *g, seqlet_t *q, layregv *regs, BitVec *rdbits, int closed){
	node_t *n1, *n2;
	reg_t *r1, *r2;
	u4i rid, beg, end;
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
			if((closed || (r1->closed == 0 && r2->closed == 0)) && (rdbits == NULL || get_bitvec(rdbits, rid))){
				if(r1->beg < r2->beg){
					if(q->dir1 ^ r1->dir){ r1 ++; r2 ++; continue; }
					beg = r1->beg; end = r2->end;
					push_layregv(regs, (lay_reg_t){rid, 0, beg, end, 0});
				} else {
					if(!(q->dir1 ^ r1->dir)){ r1 ++; r2 ++; continue; }
					beg = r2->beg; end = r1->end;
					push_layregv(regs, (lay_reg_t){rid, 1, beg, end, 0});
				}
			}
			r1 ++; r2 ++;
		}
	}
}

typedef struct {
	u4i rid, closed;
	u8i nodes[2];
} readreg_t;
define_list(readregv, readreg_t);

static inline u4i densify_seqlet_graph(Graph *g, seqlet_t *q, seqletv *qs, int minoff, readregv *rds, subnodehash *nodes, subedgev *edges, subnodev *heap, FILE *log){
	seqlet_t *q2;
	node_t *nd;
	read_t *rd;
	reg_t  *rg;
	readreg_t *rr;
	subnode_t N, *n, *n1, *n2;
	subedge_t *e;
	u8i idx, edx;
	u4i i, k, k1, k2, d, cov, flg, ret;
	int exists, off;
	if(q->len < minoff){
		push_seqletv(qs, *q);
		return 1;
	}
	clear_readregv(rds);
	for(k=0;k<2;k++){
		nd = ref_nodev(g->nodes, k? q->node2 : q->node1);
		d  = k? !q->dir2 : q->dir1;
		for(i=0;i<nd->regs.cnt;i++){
			rg = ref_regv(g->regs, nd->regs.idx + i);
			rr = next_ref_readregv(rds);
			rr->rid = rg->rid;
			rr->closed = 0;
			rr->nodes[d ^ rg->dir] = rg->node;
			rr->nodes[!d ^ rg->dir] = MAX_U8;
		}
	}
	sort_array(rds->buffer, rds->size, readreg_t, num_cmpgt(a.rid, b.rid));
	for(i=1;i<rds->size;i++){
		rr = ref_readregv(rds, i - 1);
		if(rr->rid == rds->buffer[i].rid){
			if(rds->buffer[i].nodes[0] != MAX_U8) rr->nodes[0] = rds->buffer[i].nodes[0];
			if(rds->buffer[i].nodes[1] != MAX_U8) rr->nodes[1] = rds->buffer[i].nodes[1];
			rds->buffer[i].closed = 1; // max occ of rid is 2
		}
	}
	// prepare nodes in subgraph
	clear_subnodehash(nodes);
	clear_subedgev(edges);
	next_ref_subedgev(edges);
	memset(&N, 0, sizeof(subnode_t));
	N.cov = 1;
	N.bt_nidx = MAX_BT_NIDX;
	cov = 0;
	for(i=0;i<rds->size;i++){
		rr = ref_readregv(rds, i);
		if(rr->closed) continue;
		if(rr->nodes[0] != MAX_U8 && rr->nodes[1] != MAX_U8){
			cov ++;
		}
		rd = ref_readv(g->reads, rr->rid);
		flg = (rr->nodes[0] == MAX_U8);
		idx = rd->regs.idx;
		while(idx){
			rg = ref_regv(g->regs, idx);
			idx = rg->read_link;
			if(flg == 0){
				if(rg->node == rr->nodes[0]){
					flg = 1;
				}
			}
			if(!flg) continue;
			N.node = rg->node;
			n = prepare_subnodehash(nodes, N, &exists);
			if(exists){
				n->cov ++;
			} else {
				*n = N;
			}
			if(rg->node == rr->nodes[1]){
				flg = 0;
				break;
			}
		}
	}
	// mask low cov nodes
	reset_iter_subnodehash(nodes);
	while((n = ref_iter_subnodehash(nodes))){
		if(n->cov < cov) n->closed = 1;
	}
	// build edges
	for(i=0;i<rds->size;i++){
		rr = ref_readregv(rds, i);
		if(rr->closed) continue;
		rd = ref_readv(g->reads, rr->rid);
		flg = (rr->nodes[0] == MAX_U8);
		n1 = NULL;
		k1 = 0;
		off = 0;
		idx = rd->regs.idx;
		while(idx){
			rg = ref_regv(g->regs, idx);
			idx = rg->read_link;
			if(flg == 0){
				if(rg->node == rr->nodes[0]){
					flg = 1;
				}
			}
			if(!flg) continue;
			do {
				N.node = rg->node;
				n2 = get_subnodehash(nodes, N);
				k2 = rg->dir;
				if(n2->closed) break;
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
						e->off = rg->beg - off;
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
						e->off = rg->beg - off;
						e->next = n2->edges[!k2].idx;
						n2->edges[!k2].idx = edx;
						n2->edges[!k2].cnt ++;
					}
				}
			} while(0);
			n1 = n2;
			k1 = k2;
			off = rg->end;
			if(rg->node == rr->nodes[1]){
				flg = 0;
				break;
			}
		}
	}
	// searching a most dense path from q->node1 to q->node2
	N.node = q->node1;
	n = get_subnodehash(nodes, N);
	n->visit = 0;
	n->bt_step = 0;
	n->bt_score = 0;
	n->bt_dir = !q->dir1;
	n->bt_nidx = MAX_BT_NIDX;
	// NB: below codes cannot grant to find the most dense path, but just useful in many cases
	clear_subnodev(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, subnode_t, *n, num_cmp(a.bt_score, b.bt_score));
	while(heap->size){
		N = heap->buffer[0];
		array_heap_remove(heap->buffer, heap->size, heap->cap, subnode_t, 0, num_cmp(a.bt_score, b.bt_score));
		n1 = get_subnodehash(nodes, N);
		if(n1->visit){
			if(n1->node == q->node2){
				if(n1->bt_step < N.bt_step){
					*n1 = N;
					n1->visit = 1;
				}
			}
			continue;
		} else {
			*n1 = N;
			n1->visit = 1;
		}
		k1 = !n1->bt_dir;
		idx = n1->edges[k1].idx;
		while(idx){
			e = ref_subedgev(edges, idx);
			idx = e->next;
			n2 = e->node;
			if(n2->bt_step >= n1->bt_step + 1){ // init n2->bt_step = 0
				continue;
			}
			if(n2->visit && n2->node != q->node2){
				continue;
			}
			N = *n2;
			N.bt_dir = !e->dir;
			N.bt_nidx = offset_subnodehash(nodes, n1);
			N.bt_step = n1->bt_step + 1;
			N.bt_score = n1->bt_score + g->reglen + num_max(0, e->off);
			array_heap_push(heap->buffer, heap->size, heap->cap, subnode_t, N, num_cmp(a.bt_score, b.bt_score));
		}
	}
	N.node = q->node2;
	n = get_subnodehash(nodes, N);
	if(n->visit == 0){
#if DEBUG
		if(nodes->count == 1000000){
			fprint_subgraph_dot(g, 0, nodes, edges, "1.dot");
		}
#endif
		push_seqletv(qs, *q);
		return 1;
	}
	if(n->bt_dir == q->dir2){
		push_seqletv(qs, *q);
		return 1;
	}
	n1 = NULL;
	while(1){
		if(n->bt_nidx == MAX_BT_NIDX){
			push_seqletv(qs, *q);
			return 1;
		}
		n2 = ref_subnodehash(nodes, n->bt_nidx);
		n->bt_dir = !n->bt_dir;
		n->bt_nidx = n1? offset_subnodehash(nodes, n1) : MAX_BT_NIDX;
		if(n->node == q->node1){
			break;
		}
		n1 = n;
		n = n2;
	}
	ret = 0;
	while(1){
		n2 = ref_subnodehash(nodes, n->bt_nidx);
		ret ++;
		q2 = next_ref_seqletv(qs);
		q2->node1 = n->node;
		q2->dir1 = n->bt_dir;
		q2->node2 = n2->node;
		q2->dir2 = n2->bt_dir;
		q2->off = 0;
		off = 0;
		edx = n->edges[n->bt_dir].idx;
		while(edx){
			e = ref_subedgev(edges, edx);
			edx = e->next;
			if(e->node->node == q2->node2){
				off = e->off;
				break;
			}
		}
		q2->len = 2 * g->reglen + off;
		if(n2->node == q->node2){
			break;
		}
		n = n2;
	}
	if(ret > 1 && log){
		fprintf(log, "N%llu -> N%llu edge_len=%d densified into %u\n", (u8i)q->node1, (u8i)q->node2, q->len, ret);
	}
	return ret;
}

thread_beg_def(mlay);
Graph *g;
seqletv *path;
u8i div_idx;
u8i pb, pe;
layv    *lays;
layregv *regs;
int all_regs;
FILE *log;
thread_end_def(mlay);

thread_beg_func(mlay);
subnodehash *nodes;
subnodev *heap;
subedgev *edges;
readregv *rds;
seqletv *lets;
BitVec *rdbits;
seqlet_t *_let, *let;
lay_t *lay;
u8i i;
u4i j;
lets = init_seqletv(4);
nodes = init_subnodehash(13);
edges = init_subedgev(32);
rds = init_readregv(32);
heap = init_subnodev(32);
rdbits = init_bitvec(mlay->g->reads->size);
thread_beg_loop(mlay);
clear_layv(mlay->lays);
clear_layregv(mlay->regs);
for(i=mlay->pb;i<mlay->pe;i++){
	_let = ref_seqletv(mlay->path, i);
	clear_seqletv(lets);
	densify_seqlet_graph(mlay->g, _let, lets, 12, rds, nodes, edges, heap, mlay->log);
	if(lets->size > 1){
		for(j=0;j<rds->size;j++){
			if(rds->buffer[j].closed == 0){
				one_bitvec(rdbits, rds->buffer[j].rid);
			}
		}
	}
	for(j=0;j<lets->size;j++){
		let = ref_seqletv(lets, j);
		lay = next_ref_layv(mlay->lays);
		lay->edge = *let;
		lay->roff = mlay->regs->size;
		gen_lay_regs_core_graph(mlay->g, let, mlay->regs, lets->size > 1? rdbits : NULL, mlay->all_regs);
		lay->rcnt = mlay->regs->size - lay->roff;
		sort_array(mlay->regs->buffer + lay->roff, lay->rcnt, lay_reg_t, num_cmpgt(b.end - b.beg, a.end - a.beg));
		if(lay->rcnt == 0 && mlay->log){
			thread_beg_syn(mlay);
			fprintf(mlay->log, " -- N%llu(%c) -> N%llu(%c) has no read path --\n", (u8i)let->node1, "+-"[let->dir1], (u8i)let->node2, "+-"[let->dir2]); fflush(mlay->log);
			thread_end_syn(mlay);
		}
	}
	if(lets->size > 1){
		for(j=0;j<rds->size;j++){
			if(rds->buffer[j].closed == 0){
				zero_bitvec(rdbits, rds->buffer[j].rid);
			}
		}
	}
}
thread_end_loop(mlay);
free_seqletv(lets);
free_subnodev(heap);
free_readregv(rds);
free_subedgev(edges);
free_subnodehash(nodes);
free_bitvec(rdbits);
thread_end_func(mlay);

static inline u8i print_ctgs_graph(Graph *g, u8i uid, u8i beg, u8i end, char *prefix, char *lay_suffix, u4i ncpu, FILE *log){
	FILE *o_lay;
	BufferedWriter *bw;
	layv *lays;
	layregv *regs;
	seqletv *path;
	u8v *divs;
	seqlet_t *t;
	lay_t *lay;
	lay_reg_t *reg;
	u8i i, pb, pe, d, div_idx, ret;
	u4i j, c, len, bsize, nrun;
	thread_preprocess(mlay);
	o_lay = open_file_for_write(prefix, lay_suffix, 1);
	bw = zopen_bufferedwriter(o_lay, 1024 * 1024, ncpu, 0);
	lays = init_layv(32);
	regs = init_layregv(32);
	thread_beg_init(mlay, ncpu);
	mlay->g = g;
	mlay->path = NULL;
	mlay->lays = init_layv(32);
	mlay->regs = init_layregv(32);
	mlay->pb = 0;
	mlay->pe = 0;
	mlay->div_idx = MAX_U8;
	mlay->all_regs = 1;
	mlay->log  = log;
	thread_end_init(mlay);
	ret = 0;
	bsize = 100;
	divs = init_u8v(1024);
	for(i=beg;i<end;i++){
		path = (seqletv*)get_vplist(g->ctgs, i);
		if(path->size == 0) continue;
		len = path->buffer[path->size - 1].off + path->buffer[path->size - 1].len;
		len = len * KBM_BIN_SIZE;
		//clear_and_inc_layv(lays, path->size);
		clear_layv(lays);
		clear_u8v(divs);
		clear_layregv(regs);
		div_idx = 0;
		pe = 0;
		nrun = 0;
		while(1){
			pb = pe;
			pe = num_min(pb + bsize, path->size);
			thread_wait_one(mlay);
			if(mlay->div_idx != MAX_U8){
				nrun --;
				divs->buffer[mlay->div_idx * 2 + 0] = lays->size;
				divs->buffer[mlay->div_idx * 2 + 1] = lays->size + mlay->lays->size;
				for(j=0;j<mlay->lays->size;j++){
					mlay->lays->buffer[j].roff += regs->size;
				}
				append_layv(lays, mlay->lays);
				clear_layv(mlay->lays);
				append_layregv(regs, mlay->regs);
				clear_layregv(mlay->regs);
				mlay->pb = 0;
				mlay->pe = 0;
				mlay->div_idx = MAX_U8;
			}
			if(pb < pe){
				inc_u8v(divs, 2);
				mlay->div_idx = div_idx ++;
				mlay->path = path;
				mlay->pb = pb;
				mlay->pe = pe;
				thread_wake(mlay);
				nrun ++;
			} else if(nrun == 0){
				break;
			}
		}
		//thread_apply_all(mlay, EXPR(mlay->path = path));
		uid ++;
		ret ++;
		{
			beg_bufferedwriter(bw);
			fprintf(bw->out, ">ctg%llu nodes=%llu len=%u\n", uid, (u8i)path->size + 1, len);
			if(log) fprintf(log, "OUTPUT_CTG\tctg%d -> ctg%d nodes=%llu len=%u\n", (int)i, (int)uid, (u8i)path->size + 1, len);
			for(d=0;d<div_idx;d++){
				pb = divs->buffer[d * 2 + 0];
				pe = divs->buffer[d * 2 + 1];
				flush_bufferedwriter(bw);
				for(j=pb;j<pe;j++){
					lay = ref_layv(lays, j);
					if(lay->rcnt == 0) continue;
					t = &lay->edge;
					fprintf(bw->out, "E\t%d\tN%llu\t%c\tN%llu\t%c\n", (int)t->off * KBM_BIN_SIZE, (u8i)t->node1, "+-"[t->dir1], (u8i)t->node2, "+-"[t->dir2]);
					for(c=0;c<lay->rcnt;c++){
						reg = ref_layregv(regs, lay->roff + c);
						fprintf(bw->out, "%c\t%s\t%c\t%d\t%d\t", "Ss"[reg->view], g->kbm->reads->buffer[reg->rid].tag, "+-"[reg->dir], reg->beg * KBM_BIN_SIZE, (reg->end - reg->beg) * KBM_BIN_SIZE);
						if(reg->dir){
							print_revseq_basebank(g->kbm->rdseqs, (g->kbm->reads->buffer[reg->rid].seqoff + reg->beg) * KBM_BIN_SIZE, (reg->end - reg->beg) * KBM_BIN_SIZE, bw->out);
						} else {
							print_seq_basebank(g->kbm->rdseqs, (g->kbm->reads->buffer[reg->rid].seqoff + reg->beg) * KBM_BIN_SIZE, (reg->end - reg->beg) * KBM_BIN_SIZE, bw->out);
						}
						fprintf(bw->out, "\n");
					}
				}
			}
			end_bufferedwriter(bw);
		}
	}
	close_bufferedwriter(bw);
	fclose(o_lay);
	thread_beg_close(mlay);
	free_layregv(mlay->regs);
	free_layv(mlay->lays);
	thread_end_close(mlay);
	fprintf(KBM_LOGF, "[%s] output %u contigs\n", date(), (u4i)ret);
	free_layv(lays);
	free_layregv(regs);
	free_u8v(divs);
	return uid;
}

static inline u4i print_traces_graph(Graph *g, tracev *path, FILE *out){
	String *str;
	trace_t *t1, *t2;
	node_t *n1, *n2;
	reg_t *r1, *r2;
	edge_ref_t *f;
	edge_t *e;
	int offset, fst;
	u8i beg, end;
	u4i i, rid;
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
						beg = r1->beg * KBM_BIN_SIZE; end = r2->end * KBM_BIN_SIZE;
						fprintf(out, "S\t%s\t", g->kbm->reads->buffer[rid].tag);
						fprintf(out, "+\t%d\t%d\t", (int)beg, (int)(end - beg));
						encap_string(str, end - beg);
						fwdseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[rid].seqoff * KBM_BIN_SIZE + beg, end - beg, str->string);
						fputs(str->string, out);
					} else {
						if(!(t1->dir ^ r1->dir)){ r1 ++; r2 ++; continue; }
						beg = r2->beg * KBM_BIN_SIZE; end = r1->end * KBM_BIN_SIZE;
						fprintf(out, "S\t%s\t", g->kbm->reads->buffer[rid].tag);
						fprintf(out, "-\t%d\t%d\t", (int)beg, (int)(end - beg));
						encap_string(str, end - beg);
						revseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[rid].seqoff * KBM_BIN_SIZE + beg, end - beg, str->string);
						fputs(str->string, out);
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

static inline u8i print_utgs_graph(Graph *g, char *prefix, char *utg, char *lay){
	FILE *o_seq, *o_lay, *files[4];
	tracev *path;
	String *seq;
	char *str;
	u8i i, uid, cnt, tot;
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
 * For debug in GDB
 * local_dot_node, local_dot_step, and print_local_dot_graph()
 */

static u8i local_dot_node = 1;
static u4i local_dot_step = 10;

static inline void get_subgraph_nodes_graph(Graph *g, ptrrefhash *nodes, u8v *stack, uint16_t max_step, u4i closed_val){
	node_t *n;
	edge_ref_t *f;
	edge_t *e;
	ptr_ref_t *p, *pp;
	u8i nid, idx;
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

static inline u8i print_local_dot_graph(Graph *g, char *prefix, char *suffix){
	FILE *out;
	ptrrefhash *hash;
	u8v *stack;
	ptr_ref_t *p;
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	u4i j, k, max;
	out = open_file_for_write(prefix, suffix, 1);
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
	fclose(out);
	return 0;
}

static inline u8i print_dot_full_graph(Graph *g, FILE *out){
	BufferedWriter *bw;
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	u4i j, k, max;
	bw = zopen_bufferedwriter(out, 1024 * 1024, 8, 0);
	beg_bufferedwriter(bw);
	fprintf(bw->out, "digraph {\n");
	fprintf(bw->out, "node [shape=record]\n");
	for(i=0;i<g->nodes->size;i++){
		if((i % 1000) == 0) flush_bufferedwriter(bw);
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
		fprintf(bw->out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	for(i=0;i<g->nodes->size;i++){
		if((i % 1000) == 0) flush_bufferedwriter(bw);
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
					fprintf(bw->out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1], e->closed? " style=dashed" : "");
				} else {
					fprintf(bw->out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2], e->closed? " style=dashed" : "");
				}
			}
		}
	}
	fprintf(bw->out, "}\n");
	end_bufferedwriter(bw);
	close_bufferedwriter(bw);
	return 0;
}

static inline u8i print_dot_graph(Graph *g, FILE *out){
	BufferedWriter *bw;
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t *f;
	edge_t *e;
	unsigned long long i, idx;
	u4i j, k, max;
	bw = zopen_bufferedwriter(out, 1024 * 1024, 8, 0);
	beg_bufferedwriter(bw);
	fprintf(bw->out, "digraph {\n");
	fprintf(bw->out, "node [shape=record]\n");
	for(i=0;i<g->nodes->size;i++){
		if((i % 1000) == 0) flush_bufferedwriter(bw);
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
		fprintf(bw->out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	for(i=0;i<g->nodes->size;i++){
		if((i % 1000) == 0) flush_bufferedwriter(bw);
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
					fprintf(bw->out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1]);
				} else {
					fprintf(bw->out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s]\n", i, (unsigned long long)e->node2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2]);
				}
			}
		}
	}
	fprintf(bw->out, "}\n");
	end_bufferedwriter(bw);
	close_bufferedwriter(bw);
	return 0;
}

static inline u8i print_nodes_graph(Graph *g, FILE *out){
	node_t *n;
	reg_t *r;
	unsigned long long i;
	u4i j;
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

static inline u8i print_reads_graph(Graph *g, FILE *out){
	read_t *rd;
	reg_t  *r;
	u8i idx;
	u4i i;
	for(i=0;i<g->kbm->reads->size;i++){
		rd = ref_readv(g->reads, i);
		fprintf(out, "%s\t%d\t%u", g->kbm->reads->buffer[i].tag, g->kbm->reads->buffer[i].bincnt * KBM_BIN_SIZE, rd->regs.cnt);
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

static inline u8i print_frgs_nodes_graph(Graph *g, FILE *out){
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

static inline u8i print_frgs_dot_graph(Graph *g, FILE *_out){
	BufferedWriter *bw;
	frg_t *frg;
	trace_t *t1, *t2;
	node_t *n1, *n2;
	reg_t *r1, *r2, *rr;
	edge_ref_t *f;
	lnk_t *e;
	unsigned long long i, idx;
	u4i j, k, max;
	bw = zopen_bufferedwriter(_out, 1024 * 1024, 8, 0);
	beg_bufferedwriter(bw);
	fprintf(bw->out, "digraph {\n");
	fprintf(bw->out, "node [shape=record]\n");
	for(i=0;i<g->frgs->size;i++){
		if((i % 1000) == 0) flush_bufferedwriter(bw);
		frg = ref_frgv(g->frgs, i);
		if(frg->closed){
			continue;
		}
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
		if(r1 == NULL){
			continue;
		}
		r2 = NULL; max = 0;
		for(j=0;j<n2->regs.cnt;j++){
			rr = ref_regv(g->regs, n2->regs.idx + j);
			if(g->reads->buffer[rr->rid].regs.cnt > max){
				r2 = rr;
				max = g->reads->buffer[rr->rid].regs.cnt;
			}
		}
		if(r2 == NULL){
			continue;
		}
		fprintf(bw->out, "F%llu [label=\"{F%llu %u %u/%u | { {N%llu:%c | %s | %c_%d_%d} | {N%llu:%c | %s | %c_%d_%d}}}\"]\n", i, i, frg->ty - frg->tx, frg->len, frg->length,
			t1->node, "+-"[t1->dir], g->kbm->reads->buffer[r1->rid].tag, "FR"[r1->dir], r1->beg, r1->end - r1->beg,
			t2->node, "+-"[t2->dir], g->kbm->reads->buffer[r2->rid].tag, "FR"[r2->dir], r2->beg, r2->end - r2->beg);
	}
	for(i=0;i<g->frgs->size;i++){
		if((i % 1000) == 0) flush_bufferedwriter(bw);
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
					fprintf(bw->out, "F%llu -> F%llu [label=\"%c%c:%d:%d\" color=%s style=%s]\n", i, (u8i)e->frg1, "+-"[k], "+-"[!e->dir1], e->cov, e->off, colors[k][!e->dir1], e->weak? "dashed" : "solid");
				} else {
					//if(g->frgs->buffer[e->frg2].ty - g->frgs->buffer[e->frg2].tx < (u4i)g->min_ctg_nds) continue;
					fprintf(bw->out, "F%llu -> F%llu [label=\"%c%c:%d:%d\" color=%s style=%s]\n", i, (u8i)e->frg2, "+-"[k], "+-"[e->dir2], e->cov, e->off, colors[k][e->dir2], e->weak? "dashed" : "solid");
				}
			}
		}
	}
	fprintf(bw->out, "}\n");
	end_bufferedwriter(bw);
	close_bufferedwriter(bw);
	return 0;
}

typedef u8i (*graph_print_func)(Graph *g, FILE *out);

static inline u8i generic_print_graph(Graph *g, graph_print_func func, char *prefix, char *suffix){
	FILE *out;
	char *file;
	u8i cnt;
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

#endif
