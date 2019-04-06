#ifndef __WTDBG_GRAPH_RJ_H
#define __WTDBG_GRAPH_RJ_H

#include "wtdbg.h"
#include "filewriter.h"
#include "pgzf.h"

static u8i debug_node = MAX_U8;

static inline void print_node_edges_cov_graph(Graph *g, FILE *out){
	node_t *n;
	edge_t *e;
	edge_ref_t eref;
	u8i nid;
	u4i k;
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->erefs[0].cnt + n->erefs[1].cnt == 0) continue;
		fprintf(out, "NODE_COV\tN%llu\t%u\t%u", nid, n->cov, n->regs.cnt);
		for(k=0;k<2;k++){
			beg_iter_edges_graph(n, k, &eref);
			while((e = ref_iter_edges_graph(g, &eref))){
				fprintf(out, "\t%u", get_edge_cov(e));
			}
		}
		fprintf(out, "\n");
	}
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
		fprintf(out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"%s]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg, p->idx == local_dot_node? " style=filled fillcolor=yellow" : "");
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
edge_ref_t eref;
edge_t *e;
u8i idx, nidx, hidx;
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
		beg_iter_edges_graph(n, k, &eref);
		while((e = ref_iter_edges_graph(g, &eref))){
			if(get_edge_cov(e) < g->min_edge_cov) continue;
			nidx = get_edge_didx(e, eref.flg);
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
	edge_ref_t eref;
	edge_t *e;
	n = ref_nodev(g->nodes, nid);
	beg_iter_edges_graph(n, dir, &eref);
	while((e = ref_iter_edges_graph(g, &eref))){
		fprintf(out, "N%llu\t%c\tN%llu\t%c\t%d\t%d\n", (u8i)get_edge_sidx(e, eref.flg), "+-"[get_edge_sdir(e, eref.flg)], (u8i)get_edge_didx(e, eref.flg), "+-"[get_edge_ddir(e, eref.flg)], get_edge_cov(e), get_edge_off(e));
	}
}

#define count_living_edges_graph(g, n, dir) (n)->erefs[dir].cnt

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

static inline int cal_offset_traces_graph(Graph *g, tracev *path, u8i beg, u8i end, int offset){
	trace_t *t;
	node_t *n;
	edge_t *e;
	u8i i;
	int off;
	off = offset;
	for(i=beg;i<end;i++){
		t = ref_tracev(path, i);
		t->off = off;
		n = ref_nodev(g->nodes, t->node);
		if(t->edges[t->dir].idx == 0){
			off += g->reglen;
		} else {
			e = ref_edgev(g->edges, t->edges[t->dir].idx);
			off += g->reglen + get_edge_off(e);
		}
	}
	return off;
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

static inline u4i gen_utgs_graph(Graph *g){
	utg_t *utg;
	tracev *path;
	trace_t *t;
	node_t *n, *v;
	edge_t *e;
	u8i nid;
	u4i uid, i;
	int len;
	uid = 0;
	clear_utgv(g->utgs);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		n->utg_idx  = MAX_UTG_IDX;
		n->bt_visit = 0;
		n->bt_idx = 0;
	}
	path = init_tracev(32);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		if(n->closed) continue;
		if(n->bt_visit) continue;
		n->bt_visit = 1;
		clear_tracev(path);
		t = next_ref_tracev(path);
		t->node = nid;
		t->edges[0] = t->edges[1] = EDGE_REF_NULL;
		t->dir = 0;
		exclusive_linear_trace_graph(g, t, MAX_B4, MAX_U4, 1, path, NULL);
		reverse_tracev(path);
		for(i=0;i<path->size;i++) path->buffer[i].dir = !path->buffer[i].dir;
		t = tail_tracev(path);
		exclusive_linear_trace_graph(g, t, MAX_B4, MAX_U4, 1, path, NULL);
		len = path->size * g->reglen;
		for(i=0;i<path->size;i++){
			t = ref_tracev(path, i);
			n = ref_nodev(g->nodes, t->node);
			n->utg_idx = uid;
			n->utg_dir = path->buffer[i].dir;
			n->bt_idx  = i;
			if(t->edges[t->dir].idx){
				e = ref_edgev(g->edges, t->edges[t->dir].idx);
				len += get_edge_off(e);
			}
		}
		if(len <= 0){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
		utg = next_ref_utgv(g->utgs);
		ZEROS(utg);
		utg->len = len;
		utg->cnt = path->size;
		t = tail_tracev(path);
		v = ref_nodev(g->nodes, t->node);
		utg->nodes[0].node = t->node;
		utg->nodes[0].dir  = t->dir;
		t = head_tracev(path);
		v = ref_nodev(g->nodes, t->node);
		utg->nodes[1].node = t->node;
		utg->nodes[1].dir  = !t->dir;
		uid ++;
	}
	free_tracev(path);
	return uid;
}

static inline u4i trace_utg_graph(Graph *g, u4i uid, tracev *path){
	utg_t *utg;
	trace_t *t;
	node_t *v, *w;
	edge_ref_t eref;
	edge_t *e;
	u4i off;
	utg = ref_utgv(g->utgs, uid);
	t = next_ref_tracev(path);
	t->node = utg->nodes[1].node;
	t->dir  = !utg->nodes[1].dir;
	t->edges[0] = EDGE_REF_NULL;
	t->edges[1] = EDGE_REF_NULL;
	t->off = off = 0;
	t->cov = 0;
	while(t->node != utg->nodes[0].node){
		v = ref_nodev(g->nodes, t->node);
		beg_iter_edges_graph(v, t->dir, &eref);
		while((e = ref_iter_edges_graph(g, &eref))){
			w = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
			if(w->utg_idx == uid){
				break;
			}
		}
		if(e == NULL){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		t->edges[t->dir] = (edge_ref_t){offset_edgev(g->edges, e), eref.flg, 0, 0};
		t = next_ref_tracev(path);
		t->node = get_edge_didx(e, eref.flg);
		t->dir  = get_edge_ddir(e, eref.flg);
		t->edges[!t->dir] = (edge_ref_t){offset_edgev(g->edges, e), !eref.flg, 0, 0};
		t->edges[t->dir] = EDGE_REF_NULL;
		t->cov = 0;
		off += g->reglen + get_edge_off(e);
		t->off = off;
	}
	return off + g->reglen;
}

static inline u4i revise_utg_graph(Graph *g, u4i uid, tracev *path){
	utg_t *u;
	trace_t *t;
	edge_t *e;
	node_t *v;
	u4i i, ret, cnt, len;
	ret = 0;
	u = ref_utgv(g->utgs, uid);
	cnt = len = 0;
	for(i=0;i<path->size;i++){
		t = ref_tracev(path, i);
		e = (t->edges[t->dir].idx == EDGE_REF_NULL.idx)? NULL : ref_edgev(g->edges, t->edges[t->dir].idx);
		if(u == NULL){
			u = next_ref_utgv(g->utgs);
			ret ++;
			ZEROS(u);
			cnt = 0;
			len = 0;
			u->nodes[1].node = t->node;
			u->nodes[1].dir  = !t->dir;
		}
		v = ref_nodev(g->nodes, t->node);
		v->utg_idx = offset_utgv(g->utgs, u);
		v->utg_dir = t->dir;
		cnt ++;
		len += g->reglen + (e? get_edge_off(e) : 0);
		if(e && get_edge_closed(e) == 0) continue;
		u->cnt = cnt;
		u->len = len;
		u->nodes[0].node = t->node;
		u->nodes[0].dir  = t->dir;
		u = NULL;
	}
	return ret;
}

static inline u8i print_utgs_nodes_graph(Graph *g, FILE *out){
	utg_t *utg;
	tracev *path;
	trace_t *t;
	u4i uid, len, i;
	path = init_tracev(32);
	for(uid=0;uid<g->utgs->size;uid++){
		utg = ref_utgv(g->utgs, uid);
		clear_tracev(path);
		len = trace_utg_graph(g, uid, path) * KBM_BIN_SIZE;
		fprintf(out, "UTG%u\t%u\t%u", uid, utg->cnt, len);
		for(i=0;i<path->size;i++){
			t = ref_tracev(path, i);
			fprintf(out, "\tN%llu:%c:%u", t->node, "+-"[t->dir], t->off);
		}
		fprintf(out, "\n");
	}
	free_tracev(path);
	return g->utgs->size;
}

static inline u4i count_utg_edges_graph(Graph *g, u4i uid, int dir){
	utg_t *utg;
	node_t *v;
	utg = ref_utgv(g->utgs, uid);
	v = ref_nodev(g->nodes, utg->nodes[dir].node);
	return v->erefs[utg->nodes[dir].dir].cnt;
}

static inline void filter_utgs_graph(Graph *g){
	utg_t *utg;
	u4v *lens;
	u8i i;
	lens = init_u4v(32);
	for(i=0;i<g->utgs->size;i++){
		utg = ref_utgv(g->utgs, i);
		if(utg->closed) continue;
		if(utg->len * KBM_BIN_SIZE < g->min_ctg_len || utg->cnt < g->min_ctg_nds){
			utg->closed = 1;
		} else {
			push_u4v(lens, utg->len * KBM_BIN_SIZE);
		}
	}
	fprintf(KBM_LOGF, "[%s] Estimated: ", date()); num_n50(lens, KBM_LOGF); fprintf(KBM_LOGF, "\n");
	free_u4v(lens);
}

/*
static inline void update_utg_graph(Graph *g, u4i uid, u8i visit){
	utg_t *utg;
	node_t *v;
	edge_ref_t eref;
	edge_t *e;
	u4i k, dir;
	utg = ref_utgv(g->utgs, uid);
	for(k=0;k<2;k++){
		v = ref_nodev(g->nodes, utg->nodes[k].node);
		dir = utg->nodes[k].dir;
		if(v->erefs[dir].cnt == utg->nodes[k].edge){
			continue;
		}
		utg->nodes[k].edge = v->erefs[dir].cnt;
		v->bt_visit = visit;
		while(1){
			if(v->erefs[dir].cnt != 1){
				break;
			}
			beg_iter_edges_graph(v, dir, &eref);
			e = ref_iter_edges_graph(g, &eref);
			v = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
			if(v->erefs[!get_edge_ddir(e, eref.flg)].cnt > 1){
				break;
			}
			if(v->bt_visit == visit){
				break;
			} else {
				v->bt_visit = visit;
			}
			if(v->utg_idx == uid){ // maybe a loop
				break;
			}
			utg->cnt ++;
			utg->len += g->reglen + get_edge_off(e);
			utg->nodes[k].node = get_edge_didx(e, eref.flg);
			utg->nodes[k].dir  = dir = get_edge_ddir(e, eref.flg);
			utg->nodes[k].edge = v->erefs[dir].cnt;
			g->utgs->buffer[v->utg_idx].closed = 1;
			v->utg_idx = uid;
			v->utg_dir = dir ^ k;
		}
	}
}
*/

static inline void readpath_highlight_nodes_graph(Graph *g, u8i nid, u8i visit){
	node_t *nd, *nd2;
	read_t *rd;
	reg_t  *rg, *rg2;
	u8i ri, idx;
	nd2 = ref_nodev(g->nodes, nid);
	for(ri=0;ri<nd2->regs.cnt;ri++){
		rg2 = ref_regv(g->regs, nd2->regs.idx + ri);
		rd  = ref_readv(g->reads, rg2->rid);
		idx = rd->regs.idx;
		while(idx){
			rg = ref_regv(g->regs, idx);
			idx = rg->read_link;
			if(rg->closed) continue;
			nd = ref_nodev(g->nodes, rg->node);
			if(nd->closed) continue;
			nd->bt_visit = visit;
		}
	}
}

static inline u4i readpath_transitive_reduction_graph_core(Graph *g, u8i nid, BitVec *visits){
	node_t *n, *w;
	edge_ref_t eref;
	edge_t *e;
	reg_t  *r1, *re1, *r2, *re2;
	u4i k, cnt, vst, ret;
	n = ref_nodev(g->nodes, nid);
	if(n->closed) return 0;
	if(n->erefs[0].cnt < 2 && n->erefs[1].cnt < 2) return 0;
	clear_bitvec(visits);
	encap_bitvec(visits, n->regs.cnt);
	ret = 0;
	for(k=0;k<2;k++){
		if(n->erefs[k].cnt < 2) continue;
		reg_zeros_bitvec(visits, 0, n->regs.cnt);
		beg_iter_edges_graph(n, k, &eref);
		while((e = ref_iter_edges_graph(g, &eref))){
			w = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
			if(w->regs.cnt == 0) continue;
			r1  = ref_regv(g->regs, n->regs.idx);
			re1 = ref_regv(g->regs, n->regs.idx + n->regs.cnt);
			r2  = ref_regv(g->regs, w->regs.idx);
			re2 = ref_regv(g->regs, w->regs.idx + w->regs.cnt);
			cnt = vst = 0;
			while(1){
				if(r1->rid > r2->rid){
					r2 ++;
					if(r2 > re2) break;
				} else if(r1->rid < r2->rid){
					r1 ++;
					if(r1 > re1) break;
				} else {
					if(get_edge_sdir(e, eref.flg) == ((r1->beg > r2->beg) ^ r1->dir) && get_edge_ddir(e, eref.flg) == ((r1->beg > r2->beg) ^ r2->dir)){
						cnt ++;
						if(get_bitvec(visits, offset_regv(g->regs, r1) - n->regs.idx)){
							vst ++;
						} else {
							one_bitvec(visits, offset_regv(g->regs, r1) - n->regs.idx);
						}
					}
					r1 ++;
					if(r1 > re1) break;
					r2 ++;
					if(r2 > re2) break;
				}
			}
			if(cnt >= get_edge_cov(e) + g->min_edge_cov){ // there must inter nodes
				if(vst >= g->min_edge_cov){ // sum(inter_edge_cov) is ok, but not sure any edge has enough cov. However, the edges existed
					if(get_edge_keep(e) == 0){ // not the shortest edges for each node
						cut_edge_graph(g, e);
						ret ++;
					}
				}
			}
		}
	}
	return ret;
}

static inline u8i readpath_transitive_reduction_graph(Graph *g){
	BitVec *visits;
	node_t *n;
	edge_ref_t eref;
	edge_t *e;
	u8i idx, ret;
	int k, elen;
	ret = 0;
	visits = init_bitvec(1024);
	for(idx=1;idx<g->edges->size;idx++){
		set_edge_keep(ref_edgev(g->edges, idx), 0);
	}
	for(idx=0;idx<g->nodes->size;idx++){
		n = ref_nodev(g->nodes, idx);
		for(k=0;k<2;k++){
			srt_node_edges_graph(g, idx, k);
			beg_iter_edges_graph(n, k, &eref);
			elen = MAX_B4;
			while((e = ref_iter_edges_graph(g, &eref))){
				if(get_edge_off(e) <= elen){
					elen = get_edge_off(e);
					set_edge_keep(e, 1);
				} else {
					break;
				}
			}
		}
	}
	for(idx=0;idx<g->nodes->size;idx++){
		ret += readpath_transitive_reduction_graph_core(g, idx, visits);
	}
	free_bitvec(visits);
	return ret;
}

typedef struct {
	u4i rid:30, dir:1, sel:1;
	u2i eidx[2];
} preg_t;
define_list(pregv, preg_t);

static inline u4i count_nonlinear_traces_graph(Graph *g, u8i _nid, int dir, u4i max_step, u8i visit, u8v *stack){
	node_t *v, *w;
	edge_ref_t eref;
	edge_t *e;
	u8i nid;
	u4i mstep;
	v = ref_nodev(g->nodes, _nid);
	v->bt_visit = visit;
	v->bt_dir   = dir;
	v->bt_idx   = mstep = 1;
	clear_u8v(stack);
	push_u8v(stack, _nid);
	while(mstep < max_step && pop_u8v(stack, &nid)){
		v = ref_nodev(g->nodes, nid);
		beg_iter_edges_graph(v, v->bt_dir, &eref);
		while((e = ref_iter_edges_graph(g, &eref))){
			w = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
			if(w->bt_visit == visit) continue;
			w->bt_visit = visit;
			w->bt_dir   = get_edge_ddir(e, eref.flg);
			w->bt_idx   = v->bt_idx + 1;
			if(w->bt_idx >= mstep){
				mstep = w->bt_idx;
			}
			push_u8v(stack, offset_nodev(g->nodes, w));
		}
	}
	return mstep;
}

static inline u4i split_utg_graph(Graph *g, u4i uid, edge_ref_t *erefs[2], tracev *path, BitVec *rdbits){
	utg_t *utg, *utg2;
	trace_t *t;
	node_t *v, *w;
	edge_ref_t eref;
	edge_t *e, *f;
	read_t *rd;
	reg_t *r, *s;
	u8i idx, ndoff;
	u4i k, ti, ri, uid2, eoff, ecov, cuts[2], i, j;
	encap_utgv(g->utgs, 1);
	utg  = ref_utgv(g->utgs, uid);
	// check whether need to split
	if(count_utg_edges_graph(g, uid, 0) <= 1 && count_utg_edges_graph(g, uid, 1) <= 1){
		return uid;
	}
	uid2 = g->utgs->size;
	utg2 = next_ref_utgv(g->utgs);
	*utg2 = *utg;
	encap_nodev(g->nodes, path->size);
	encap_edgev(g->edges, path->size - 1);
	ndoff = g->nodes->size;
	cuts[0] = cuts[1] = 0;
	eoff = 0;
	for(ti=0;ti<path->size;ti++){
		t = ref_tracev(path, ti);
		v = ref_nodev(g->nodes, t->node);
		w = next_ref_nodev(g->nodes);
		ZEROS(w);
		w->utg_idx = uid2;
		w->utg_dir = v->utg_dir;
		w->regs = v->regs;
		v->bt_idx = eoff; // for recovery
		w->bt_idx = eoff;
		e = ref_edgev(g->edges, t->edges[t->dir].idx);
		eoff += g->reglen + get_edge_off(e);
		for(ri=0;ri<v->regs.cnt;ri++){
			r = ref_regv(g->regs, v->regs.idx + ri);
			if(get_bitvec(rdbits, r->rid)){
				r->node = offset_nodev(g->nodes, w);
			}
		}
		sort_array(g->regs->buffer + v->regs.idx, v->regs.cnt, reg_t, num_cmpgtx(a.node, b.node, a.rid, b.rid));
		for(v->regs.cnt=0;v->regs.cnt<w->regs.cnt;v->regs.cnt++){
			r = ref_regv(g->regs, v->regs.idx + v->regs.cnt);
			if(r->node != t->node) break;
		}
		v->cov = v->regs.cnt;
		w->cov = w->regs.cnt - v->regs.cnt;
		// revise read_link
		for(ri=0;ri<w->regs.cnt;ri++){
			r = ref_regv(g->regs, w->regs.idx + ri);
			rd = ref_readv(g->reads, r->rid);
			idx = rd->regs.idx;
			s = NULL;
			while(idx){
				if(idx >= w->regs.idx && idx < w->regs.idx + w->regs.cnt){
					if(s){
						s->read_link = w->regs.idx + ri;
					} else {
						rd->regs.idx = w->regs.idx + ri;
					}
					break;
				}
				s = ref_regv(g->regs, idx);
				idx = s->read_link;
			}
		}
		w->regs.idx = v->regs.idx + v->regs.cnt;
		w->regs.cnt -= v->regs.cnt;
		if(ti){ // add edge between new vertices
			e = ref_edgev(g->edges, t->edges[!t->dir].idx);
			f = next_ref_edgev(g->edges);
			*f = *e;
			set_edge_sidx(f, t->edges[!t->dir].flg, ndoff + ti);
			set_edge_didx(f, t->edges[!t->dir].flg, ndoff + ti - 1);
			ecov = cal_edge_cov_graph(g, e);
			set_edge_cov(e, ecov);
			if(ecov == 0 || (ecov < g->min_edge_cov && !get_edge_hard(e))){
				cut_edge_graph(g, e);
				cuts[0] ++;
			}
			ecov = cal_edge_cov_graph(g, f);
			set_edge_cov(f, ecov);
			add_edge_graph(g, f);
			if(ecov == 0 || (ecov < g->min_edge_cov && !get_edge_hard(f))){
				cut_edge_graph(g, f);
				cuts[1] ++;
			}
		}
	}
	// recovery splited utgs
	for(k=0;k<2;k++){
		if(cuts[k] == 0) continue;
		for(i=0;i+2<path->size;i++){
			v = ref_nodev(g->nodes, k? ndoff + i : path->buffer[i].node);
			for(j=i+2;j<path->size;j++){
				w = ref_nodev(g->nodes, k? ndoff + j : path->buffer[j].node);
				encap_edgev(g->edges, 1);
				e = edge_node2node_graph(g, offset_nodev(g->nodes, v), v->utg_dir, offset_nodev(g->nodes, w), w->utg_dir, NULL);
				if(e){
					if(is_edge_closed(e)){
						ecov = cal_edge_cov_graph(g, e);
						if(ecov >= g->min_edge_cov || (ecov && get_edge_hard(e))){
							revive_edge_graph(g, e);
						}
					}
				} else {
					e = head_edgev(g->edges);
					set_edge_sidx(e, 0, offset_nodev(g->nodes, v));
					set_edge_didx(e, 0, offset_nodev(g->nodes, w));
					set_edge_sdir(e, 0, v->utg_dir);
					set_edge_ddir(e, 0, w->utg_dir);
					ecov = cal_edge_cov_graph(g, e);
					set_edge_cov(e, ecov);
					set_edge_off(e, w->bt_idx - v->bt_idx - g->reglen);
					if(ecov >= g->min_edge_cov){
						f = next_ref_edgev(g->edges);
						*f = *e;
						set_edge_closed(f, 0);
						add_edge_graph(g, f);
					}
				}
			}
		}
		if(k){
			utg2->closed = 1;
		} else {
			utg->closed = 1;
		}
	}
	utg2->nodes[1].node = ndoff;
	utg2->nodes[0].node = ndoff + path->size - 1;
	// add edges to new utg's end
	for(k=0;k<2;k++){
		if(erefs[k] == NULL || erefs[k]->idx == 0) continue;
		encap_edgev(g->edges, 1);
		e = ref_edgev(g->edges, erefs[k]->idx);
		f = next_ref_edgev(g->edges);
		*f = *e;
		set_edge_sidx(f, erefs[k]->flg, utg2->nodes[k].node);
		add_edge_graph(g, f);
		cut_edge_graph(g, e);
		ecov = cal_edge_cov_graph(g, f);
		set_edge_cov(f, ecov);
		if(ecov == 0 || (ecov < g->min_edge_cov && !get_edge_hard(f))){
			cut_edge_graph(g, f);
		}
	}
	// check edge covs at old utg's end
	for(k=0;k<2;k++){
		t = ref_tracev(path, k? 0 : path->size - 1);
		beg_iter_edges_graph(ref_nodev(g->nodes, t->node), t->dir ^ k, &eref);
		while((e = ref_iter_edges_graph(g, &eref))){
			ecov = cal_edge_cov_graph(g, e);
			set_edge_cov(e, ecov);
			if(ecov == 0 || (ecov < g->min_edge_cov && !get_edge_hard(e))){
				cut_edge_graph(g, e);
			}
		}
	}
	return uid2;
}

static inline u4i _readpath_split_utg_graph(Graph *g, u4i uid, u4i eidxs[2], u4i *cnts, u8v *matrix, tracev *path, BitVec *rdbits, pregv *pregs){
	preg_t *pr;
	edge_ref_t *erefs[2], FWD_EF, REV_EF;
	u4i i, ret;
	if(g->utgs->buffer[uid].closed) return 0;
	for(i=0;i<pregs->size;i++){
		pr = ref_pregv(pregs, i);
		if((pr->sel = (pr->eidx[0] == eidxs[0] && pr->eidx[1] == eidxs[1]))){
			one_bitvec(rdbits, pr->rid);
		}
	}
	if(eidxs[0] == cnts[0]){
		FWD_EF = EDGE_REF_NULL;
	} else {
		FWD_EF = (edge_ref_t){matrix->buffer[eidxs[0]] >> 1, matrix->buffer[eidxs[0]] & 0x01, 0, 0};
	}
	if(eidxs[1] == cnts[1]){
		REV_EF = EDGE_REF_NULL;
	} else {
		REV_EF = (edge_ref_t){matrix->buffer[eidxs[1] + cnts[0]] >> 1, matrix->buffer[eidxs[1] + cnts[0]] & 0x01, 0, 0};
	}
	erefs[0] = &FWD_EF;
	erefs[1] = &REV_EF;
	ret = (split_utg_graph(g, uid, erefs, path, rdbits) != uid);
	for(i=0;i<pregs->size;i++){
		pr = ref_pregv(pregs, i);
		if(pr->sel){
			zero_bitvec(rdbits, pr->rid);
			pr->sel = 0;
		}
	}
	return ret;
}

static inline u8i get_utg_readpaths_graph(Graph *g, u4i uid, tracev *path, int ends, u4i len2end, BitVec *rdbits, pregv *pregs){
	trace_t *t;
	node_t *v;
	reg_t *r;
	preg_t *pr;
	u4i i, j;
	UNUSED(uid);
	//clear_tracev(path);
	//trace_utg_graph(g, uid, path);
	clear_pregv(pregs);
	for(i=0;i<path->size;i++){
		if(!((i < len2end && (ends & 2)) || (i + len2end >= path->size && (ends & 1)))) continue;
		t = ref_tracev(path, i);
		v = ref_nodev(g->nodes, t->node);
		for(j=0;j<v->regs.cnt;j++){
			r = ref_regv(g->regs, v->regs.idx + j);
			if(r->closed) continue;
			if(get_bitvec(rdbits, r->rid)) continue;
			one_bitvec(rdbits, r->rid);
			pr = next_ref_pregv(pregs);
			pr->rid = r->rid;
			pr->dir = v->utg_dir ^ r->dir;
			pr->sel = 0;
			pr->eidx[0] = MAX_U2;
			pr->eidx[1] = MAX_U2;
		}
	}
	return pregs->size;
}

static inline u4i readpath_utg_rescue_edges_graph_core(Graph *g, u4i uid, int ends, u4i mercy_edge_cov, tracev *path, BitVec *rdbits, pregv *pregs, edgev *edges, edgehash *ehash, edgeoffv *offs, u8v *rgidxs, UUhash *hash){
	preg_t *pr;
	read_t *rd;
	reg_t *r1, *r2;
	node_t *n1, *n2;
	edge_t *E, *e, *f;
	UUhash_t *U;
	u8i idx, lst, *u;
	u4i pridx, i, j, k, cnt;
	int exists, eoff;
	clear_tracev(path);
	trace_utg_graph(g, uid, path);
#if DEBUG
	int print_rst;
	print_rst = 0;
	if(g->utgs->buffer[uid].nodes[0].node == debug_node || g->utgs->buffer[uid].nodes[1].node == debug_node){
		print_rst = 1;
		print_local_dot_graph(g, path->buffer[0].node, 10, "1.dot");
	}
#endif
	get_utg_readpaths_graph(g, uid, path, ends, g->tip_step, rdbits, pregs);
	clear_edgev(edges);
	E = next_ref_edgev(edges);
	ZEROS(E);
	set_edge_cov(E, 1);
	set_edge_closed(E, WT_EDGE_CLOSED_MASK);
	clear_edgehash(ehash);
	set_userdata_edgehash(ehash, edges);
	clear_edgeoffv(offs);
	for(pridx=0;pridx<pregs->size;pridx++){
		pr = ref_pregv(pregs, pridx);
		zero_bitvec(rdbits, pr->rid);
		rd = ref_readv(g->reads, pr->rid);
		idx = rd->regs.idx;
		clear_u8v(rgidxs);
		while(idx){
			r1 = ref_regv(g->regs, idx);
			if(r1->closed == 0 && g->nodes->buffer[r1->node].closed == 0 && g->nodes->buffer[r1->node].init_end){
				push_u8v(rgidxs, idx);
			}
			idx = r1->read_link;
		}
		if(print_rst && (strcmp(g->kbm->reads->buffer[pr->rid].tag, "S0000216356") == 0 || strcmp(g->kbm->reads->buffer[pr->rid].tag, "S0000114166") == 0)){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
		for(i=0;i+1<rgidxs->size;i++){
			r1 = ref_regv(g->regs, rgidxs->buffer[i]);
			n1 = ref_nodev(g->nodes, r1->node);
			for(j=i+1;j<rgidxs->size;j++){
				r2 = ref_regv(g->regs, rgidxs->buffer[j]);
				n2 = ref_nodev(g->nodes, r2->node);
				if(n1->utg_idx == n2->utg_idx) continue;
				encap_edgev(edges, 1);
				E = head_edgev(edges);
				if(print_rst && ((r1->node == 11739 && r2->node == 11394) || (r2->node == 11739 && r1->node == 11394))){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				}
				k = r1->node > r2->node;
				set_edge_sidx(E, k, r1->node);
				set_edge_didx(E, k, r2->node);
				set_edge_sdir(E, k, r1->dir);
				set_edge_ddir(E, k, r2->dir);
				eoff = Int(r2->beg) - Int(r1->end);
				set_edge_off(E, eoff);
				u = prepare_edgehash(ehash, 0, &exists);
				if(exists){
					e = ref_edgev(edges, *u);
					if(get_edge_cov(e) < WT_MAX_EDGE_COV) set_edge_cov(e, get_edge_cov(e) + 1);
				} else {
					*u = edges->size;
					e = next_ref_edgev(edges);
					*e = *E;
				}
				push_edgeoffv(offs, (edge_off_t){*u, eoff});
			}
		}
	}
	sort_array(offs->buffer, offs->size, edge_off_t, num_cmpgtx(a.idx, b.idx, a.off, b.off));
	lst = 0;
	for(idx=1;idx<=offs->size;idx++){
		if(idx < offs->size && offs->buffer[idx].idx == offs->buffer[lst].idx) continue;
		set_edge_off(ref_edgev(edges, offs->buffer[lst].idx), offs->buffer[(lst+idx)/2].off);
		set_edge_cov(ref_edgev(edges, offs->buffer[lst].idx), idx - lst);
		lst = idx;
	}
	// check whether out-links point to one other utg
	clear_UUhash(hash);
	for(idx=1;idx<edges->size;idx++){
		u8i key, dis;
		e = ref_edgev(edges, idx);
		if(get_edge_cov(e) < mercy_edge_cov) continue;
		f = edge_node2node_graph(g, get_edge_sidx(e, 0), get_edge_sdir(e, 0), get_edge_didx(e, 0), get_edge_ddir(e, 0), NULL);
		if(f) continue;
		k = (ref_nodev(g->nodes, get_edge_sidx(e, 1))->utg_idx == uid);
		n1 = ref_nodev(g->nodes, get_edge_sidx(e, k));
		if(((n1->init_end >> (get_edge_sdir(e, k) ^ n1->utg_dir)) & 0x01) == 0){
			continue;
		}
		n2 = ref_nodev(g->nodes, get_edge_didx(e, k));
		if(((n2->init_end >> ((!get_edge_ddir(e, k)) ^ n2->utg_dir)) & 0x01) == 0){
			continue;
		}
		dis  = (get_edge_sdir(e, k) ^ n1->utg_dir)? g->utgs->buffer[n1->utg_idx].cnt - 1 - n1->bt_idx : n1->bt_idx;
		dis += (get_edge_ddir(e, k) ^ n2->utg_dir)? g->utgs->buffer[n2->utg_idx].cnt - 1 - n2->bt_idx : n2->bt_idx;
		key = (((u8i)n2->utg_idx) << 2) | ((get_edge_ddir(e, k) ^ n2->utg_dir) << 1) | (get_edge_sdir(e, k) ^ n1->utg_dir);
		U = prepare_UUhash(hash, key, &exists);
		if(exists){
			if((U->val >> 10) > dis){
				U->val = (dis << 10) | idx;
			}
		} else {
			U->key = key;
			U->val = (dis << 10) | idx;
		}
	}
	cnt = 0;
	reset_iter_UUhash(hash);
	while((U = ref_iter_UUhash(hash))){
		e = ref_edgev(edges, U->val & 0x3FF);
		set_edge_closed(e, 0);
		set_edge_hard(e, 1); // this edge won't be cut when be detected to be less than g->min_edge_cov
		f = next_ref_edgev(g->edges);
		*f = *e;
		add_edge_graph(g, f);
		cnt ++;
#if DEBUG
		if(print_rst){
			fprintf(stdout, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", (u8i)get_edge_sidx(e, 0), (u8i)get_edge_didx(e, 0),
				"+-"[get_edge_sdir(e, 0)], "+-"[get_edge_ddir(e, 0)], get_edge_cov(e), get_edge_off(e), colors[get_edge_sdir(e, 0)][get_edge_ddir(e, 0)], is_edge_closed(e)? " style=dashed" : "");
		}
#endif
	}
	return cnt;
}

static inline u8i readpath_utg_rescue_edges_graph(Graph *g){
	BitVec *rdbits;
	tracev *path;
	pregv *pregs;
	edgev *edges;
	edgehash *ehash;
	edgeoffv *offs;
	u8v *rgidxs;
	UUhash *hash;
	u8i ret, idx;
	u4i uid, cnt, i, k, mercy_edge_cov;
	ret = 0;
	mercy_edge_cov = num_min(g->min_edge_cov, 1);
	rdbits = init_bitvec(g->reads->size);
	path = init_tracev(32);
	pregs = init_pregv(32);
	edges = init_edgev(32);
	ehash = init_edgehash(13);
	offs = init_edgeoffv(32);
	rgidxs = init_u8v(8);
	hash = init_UUhash(13);
	gen_utgs_graph(g);
	for(idx=0;idx<g->nodes->size;idx++){
		g->nodes->buffer[idx].init_end = 0;
	}
	for(uid=0;uid<g->utgs->size;uid++){
		if(g->utgs->buffer[uid].cnt < 2 * g->tip_step) continue;
		k = (count_utg_edges_graph(g, uid, 0)? 0 : 1) | (count_utg_edges_graph(g, uid, 1)? 0 : 2);
		if(k == 0) continue;
		clear_tracev(path);
		trace_utg_graph(g, uid, path);
		if(k & 0x02){
			for(i=0;i<g->tip_step;i++){
				g->nodes->buffer[path->buffer[i].node].init_end = 2;
			}
		}
		if(k & 0x01){
			for(i=0;i<g->tip_step;i++){
				g->nodes->buffer[path->buffer[path->size - 1 - i].node].init_end = 1;
			}
		}
	}
	for(uid=0;uid<g->utgs->size;uid++){
		if(g->utgs->buffer[uid].cnt < 2 * g->tip_step) continue;
		k = (count_utg_edges_graph(g, uid, 0)? 0 : 1) | (count_utg_edges_graph(g, uid, 1)? 0 : 2);
		if(k == 0) continue;
		cnt = readpath_utg_rescue_edges_graph_core(g, uid, k, mercy_edge_cov, path, rdbits, pregs, edges, ehash, offs, rgidxs, hash);
		ret += cnt;
	}
	free_bitvec(rdbits);
	free_tracev(path);
	free_edgev(edges);
	free_edgehash(ehash);
	free_edgeoffv(offs);
	free_u8v(rgidxs);
	free_UUhash(hash);
	return ret;
}

static inline u4i readpath_detach_branches_graph_core(Graph *g, u4i uid, tracev *path, BitVec *rdbits, pregv *pregs, uuhash *nhash, u8v *matrix, u8v *rgidxs, u8i *visit, u8v *heap){
	utg_t *utg, *wtg;
	trace_t *t;
	node_t *v, *w;
	edge_ref_t eref;
	edge_t *e, *f;
	reg_t *r, *r1, *r2;
	preg_t *pr;
	read_t *rd;
	uuhash_t *nh;
	u8i *mtx;
	u8i idx;
	u4i i, j, beg, end, dir, k, repidx, flg, rid, nid, len, ret, cnts[3];
	u4i m, n, sum, eidxs[2];
	int exists;
#if DEBUG
	int print_rst;
	print_rst = 0;
	if(g->utgs->buffer[uid].nodes[0].node == debug_node || g->utgs->buffer[uid].nodes[1].node == debug_node){
		print_rst = 1;
		print_local_dot_graph(g, path->buffer[0].node, 10, "1.dot");
	}
#endif
	clear_uuhash(nhash);
	clear_pregv(pregs);
	ret = 0;
	clear_u8v(matrix);
	for(k=0;k<2;k++){
		if(k == 0){
			t = tail_tracev(path);
			v = ref_nodev(g->nodes, t->node);
			beg_iter_edges_graph(v, t->dir, &eref);
		} else {
			t = head_tracev(path);
			v = ref_nodev(g->nodes, t->node);
			beg_iter_edges_graph(v, !t->dir, &eref);
		}
		while((e = ref_iter_edges_graph(g, &eref))){
			w = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
			nh = prepare_uuhash(nhash, w->utg_idx, &exists);
			if(exists){
			} else {
				nh->key = w->utg_idx;
				nh->val = ((eref.cnt - 1) << 1) | k;
				push_u8v(matrix, (offset_edgev(g->edges, e) << 1) | eref.flg);
			}
		}
		cnts[k] = eref.cnt;
	}
	if(cnts[0] == 0 && cnts[1] == 0) return 0;
	cnts[2] = cnts[0] + 2 + cnts[1] + 2;
	encap_and_zeros_u8v(matrix, cnts[2] * cnts[2]);
	mtx = matrix->buffer + cnts[0] + cnts[1];
	get_utg_readpaths_graph(g, uid, path, 3, 2 * g->tip_step, rdbits, pregs);
	for(i=0;i<pregs->size;i++){
		pr = ref_pregv(pregs, i);
		zero_bitvec(rdbits, pr->rid);
		rid = pr->rid;
		rd = ref_readv(g->reads, rid);
		r1 = NULL;
		idx = rd->regs.idx;
		clear_u8v(rgidxs);
		flg = MAX_U4;
		n = 0;
		while(idx){
			r2 = ref_regv(g->regs, idx);
			idx = r2->read_link;
			if(r2->closed) continue;
			if(ref_nodev(g->nodes, r2->node)->closed) continue;
			repidx = ref_nodev(g->nodes, r2->node)->utg_idx;
			if(repidx == uid){
				n ++;
			}
			if(r1 && repidx == ref_nodev(g->nodes, r1->node)->utg_idx){
				continue;
			}
			r1 = r2;
			if(repidx == uid){
				flg = rgidxs->size;
			} else if((nh = get_uuhash(nhash, repidx))){
				pr->eidx[nh->val & 0x01] = nh->val >> 1;
			}
			push_u8v(rgidxs, offset_regv(g->regs, r2));
		}
		if(flg == MAX_U4){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			continue;
		}
		for(k=0;k<2;k++){
			if(pr->eidx[k] != MAX_U2) continue;
			dir = pr->dir ^ k;
			if(dir){
				if(flg){
					j = flg - 1;
				} else {
					pr->eidx[k] = cnts[k];
					continue;
				}
			} else {
				if(flg + 1 == rgidxs->size){
					pr->eidx[k] = cnts[k];
					continue;
				} else {
					j = flg + 1;
				}
			}
			r = ref_regv(g->regs, rgidxs->buffer[j]);
			len = num_diff(r->beg, g->regs->buffer[rgidxs->buffer[flg]].beg);
			len = num_min((len + 1) * 2, len + 8);
			dir = dir ^ r->dir;
			// trace to a node existing in nhash
			v = ref_nodev(g->nodes, r->node);
			(*visit) ++;
			v->bt_visit = *visit;
			v->unvisit = dir;
			v->bt_idx = 0; // distance
			clear_u8v(heap);
			push_u8v(heap, r->node);
			while(heap->size){
				nid = heap->buffer[0];
				array_heap_remove(heap->buffer, heap->size, heap->cap, u8i, 0, num_cmp(g->nodes->buffer[a].bt_idx, g->nodes->buffer[b].bt_idx));
				v = ref_nodev(g->nodes, nid);
				if((nh = get_uuhash(nhash, v->utg_idx))){
					if((nh->val & 0x01) == k){
						if(pr->eidx[k] == MAX_U2){
							pr->eidx[k] = nh->val >> 1;
						} else if(pr->eidx[k] != (nh->val >> 1)){
							pr->eidx[k] = MAX_U2;
							break;
						}
					}
				}
				beg_iter_edges_graph(v, v->unvisit, &eref);
				while((e = ref_iter_edges_graph(g, &eref))){
					w = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
					if(w->bt_visit == (*visit)) continue;
					w->bt_visit = (*visit);
					w->unvisit  = get_edge_ddir(e, eref.flg);
					w->bt_idx   = v->bt_idx + g->reglen + get_edge_off(e);
					if(w->bt_idx > len) continue;
					array_heap_push(heap->buffer, heap->size, heap->cap, u8i, get_edge_didx(e, eref.flg), num_cmp(g->nodes->buffer[a].bt_idx, g->nodes->buffer[b].bt_idx));
				}
			}
			if(pr->eidx[k] == MAX_U2){
				pr->eidx[k] = cnts[k] + 1; // disconnect subgraph
			}
		}
		if(pr->eidx[0] == cnts[0] && pr->eidx[1] == cnts[1]){
			if(n < 2){ // not a real inner read path
				pr->eidx[0] = cnts[0] + 1;
				pr->eidx[1] = cnts[1] + 1;
			}
		}
		mtx[pr->eidx[0] * cnts[2] + (pr->eidx[1] + cnts[0] + 2)] ++;
		mtx[(pr->eidx[1] + cnts[0] + 2) * cnts[2] + pr->eidx[0]] ++;
	}
#if DEBUG
	if(print_rst){
		u8v *bs[2];
		printf("PATH:");
		for(i=0;i<path->size;i++){
			t = ref_tracev(path, i);
			printf(" N%llu:%c", t->node, "+-"[t->dir]);
		}
		printf("\n");
		bs[0] = init_u8v(4);
		bs[1] = init_u8v(4);
		for(k=0;k<2;k++){
			if(k == 0){
				t = tail_tracev(path);
				v = ref_nodev(g->nodes, t->node);
				beg_iter_edges_graph(v, t->dir, &eref);
			} else {
				t = head_tracev(path);
				v = ref_nodev(g->nodes, t->node);
				beg_iter_edges_graph(v, !t->dir, &eref);
			}
			while((e = ref_iter_edges_graph(g, &eref))){
				idx = get_edge_didx(e, eref.flg);
				push_u8v(bs[k], idx);
				v = ref_nodev(g->nodes, idx);
				utg = ref_utgv(g->utgs, v->utg_idx);
				printf("EDGE[%d][%d] N%llu UTG%llu:%c:%d:[%d, %d]\n", k, eref.cnt - 1, idx, (u8i)v->utg_idx, "+-"[get_edge_ddir(e, eref.flg) ^ v->utg_dir], utg->cnt, count_utg_edges_graph(g, v->utg_idx, 0), count_utg_edges_graph(g, v->utg_idx, 1));
			}
		}
		printf("MATRIX");
		for(j=0;j<cnts[1]+2;j++){
			if(j < cnts[1]){
				printf("\tN%llu", bs[1]->buffer[j]);
			} else if(j == cnts[1]){
				printf("\tINNER");
			} else {
				printf("\tUNKNOWN");
			}
		}
		printf("\n");
		for(i=0;i<cnts[0]+2;i++){
			if(i < cnts[0]){
				printf("N%llu", bs[0]->buffer[i]);
			} else if(i == cnts[0]){
				printf("INNER");
			} else {
				printf("UNKNOWN");
			}
			for(j=0;j<cnts[1]+2;j++){
				printf("\t%llu", mtx[i * cnts[2] + j + cnts[0] + 2]);
			}
			printf("\n");
		}
		free_u8v(bs[0]);
		free_u8v(bs[1]);
	}
#endif
	// detach
	clear_uuhash(nhash);
	for(i=0;i<=cnts[2]-2;i++){
		if(i == cnts[0] || i == cnts[0] + 1) continue;
		k = i > cnts[0];
		beg = k? 0 : cnts[0] + 2;
		end = k? cnts[0] : cnts[2] - 2;
		m = cnts[2];
		//sum = 0;
		for(j=beg;j<end;j++){
			n = mtx[i * cnts[2] + j];
			//sum += n;
			if(n >= g->min_edge_cov){
				if(m == cnts[2]){
					m = j;
				} else if(m < cnts[2]){
					m = cnts[2] + 1;
				}
			}
		}
		//sum += mtx[i * cnts[2] + end] + mtx[i * cnts[2] + end + 1];
		if(m == cnts[2] + 1){
			// multiple partners
			continue;
		} else if(m < cnts[2]){
			// maybe unique partner, need further check whether both unique
			put_uuhash(nhash, (uuhash_t){i, m});
		} else {
			// Remove tips
			//n = mtx[i * cnts[2] + end + 1];
			//if(n < g->min_edge_cov){
			{
				// maybe tip
				idx = matrix->buffer[i - 2 * k];
				e = ref_edgev(g->edges, idx >> 1);
				if(is_edge_closed(e)){ // self-loop, e has been cut in previous i
					continue;
				}
				flg = idx & 0x01;
				w = ref_nodev(g->nodes, get_edge_didx(e, flg));
				utg = ref_utgv(g->utgs, uid);
				wtg = ref_utgv(g->utgs, w->utg_idx);
				//if(wtg->closed) continue;
				dir = get_edge_ddir(e, flg) ^ w->utg_dir;
				if(wtg->cnt <= g->tip_step){
					if(count_utg_edges_graph(g, w->utg_idx, dir) == 0){
						// found tip candidate, need to check whether shorter than other branches
						v = ref_nodev(g->nodes, get_edge_sidx(e, flg));
						beg_iter_edges_graph(v, get_edge_sdir(e, flg), &eref);
						while((f = ref_iter_edges_graph(g, &eref))){
							if(f == e) continue;
							(*visit) ++;
							n = count_nonlinear_traces_graph(g, get_edge_didx(f, eref.flg), get_edge_ddir(f, eref.flg), wtg->cnt + 1, *visit, heap);
							if(n > wtg->cnt || (n == wtg->cnt && get_edge_cov(e) <= get_edge_cov(f))){
								ret ++;
								cut_edge_graph(g, e);
#if DEBUG
								if(print_rst){
									printf("TIP\tN%llu%c\n", get_edge_didx(e, flg), "+-"[!get_edge_ddir(e, flg)]);
								}
#endif
								break;
							}
						}
					}
				}
			}
		}
	}
	if(ret == 0){
		reset_iter_uuhash(nhash);
		while((nh = ref_iter_uuhash(nhash))){
			if(nh->key > nh->val) continue;
			if(!exists_uuhash(nhash, nh->val)){
				continue;
			}
			// Detach repeats
			eidxs[0] = nh->key;
			eidxs[1] = nh->val - 2 - cnts[0];
			if(_readpath_split_utg_graph(g, uid, eidxs, cnts, matrix, path, rdbits, pregs)){
				ret ++;
#if DEBUG
				if(print_rst){
					printf("REP\t%d\t%d\tYES\n", eidxs[0], eidxs[1]);
				}
#endif
			} else {
#if DEBUG
				if(print_rst){
					printf("REP\t%d\t%d\tNO\n", eidxs[0], eidxs[1]);
				}
#endif
			}
		}
	}
	while(ret == 0 && g->utgs->buffer[uid].cnt <= g->tip_step){
		// Split tip-like utg
		if(cnts[0] == 0){
			k = 0;
		} else if(cnts[1] == 0){
			k = 1;
		} else {
			break;
		}
		i = k? cnts[0] + 2 : 0;
		beg = k? 0 : cnts[0] + 2;
		end = k? cnts[0] : cnts[2] - 2;
		sum = 0;
		for(j=beg;j<end;j++){
			sum += mtx[i * cnts[2] + j];
		}
		{
			eidxs[0] = 0;
			eidxs[1] = 0;
			if(_readpath_split_utg_graph(g, uid, eidxs, cnts, matrix, path, rdbits, pregs)){
				ret ++;
#if DEBUG
				if(print_rst){
					printf("SPLIT\t%d\t%d\tYES\n", eidxs[0], eidxs[1]);
				}
#endif
			}
		}
		break;
	}
	if(ret == 0 && cnts[0] > 1 && cnts[1] > 1){
		// Split longest tip-like branch
		u4i tlen;
		tlen = 0;
		for(k=0;k<2;k++){
			if(k == 0){
				t = tail_tracev(path);
				v = ref_nodev(g->nodes, t->node);
				beg_iter_edges_graph(v, t->dir, &eref);
			} else {
				t = head_tracev(path);
				v = ref_nodev(g->nodes, t->node);
				beg_iter_edges_graph(v, !t->dir, &eref);
			}
			while((e = ref_iter_edges_graph(g, &eref))){
				w = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
				wtg = ref_utgv(g->utgs, w->utg_idx);
				if(wtg->cnt > g->tip_step) continue;
				dir = get_edge_ddir(e, eref.flg) ^ w->utg_dir;
				if(count_utg_edges_graph(g, w->utg_idx, dir) == 0){
					if(wtg->cnt > tlen){
						tlen = wtg->cnt;
						eidxs[k]  = eref.cnt - 1;
						eidxs[!k] = cnts[!k];
					}
				}
			}
		}
		if(tlen){
			if(_readpath_split_utg_graph(g, uid, eidxs, cnts, matrix, path, rdbits, pregs)){
				ret ++;
#if DEBUG
				if(print_rst){
					printf("SPLIT\t%d\t%d\tYES\tTCNT=%d\n", eidxs[0], eidxs[1], tlen);
				}
#endif
			}
		}
	}
	return ret;
}

static inline u8i readpath_detach_branches_graph(Graph *g){
	BitVec *rdbits;
	pregv *pregs;
	uuhash *nhash;
	u8v *matrix, *rgidxs;
	u8v *heap;
	tracev *path;
	utg_t *utg;
	node_t *n;
	u8i nid, vst, ret;
	u4i uid;
	ret = 0;
	rdbits = init_bitvec(g->reads->size);
	pregs = init_pregv(1024);
	nhash = init_uuhash(13);
	matrix = init_u8v(1024);
	rgidxs = init_u8v(4);
	heap = init_u8v(32);
	path = init_tracev(32);
	gen_utgs_graph(g);
	for(nid=0;nid<g->nodes->size;nid++){
		n = ref_nodev(g->nodes, nid);
		n->bt_visit = 0;
	}
	vst = 1;
	for(uid=0;uid<g->utgs->size;uid++){
		utg = ref_utgv(g->utgs, uid);
		if(utg->closed) continue;
		if(count_utg_edges_graph(g, uid, 0) < 2 && count_utg_edges_graph(g, uid, 1) < 2){
			continue;
		}
		clear_tracev(path);
		trace_utg_graph(g, uid, path);
		ret += readpath_detach_branches_graph_core(g, uid, path, rdbits, pregs, nhash, matrix, rgidxs, &vst, heap);
	}
	free_bitvec(rdbits);
	free_pregv(pregs);
	free_uuhash(nhash);
	free_u8v(matrix);
	free_u8v(rgidxs);
	free_u8v(heap);
	free_tracev(path);
	return ret;
}

// find a bubble, label visited utgs with *vst
static inline u8i fwd_search_utg_bubble_graph(Graph *g, u4i uid, int dir, u8i *vst, u4i max_len, u4v *heap){
	utg_t *s, *t;
	node_t *v, *w;
	edge_ref_t eref;
	edge_t *e;
	u4i idx, k;
	int ms;
	(*vst) ++;
	s = ref_utgv(g->utgs, uid);
	s->bt_vst = *vst;
	s->bt_dir = dir;
	s->bt_idx = 0;
	ms = max_len + s->len;
	clear_u4v(heap);
	array_heap_push(heap->buffer, heap->size, heap->cap, u4i, uid, num_cmp(g->utgs->buffer[a].bt_idx, g->utgs->buffer[b].bt_idx));
	while(heap->size){
		idx = heap->buffer[0];
		array_heap_remove(heap->buffer, heap->size, heap->cap, u4i, 0, num_cmp(g->utgs->buffer[a].bt_idx, g->utgs->buffer[b].bt_idx));
		s = ref_utgv(g->utgs, idx);
		v = ref_nodev(g->nodes, s->nodes[s->bt_dir].node);
		k = s->nodes[s->bt_dir].dir;
		beg_iter_edges_graph(v, k, &eref);
		while((e = ref_iter_edges_graph(g, &eref))){
			if(Int(s->bt_idx + s->len) + get_edge_off(e) > ms) continue;
			w = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
			t = ref_utgv(g->utgs, w->utg_idx);
			if(s == t){
				continue;
			}
			if(t->bt_vst == *vst){
				if(get_edge_ddir(e, eref.flg) ^ w->utg_dir ^ t->bt_dir){
					// loop
					return MAX_U8;
				} else {
					// bubble
					return (((u8i)w->utg_idx) << 1) | (get_edge_ddir(e, eref.flg) ^ w->utg_dir);
				}
			} else {
				t->bt_vst = *vst;
				t->bt_dir = get_edge_ddir(e, eref.flg) ^ w->utg_dir;
				t->bt_idx = s->bt_idx + s->len + get_edge_off(e);
				array_heap_push(heap->buffer, heap->size, heap->cap, u4i, w->utg_idx, num_cmp(g->utgs->buffer[a].bt_idx, g->utgs->buffer[b].bt_idx));
			}
		}
	}
	return MAX_U8;
}

// trace back along with utgs with *vst, store all utgs within bubble into uids
static inline u4i rev_search_utg_bubble_graph(Graph *g, u4i uid, int dir, u8i *vst, u4v *stack, u4v *uids){
	utg_t *s, *t;
	node_t *v, *w;
	edge_ref_t eref;
	edge_t *e;
	u8i vts;
	u4i idx, k;
	vts = *vst;
	(*vst) ++;
	clear_u4v(stack);
	clear_u4v(uids);
	s = ref_utgv(g->utgs, uid);
	s->bt_vst = *vst;
	s->bt_dir = dir;
	push_u4v(stack, uid);
	while(pop_u4v(stack, &idx)){
		s = ref_utgv(g->utgs, idx);
		push_u4v(uids, idx);
		v = ref_nodev(g->nodes, s->nodes[s->bt_dir].node);
		k = s->nodes[s->bt_dir].dir;
		beg_iter_edges_graph(v, k, &eref);
		while((e = ref_iter_edges_graph(g, &eref))){
			w = ref_nodev(g->nodes, get_edge_didx(e, eref.flg));
			t = ref_utgv(g->utgs, w->utg_idx);
			if(t->bt_vst != vts){
				continue;
			}
			t->bt_vst = *vst;
			t->bt_dir = get_edge_ddir(e, eref.flg) ^ w->utg_dir;
			push_u4v(stack, w->utg_idx);
		}
	}
	return uids->size;
}

static inline int readpath_merge_bubble_graph_core(Graph *g, u4i bidxs[2], u4v *uids, tracev *path, BitVec *rdbits[3], u8v *rdidxs[3]){
	utg_t *ub, *ue;
	trace_t *t;
	node_t *v;
	reg_t *r;
	edge_ref_t eref;
	edge_t *e;
	u4i uid, i, j, l, k, cnt, ecov, ret, inner, ends[2];
#if DEBUG
	int print_rst;
	if(g->utgs->buffer[bidxs[0]].nodes[0].node == debug_node || g->utgs->buffer[bidxs[0]].nodes[1].node == debug_node){
		print_rst = 1;
	} else {
		print_rst = 0;
	}
#endif
	clear_u8v(rdidxs[0]);
	clear_u8v(rdidxs[1]);
	clear_u8v(rdidxs[2]);
	ret = 0;
	for(k=0;k<2;k++){
		clear_tracev(path);
		trace_utg_graph(g, bidxs[k], path);
		for(j=0;j<path->size;j++){
			t = ref_tracev(path, j);
			v = ref_nodev(g->nodes, t->node);
			for(l=0;l<v->regs.cnt;l++){
				r = ref_regv(g->regs, v->regs.idx + l);
				if(r->closed) continue;
				if(get_bitvec(rdbits[k], r->rid) == 0){
					one_bitvec(rdbits[k], r->rid);
					push_u8v(rdidxs[k], r->rid);
				}
			}
		}
	}
	inner = 0;
	for(i=0;i<uids->size;i++){
		uid = uids->buffer[i];
		if(uid == bidxs[0] || uid == bidxs[1]){
			continue;
		}
		inner ++;
		clear_tracev(path);
		trace_utg_graph(g, uid, path);
		for(j=0;j<path->size;j++){
			t = ref_tracev(path, j);
			v = ref_nodev(g->nodes, t->node);
			for(l=0;l<v->regs.cnt;l++){
				r = ref_regv(g->regs, v->regs.idx + l);
				if(r->closed) continue;
				// find reads ending within bubble
				if(get_bitvec(rdbits[2], r->rid) == 0){
					k = (get_bitvec(rdbits[1], r->rid) << 1) | get_bitvec(rdbits[0], r->rid);
					if(k == 0 || k == 3) continue;
					k = k - 1; // 1 -> 0, 2 -> 1
					read_t *rd;
					reg_t *rg;
					node_t *w;
					u8i idx;
					rd = ref_readv(g->reads, r->rid);
					idx = rd->regs.idx;
					cnt = 0;
					ecov = 0;
					ends[0] = ends[1] = MAX_U4;
					while(idx){
						rg = ref_regv(g->regs, idx);
						idx = rg->read_link;
						if(rg->closed) continue;
						w = ref_nodev(g->nodes, rg->node);
						if(w->closed) continue;
						ends[cnt > 0] = w->utg_idx;
						cnt ++;
						if(w->utg_idx == uid) continue;
						if(w->utg_idx == bidxs[k]) continue;
						ecov ++;
						break;
					}
					if(cnt && ends[0] != bidxs[k] && ends[1] != bidxs[k]){
#if DEBUG
						if(print_rst){
							fprintf(stderr, " -- Keep %s in %s -- %s:%d --\n", g->kbm->reads->buffer[r->rid].tag, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						}
#endif
						continue;
					}
					one_bitvec(rdbits[2], r->rid);
					push_u8v(rdidxs[2], r->rid);
				}
			}
		}
	}
	for(i=0;i<uids->size;i++){
		uid = uids->buffer[i];
		if(uid == bidxs[0] || uid == bidxs[1]){
			continue;
		}
		clear_tracev(path);
		trace_utg_graph(g, uid, path);
		for(j=0;j<path->size;j++){
			t = ref_tracev(path, j);
			v = ref_nodev(g->nodes, t->node);
			cnt = 0;
			for(l=0;l<v->regs.cnt;l++){
				r = ref_regv(g->regs, v->regs.idx + l);
				if(r->closed) continue;
				 k = (get_bitvec(rdbits[1], r->rid) << 1) | get_bitvec(rdbits[0], r->rid);
				 if(k == 3){
					r->closed = 1;
					cnt ++;
				 } else if(k && get_bitvec(rdbits[2], r->rid)){
					r->closed = 1;
					cnt ++;
				 }
			}
			if(cnt == 0) continue;
			cnt = 0;
			for(k=0;k<2;k++){
				beg_iter_edges_graph(v, k, &eref);
				while((e = ref_iter_edges_graph(g, &eref))){
					ecov = cal_edge_cov_graph(g, e);
					set_edge_cov(e, ecov);
					if(ecov == 0 || (ecov < g->min_edge_cov && !get_edge_hard(e))){
						cut_edge_graph(g, e);
						cnt ++;
					}
				}
			}
			if(cnt){
				ret ++;
			}
#if DEBUG
			if(cnt == 10000){ // for debug call
				print_local_dot_graph(g, t->node, 10, "1.dot");
			}
#endif
			if(cnt == 0 || path->size < 2) continue;
			// revise utg
			revise_utg_graph(g, uid, path);
		}
	}
	for(k=0;k<3;k++){
		for(i=0;i<rdidxs[k]->size;i++){
			zero_bitvec(rdbits[k],  rdidxs[k]->buffer[i]);
		}
	}
	if(inner == 0){
#if DEBUG
		if(bidxs[0] != bidxs[1]){ // not a self-loop
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
#endif
		return 0;
	}
	if(ret == 0){
#if DEBUG
		//fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		//abort();
#endif
		return 0;
	}
	// try add edge between beg and end of this bubble
	ub = ref_utgv(g->utgs, bidxs[0]);
	ue = ref_utgv(g->utgs, bidxs[1]);
	e = edge_node2node_graph(g, ub->nodes[!ub->bt_dir].node, ub->nodes[!ub->bt_dir].dir, ue->nodes[ue->bt_dir].node, !ue->nodes[ue->bt_dir].dir, NULL);
	if(e){
		ecov = cal_edge_cov_graph(g, e);
		set_edge_cov(e, ecov);
		if(ecov == 0 || (ecov < g->min_edge_cov && !get_edge_hard(e))){
			cut_edge_graph(g, e);
		} else {
			revive_edge_graph(g, e);
		}
	} else {
		e = next_ref_edgev(g->edges);
		ZEROS(e);
		set_edge_sidx(e, 0, ub->nodes[!ub->bt_dir].node);
		set_edge_didx(e, 0, ue->nodes[ue->bt_dir].node);
		set_edge_sdir(e, 0, ub->nodes[!ub->bt_dir].dir);
		set_edge_ddir(e, 0, !ue->nodes[ue->bt_dir].dir);
		set_edge_off(e, (Int(ue->bt_idx - ub->bt_idx) - ub->len));
		ecov = cal_edge_cov_graph(g, e);
		set_edge_cov(e, ecov);
		if(ecov >= g->min_edge_cov){
			add_edge_graph(g, e);
		}
	}
	return 1;
}

static inline u4i readpath_merge_bubbles_graph(Graph *g, u4i max_len){
	BitVec *rdbits[3];
	tracev *path;
	u8v *rdidxs[3];
	u4v *heap, *uids;
	utg_t *utg;
	u8i vst, udi;
	u4i ret, uid, usize, mgr, k, bidxs[2];
	vst = 0;
	rdbits[0] = init_bitvec(g->reads->size);
	rdbits[1] = init_bitvec(g->reads->size);
	rdbits[2] = init_bitvec(g->reads->size);
	path = init_tracev(32);
	rdidxs[0] = init_u8v(32);
	rdidxs[1] = init_u8v(32);
	rdidxs[2] = init_u8v(32);
	heap = init_u4v(32);
	uids = init_u4v(32);
	gen_utgs_graph(g);
	ret = 0;
	usize = g->utgs->size; // utgs will increase, leave the new utg away
	for(uid=0;uid<usize;uid++){
		utg = ref_utgv(g->utgs, uid);
		for(k=0;k<2;k++){
#if DEBUG
			if(utg->nodes[k].node == debug_node){
				fprintf(stderr, " -- N%llu in %s -- %s:%d --\n", debug_node, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
#endif
			if(count_utg_edges_graph(g, uid, k) < 2) continue;
			udi = fwd_search_utg_bubble_graph(g, uid, k, &vst, max_len, heap);
			if(udi == MAX_U8) continue;
			if(uid == (udi >> 1)) continue;
			rev_search_utg_bubble_graph(g, udi >> 1, !(udi & 0x01), &vst, heap, uids);
			bidxs[0] = uid;
			bidxs[1] = udi >> 1;
			mgr = readpath_merge_bubble_graph_core(g, bidxs, uids, path, rdbits, rdidxs);
			if(mgr){
				ret ++;
			}
		}
	}
	free_bitvec(rdbits[0]);
	free_bitvec(rdbits[1]);
	free_bitvec(rdbits[2]);
	free_tracev(path);
	free_u8v(rdidxs[0]);
	free_u8v(rdidxs[1]);
	free_u8v(rdidxs[2]);
	free_u4v(heap);
	free_u4v(uids);
	return ret;
}

static inline u8i mask_all_branching_nodes_graph(Graph *g){
	node_t *n;
	u8i node, ret;
	ret = 0;
	for(node=0;node<g->nodes->size;node++){
		n = ref_nodev(g->nodes, node);
		if(n->closed) continue;
		if(n->erefs[0].cnt > 1 || n->erefs[1].cnt > 1){
			n->unvisit = 1;
			ret ++;
		} else {
			n->unvisit = 0;
		}
	}
	for(node=0;node<g->nodes->size;node++){
		n = ref_nodev(g->nodes, node);
		if(n->closed) continue;
		if(n->unvisit == 0) continue;
		del_node_graph(g, n);
	}
	return ret;
}

static inline u4i readpath_densify_nodes_graph(Graph *g, u8i eidx){
}

typedef struct {
	u8i rid:30, dir:1, beg:16, end:16, view:1;
} lay_reg_t;
define_list(layregv, lay_reg_t);

typedef struct {
	u4i tidx;
	u4i seqoff;
	u8i roff:48, rcnt:16;
} lay_t;
define_list(layv, lay_t);

static inline void traces2seqlets_graph(Graph *g, tracev *path, seqletv *sqs){
	trace_t *s, *t;
	seqlet_t *r;
	u4i i;
	s = NULL;
	for(i=0;i<path->size;i++){
		t = ref_tracev(path, i);
		if(s){
			r = next_ref_seqletv(sqs);
			r->node1 = s->node;
			r->dir1  = s->dir;
			r->node2 = t->node;
			r->dir2  = t->dir;
			r->off   = s->off;
			r->len   = g->reglen + t->off - s->off;
		}
		s = t;
	}
}

static inline void gen_lay_regs_core_graph(Graph *g, seqlet_t *q, layregv *regs){
	node_t *n1, *n2;
	reg_t *r1, *r2, *re1, *re2;
	u4i rid, beg, end;
	n1 = ref_nodev(g->nodes, q->node1);
	n2 = ref_nodev(g->nodes, q->node2);
	if(n1->regs.cnt == 0 || n2->regs.cnt == 0) return;
	r1 = ref_regv(g->regs, n1->regs.idx);
	re1 = r1 + n1->regs.cnt;
	r2 = ref_regv(g->regs, n2->regs.idx);
	re2 = r2 + n2->regs.cnt;
	while(1){
		if(r1->rid > r2->rid){
			r2 ++;
			if(r2 >= re2) break;
		} else if(r1->rid < r2->rid){
			r1 ++;
			if(r1 >= re1) break;
		} else {
			rid = r1->rid;
			if(r1->beg < r2->beg){
				if(q->dir1 ^ r1->dir){ r1 ++; r2 ++; continue; }
				beg = r1->beg; end = r2->end;
				if(beg > end){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				} else {
					push_layregv(regs, (lay_reg_t){rid, 0, beg, end, 0});
				}
			} else {
				if(!(q->dir1 ^ r1->dir)){ r1 ++; r2 ++; continue; }
				beg = r2->beg; end = r1->end;
				if(beg > end){
					fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				} else {
					push_layregv(regs, (lay_reg_t){rid, 1, beg, end, 0});
				}
			}
			r1 ++;
			if(r1 >= re1) break;
			r2 ++;
			if(r2 >= re2) break;
		}
	}
}

thread_beg_def(mlay);
Graph *g;
seqletv *path;
layv    *lays;
layregv *regs;
String  *seqs;
FILE *log;
thread_end_def(mlay);

thread_beg_func(mlay);
seqlet_t *let;
lay_t *lay;
lay_reg_t *reg;
Graph *g;
String *seqs;
u8i i;
u4i j;
thread_beg_loop(mlay);
clear_layregv(mlay->regs);
clear_string(mlay->seqs);
seqs = mlay->seqs;
g = mlay->g;
for(i=mlay->t_idx;i<mlay->path->size;i+=mlay->n_cpu){
	let = ref_seqletv(mlay->path, i);
	lay = ref_layv(mlay->lays, i);
	lay->tidx = mlay->t_idx;
	lay->roff = mlay->regs->size;
	lay->seqoff = mlay->seqs->size;
	gen_lay_regs_core_graph(mlay->g, let, mlay->regs);
	lay->rcnt = mlay->regs->size - lay->roff;
	sort_array(mlay->regs->buffer + lay->roff, lay->rcnt, lay_reg_t, num_cmpgt(b.end - b.beg, a.end - a.beg));
	if(lay->rcnt == 0 && mlay->log){
		thread_beg_syn(mlay);
		fprintf(mlay->log, " -- N%llu(%c) -> N%llu(%c) has no read path --\n", (u8i)let->node1, "+-"[let->dir1], (u8i)let->node2, "+-"[let->dir2]); fflush(mlay->log);
		thread_end_syn(mlay);
	}
	for(j=0;j<lay->rcnt;j++){
		reg = ref_layregv(mlay->regs, lay->roff + j);
		encap_string(seqs, (reg->end - reg->beg) * KBM_BIN_SIZE + 1);
		if(reg->dir){
			fwdseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[reg->rid].rdoff + reg->beg * KBM_BIN_SIZE, (reg->end - reg->beg) * KBM_BIN_SIZE, seqs->string + seqs->size);
		} else {
			fwdseq_basebank(g->kbm->rdseqs, g->kbm->reads->buffer[reg->rid].rdoff + reg->beg * KBM_BIN_SIZE, (reg->end - reg->beg) * KBM_BIN_SIZE, seqs->string + seqs->size);
		}
		seqs->string[seqs->size + (reg->end - reg->beg) * KBM_BIN_SIZE] = 0;
		seqs->size += (reg->end - reg->beg) * KBM_BIN_SIZE + 1;
	}
}
thread_end_loop(mlay);
thread_end_func(mlay);

static inline u4i print_utgs_graph(Graph *g, char *prefix, char *lay_suffix, u4i ncpu, FILE *log){
	FILE *o_lay;
	BufferedWriter *bw;
	utg_t *utg;
	layv *lays;
	layregv *regs;
	String *seqs;
	tracev *tpath;
	seqletv *path;
	seqlet_t *t;
	lay_t *lay;
	lay_reg_t *reg;
	u4i uid, ret, beg, end, j, c, len, soff;
	thread_preprocess(mlay);
	o_lay = open_file_for_write(prefix, lay_suffix, 1);
	bw = zopen_bufferedwriter(o_lay, 1024 * 1024, ncpu, 0);
	lays = init_layv(32);
	tpath = init_tracev(32);
	path = init_seqletv(32);
	thread_beg_init(mlay, ncpu);
	mlay->g = g;
	mlay->path = NULL;
	mlay->lays = lays;
	mlay->regs = init_layregv(32);
	mlay->seqs = init_string(1024);
	mlay->log  = log;
	thread_end_init(mlay);
	uid = 0;
	beg = 0;
	end = g->utgs->size;
	ret = 0;
	for(uid=beg;uid<end;uid++){
		utg = ref_utgv(g->utgs, uid);
		if(utg->closed) continue;
		clear_tracev(tpath);
		trace_utg_graph(g, uid, tpath);
		if(tpath->size < 2) continue;
		clear_seqletv(path);
		traces2seqlets_graph(g, tpath, path);
		clear_and_inc_layv(lays, path->size);
		thread_apply_all(mlay, EXPR(mlay->path = path));
		uid ++;
		ret ++;
		len = path->buffer[path->size - 1].off + path->buffer[path->size - 1].len;
		len = len * KBM_BIN_SIZE;
		{
			beg_bufferedwriter(bw);
			fprintf(bw->out, ">ctg%u nodes=%llu len=%u\n", uid, (u8i)path->size + 1, len);
			for(j=0;j<lays->size;j++){
				if((j % 100) == 0){
					flush_bufferedwriter(bw);
				}
				lay = ref_layv(lays, j);
				if(lay->rcnt == 0) continue;
				t = ref_seqletv(path, j);
				fprintf(bw->out, "E\t%d\tN%llu\t%c\tN%llu\t%c\n", (int)t->off * KBM_BIN_SIZE, (u8i)t->node1, "+-"[t->dir1], (u8i)t->node2, "+-"[t->dir2]);
				regs = thread_access(mlay, (j % ncpu))->regs;
				seqs = thread_access(mlay, (j % ncpu))->seqs;
				soff = lay->seqoff;
				for(c=0;c<lay->rcnt;c++){
					reg = ref_layregv(regs, lay->roff + c);
					fprintf(bw->out, "%c\t%s\t%c\t%d\t%d\t", "Ss"[reg->view], g->kbm->reads->buffer[reg->rid].tag, "+-"[reg->dir], reg->beg * KBM_BIN_SIZE, (reg->end - reg->beg) * KBM_BIN_SIZE);
					fprintf(bw->out, "%s\n", seqs->string + soff);
					soff += (reg->end - reg->beg) * KBM_BIN_SIZE + 1;
				}
			}
			end_bufferedwriter(bw);
		}
	}
	close_bufferedwriter(bw);
	fclose(o_lay);
	thread_beg_close(mlay);
	free_layregv(mlay->regs);
	free_string(mlay->seqs);
	thread_end_close(mlay);
	fprintf(KBM_LOGF, "[%s] output %u contigs\n", date(), (u4i)ret);
	free_tracev(tpath);
	free_seqletv(path);
	free_layv(lays);
	return uid;
}

static inline u8i print_dot_full_graph(Graph *g, FILE *out){
	BufferedWriter *bw;
	node_t *n;
	reg_t *r, *rr;
	edge_ref_t eref;
	edge_t *e;
	unsigned long long i;
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
		if(g->rdmaps && g->rdmaps->buffer[r->rid].mat){
			read_map_t *map;
			u4i refoff;
			map = ref_readmapv(g->rdmaps, r->rid);
			if(r->beg >= map->qb && r->end <= map->qe){
				if(map->refdir){
					refoff = map->tb + (map->qe - r->end) * KBM_BIN_SIZE;
				} else {
					refoff = map->tb + (r->beg - map->qb) * KBM_BIN_SIZE;
				}
				fprintf(bw->out, "N%llu [label=\"{N%llu %d | %s | %s_%c_%u_%u}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, g->reftags->buffer[map->refidx], "FR"[map->refdir ^ r->dir], refoff, (r->end - r->beg) * KBM_BIN_SIZE);
				continue;
			}
		}
		fprintf(bw->out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	for(i=0;i<g->nodes->size;i++){
		if((i % 1000) == 0) flush_bufferedwriter(bw);
		n = ref_nodev(g->nodes, i);
		//if(n->closed) continue;
		for(k=0;k<2;k++){
			beg_iter_edges_graph(n, k, &eref);
			while((e = ref_iter_edges_graph_core(g->edges, &eref, 1))){
				fprintf(bw->out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", (u8i)get_edge_sidx(e, eref.flg), (u8i)get_edge_didx(e, eref.flg), "+-"[get_edge_sdir(e, eref.flg)], "+-"[get_edge_ddir(e, eref.flg)], get_edge_cov(e), get_edge_off(e), colors[get_edge_sdir(e, eref.flg)][get_edge_ddir(e, eref.flg)], is_edge_closed(e)? " style=dashed" : "");
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
	edge_ref_t eref;
	edge_t *e;
	unsigned long long i;
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
		if(g->rdmaps && g->rdmaps->buffer[r->rid].mat){
			read_map_t *map;
			u4i refoff;
			map = ref_readmapv(g->rdmaps, r->rid);
			if(r->beg >= map->qb && r->end <= map->qe){
				if(map->refdir){
					refoff = map->tb + (map->qe - r->end) * KBM_BIN_SIZE;
				} else {
					refoff = map->tb + (r->beg - map->qb) * KBM_BIN_SIZE;
				}
				fprintf(bw->out, "N%llu [label=\"{N%llu %d | %s | %s_%c_%u_%u}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, g->reftags->buffer[map->refidx], "FR"[map->refdir ^ r->dir], refoff, (r->end - r->beg) * KBM_BIN_SIZE);
				continue;
			}
		}
		fprintf(bw->out, "N%llu [label=\"{N%llu %d | %s | %c_%d_%d}\"]\n", i, i, n->regs.cnt, g->kbm->reads->buffer[r->rid].tag, "FR"[r->dir], r->beg, r->end - r->beg);
	}
	for(i=0;i<g->nodes->size;i++){
		if((i % 1000) == 0) flush_bufferedwriter(bw);
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		for(k=0;k<2;k++){
			beg_iter_edges_graph(n, k, &eref);
			while((e = ref_iter_edges_graph(g, &eref))){
				fprintf(bw->out, "N%llu -> N%llu [label=\"%c%c:%d:%d\" color=%s%s]\n", (u8i)get_edge_sidx(e, eref.flg), (u8i)get_edge_didx(e, eref.flg), "+-"[get_edge_sdir(e, eref.flg)], "+-"[get_edge_ddir(e, eref.flg)], get_edge_cov(e), get_edge_off(e), colors[get_edge_sdir(e, eref.flg)][get_edge_ddir(e, eref.flg)], is_edge_closed(e)? " style=dashed" : "");
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
