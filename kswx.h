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

#ifndef __KSW_EXT_RJ_H
#define __KSW_EXT_RJ_H

#include "ksw.h"
#include "dna.h"
#include "chararray.h"
#include "list.h"
#include <math.h>
#include <float.h>

typedef struct {
	int score;
	int tb, te, qb, qe;
	int aln, mat, mis, ins, del;
} kswx_t;

static const kswr_t KSWR_NULL = (kswr_t){0, 0, 0, 0, 0, 0, 0};
static const kswx_t KSWX_NULL = (kswx_t){0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static inline kswr_t kswx2kswr(kswx_t x){
	kswr_t r;
	r.score = x.score;
	r.tb = x.tb;
	r.te = x.te;
	r.qb = x.qb;
	r.qe = x.qe;
	r.score2 = 0;
	r.te2 = 0;
	return r;
}

static inline void kswx_push_cigar(u32list *cigars, uint32_t op, uint32_t len){
	if(len == 0) return;
	if(cigars->size && (cigars->buffer[cigars->size-1] & 0xF) == op){
		cigars->buffer[cigars->size-1] += len << 4;
	} else push_u32list(cigars, (len << 4) | op);
}

static inline void kswx_push_cigars(u32list *cigars, uint32_t *cigar, size_t size){
	if(size == 0) return;
	if(cigars->size && (cigars->buffer[cigars->size-1] & 0x0F) == (cigar[0] & 0x0FU)){
		cigars->buffer[cigars->size-1] += cigar[0] & 0xFFFFFFF0U;
		append_array_u32list(cigars, cigar + 1, size - 1);
	} else append_array_u32list(cigars, cigar, size);
}

static inline void kswx_print_cigars(uint32_t *cigar, size_t size, FILE *out){
	size_t i;
	uint32_t op, len;
	for(i=0;i<size;i++){
		op  = cigar[i] & 0x0F;
		len = cigar[i] >> 4;
		fprintf(out, "%d%c", len, "MIDX"[op]);
	}
	fflush(out);
}

static inline kswx_t kswx_stat_cigars(uint32_t *cigar, size_t size){
	kswx_t x;
	size_t i;
	uint32_t op, len;
	x = KSWX_NULL;
	for(i=0;i<size;i++){
		op  = cigar[i] & 0x0F;
		len = cigar[i] >> 4;
		x.aln += len;
		switch(op){
			case 0: x.mat += len; x.te += len; x.qe += len; break;
			case 1: x.ins += len; x.qe += len; break;
			case 2: x.del += len; x.te += len; break;
		}
	}
	return x;
}

static inline void revseq_bytes(uint8_t *ary, int size){
	int i, v;
	for(i=0;i<size>>1;i++){
		v = ary[i]; ary[i] = ary[size - 1 - i]; ary[size - 1 - i] = v;
	}
}

static inline void revseq_4bytes(uint32_t *ary, int size){
	int i;
	uint32_t v;
	for(i=0;i<size>>1;i++){
		v = ary[i]; ary[i] = ary[size - 1 - i]; ary[size - 1 - i] = v;
	}
}

#define kswx_roundup8x(n) (((n) + 0x7LLU) & 0xFFFFFFFFFFFFFFF8LLU)

static inline void kswx_overlap_align_core(kswx_t *xs[2], u32list *cigars[2], int qlen, uint8_t *query, int tlen, uint8_t *target, int strand, int M, int X, int I, int D, int E, u8list *mem_cache){
	kswx_t x;
	int *rh, *re;
	uint8_t *z, *zi, d;
	int ql, tl, n_col;
	int i, j, k, h1, h, m, e, f, t;
	int imax, mj2, gmax, gi;
	if(xs[0]) *xs[0] = KSWX_NULL;
	if(xs[1]) *xs[1] = KSWX_NULL;
	if(qlen <= 0 || tlen <= 0){
		if(xs[0]) clear_u32list(cigars[0]);
		if(xs[1]) clear_u32list(cigars[1]);
		return;
	}
	ql = qlen; tl = tlen;
	n_col = tl;
	encap_u8list(mem_cache, kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x(((long long)ql) * n_col));
	rh = (int*)(mem_cache->buffer + mem_cache->size);
	re = (int*)(mem_cache->buffer + mem_cache->size + kswx_roundup8x((tl + 2) * sizeof(int)));
	z  =       (mem_cache->buffer + mem_cache->size + kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x((tl + 2) * sizeof(int)));
	for(j=0;j<=tl;j++) rh[j] = 0;
	for(j=0;j<=tl;j++) re[j] = -10000;
	gmax = 0; gi = -1;
	imax = 0; mj2 = -1;
	for(i=0;i<ql;i++){
		h1 = 0; f = -10000;
		zi = &z[i * n_col];
		for(j=0;j<tl;j++){
			m = rh[j] + ((query[i * strand] == target[j * strand])? M : X);
			rh[j] = h1;
			e = re[j];
			d = m >= e? 0 : 1;
			h = m >= e? m : e;
			d = h >= f? d : 2;
			h = h >= f? h : f;
			h1 = h;
			if(i + 1 == ql){
				imax = imax > h? imax : h;
				mj2  = imax > h? mj2  : j;
			}
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
		rh[j] = h1; re[j] = -10000;
		if(gmax < h1){ gmax = h1; gi = i; }
	}
	for(k=0;k<2;k++){
		if(xs[k] == NULL) continue;
		x = KSWX_NULL;
		if(k){
			x.score = gmax; x.qe = gi; x.te = tl - 1;
		} else {
			x.score = imax; x.qe = ql - 1; x.te = mj2;
		}
		i = x.qe; j = x.te;
		d = 0;
		if(cigars[k]) clear_u32list(cigars[k]);
		while(i >= 0 && j >= 0){
			d = (z[i * n_col + j] >> (d << 1)) & 0x03;
			if(d == 0){
				if(query[i * strand] == target[j * strand]){ x.mat ++; } else { x.mis ++; }  i --; j --; 
			} else if(d == 1){
				i --; x.ins ++;
			} else {
				j --; x.del ++;
			}
			if(cigars[k]) kswx_push_cigar(cigars[k], d, 1);
		}
		switch(d){
			case 1: x.qb = i + 1; x.tb = j; break;
			case 2: x.tb = i; x.tb = j + 1; break;
			default: x.qb = i + 1; x.tb = j + 1; break;
		}
		if(cigars[k]) reverse_u32list(cigars[k]);
		x.aln = x.mat + x.mis + x.ins + x.del;
		x.qe ++; x.te ++;
		*xs[k] = x;
	}
}

// diagonal will shift 1 bp to max score of each row
static inline kswx_t kswx_extend_align_shift_core(int qlen, uint8_t *query, int tlen, uint8_t *target, int strand, int init_score, int W, int M, int X, int I, int D, int E, int T, u8list *mem_cache, u32list *cigars){
	kswx_t x;
	int *rh, *re, *zb;
	uint8_t *z, *zi;
	int8_t  *qp, *qpp;
	register uint8_t d;
	int ql, tl, n_col;
	int i, j, jb, je, h1, c;
	register int h, m, e, f, t;
	int max_gap, max, mi, mj, imax, mj2, gmax, gi, gj;
	clear_u32list(cigars);
	x = KSWX_NULL;
	if(init_score < 0) init_score = 0;
	if(qlen <= 0 || tlen <= 0){ x.score = init_score; return x; }
	if(W > 0){
		max = ((qlen < tlen)? qlen : tlen) * M + init_score + (- T);
		max_gap = (max + (((I > D)? I : D))) / (- E) + 1;
		if(max_gap < 1) max_gap = 1;
		if(W > max_gap) W = max_gap;
	} else W = -W;
	W = num_min(W, num_max(qlen, tlen));
	if(qlen < tlen){
		if(qlen + W < tlen){ ql = qlen; tl = qlen + W; }
		else { ql = qlen; tl = tlen; }
	} else {
		if(tlen + W < qlen){ tl = tlen; ql = tlen + W; }
		else { tl = tlen; ql = qlen; }
	}
	n_col = tl < 2 * W + 1? tl : 2 * W + 1;
	encap_u8list(mem_cache, kswx_roundup8x(mem_cache->size) - mem_cache->size
		+ kswx_roundup8x((tl + 2) * sizeof(int)) * 2
		+ kswx_roundup8x((ql + 2) * sizeof(int))
		+ kswx_roundup8x(((long long)ql) * n_col)
		+ kswx_roundup8x(tl + 2) * sizeof(int8_t) * 4); // assume only A,C,G,T
	rh = (int*)(mem_cache->buffer + kswx_roundup8x(mem_cache->size));
	re = (int*)(mem_cache->buffer + kswx_roundup8x(mem_cache->size) + kswx_roundup8x((tl + 2) * sizeof(int)));
	zb = (int*)(mem_cache->buffer + kswx_roundup8x(mem_cache->size) + kswx_roundup8x((tl + 2) * sizeof(int)) * 2);
	z  =       (mem_cache->buffer + kswx_roundup8x(mem_cache->size) + kswx_roundup8x((tl + 2) * sizeof(int)) * 2 + kswx_roundup8x((ql + 2) * sizeof(int)));
	qp = (int8_t*)(mem_cache->buffer + kswx_roundup8x(mem_cache->size) + kswx_roundup8x((tl + 2) * sizeof(int)) * 2 + kswx_roundup8x((ql + 2) * sizeof(int)) + kswx_roundup8x(((long long)ql) * n_col));
	for(c=i=0;c<4;c++){
		for(j=0;j<tl;j++) qp[i++] = (c == target[j * strand])? M : X;
	}
	rh[0] = init_score; rh[1] = init_score + D + E;
	for(j=2;j<=tl;j++) rh[j] = rh[j - 1] + E;
	for(j=0;j<=tl;j++) re[j] = -10000;
	max = init_score; mi = -1; mj = -1; gmax = 0; gi = -1; gj = -1;
	jb = 0; je = tl;
	for(i=c=0;i<ql;i++){
		if(jb < c - W) jb = c - W;
		if(je > c + W + 1) je = c + W + 1;
		if(je > tl) je = tl;
		if(jb == 0) h1 = init_score + I + E * (i + 1);
		else h1 = -10000;
		zi = z + (i * n_col);
		zb[i] = jb;
		imax = 0; mj2 = -1; f = -10000;
		qpp = qp + query[i * strand] * tl;
		for(j=jb;j<je;j++){
			//m = rh[j] + ((query[i * strand] == target[j * strand])? M : X);
			m = rh[j] + qpp[j];
			rh[j] = h1;
			e = re[j];
			if(m >= e){ d = 0; h = m; }
			else      { d = 1; h = e; }
			if(h < f){ d = 2; h = f; }
			h1 = h;
			if(h > imax){ imax = h; mj2 = j; }
			t = m + I + E;
			e = e + E;
			if(e > t){ d |= 1 << 2; }
			else     { e = t; }
			re[j] = e;

			t = m + D + E;
			f = f + E;
			if(f > t){ d |= 2 << 4; }
			else     { f = t; }
			zi[j-jb] = d;
		}
		rh[j] = h1; re[j] = -10000;
		if(j == tlen && gmax < h1){ gmax = h1; gi = i; gj = j - 1; }
		if(i + 1 == qlen && gmax < imax){ gmax = imax; gi = i; gj = mj2; }
		if(imax > max){
			max = imax; mi = i; mj = mj2;
		} else if(imax <= 0) break;
		c ++;
		if(c < mj2){
			c ++;
			if(je < tl){
				rh[je + 1] = -10000;
				re[je + 1] = -10000;
			}
		} else if(c > mj2) {
			c --;
			if(jb){
				rh[jb-1] = -10000;
				re[jb-1] = -10000;
			}
		}
		jb = 0; je = tl;
	}
	if(gmax > 0 && gmax >= max + T){
		x.score = gmax; x.qe = gi; x.te = gj;
	} else {
		x.score =  max; x.qe = mi; x.te = mj;
	}
	i = x.qe; j = x.te;
	d = 0;
	while(i >= 0 && j >= 0){
		d = (z[i * n_col + (j - zb[i])] >> (d << 1)) & 0x03;
		if(d == 0){
			if(query[i * strand] == target[j * strand]){ x.mat ++; } else { x.mis ++; }  i --; j --; 
		} else if(d == 1){
			i --; x.ins ++;
		} else {
			j --; x.del ++;
		}
		if(cigars) kswx_push_cigar(cigars, d, 1);
	}
	if(i >= 0){
		x.ins += i + 1;
		if(cigars) kswx_push_cigar(cigars, 1, i + 1);
	}
	if(j >= 0){
		x.del += j + 1;
		if(cigars) kswx_push_cigar(cigars, 2, j + 1);
	}
	if(cigars) reverse_u32list(cigars);
	x.aln = x.mat + x.mis + x.ins + x.del;
	x.qe ++; x.te ++;
	return x;
}

static inline kswx_t kswx_extend_align_core(int qlen, uint8_t *query, int tlen, uint8_t *target, int strand, int init_score, int W, int M, int X, int I, int D, int E, int T, u8list *mem_cache, u32list *cigars){
	kswx_t x;
	int *rh, *re;
	uint8_t *z, *zi, d;
	int ql, tl, n_col;
	int i, j, jb, je, h1, h, m, e, f, t;
	int max_gap, max, mi, mj, imax, mj2, gmax, gi, gj;
	x = KSWX_NULL;
	if(init_score < 0) init_score = 0;
	if(qlen <= 0 || tlen <= 0){ x.score = init_score; return x; }
	if(W > 0){
		max = ((qlen < tlen)? qlen : tlen) * M + init_score + (- T);
		max_gap = (max + (((I > D)? I : D))) / (- E) + 1;
		if(max_gap < 1) max_gap = 1;
		if(W > max_gap) W = max_gap;
	} else W = -W;
	W = num_min(W, num_max(qlen, tlen));
	if(qlen < tlen){
		if(qlen + W < tlen){ ql = qlen; tl = qlen + W; }
		else { ql = qlen; tl = tlen; }
	} else {
		if(tlen + W < qlen){ tl = tlen; ql = tlen + W; }
		else { tl = tlen; ql = qlen; }
	}
	n_col = tl < 2 * W + 1? tl : 2 * W + 1;
	encap_u8list(mem_cache, kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x(((long long)ql) * n_col));
	rh = (int*)(mem_cache->buffer + mem_cache->size);
	re = (int*)(mem_cache->buffer + mem_cache->size + kswx_roundup8x((tl + 2) * sizeof(int)));
	z  =       (mem_cache->buffer + mem_cache->size + kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x((tl + 2) * sizeof(int)));
	rh[0] = init_score; rh[1] = init_score + D + E;
	for(j=2;j<=tl;j++) rh[j] = rh[j - 1] + E;
	for(j=0;j<=tl;j++) re[j] = -10000;
	max = init_score; mi = -1; mj = -1; gmax = 0; gi = -1; gj = -1;
	jb = 0; je = tl;
	for(i=0;i<ql;i++){
		jb = i - W; if(jb < 0) jb = 0;
		je = i + W + 1; if(je > tl) je = tl;
		if(jb == 0) h1 = init_score + I + E * (i + 1);
		else h1 = -10000;
		zi = &z[i * n_col];
		imax = 0; mj2 = -1; f = -10000;
		for(j=jb;j<je;j++){
			m = rh[j] + ((query[i * strand] == target[j * strand])? M : X);
			rh[j] = h1;
			e = re[j];
			d = m >= e? 0 : 1;
			h = m >= e? m : e;
			d = h >= f? d : 2;
			h = h >= f? h : f;
			h1 = h;
			imax = imax > h? imax : h;
			mj2  = imax > h? mj2  : j;
			t = m + I + E;
			e = e + E;
			d |= e > t? 1<<2 : 0;
			e = e > t? e : t;
			re[j] = e;
			t = m + D + E;
			f = f + E;
			d |= f > t? 2<<4 : 0;
			f = f > t? f : t;
			zi[j-jb] = d;
		}
		rh[j] = h1; re[j] = -10000;
		if(j == tlen && gmax < h1){ gmax = h1; gi = i; gj = j - 1; }
		if(i + 1 == qlen && gmax < imax){ gmax = imax; gi = i; gj = mj2; }
		if(imax > max){
			max = imax; mi = i; mj = mj2;
		} else if(imax <= 0) break;
	}
	if(gmax > 0 && gmax >= max + T){
		x.score = gmax; x.qe = gi; x.te = gj;
	} else {
		x.score =  max; x.qe = mi; x.te = mj;
	}
	i = x.qe; j = x.te;
	d = 0;
	clear_u32list(cigars);
	while(i >= 0 && j >= 0){
		d = (z[i * n_col + j - (i > W? i - W : 0)] >> (d << 1)) & 0x03;
		if(d == 0){
			if(query[i * strand] == target[j * strand]){ x.mat ++; } else { x.mis ++; }  i --; j --; 
		} else if(d == 1){
			i --; x.ins ++;
		} else {
			j --; x.del ++;
		}
		if(cigars) kswx_push_cigar(cigars, d, 1);
	}
	if(i >= 0){
		x.ins += i + 1;
		if(cigars) kswx_push_cigar(cigars, 1, i + 1);
	}
	if(j >= 0){
		x.del += j + 1;
		if(cigars) kswx_push_cigar(cigars, 2, j + 1);
	}
	if(cigars) reverse_u32list(cigars);
	x.aln = x.mat + x.mis + x.ins + x.del;
	x.qe ++; x.te ++;
	return x;
}

static inline kswx_t kswx_mismatch_free_extend_align_core(int qlen, uint8_t *query, int tlen, uint8_t *target, int strand, int init_score, int W, int M, int I, int D, int E, int T, u8list *mem_cache, u32list *cigars){
	kswx_t x;
	int *rh, *re;
	uint8_t *z, *zi, d;
	int ql, tl, n_col;
	int i, j, jb, je, h1, h, m, e, f;
	int max_gap, max, mi, mj, imax, mj2, gmax, gi, gj;
	x = KSWX_NULL;
	if(init_score < 0) init_score = 0;
	if(qlen <= 0 || tlen <= 0){ x.score = init_score; return x; }
	if(W > 0){
		max = ((qlen < tlen)? qlen : tlen) * M + init_score + (- T);
		max_gap = (max + (((I > D)? I : D))) / (- E) + 1;
		if(max_gap < 1) max_gap = 1;
		if(W > max_gap) W = max_gap;
	} else W = -W;
	W = num_min(W, num_max(qlen, tlen));
	if(qlen < tlen){
		if(qlen + W < tlen){ ql = qlen; tl = qlen + W; }
		else { ql = qlen; tl = tlen; }
	} else {
		if(tlen + W < qlen){ tl = tlen; ql = tlen + W; }
		else { tl = tlen; ql = qlen; }
	}
	n_col = tl < 2 * W + 1? tl : 2 * W + 1;
	encap_u8list(mem_cache, kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x(((long long)ql) * n_col));
	rh = (int*)(mem_cache->buffer + mem_cache->size);
	re = (int*)(mem_cache->buffer + mem_cache->size + kswx_roundup8x((tl + 2) * sizeof(int)));
	z  =       (mem_cache->buffer + mem_cache->size + kswx_roundup8x((tl + 2) * sizeof(int)) + kswx_roundup8x((tl + 2) * sizeof(int)));
	rh[0] = init_score; rh[1] = init_score + D + E;
	for(j=2;j<=tl;j++) rh[j] = rh[j - 1] + E;
	for(j=0;j<=tl;j++) re[j] = -10000;
	max = init_score; mi = -1; mj = -1; gmax = 0; gi = -1; gj = -1;
	jb = 0; je = tl;
	for(i=0;i<ql;i++){
		jb = i - W; if(jb < 0) jb = 0;
		je = i + W + 1; if(je > tl) je = tl;
		if(jb == 0) h1 = init_score + I + E * (i + 1);
		else h1 = -10000;
		zi = &z[i * n_col];
		imax = 0; mj2 = -1; f = h1 + D + E;
		for(j=jb;j<je;j++){
			e = re[j];
			if(query[i * strand] == target[j * strand]){
				m = rh[j] + M;
				rh[j] = h1;
				if(m >= e){
					if(m >= f){
						d = 0;
						h = m;
						f = m + D + E;
						e = m + I + E;
					} else {
						d = 2;
						h = f;
						f = f + E;
						e = f + I + E;
					}
				} else {
					if(e >= f){
						d = 1;
						h = e;
						f = e + D + E;
						e = e + E;
					} else {
						d = 2;
						h = f;
						f = f + E;
						e = f + I + E;
					}
				}
			} else {
				rh[j] = h1;
				{
					if(e >= f){
						d = 1;
						h = e;
						f = e + D + E;
						e = e + E;
					} else {
						d = 2;
						h = f;
						f = f + E;
						e = f + I + E;
					}
				}
			}
			h1 = h;
			if(imax < h){ imax = h; mj2 = j; }
			re[j] = e;
			zi[j-jb] = d;
		}
		rh[j] = h1; re[j] = -10000;
		if(j == tlen && gmax < h1){ gmax = h1; gi = i; gj = j - 1; }
		if(i + 1 == qlen && gmax < imax){ gmax = imax; gi = i; gj = mj2; }
		if(imax > max){
			max = imax; mi = i; mj = mj2;
		} else if(imax <= 0) break;
	}
	if(gmax > 0 && gmax >= max + T){
		x.score = gmax; x.qe = gi; x.te = gj;
	} else {
		x.score =  max; x.qe = mi; x.te = mj;
	}
	i = x.qe; j = x.te;
	d = 0;
	clear_u32list(cigars);
	while(i >= 0 && j >= 0){
		d = (z[i * n_col + j - (i > W? i - W : 0)]) & 0x03;
		if(d == 0){
			x.mat ++; i --; j --; 
		} else if(d == 1){
			i --; x.ins ++;
		} else {
			j --; x.del ++;
		}
		if(cigars) kswx_push_cigar(cigars, d, 1);
	}
	if(i >= 0){
		x.ins += i + 1;
		if(cigars) kswx_push_cigar(cigars, 1, i + 1);
	}
	if(j >= 0){
		x.del += j + 1;
		if(cigars) kswx_push_cigar(cigars, 2, j + 1);
	}
	if(cigars) reverse_u32list(cigars);
	x.aln = x.mat + x.mis + x.ins + x.del;
	x.qe ++; x.te ++;
	return x;
}

static inline kswx_t kswx_extend_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int strand, int init_score, int W, int M, int X, int I, int D, int E, int T, int *n_cigar, uint32_t **cigar){
	u8list  *caches;
	u32list *cigars;
	kswx_t x;
	caches = init_u8list(256);
	cigars = init_u32list(64);
	x = kswx_extend_align_shift_core(qlen, query, tlen, target, strand, init_score, W, M, X, I, D, E, T, caches, cigars);
	free_u8list(caches);
	*n_cigar = cigars->size;
	*cigar = cigars->buffer;
	free(cigars);
	return x;
}

static inline kswx_t kswx_flexible_banded_affine_alignment(uint8_t *query, int ql, uint8_t *target, int tl, int M, int X, int I, int D, int E, int T, int is_local, int *bands[2], u8list *mem, u32list *cigars){
	kswx_t x;
	int *rh, *re;
	uint8_t *z, *zi;
	long long tot_z;
	int8_t *qp, *qpp;
	register uint8_t d;
	register int h, h1, m, e, f, t, imax, mj2;
	int i, j, c, max, mi, mj, gmax, gi, gj;
	long long ii;
	clear_u32list(cigars);
	tot_z = 0;
	for(i=0;i<ql;i++){ tot_z += bands[1][i] - bands[0][i]; }
	clear_and_encap_u8list(mem, kswx_roundup8x((tl + 2) * sizeof(int)) * 2 // rh, re
		+ kswx_roundup8x(tot_z) // z
		+ kswx_roundup8x(tl + 2) * sizeof(int8_t) * 4); // qp
	rh = (int*)(mem->buffer);
	re = (int*)(((void*)rh) + kswx_roundup8x((tl + 2) * sizeof(int)));
	z  = (uint8_t*)(((void*)re) + kswx_roundup8x((tl + 2) * sizeof(int)));
	qp = (int8_t*)(((void*)z) + kswx_roundup8x(tot_z));
	for(c=ii=0;c<4;c++){
		for(j=0;j<tl;j++) qp[ii++] = (c == target[j])? M : X;
	}
	rh[0] = 0; rh[1] = D + E;
	for(i=2;i<=tl;i++) rh[i] = rh[i-1] + E;
	for(i=0;i<=tl;i++) re[i] = -10000;
	max = 0; mi = 0; mj = 0;
	gmax = 0; gi = 0; gj = 0;
	zi = z;
	for(i=0;i<ql;i++){
		qpp = qp + query[i] * tl;
		h1 = (bands[0][i] == 0)? I + E * (i + 1) : -10000;
		f = -10000;
		imax = -10000; mj2 = -1;
		for(j=bands[0][i];j<bands[1][i];j++){
			m = rh[j] + qpp[j];
			rh[j] = h1;
			e = re[j];
			if(m >= e){ d = 0; h = m; }
			else      { d = 1; h = e; }
			if(h < f){ d = 2; h = f; }
			if(is_local && h < 0){ d = 3; h = 0; }
			h1 = h;
			if(imax < h){ imax = h; mj2 = j; }
			t = m + I + E;
			e = e + E;
			if(e > t){ d |= 1 << 2; }
			else     { e = t; }
			re[j] = e;
			t = m + D + E;
			f = f + E;
			if(f > t){ d |= 2 << 4; }
			else     { f = t; }
			*zi = d;
			zi ++;
		}
		rh[j] = h1; re[j] = -10000;
		if(j == tl && gmax < h1){ gmax = h1; gi = i; gj = j - 1; }
		if(i + 1 == ql && gmax < imax){ gmax = imax; gi = i; gj = mj2; }
		if(imax > max){ max = imax; mi = i; mj = mj2; }
	}
	if(gmax >= max + T){
		x.score = gmax; x.qe = gi; x.te = gj;
	} else {
		x.score = max; x.qe = mi; x.te = mj;
	}
	i = x.qe; j = x.te;
	d = 0;
	for(i=ql-1;i>=x.qe;i--){ zi -= bands[1][i] - bands[0][i]; }
	clear_u32list(cigars);
	while(i >= 0 && j >= 0){
		d = (zi[j - bands[0][i]] >> (d << 1)) & 0x03;
		if(d == 3) break;
		if(d == 0){
			x.mat ++; i --; zi -= bands[1][i] - bands[0][i]; j --; 
		} else if(d == 1){
			i --; zi -= bands[1][i] - bands[0][i]; x.ins ++;
		} else {
			j --; x.del ++;
		}
		if(cigars) kswx_push_cigar(cigars, d, 1);
	}
	if(d != 3){
		if(i >= 0){
			x.ins += i + 1;
			if(cigars) kswx_push_cigar(cigars, 1, i + 1);
		}
		if(j >= 0){
			x.del += j + 1;
			if(cigars) kswx_push_cigar(cigars, 2, j + 1);
		}
	}
	x.qb = i; x.tb = j;
	if(cigars) reverse_u32list(cigars);
	x.aln = x.mat + x.mis + x.ins + x.del;
	x.qe ++; x.te ++;
	return x;
}

static inline kswx_t kswx_refine_alignment(uint8_t *query, int qb, uint8_t *target, int tb, int W, int M, int X, int I, int D, int E, u32list *cigars, u8list *mem_cache, u32list *cigars2){
	kswx_t y;
	int *rh, *re;
	int *zb, *ze, *zw;
	int ql, tl, _zb, _ze, tx, qx;
	int qe, te;
	uint8_t *z, *zi;
	int8_t  *qp, *qpp;
	register uint8_t d;
	register int h, m, e, f, t;
	long long i, j;
	int op, len, h1, c;
	clear_u32list(cigars2);
	qe = qb; te = tb;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0: qe += len; te += len; break;
			case 1: qe += len; break;
			default: te += len;
		}
	}
	ql = qe - qb;
	tl = te - tb;
	if(ql == 0 || tl == 0) return KSWX_NULL;
	clear_and_encap_u8list(mem_cache, 
		  kswx_roundup8x((tl + 2) * sizeof(int)) * 2 // rh, re
		+ kswx_roundup8x((ql + 2) * sizeof(int)) * 3 // zb, ze, zw
		+ kswx_roundup8x(((long long)ql) * tl) // z
		+ kswx_roundup8x(tl + 2) * sizeof(int8_t) * 4); // qp
	rh = (int*)(mem_cache->buffer);
	re = (int*)(((void*)rh) + kswx_roundup8x((tl + 2) * sizeof(int)));
	zb = (int*)(((void*)re) + kswx_roundup8x((tl + 2) * sizeof(int)));
	ze = (int*)(((void*)zb) + kswx_roundup8x((ql + 2) * sizeof(int)));
	zw = (int*)(((void*)ze) + kswx_roundup8x((ql + 2) * sizeof(int)));
	z  = (uint8_t*)(((void*)zw) + kswx_roundup8x((ql + 2) * sizeof(int)));
	qp = (int8_t*)(((void*)z) + kswx_roundup8x(((long long)ql) * tl));
	for(c=i=0;c<4;c++){
		for(j=0;j<tl;j++) qp[i++] = (c == target[j + tb])? M : X;
	}
	// calculate alignment band
	// basic band
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			for(j=0;j<len;j++) zw[qx++] = W;
			break;
			case 1:
			for(j=0;j<len;j++) zw[qx++] = W + len;
			break;
			default:
			zw[qx] += len;
		}
	}
	// extending bandwidth around indels
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			qx += len;
			break;
			case 1:
			for(j=1;j<len&&j<qx;j++) zw[qx-j] += len - j;
			qx += len - 1;
			for(j=1;j<len&&j+qx<ql;j++) zw[qx+j] += len - j;
			qx ++;
			break;
			default:
			for(j=1;j<len&&j<qx;j++) zw[qx-j] += len - j;
			for(j=1;j<len&&j+qx<ql;j++) zw[qx+j] += len - j;
		}
	}
	// generate band region
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			for(j=0;j<len;j++){
				_zb = tx - zw[qx]; if(_zb < 0) _zb = 0;
				_ze = tx + 1 + zw[qx]; if(_ze > tl) _ze = tl;
				zb[qx] = _zb;
				ze[qx] = _ze;
				tx ++;
				qx ++;
			}
			break;
			case 1:
			for(j=0;j<len;j++){
				_zb = tx - zw[qx]; if(_zb < 0) _zb = 0;
				_ze = tx + 1 + zw[qx]; if(_ze > tl) _ze = tl;
				zb[qx] = _zb;
				ze[qx] = _ze;
				qx ++;
			}
			break;
			default:
			tx += len;
		}
	}
	// trim band beg
	_zb = 0;
	for(i=0;i<ql;i++){
		if(zb[i] < _zb) zb[i] = _zb;
		else if(zb[i] > _zb) _zb = zb[i];
	}
	// trim band end
	_ze = tl;
	for(i=ql-1;i>=0;i--){
		if(ze[i] > _ze) ze[i] = _ze;
		else if(ze[i] < _ze) _ze = ze[i];
	}
	// perform alignment
	rh[0] = 0;
	for(i=1;i<=tl;i++) rh[i] = -10000;
	for(i=0;i<=tl;i++) re[i] = -10000;
	for(i=0;i<ql;i++){
		qpp = qp + query[i + qb] * tl;
		h1 = f = -10000;
		zi = z + (((long long)i) * tl);
		//fprintf(stdout, "ZW\t%d\t%d\t%d\t%d\n", (int)i, zb[i], ze[i], ze[i] - zb[i]);
		for(j=zb[i];j<ze[i];j++){
			m = rh[j] + qpp[j];
			rh[j] = h1;
			e = re[j];
			if(m >= e){ d = 0; h = m; }
			else      { d = 1; h = e; }
			if(h < f){ d = 2; h = f; }
			h1 = h;
			t = m + I + E;
			e = e + E;
			if(e > t){ d |= 1 << 2; }
			else     { e = t; }
			re[j] = e;
			t = m + D + E;
			f = f + E;
			if(f > t){ d |= 2 << 4; }
			else     { f = t; }
			zi[j-zb[i]] = d;
		}
		rh[j] = h1; re[j] = -10000;
	}
	y.qb = qb; y.qe = qe; y.tb = tb; y.te = te;
	y.score = rh[tl];
	y.mat = y.mis = y.ins = y.del = y.aln = 0;
	d = 0;
	i = ql - 1; j = tl - 1;
	while(i >= 0 && j >= 0){
		d = (z[i * tl + (j - zb[i])] >> (d << 1)) & 0x03;
		if(d == 0){
			if(query[i + y.qb] == target[j + y.tb]){ y.mat ++; } else { y.mis ++; }  i --; j --; 
		} else if(d == 1){
			i --; y.ins ++;
		} else {
			j --; y.del ++;
		}
		kswx_push_cigar(cigars2, d, 1);
	}
	if(i >= 0){
		y.ins += i + 1;
		kswx_push_cigar(cigars2, 1, i + 1);
	}
	if(j >= 0){
		y.del += j + 1;
		kswx_push_cigar(cigars2, 2, j + 1);
	}
	reverse_u32list(cigars2);
	y.aln = y.mat + y.mis + y.ins + y.del;
	return y;
}

// qvs: 0:errQV, 1:misQV, 2:insQV, 3:delQV, 4:mrgQV, 5, misTag, 6, delTag
static inline kswx_t kswx_refine_alignment_5q(uint8_t *query, int qb, uint8_t *target, int tb, int W, uint8_t QCLP, uint8_t QMIS, uint8_t QDEL, uint8_t *qvs[7], u32list *cigars, u8list *mem_cache, u32list *cigars2){
	kswx_t y;
	uint32_t *rh, SCORE_INF;
	int *zb, *ze, *zw;
	int ql, tl, rl, _zb, _ze, tx, qx;
	int qe, te;
	uint8_t *z, *zi;
	uint8_t d, qmat[4], qins, qdel[4];
	register uint32_t h1, h, m, e, f;
	long long i, j;
	int op, len;
	SCORE_INF = 0xFFFFFFFFU >> 1;
	clear_u32list(cigars2);
	qe = qb; te = tb;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0: qe += len; te += len; break;
			case 1: qe += len; break;
			default: te += len;
		}
	}
	ql = qe - qb;
	tl = te - tb;
	if(ql == 0 || tl == 0) return KSWX_NULL;
	rl = tl;
	clear_and_encap_u8list(mem_cache, 
		  kswx_roundup8x((tl + 2) * sizeof(uint32_t)) // rh
		+ kswx_roundup8x((ql + 2) * sizeof(int)) * 3 // zb, ze, zw
		+ kswx_roundup8x(((long long)ql) * rl)); // z
	rh = (uint32_t*)(mem_cache->buffer);
	zb = (int*)(((void*)rh) + kswx_roundup8x((tl + 2) * sizeof(uint32_t)));
	ze = (int*)(((void*)zb) + kswx_roundup8x((ql + 2) * sizeof(int)));
	zw = (int*)(((void*)ze) + kswx_roundup8x((ql + 2) * sizeof(int)));
	z  = (uint8_t*)(((void*)zw) + kswx_roundup8x((ql + 2) * sizeof(int)));
	// calculate alignment band
	// basic band
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			for(j=0;j<len;j++) zw[qx++] = W;
			break;
			case 1:
			for(j=0;j<len;j++) zw[qx++] = W + len;
			break;
			default:
			zw[qx] += len;
		}
	}
	// extending bandwidth around indels
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			qx += len;
			break;
			case 1:
			for(j=1;j<len&&j<qx;j++) zw[qx-j] += len - j;
			qx += len - 1;
			for(j=1;j<len&&j+qx<ql;j++) zw[qx+j] += len - j;
			qx ++;
			break;
			default:
			for(j=1;j<len&&j<qx;j++) zw[qx-j] += len - j;
			for(j=1;j<len&&j+qx<ql;j++) zw[qx+j] += len - j;
		}
	}
	// generate band region
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			for(j=0;j<len;j++){
				_zb = tx - zw[qx]; if(_zb < 0) _zb = 0;
				_ze = tx + 1 + zw[qx]; if(_ze > tl) _ze = tl;
				zb[qx] = _zb;
				ze[qx] = _ze;
				tx ++;
				qx ++;
			}
			break;
			case 1:
			for(j=0;j<len;j++){
				_zb = tx - zw[qx]; if(_zb < 0) _zb = 0;
				_ze = tx + 1 + zw[qx]; if(_ze > tl) _ze = tl;
				zb[qx] = _zb;
				ze[qx] = _ze;
				qx ++;
			}
			break;
			default:
			tx += len;
		}
	}
	// trim band beg
	_zb = 0;
	for(i=0;i<ql;i++){
		if(zb[i] < _zb) zb[i] = _zb;
		else if(zb[i] > _zb) _zb = zb[i];
	}
	// trim band end
	_ze = tl;
	for(i=ql-1;i>=0;i--){
		if(ze[i] > _ze) ze[i] = _ze;
		else if(ze[i] < _ze) _ze = ze[i];
	}
	// perform alignment
	for(i=0;i<=tl;i++) rh[i] = SCORE_INF;
	for(i=0;i<ql;i++){
		if(i == 0){
			h1 = i * QCLP;
			f = (i + 1) * QCLP;
		} else {
			h1 = f = SCORE_INF;
		}
		zi = z + (i * rl);
		for(j=0;j<4;j++){
			if(j == query[i + qb]) qmat[j] = 0;
			else if(j == qvs[5][i + qb]) qmat[j] = qvs[1][i + qb];
			else qmat[j] = QMIS;
		}
		qins = qvs[2][i + qb];
		if(i + 1 == ql) qdel[0] = qdel[1] = qdel[2] = qdel[3] = QCLP;
		else {
			for(j=0;j<4;j++){
				if(j == qvs[6][i + 1 + qb]){
					qdel[j] = qvs[3][i + 1 + qb];
				} else qdel[j] = QDEL;
			}
		}
		for(j=zb[i];j<ze[i];j++){
			m = h1 + qmat[target[j + tb]];
			f = f  + qdel[target[j + tb]];
			if(m <= f){ d = 0; h = m; }
			else      { d = 2; h = f; }
			e = rh[j] + qins;
			if(e < h){ d = 1; h = e; }
			h1 = rh[j];
			rh[j] = f = h;
			zi[j-zb[i]] = d;
		}
	}
	y.qb = qb; y.qe = qe; y.tb = tb; y.te = te;
	y.score = - (int)rh[tl-1];
	y.mat = y.mis = y.ins = y.del = y.aln = 0;
	d = 0;
	i = ql - 1; j = tl - 1;
	//String *qs = init_string(1024);
	//String *ts = init_string(1024);
	//String *ms = init_string(1024);
	while(i >= 0 && j >= 0){
		d = z[i * rl + (j - zb[i])] & 0x03;
		if(d == 0){
			//add_char_string(qs, bit_base_table[query[i + y.qb]]);
			//add_char_string(ts, bit_base_table[target[j + y.tb]]);
			if(query[i + y.qb] == target[j + y.tb]){
				//add_char_string(ms, '|');
				y.mat ++;
			} else {
				//add_char_string(ms, '*');
				y.mis ++;
			}
			i --; j --; 
		} else if(d == 1){
			//add_char_string(qs, bit_base_table[query[i + y.qb]]);
			//add_char_string(ts, '-');
			//add_char_string(ms, '-');
			i --; y.ins ++;
		} else {
			//add_char_string(qs, '-');
			//add_char_string(ts, bit_base_table[target[j + y.tb]]);
			//add_char_string(ms, '-');
			j --; y.del ++;
		}
		kswx_push_cigar(cigars2, d, 1);
	}
	//reverse_string(qs);
	//reverse_string(ts);
	//reverse_string(ms);
	//printf("%s\n%s\n%s\n", qs->string, ms->string, ts->string);
	//free_string(qs);
	//free_string(ts);
	//free_string(ms);
	if(1){
		if(i >= 0){
			y.ins += i + 1;
			kswx_push_cigar(cigars2, 1, i + 1);
		}
		if(j >= 0){
			y.del += j + 1;
			kswx_push_cigar(cigars2, 2, j + 1);
		}
	} else {
		y.qb += i + 1;
		y.tb += j + 1;
	}
	reverse_u32list(cigars2);
	y.aln = y.mat + y.mis + y.ins + y.del;
	return y;
}

static inline kswx_t kswx_refine_affine_alignment_5q(uint8_t *query, int qb, uint8_t *target, int tb, int W, uint8_t QCLP, uint8_t QMIS, uint8_t QDEL, uint8_t QEXT, uint8_t *qvs[7], u32list *cigars, u8list *mem_cache, u32list *cigars2){
	kswx_t y;
	uint32_t *rh, *re, SCORE_INF;
	int *zb, *ze, *zw;
	int ql, tl, rl, _zb, _ze, tx, qx;
	int qe, te;
	uint8_t *z, *zi;
	uint8_t d, qmat[4], qins, qdel[4];
	register uint32_t h1, h, m, e, f, t;
	long long i, j;
	int op, len;
	SCORE_INF = 0xFFFFFFFFU >> 1;
	clear_u32list(cigars2);
	qe = qb; te = tb;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0: qe += len; te += len; break;
			case 1: qe += len; break;
			default: te += len;
		}
	}
	ql = qe - qb;
	tl = te - tb;
	if(ql == 0 || tl == 0) return KSWX_NULL;
	rl = tl;
	clear_and_encap_u8list(mem_cache, 
		  kswx_roundup8x((tl + 2) * sizeof(uint32_t)) * 2 // rh, re
		+ kswx_roundup8x((ql + 2) * sizeof(int)) * 3 // zb, ze, zw
		+ kswx_roundup8x(((long long)ql) * rl)); // z
	rh = (uint32_t*)(mem_cache->buffer);
	re = (uint32_t*)(((void*)rh) + kswx_roundup8x((tl + 2) * sizeof(uint32_t)));
	zb = (int*)(((void*)re) + kswx_roundup8x((tl + 2) * sizeof(uint32_t)));
	ze = (int*)(((void*)zb) + kswx_roundup8x((ql + 2) * sizeof(int)));
	zw = (int*)(((void*)ze) + kswx_roundup8x((ql + 2) * sizeof(int)));
	z  = (uint8_t*)(((void*)zw) + kswx_roundup8x((ql + 2) * sizeof(int)));
	// calculate alignment band
	// basic band
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			for(j=0;j<len;j++) zw[qx++] = W;
			break;
			case 1:
			for(j=0;j<len;j++) zw[qx++] = W + len;
			break;
			default:
			zw[qx] += len;
		}
	}
	// extending bandwidth around indels
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			qx += len;
			break;
			case 1:
			for(j=1;j<len&&j<qx;j++) zw[qx-j] += len - j;
			qx += len - 1;
			for(j=1;j<len&&j+qx<ql;j++) zw[qx+j] += len - j;
			qx ++;
			break;
			default:
			for(j=1;j<len&&j<qx;j++) zw[qx-j] += len - j;
			for(j=1;j<len&&j+qx<ql;j++) zw[qx+j] += len - j;
		}
	}
	// generate band region
	tx = qx = 0;
	for(i=0;i<(int)cigars->size;i++){
		op  = cigars->buffer[i] & 0xFU;
		len = cigars->buffer[i] >> 4;
		switch(op){
			case 0:
			for(j=0;j<len;j++){
				_zb = tx - zw[qx]; if(_zb < 0) _zb = 0;
				_ze = tx + 1 + zw[qx]; if(_ze > tl) _ze = tl;
				zb[qx] = _zb;
				ze[qx] = _ze;
				tx ++;
				qx ++;
			}
			break;
			case 1:
			for(j=0;j<len;j++){
				_zb = tx - zw[qx]; if(_zb < 0) _zb = 0;
				_ze = tx + 1 + zw[qx]; if(_ze > tl) _ze = tl;
				zb[qx] = _zb;
				ze[qx] = _ze;
				qx ++;
			}
			break;
			default:
			tx += len;
		}
	}
	// trim band beg
	_zb = 0;
	for(i=0;i<ql;i++){
		if(zb[i] < _zb) zb[i] = _zb;
		else if(zb[i] > _zb) _zb = zb[i];
	}
	// trim band end
	_ze = tl;
	for(i=ql-1;i>=0;i--){
		if(ze[i] > _ze) ze[i] = _ze;
		else if(ze[i] < _ze) _ze = ze[i];
	}
	// perform alignment
	rh[0] = 0;
	for(i=1;i<=tl;i++) rh[i] = rh[i - 1] + QCLP;
	for(i=0;i<=tl;i++) re[i] = SCORE_INF;
	for(i=0;i<ql;i++){
		if(zb[i] == 0){
			h1 = i * QCLP;
			f = h1 + QCLP;
		} else {
			h1 = f = SCORE_INF;
		}
		zi = z + (i * rl);
		for(j=0;j<4;j++){
			if(j == query[i + qb]) qmat[j] = 0;
			else if(j == qvs[5][i + qb]) qmat[j] = qvs[1][i + qb];
			else qmat[j] = QMIS;
		}
		qins = i + 1 == ql? QCLP : qvs[2][i + 1 + qb];
		if(i + 1 == ql) qdel[0] = qdel[1] = qdel[2] = qdel[3] = QCLP;
		else {
			for(j=0;j<4;j++){
				if(j == qvs[6][i + 1 + qb]){
					qdel[j] = qvs[3][i + 1 + qb];
				} else qdel[j] = QDEL;
			}
		}
		for(j=zb[i];j<ze[i];j++){
			m = rh[j] + qmat[target[j + tb]];
			rh[j] = h1;
			e = re[j];
			if(m <= e){ d = 0; h = m; }
			else      { d = 1; h = e; }
			if(h > f){ d = 2; h = f; }
			h1 = h;
			t = m + qins;
			//e = e + QEXT;
			e = e + qins;
			if(e < t){ d |= 1 << 2; }
			else     { e = t; }
			re[j] = e;
			t = m + qdel[target[j + tb]];
			f = f + QEXT;
			if(f < t){ d |= 2 << 4; }
			else     { f = t; }
			zi[j-zb[i]] = d;
		}
		rh[j] = h1; re[j] = SCORE_INF;
	}
	y.qb = qb; y.qe = qe; y.tb = tb; y.te = te;
	y.score = - (int)rh[tl];
	y.mat = y.mis = y.ins = y.del = y.aln = 0;
	d = 0;
	i = ql - 1; j = tl - 1;
	//String *qs = init_string(1024);
	//String *ts = init_string(1024);
	//String *ms = init_string(1024);
	while(i >= 0 && j >= 0){
		d = (z[i * rl + (j - zb[i])] >> (d << 1)) & 0x03;
		if(d == 0){
			//add_char_string(qs, bit_base_table[query[i + y.qb]]);
			//add_char_string(ts, bit_base_table[target[j + y.tb]]);
			if(query[i + y.qb] == target[j + y.tb]){
				//add_char_string(ms, '|');
				y.mat ++;
			} else {
				//add_char_string(ms, '*');
				y.mis ++;
			}
			i --; j --; 
		} else if(d == 1){
			//add_char_string(qs, bit_base_table[query[i + y.qb]]);
			//add_char_string(ts, '-');
			//add_char_string(ms, '-');
			i --; y.ins ++;
		} else {
			//add_char_string(qs, '-');
			//add_char_string(ts, bit_base_table[target[j + y.tb]]);
			//add_char_string(ms, '-');
			j --; y.del ++;
		}
		kswx_push_cigar(cigars2, d, 1);
	}
	//reverse_string(qs);
	//reverse_string(ts);
	//reverse_string(ms);
	//printf("%s\n%s\n%s\n", qs->string, ms->string, ts->string);
	//free_string(qs);
	//free_string(ts);
	//free_string(ms);
	if(1){
		if(i >= 0){
			y.ins += i + 1;
			kswx_push_cigar(cigars2, 1, i + 1);
		}
		if(j >= 0){
			y.del += j + 1;
			kswx_push_cigar(cigars2, 2, j + 1);
		}
	} else {
		y.qb += i + 1;
		y.tb += j + 1;
	}
	reverse_u32list(cigars2);
	y.aln = y.mat + y.mis + y.ins + y.del;
	return y;
}

static inline void kswx_cigar2string(String *str, int n, uint32_t *cigar){
	int i, j, c, op, len;
	char ch;
	for(i=0;i<n;i++){
		op = cigar[i] & 0xFU;
		len = cigar[i] >> 4;
		if(len == 0) continue;
		if(op > 2){
			fprintf(stderr, " -- CIGAR only support M(0),I(1),D(2) cigar, but met ?(%d) in %s -- %s:%d --\n", op, __FUNCTION__, __FILE__, __LINE__);
			exit(1);
		}
		c = 0;
		encap_string(str, 20);
		while(len){
			str->string[str->size + c] = '0' + (len % 10);
			len = len / 10;
			c ++;
		}
		for(j=0;j<c>>1;j++){
			ch = str->string[str->size + j];
			str->string[str->size + j] = str->string[str->size + c - 1 - j];
			str->string[str->size + c - 1 - j] = ch;
		}
		str->size += c;
		str->string[str->size++] = "MIDX"[op];
	}
	str->string[str->size] = '\0';
}

static inline uint32_t kswx_string2cigar(u32list *cigar, char *str){
	uint32_t c, len, op, n;
	char *p;
	p = str;
	len = 0;
	n = 0;
	while(*p){
		if(*p >= '0' && *p <= '9'){
			len = len * 10 + (*p) - '0';
		} else {
			switch(*p){
				case 'M': op = 0; break;
				case 'I': op = 1; break;
				case 'D': op = 2; break;
				default: op = 3;
			}
			c = (len << 4) | op;
			if(cigar->size && (cigar->buffer[cigar->size - 1] & 0xFU) == op){
				cigar->buffer[cigar->size - 1] += len << 4;
			} else push_u32list(cigar, c);
			n ++;
			len = 0;
		}
		p ++;
	}
	return n;
}

static inline void kswx_cigar2pairwise(String *alns[2], uint8_t *seqs[2], int n, uint32_t *cigar){
	char *p1, *p2;
	int i, j, op, len, off[2];
	off[0] = off[1] = 0;
	for(i=0;i<n;i++){
		op = cigar[i] & 0xFU;
		len = cigar[i] >> 4;
		if(len == 0) continue;
		encap_string(alns[0], len);
		encap_string(alns[1], len);
		p1 = alns[0]->string + alns[0]->size;
		p2 = alns[1]->string + alns[1]->size;
		switch(op){
			case 0:
			for(j=0;j<len;j++){
				p1[j] = bit_base_table[(int)seqs[0][off[0] + j]];
				p2[j] = bit_base_table[(int)seqs[1][off[1] + j]];
			}
			off[0] += len;
			off[1] += len;
			break;
			case 1:
			for(j=0;j<len;j++){
				p1[j] = bit_base_table[(int)seqs[0][off[0] + j]];
				p2[j] = '-';
			}
			off[0] += len;
			break;
			case 2:
			for(j=0;j<len;j++){
				p1[j] = '-';
				p2[j] = bit_base_table[(int)seqs[1][off[1] + j]];
			}
			off[1] += len;
			break;
			default:
			fprintf(stderr, " -- CIGAR only support M(0),I(1),D(2) cigar, but met ?(%d) in %s -- %s:%d --\n", op, __FUNCTION__, __FILE__, __LINE__);
			exit(1);
		}
		alns[0]->size += len;
		alns[1]->size += len;
	}
	alns[0]->string[alns[0]->size] = 0;
	alns[1]->string[alns[1]->size] = 0;
}

static inline void kswx_polish_pairwise(char *alns[2], int len){
	int i, j, m, c, gaps[2];
	while(1){
		c = 0;
		gaps[0] = gaps[1] = 0;
		for(i=0;i<len;i++){
			for(j=0;j<2;j++){
				if(alns[j][i] == '-'){ gaps[j] ++; continue; }
				if(gaps[j] == 0) continue;
				for(m=i-gaps[j];m<i;m++){
					if(alns[!j][m] == alns[j][i]){
						alns[j][m] = alns[j][i];
						alns[j][i] = '-';
						c ++;
						break;
					}
				}
				gaps[j] = i - m;
			}
		}
		if(c == 0) break;
	}
}

static inline uint32_t kswx_pairwise2cigar(u32list *cigar, char *alns[2], uint32_t size){
	uint32_t i, op, len, n, s;
	op = len = n = 0;
	for(i=0;i<=size;i++){
		if(i == n) s = 10;
		if(alns[0][i] == '-') s = 2;
		else if(alns[1][i] == '-') s = 1;
		else s = 0;
		if(op == s) len ++;
		else {
			push_u32list(cigar, (len << 4) | op);
			n ++;
			op = s;
			len = 1;
		}
	}
	return n;
}

#define KSWX_ZIGAR_ZERO	32
#define KSWX_ZIGAR_MAX_LEN	30

static inline void kswx_cigar2zigar(String *str, int n, uint32_t *cigar){
	int i, op, len, t, s, ss;
	t = s = 0;
	for(i=0;i<=n;i++){
		if(i < n){
			op = cigar[i] & 0xFU;
			if(op > 2){
				fprintf(stderr, " -- ZIGAR only support M(0),I(1),D(2) cigar, but met ?(%d) in %s -- %s:%d --\n", op, __FUNCTION__, __FILE__, __LINE__);
				exit(1);
			}
			len = cigar[i] >> 4;
			if(len == 0) continue;
		} else { op = 0x1FU; len = 0; }
		if(op == t) s += len;
		else {
			while(s){
				ss = (s <= KSWX_ZIGAR_MAX_LEN)? s : KSWX_ZIGAR_MAX_LEN;
				s -= ss;
				ss = ss + t * KSWX_ZIGAR_MAX_LEN + KSWX_ZIGAR_ZERO;
				add_char_string(str, ss);
			}
			t = op;
			s = len;
		}
	}
}

static inline uint32_t kswx_zigar2cigar(u32list *cigar, char *str){
	uint32_t c, len, op, s, t;
	char *p;
	p = str;
	op = 0x1FU;
	len = 0;
	c = 0;
	while(1){
		if(*p){
			s = (*p) - KSWX_ZIGAR_ZERO;
			t = s / KSWX_ZIGAR_MAX_LEN;
			s = s % KSWX_ZIGAR_MAX_LEN;
		} else {
			t = 0x1FU;
			s = 0;
		}
		if(op == t) len += s;
		else {
			if(len){
				push_u32list(cigar, (len << 4) | op);
				c ++;
			}
			op = t;
			len = s;
		}
		if(*p == 0) break;
		p ++;
	}
	return c;
}

static inline void kswx_zigar2string(String *str, char *zigar){
	int j, c, op, len, s, t;
	char ch, *p;
	p = zigar;
	op = 0x1FU;
	len = 0;
	while(1){
		if(*p){
			s = (*p) - KSWX_ZIGAR_ZERO;
			t = s / KSWX_ZIGAR_MAX_LEN;
			s = s % KSWX_ZIGAR_MAX_LEN;
		} else {
			t = 0x1FU;
			s = 0;
		}
		if(op == t) len += s;
		else {
			if(len){
				c = 0;
				encap_string(str, 20);
				while(len){
					str->string[str->size + c] = '0' + (len % 10);
					len = len / 10;
					c ++;
				}
				for(j=0;j<c>>1;j++){
					ch = str->string[str->size + j];
					str->string[str->size + j] = str->string[str->size + c - 1 - j];
					str->string[str->size + c - 1 - j] = ch;
				}
				str->size += c;
				str->string[str->size++] = "MIDX"[op];
			}
			op = t;
			len = s;
		}
		if(*p == '\0') break;
		p ++;
	}
	str->string[str->size] = '\0';
}

static inline int kswx_string2zigar(String *zigar, char *str){
	uint32_t len, op, n, t, s, ss;
	char *p;
	p = str;
	op = 0;
	len = 0;
	n = 0;
	t = s = 0;
	while(1){
		if(*p >= '0' && *p <= '9'){
			len = len * 10 + (*p) - '0';
			p ++;
			continue;
		}
		if(*p){
			switch(*p){
				case 'M': op = 0; break;
				case 'I': op = 1; break;
				case 'D': op = 2; break;
				default: op = 3;
					fprintf(stderr, " -- ZIGAR only support M(0),I(1),D(2) cigar, but met %c(?) in %s -- %s:%d --\n", *p, __FUNCTION__, __FILE__, __LINE__);
					fprintf(stderr, "%s\n", str); fflush(stderr);
					exit(1);
			}
		} else { op = 0x1FU; len = 0; }
		if(op == t){ s += len; len = 0; }
		else {
			while(s){
				ss = (s <= KSWX_ZIGAR_MAX_LEN)? s : KSWX_ZIGAR_MAX_LEN;
				s -= ss;
				ss = ss + t * KSWX_ZIGAR_MAX_LEN + KSWX_ZIGAR_ZERO;
				add_char_string(zigar, ss);
				n ++;
			}
			t = op;
			s = len;
			len = 0;
		}
		if(*p == '\0') break;
		p ++;
	};
	return n;
}

static inline kswr_t kswx_extend_core(int qlen, uint8_t *query, int tlen, uint8_t *target, kswr_t r, int m, const int8_t *matrix, int w, int I, int D, int E, int T){
	int y1, y2, x1, x2, x3, score, gscore, max_off;
	if(T < 0){
		// try extend left
		do {
			if(r.qb == 0 || r.tb == 0) break;
			if(r.tb >= r.qb){
				y1 = (r.qb + w > r.tb)? r.tb : r.qb + w;
				y2 = r.qb;
				revseq_bytes(target + r.tb - y1, y1); revseq_bytes(query  + r.qb - y2, y2);
				score = ksw_extend2(y2, query + r.qb - y2, y1, target + r.tb - y1, m, matrix, - D, - E, - I, - E, w, - T, -1, r.score, &x2, &x1, &x3, &gscore, &max_off);
				revseq_bytes(target + r.tb - y1, y1); revseq_bytes(query  + r.qb - y2, y2);
				if(gscore <= 0 || gscore <= score + T){
					r.tb = r.tb - x1; r.qb = r.qb - x2; r.score = score;
				} else {
					r.tb = r.tb - x3; r.qb = 0; r.score = gscore;
				}
			} else {
				y1 = r.tb;
				y2 = (r.tb + w > r.qb)? r.qb : r.tb + w;
				revseq_bytes(target + r.tb - y1, y1); revseq_bytes(query  + r.qb - y2, y2);
				score = ksw_extend2(y1, target + r.tb - y1, y2, query + r.qb - y2, m, matrix, - I, - E, - D, - E, w, - T, -1, r.score, &x1, &x2, &x3, &gscore, &max_off);
				revseq_bytes(target + r.tb - y1, y1); revseq_bytes(query  + r.qb - y2, y2);
				if(gscore <= 0 || gscore <= score + T){
					r.tb = r.tb - x1; r.qb = r.qb - x2; r.score = score;
				} else {
					r.qb = r.qb - x3; r.tb = 0; r.score = gscore;
				}
			}
		} while(0);
		// try extend right
		do {
			if(r.qe == qlen || r.te == tlen) break;
			if(tlen - r.te >= qlen - r.qe){
				y1 = (qlen - r.qe + w > tlen - r.te)? tlen - r.te : qlen - r.qe + w;
				y2 = qlen - r.qe;
				score = ksw_extend2(y2, query + r.qe, y1, target + r.te, m, matrix, - D, - E, - I, - E, w, - T, -1, r.score, &x2, &x1, &x3, &gscore, &max_off);
				if(gscore <= 0 || gscore <= score + T){
					r.te += x1; r.qe += x2; r.score = score;
				} else {
					r.te += x3; r.qe = qlen; r.score = gscore;
				}
			} else {
				y1 = tlen - r.te;
				y2 = (tlen - r.te + w > qlen - r.qe)? qlen - r.qe : tlen - r.te + w;
				score = ksw_extend2(y1, target + r.te, y2, query + r.qe, m, matrix, - I, - E, - D, - E, w, - T, -1, r.score, &x1, &x2, &x3, &gscore, &max_off);
				if(gscore <= 0 || gscore <= score + T){
					r.te += x1; r.qe += x2; r.score = score;
				} else {
					r.te = tlen; r.qe += x3; r.score = gscore;
				}
			}
		} while(0);
	}
	return r;
}

static inline kswx_t kswx_gen_cigar_core2(int qlen, uint8_t *query, int tlen, uint8_t *target, kswr_t r, int m, const int8_t *matrix, int w, int I, int D, int E, int *_n_cigar, uint32_t **_cigar){
	kswx_t x;
	uint32_t *cigar;
	int i, j, n_cigar, w2, x1, x2, op, len;
	UNUSED(qlen);
	UNUSED(tlen);
	n_cigar = *_n_cigar;
	cigar   = *_cigar;
	if(w >= 0){
		x1 = r.te - r.tb;
		x2 = r.qe - r.qb;
		x1 -= x2;
		if(x1 < 0) x1 = -x1;
		if(w < x1 * 1.5) w = x1 * 1.5;
		w2 = (((r.te - r.tb < r.qe - r.qb)? r.te - r.tb : r.qe - r.qb) * matrix[0] - r.score) / (-E) + 1;
		if(w2 < (w / 4) + 1) w2 = (w / 4) + 1;
		if(w > w2) w = w2;
	} else w = - w;
	x.score = ksw_global2(r.qe - r.qb, query + r.qb, r.te - r.tb, target + r.tb, m, matrix, - D, - E, - I, - E, w, &n_cigar, &cigar);
	x.tb = r.tb;
	x.te = r.te;
	x.qb = r.qb;
	x.qe = r.qe;
	x.aln = x.mat = x.mis = x.ins = x.del = 0;
	for(i=x1=x2=0;i<n_cigar;i++){
		op = cigar[i] & 0xF;
		len = cigar[i] >> 4;
		x.aln += len;
		switch(op){
			case 0: for(j=0;j<len;j++){ if(query[r.qb + x1 + j] == target[r.tb + x2 + j]) x.mat ++; else x.mis ++; } x1 += len; x2 += len; break;
			case 1: x1 += len; x.del += len; break;
			case 2: x2 += len; x.ins += len; break;
		}
	}
	*_n_cigar = n_cigar;
	*_cigar   = cigar;
	return x;
}

static inline kswx_t kswx_gen_cigar_core(int qlen, uint8_t *query, int tlen, uint8_t *target, kswr_t r, int m, const int8_t *matrix, int w, int I, int D, int E){
	kswx_t x;
	uint32_t *cigar;
	int n_cigar;
	n_cigar = 0;
	cigar = NULL;
	x = kswx_gen_cigar_core2(qlen, query, tlen, target, r, m, matrix, w, I, D, E, &n_cigar, &cigar);
	if(cigar) free(cigar);
	return x;
}

static inline kswx_t kswx_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *matrix, int w, int I, int D, int E, int T){
	kswr_t r;
	r = ksw_align2(qlen, query, tlen, target, m, matrix, - D, - E, - I, - E, KSW_XSTART, NULL);
	if(r.qb <= -1 || r.tb <= -1 || r.qe <= -1 || r.te <= -1) return KSWX_NULL;
	r.qe ++; r.te ++;
	r = kswx_extend_core(qlen, query, tlen, target, r, m, matrix, w, I, D, E, T);
	return kswx_gen_cigar_core(qlen, query, tlen, target, r, m, matrix, w, I, D, E);
}

static inline kswr_t kswx_align_no_stat(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *matrix, int w, int I, int D, int E, int T){
	kswr_t r;
	r = ksw_align2(qlen, query, tlen, target, m, matrix, - D, - E, - I, - E, KSW_XSTART, NULL);
	if(r.qb <= -1 || r.tb <= -1 || r.qe <= -1 || r.te <= -1) return KSWR_NULL;
	r.qe ++; r.te ++;
	r = kswx_extend_core(qlen, query, tlen, target, r, m, matrix, w, I, D, E, T);
	return r;
}

static inline kswx_t kswx_align_with_cigar(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *matrix, int w, int I, int D, int E, int T, int *n_cigar, uint32_t **cigar){
	kswr_t r;
	r = ksw_align2(qlen, query, tlen, target, m, matrix, - D, - E, - I, - E, KSW_XSTART, NULL);
	if(r.qb <= -1 || r.tb <= -1 || r.qe <= -1 || r.te <= -1) return KSWX_NULL;
	r.qe ++; r.te ++;
	r = kswx_extend_core(qlen, query, tlen, target, r, m, matrix, w, I, D, E, T);
	return kswx_gen_cigar_core2(qlen, query, tlen, target, r, m, matrix, w, I, D, E, n_cigar, cigar);
}

// qb, qe, tb, te are the positions of seed match region, will be extended to whole reads
static inline kswx_t kswx_fast_align(int qlen, uint8_t *query, int qb, int qe, int tlen, uint8_t *target, int tb, int te, int W, int M, int X, int I, int D, int E, int T, float min_sm, String *cigar_str){
	kswr_t r;
	kswx_t x, x1;
	int8_t i, m, matrix[16];
	int ol, n_cigar[3], w;
	uint32_t *cigar[3];
	m = 4;
	for(i=0;i<16;i++) matrix[i] = ((i / 4) == (i % 4))? M : X;
	r = ksw_align2(qe - qb, query + qb, te - tb, target + tb, m, matrix, - D, - E, - I, - E, KSW_XSTART, NULL);
	if(r.qb <= -1 || r.tb <= -1 || r.qe <= -1 || r.te <= -1) return KSWX_NULL;
	r.qe ++; r.te ++;
	// check whether the seed region was well aligned
	ol = (r.qe - r.qb < r.te - r.tb)? r.qe - r.qb : r.te - r.tb;
	if(ol < 30) return KSWX_NULL;
	n_cigar[0] = n_cigar[1] = n_cigar[2] = 0;
	cigar[0] = cigar[1] = cigar[2] = NULL;
	w = (r.qe - r.qb < r.te - r.tb)? ((r.te - r.tb) - (r.qe - r.qb)) : ((r.qe - r.qb) - (r.te - r.tb));
	if(w < W) w = W;
	if(cigar_str) x = kswx_gen_cigar_core2(qe - qb, query + qb, te - tb, target + tb, r, m, matrix, w, I, D, E, n_cigar + 1, cigar + 1);
	else          x = kswx_gen_cigar_core2(qe - qb, query + qb, te - tb, target + tb, r, m, matrix, w, I, D, E, NULL, NULL);
	if(x.mat < x.aln * min_sm){
		if(cigar[1]) free(cigar[1]);
		return KSWX_NULL;
	}
	x.qb += qb; x.qe += qb;
	x.tb += tb; x.te += tb;
	if(cigar_str) x1 = kswx_extend_align(qlen - x.qe, query + x.qe, tlen - x.te, target + x.te, 1, x.score, W, M, X, I, D, E, T, n_cigar + 2, cigar + 2);
	else          x1 = kswx_extend_align(qlen - x.qe, query + x.qe, tlen - x.te, target + x.te, 1, x.score, W, M, X, I, D, E, T, NULL, NULL);
	x.score = x1.score;
	x.aln   += x1.aln;
	x.mat   += x1.mat;
	x.mis   += x1.mis;
	x.ins   += x1.ins;
	x.del   += x1.del;
	x.qe    += x1.qe;
	x.te    += x1.te;
	if(cigar_str) x1 = kswx_extend_align(x.qb, query + x.qb - 1, x.tb, target + x.tb - 1, -1, x.score, W, M, X, I, D, E, T, n_cigar + 0, cigar + 0);
	else          x1 = kswx_extend_align(x.qb, query + x.qb - 1, x.tb, target + x.tb - 1, -1, x.score, W, M, X, I, D, E, T, NULL, NULL);
	x.score = x1.score;
	x.aln   += x1.aln;
	x.mat   += x1.mat;
	x.mis   += x1.mis;
	x.ins   += x1.ins;
	x.del   += x1.del;
	x.qb    -= x1.qe;
	x.tb    -= x1.te;
	if(cigar_str){
		clear_string(cigar_str);
		revseq_4bytes(cigar[0], n_cigar[0]);
		kswx_cigar2string(cigar_str, n_cigar[0], cigar[0]);
		kswx_cigar2string(cigar_str, n_cigar[1], cigar[1]);
		kswx_cigar2string(cigar_str, n_cigar[2], cigar[2]);
		if(cigar[0]) free(cigar[0]);
		if(cigar[1]) free(cigar[1]);
		if(cigar[2]) free(cigar[2]);
	}
	return x;
}

#endif
