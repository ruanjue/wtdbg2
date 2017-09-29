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

#include "list.h"
#include "hashset.h"
#include "dna.h"
#include "file_reader.h"
#include "bitvec.h"
#include "thread.h"

#ifndef __DOT_MATRIX_OVERLAPPER_RJ_H
#define __DOT_MATRIX_OVERLAPPER_RJ_H

#define hzm_debug_out stderr

#define HZMH_KMER_MOD	1024
#define hzmh_kmer_smear(K) ((K) ^ ((K) >> 4) ^ ((K) >> 7) ^ ((K) >> 12))

#define DMO_MAX_KSIZE	32
#define DMO_MAX_KLEN	64
#define DMO_MAX_RD_LEN	0x0007FFFF

#define hzm_debug_out stderr
static int hzm_debug = 0;

typedef struct {
	uint32_t rd_id;
	uint32_t dir:1, off:19, len:6, extra:6;
} hzm_t;
define_list(hzmv, hzm_t);

typedef struct {
	uint32_t mer;
	uint32_t dir:1, off:19, len:6, extra:6;
} hzmr_t;
define_list(hzmrv, hzmr_t);

typedef struct {
	uint64_t mer;
	uint64_t off:47, cnt:16, flt:1;
} hzmh_t;
define_list(hzmhv, hzmh_t);
#define hzmh_hashcode(E) u64hashcode((E).mer)
#define hzmh_hashequals(E1, E2) ((E1).mer == (E2).mer)
define_hashset(hzmhash, hzmh_t, hzmh_hashcode, hzmh_hashequals);

typedef struct {
	hzm_t p1;
	hzm_t *p2;
	uint64_t beg, end;
} hz_ref_t;
define_list(hzrefv, hz_ref_t);
#define tmp_util1(r) ((r)->p2->rd_id)
#define heap_cmp_hz_ref_macro(refs, a, b)	\
({	\
	hz_ref_t *r1, *r2;	\
	r1 = ref_hzrefv(refs, a);	\
	r2 = ref_hzrefv(refs, b);	\
	(tmp_util1(r1) > tmp_util1(r2))? 1 : ((tmp_util1(r1) < tmp_util1(r2))? -1 : ((r1->p1.off > r2->p1.off)? 1 : ((r1->p1.off < r2->p1.off)? -1 : 0)));	\
})

typedef struct {
	uint32_t dir1:1, off1:31;
	uint32_t dir2:1, off2:31;
	uint32_t len1:16, len2:16;
	uint32_t gid; // group id
} hzmp_t;
define_list(hzmpv, hzmp_t);

typedef struct {
	uint32_t pb2;
	uint32_t ovl:29, dir:1, closed:2;
	int beg[2], end[2];
	uint32_t anchors[2];
	//uint32_t kmat;
} wt_seed_t;
define_list(wtseedv, wt_seed_t);

typedef struct {
	int offset;
	uint32_t off, cnt;
} diag_t;
define_list(diagv, diag_t);

typedef struct {
	uint32_t off1, off2;
} kigar_t;
define_list(kigarv, kigar_t);

typedef struct {
	uint32_t pb1:31, dir1:1, pb2:31, dir2:1;
	int qb, qe, tb, te;
	int score, mat, aln;
	uint32_t kigar_off, kigar_len;
} wt_ovl_t;
define_list(wtovlv, wt_ovl_t);
#define WT_OVL_NULL_SCORE	19830203

#define _ovl_uniq_long_id(id1, id2, dir) ((((uint64_t)(id1)) << 33) | (((uint64_t)(id2)) << 1) | (dir))
#define ovl_uniq_long_id(id1, id2, dir) (((id1) < (id2))? _ovl_uniq_long_id(id1, id2, dir) : _ovl_uniq_long_id(id2, id1, dir))

typedef struct {
	uint64_t rdoff:44, rdlen:19, closed:1;
	u4i rdidx, gpidx; // rdidx records the original index, gpidx records the read group index
	u8i ttr_beg:19, ttr_end:19, ttr_idx:26;
	char *tag;
} pbread_t;
define_list(pbreadv, pbread_t);

typedef struct {
	int hk;
	uint32_t hzmh_kmer_mod, hzmh_kmer_win;
	uint32_t ksize, kmax, kmin, kgap; // kgap: max skip bases in the middle of kmer
	int strand_mask; // 1: forward, 2: reverse, 3: both. Now only supported by query_reg_dmo
	float    ktop; // when kmax == 0, mask top <ktop> high frequency kmers as repetitive
	uint32_t max_hit;
	int      min_aln, min_mat;
	float    min_sm, aln_var;
	int      xvar, yvar, zmin, max_overhang;
	float    deviation_penalty, gap_penalty;
	int      ttr_len, ttr_cnt, ttr_win;
	float    ttr_drate, ttr_crate, ttr_xrate;
} DMOPar;

typedef struct {
	uint32_t n_rd;
	int len_order, index_ttr;
	u8i nbp;
	BaseBank *rdseqs;
	pbreadv  *reads;
	cuhash   *tag2idx;
	hzmv     *seeds;
	hzmhash  *hashs[HZMH_KMER_MOD];
	DMOPar   *par;
	//TTR data
	hzmhash  *ttrhs[HZMH_KMER_MOD];
	u4v      *ttrvs;
} DMO;

typedef struct {
	hzrefv *refs;
	u32list *heap, *hzoff;
	hzmpv *rs, *dst[2];
	wtseedv *regs[4];
	diagv *diags;
	u4v *block, *grps;
	u8list *mem;
	u8i trials;
} DMOAux;

static inline DMOPar* init_dmopar(uint32_t ksize, uint32_t kgap, int ksave, int wsize, int xvar, int yvar, int zmin, int max_overhang, float deviation_penalty, float gap_penalty){
	DMOPar *par;
	par = malloc(sizeof(DMOPar));
	par->hzmh_kmer_mod = HZMH_KMER_MOD * (ksave < 1? 1 : ksave);
	par->hzmh_kmer_win = wsize < 1? 1 : wsize;
	par->ksize = ksize;
	par->kmax  = 0;
	par->ktop  = 0.05;
	par->kmin  = 2;
	par->kgap  = kgap;
	par->hk    = 1;
	par->max_hit = 0;
	par->min_aln = 1000;
	par->min_mat = 50;
	par->min_sm  = 0.05;
	par->aln_var = 0.5;
	par->xvar = xvar;
	par->yvar = yvar;
	par->zmin = zmin;
	par->strand_mask = 3;
	par->max_overhang = max_overhang;
	par->deviation_penalty = deviation_penalty;
	par->gap_penalty = gap_penalty;
	par->ttr_len = 5000;
	par->ttr_cnt = 2;
	par->ttr_win = 100;
	par->ttr_drate = 0.05;
	par->ttr_crate = 0.2;
	par->ttr_xrate = 0.8;
	return par;
}

static inline void free_dmopar(DMOPar *par){ free(par); }

static inline DMO* init_dmo(DMOPar *par){
	DMO *wt;
	int i;
	wt = malloc(sizeof(DMO));
	wt->rdseqs = init_basebank();
	wt->reads  = init_pbreadv(64);
	wt->tag2idx = init_cuhash(1023);
	wt->seeds = init_hzmv(64);
	for(i=0;i<HZMH_KMER_MOD;i++){
		wt->hashs[i] = init_hzmhash(1023);
	}
	wt->ttrvs = init_u4v(64);
	for(i=0;i<HZMH_KMER_MOD;i++){
		wt->ttrhs[i] = init_hzmhash(1023);
	}
	wt->n_rd = 0;
	wt->nbp = 0;
	wt->len_order = 1;
	wt->index_ttr = 0;
	wt->par = par;
	return wt;
}

static inline void reset_index_dmo(DMO *wt){
	uint32_t i;
	for(i=0;i<HZMH_KMER_MOD;i++){
		free_hzmhash(wt->hashs[i]);
		wt->hashs[i] = init_hzmhash(1023);
	}
	free_hzmv(wt->seeds);
	wt->seeds = init_hzmv(64);
}

static inline void clear_index_dmo(DMO *wt){
	uint32_t i;
	for(i=0;i<HZMH_KMER_MOD;i++){
		clear_hzmhash(wt->hashs[i]);
	}
	clear_hzmv(wt->seeds);
}

static inline void free_dmo(DMO *wt){
	uint32_t i;
	free_basebank(wt->rdseqs);
	for(i=0;i<wt->n_rd;i++){
		if(wt->reads->buffer[i].tag) free(wt->reads->buffer[i].tag);
	}
	free_pbreadv(wt->reads);
	free_cuhash(wt->tag2idx);
	for(i=0;i<HZMH_KMER_MOD;i++){
		free_hzmhash(wt->hashs[i]);
	}
	free_hzmv(wt->seeds);
	for(i=0;i<HZMH_KMER_MOD;i++){
		free_hzmhash(wt->ttrhs[i]);
	}
	free_u4v(wt->ttrvs);
	free(wt);
}

static inline void push_long_read_adv_dmo(DMO *wt, u4i gpidx, char *name, int name_len, char *seq, int seq_len){
	char *ptr;
	if(name_len){
		ptr = malloc(name_len + 1);
		memcpy(ptr, name, name_len);
		ptr[name_len] = 0;
	} else {
		ptr = NULL;
	}
	push_pbreadv(wt->reads, (pbread_t){wt->rdseqs->size, seq_len, 0, wt->n_rd, gpidx, 0, 0, 0, ptr});
	seq2basebank(wt->rdseqs, seq, seq_len);
	wt->n_rd ++;
	wt->nbp += seq_len;
}

static inline void push_long_read_dmo(DMO *wt, char *name, int name_len, char *seq, int seq_len){
	push_long_read_adv_dmo(wt, wt->n_rd, name, name_len, seq, seq_len);
}

void ready_dmo(DMO *wt){
	uint32_t i;
	if(wt->len_order) sort_array(wt->reads->buffer, wt->reads->size, pbread_t, b.rdlen > a.rdlen);
	for(i=0;i<wt->reads->size;i++){
		if(wt->reads->buffer[i].tag) put_cuhash(wt->tag2idx, (cuhash_t){wt->reads->buffer[i].tag, i});
	}
}

//TODO: unfinished
/*
typedef struct {
	u4i beg, end, den, cov, flag, mask;
} ttr_reg_t;
define_list(ttrregv, ttr_reg_t);

int detect_TTR_dmo(DMO *wt, u4i rid, hzmhv *kmers, u4v *dens, u4v *covs, u4v *brks, ttrregv *regs){
}
*/

typedef struct {
	u8i kval;
	u4i dir:1, kidx:31;
	int off;
} kmin_t;
define_list(kminv, kmin_t);

typedef struct {
	uint64_t mer:54, kidx:10;
} midx_t;
define_list(midxv, midx_t);

typedef struct {
	midx_t m;
	hzm_t  h;
} hzmm_t;
define_list(hzmmv, hzmm_t);

thread_beg_def(midx);
DMO *wt;
uint32_t beg, end;
midxv *mcache;
hzmmv *hcache;
uint64_t ktot, nrem, Nrem, none, nflt, offset;
int filter_rep_kmer;
int task;
thread_end_def(midx);

thread_beg_func(midx);
DMO *wt;
hzmh_t *u, U;
pbread_t *rd;
midxv *mcache;
hzmmv *hcache;
midx_t *m;
hzmm_t *h;
kminv *kmins;
uint64_t kmer, krev, kmin, ktmp, kmask, off, kidx;
uint32_t beg, end, pbid, pblen, i, j, ncpu, tidx, w;
uint8_t b, c, dir;
int exists;
wt = midx->wt;
ncpu = midx->n_cpu;
tidx = midx->t_idx;
kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - wt->par->ksize) << 1);
memset(&U, 0, sizeof(hzmh_t));
kmins = init_kminv(8);
thread_beg_loop(midx);
beg = midx->beg;
end = midx->end;
if(midx->task == 1){
	fprintf(stderr, " -- Obsolete code in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	for(pbid=beg;pbid<end;pbid++){
		if(hzm_debug == 0 && tidx == 0 && ((pbid - beg) % 1000) == 0){ fprintf(hzm_debug_out, "\r%u", pbid - beg); fflush(hzm_debug_out); }
		rd = ref_pbreadv(wt->reads, pbid);
		if(rd->closed) continue;
		pblen = rd->rdlen;
		b = 4;
		kmer = 0;
		kmin = kmask;
		w = 0;
		off = rd->rdoff;
		for(i=0;i<pblen;i++){
			c = bits2bit(wt->rdseqs->bits, off); off ++;
			if(wt->par->hk && c == b) continue;
			b = c;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < wt->par->ksize) continue;
			krev = dna_rev_seq(kmer, wt->par->ksize);
			if(krev == kmer) continue;
			dir  = krev > kmer? 0 : 1;
			krev = krev > kmer? kmer : krev;
			kidx = hzmh_kmer_smear(krev) % wt->par->hzmh_kmer_mod;
			if(kidx >= HZMH_KMER_MOD) continue;
			ktmp = invertible_hashcode(krev, wt->par->ksize * 2);
			if(ktmp < kmin){
				clear_kminv(kmins);
				push_kminv(kmins, (kmin_t){krev, dir, kidx, i});
			} else if(ktmp == kmin){
				push_kminv(kmins, (kmin_t){krev, dir, kidx, i});
			}
			w ++;
			if(w < wt->par->hzmh_kmer_win && i + 1 < pblen) continue;
			w = 0; kmin = kmask;
			for(j=0;j<kmins->size;j++){
				krev = kmins->buffer[j].kval;
				kidx = kmins->buffer[j].kidx;
				if((kidx % ncpu) != tidx) continue;
				U.mer = krev;
				u = prepare_hzmhash(wt->hashs[kidx], U, &exists);
				if(exists){
					if(u->cnt < 0xFFFFU) u->cnt ++;
				} else {
					u->mer = krev;
					u->off = 0;
					u->cnt = 1;
					u->flt = 0;
				}
			}
			clear_kminv(kmins);
		}
	}
	if(hzm_debug == 0 && tidx == 0){ fprintf(hzm_debug_out, "\r%u reads\n", end - beg); fflush(hzm_debug_out); }
} else if(midx->task == 2){
	fprintf(stderr, " -- Obsolete code in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	for(pbid=beg;pbid<end;pbid++){
		if(hzm_debug == 0 && tidx == 0 && ((pbid - beg) % 1000) == 0){ fprintf(hzm_debug_out, "\r%u", pbid - beg); fflush(hzm_debug_out); }
		rd = ref_pbreadv(wt->reads, pbid);
		if(rd->closed) continue;
		pblen = rd->rdlen;
		b = 4;
		kmer = 0;
		kmin = kmask;
		w = 0;
		off = rd->rdoff;
		for(i=0;i<pblen;i++){
			c = bits2bit(wt->rdseqs->bits, off); off ++;
			if(wt->par->hk && c == b) continue;
			b = c;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < wt->par->ksize) continue;
			krev = dna_rev_seq(kmer, wt->par->ksize);
			if(krev == kmer) continue;
			dir  = krev > kmer? 0 : 1;
			krev = krev > kmer? kmer : krev;
			kidx = hzmh_kmer_smear(krev) % wt->par->hzmh_kmer_mod;
			if(kidx >= HZMH_KMER_MOD) continue;
			ktmp = invertible_hashcode(krev, wt->par->ksize * 2);
			if(ktmp < kmin){
				clear_kminv(kmins);
				push_kminv(kmins, (kmin_t){krev, dir, kidx, i});
			} else if(ktmp == kmin){
				push_kminv(kmins, (kmin_t){krev, dir, kidx, i});
			}
			w ++;
			if(w < wt->par->hzmh_kmer_win && i + 1 < pblen) continue;
			w = 0; kmin = kmask;
			for(j=0;j<kmins->size;j++){
				krev = kmins->buffer[j].kval;
				kidx = kmins->buffer[j].kidx;
				if((kidx % ncpu) != tidx) continue;
				U.mer = krev;
				u = get_hzmhash(wt->hashs[kidx], U);
				if(u == NULL || u->flt || u->cnt >= wt->par->kmax) continue; // too high or singleton
				wt->seeds->buffer[u->off + u->cnt] = (hzm_t){pbid, kmins->buffer[j].dir, kmins->buffer[j].off - wt->par->ksize, wt->par->ksize, 0};
				u->cnt ++;
			}
			clear_kminv(kmins);
		}
	}
	if(hzm_debug == 0 && tidx == 0) fprintf(hzm_debug_out, "\r%u reads\n", pbid - beg); fflush(hzm_debug_out);
} else if(midx->task == 3){
	for(pbid=beg;pbid<end;pbid++){
		rd = ref_pbreadv(wt->reads, pbid);
		if(rd->closed) continue;
		pblen = rd->rdlen;
		b = 4;
		kmer = 0;
		kmin = kmask;
		w = 0;
		off = rd->rdoff;
		for(i=0;i<pblen;i++){
			c = bits2bit(wt->rdseqs->bits, off); off ++;
			if(wt->par->hk && c == b) continue;
			b = c;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < wt->par->ksize) continue;
			krev = dna_rev_seq(kmer, wt->par->ksize);
			if(krev == kmer) continue;
			dir  = krev > kmer? 0 : 1;
			krev = krev > kmer? kmer : krev;
			kidx = hzmh_kmer_smear(krev) % wt->par->hzmh_kmer_mod;
			if(kidx >= HZMH_KMER_MOD) continue;
			ktmp = invertible_hashcode(krev, wt->par->ksize * 2);
			if(ktmp < kmin){
				clear_kminv(kmins);
				push_kminv(kmins, (kmin_t){krev, dir, kidx, i});
			} else if(ktmp == kmin){
				push_kminv(kmins, (kmin_t){krev, dir, kidx, i});
			}
			w ++;
			if(w < wt->par->hzmh_kmer_win && i + 1 < pblen) continue;
			w = 0; kmin = kmask;
			for(j=0;j<kmins->size;j++){
				krev = kmins->buffer[j].kval;
				kidx = kmins->buffer[j].kidx;
				push_midxv(midx->mcache, (midx_t){krev, kidx});
			}
			clear_kminv(kmins);
		}
	}
	sort_array(midx->mcache->buffer, midx->mcache->size, midx_t, num_cmpgt(a.mer, b.mer));
} else if(midx->task == 4){
	for(pbid=beg;pbid<end;pbid++){
		rd = ref_pbreadv(wt->reads, pbid);
		if(rd->closed) continue;
		pblen = rd->rdlen;
		b = 4;
		kmer = 0;
		kmin = kmask;
		w = 0;
		off = rd->rdoff;
		for(i=0;i<pblen;i++){
			c = bits2bit(wt->rdseqs->bits, off); off ++;
			if(wt->par->hk && c == b) continue;
			b = c;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < wt->par->ksize) continue;
			krev = dna_rev_seq(kmer, wt->par->ksize);
			if(krev == kmer) continue;
			dir  = krev > kmer? 0 : 1;
			krev = krev > kmer? kmer : krev;
			kidx = hzmh_kmer_smear(krev) % wt->par->hzmh_kmer_mod;
			if(kidx >= HZMH_KMER_MOD) continue;
			ktmp = invertible_hashcode(krev, wt->par->ksize * 2);
			if(ktmp < kmin){
				clear_kminv(kmins);
				push_kminv(kmins, (kmin_t){krev, dir, kidx, i});
			} else if(ktmp == kmin){
				push_kminv(kmins, (kmin_t){krev, dir, kidx, i});
			}
			w ++;
			if(w < wt->par->hzmh_kmer_win && i + 1 < pblen) continue;
			w = 0; kmin = kmask;
			for(j=0;j<kmins->size;j++){
				krev = kmins->buffer[j].kval;
				kidx = kmins->buffer[j].kidx;
				if(kidx >= HZMH_KMER_MOD) continue;
				push_hzmmv(midx->hcache, (hzmm_t){(midx_t){krev, kidx}, (hzm_t){pbid, kmins->buffer[j].dir, kmins->buffer[j].off - wt->par->ksize, wt->par->ksize, 0}});
			}
			clear_kminv(kmins);
		}
	}
	sort_array(midx->hcache->buffer, midx->hcache->size, hzmm_t, num_cmpgt(a.m.mer, b.m.mer));
} else if(midx->task == 5){
	for(i=0;i<ncpu;i++){
		mcache = midx->midx_params[i].mcache;
		u = NULL;
		for(off=0;off<mcache->size;off++){
			m = ref_midxv(mcache, off);
			if((m->kidx % ncpu) != tidx) continue;
			if(u && m->mer == u->mer){
				if(midx->filter_rep_kmer) continue;
				if(u->cnt < 0xFFFFU) u->cnt ++;
				continue;
			}
			U.mer = m->mer;
			u = prepare_hzmhash(wt->hashs[m->kidx], U, &exists);
			if(exists){
				if(u->cnt < 0xFFFFU) u->cnt ++;
			} else {
				u->mer = m->mer;
				u->off = 0;
				u->cnt = 1;
				u->flt = 0;
			}
		}
	}
} else if(midx->task == 6){
	for(i=0;i<ncpu;i++){
		hcache = midx->midx_params[i].hcache;
		u = NULL;
		for(off=0;off<hcache->size;off++){
			h = ref_hzmmv(hcache, off);
			if((h->m.kidx % ncpu) != tidx) continue;
			if(u && h->m.mer == u->mer){
				if(midx->filter_rep_kmer) continue;
			} else {
				U.mer = h->m.mer;
				u = get_hzmhash(wt->hashs[h->m.kidx], U);
				if(u == NULL) continue;
			}
			if(u->flt || u->cnt >= wt->par->kmax) continue; // too high or singleton
			wt->seeds->buffer[u->off + u->cnt] = h->h;
			u->cnt ++;
		}
	}
} else if(midx->task == 7){
	midx->offset = 0;
	midx->ktot = midx->nrem = midx->Nrem = midx->none = midx->nflt = 0;
	for(i=tidx;i<HZMH_KMER_MOD;i+=ncpu){
		reset_iter_hzmhash(wt->hashs[i]);
		while((u = ref_iter_hzmhash(wt->hashs[i]))){
			midx->ktot += u->cnt;
			if(u->cnt < wt->par->kmin) midx->none ++;
			else if(u->cnt > wt->par->kmax){ u->cnt = 0; u->flt = 1; midx->nflt ++; }
			{
				u->off = midx->offset;
				midx->offset += u->cnt;
				if(u->cnt < wt->par->kmin) u->flt = 1;
				else { midx->nrem ++; midx->Nrem += u->cnt; }
				u->cnt = 0;
			}
		}
	}
} else if(midx->task == 8){
	for(i=tidx;i<HZMH_KMER_MOD;i+=ncpu){
		reset_iter_hzmhash(wt->hashs[i]);
		while((u = ref_iter_hzmhash(wt->hashs[i]))){ u->off += midx->offset; }
	}
}
thread_end_loop(midx);
free_kminv(kmins);
thread_end_func(midx);

thread_beg_def(msrt);
DMO *wt;
thread_end_def(msrt);

thread_beg_func(msrt);
DMO *wt;
hzmhash *hash;
hzmh_t *e;
hzm_t  *hs;
uint32_t ncpu, tidx, i;
wt = msrt->wt;
ncpu = msrt->n_cpu;
tidx = msrt->t_idx;
thread_beg_loop(msrt);
for(i=tidx;i<HZMH_KMER_MOD;i+=ncpu){
	hash = wt->hashs[i];
	reset_iter_hzmhash(hash);
	while((e = ref_iter_hzmhash(hash))){
		if(e->cnt <= 1) continue;
		hs = wt->seeds->buffer + e->off;
		//sort_array(hs, e->cnt, hzm_t, num_cmpgt(a.rd_id, b.rd_id));
		sort_array(hs, e->cnt, hzm_t, num_cmpgtx(a.rd_id, b.rd_id, a.off, b.off));
	}
}
thread_end_loop(msrt);
thread_end_func(msrt);

static inline void index_dmo(DMO *wt, uint32_t beg, uint32_t end, int filter_rep_kmer, uint32_t ncpu){
	hzmh_t  *u;
	unsigned long long ktyp, nflt, nrem, Nrem, none, ktot, off, cnt, totlen, *kcnts, MAX;
	uint32_t kavg, i, b, e, batch_size;
	thread_preprocess(midx);
	thread_preprocess(msrt);
	batch_size = 100;
	totlen = 0;
	clear_hzmv(wt->seeds);
	for(i=0;i<HZMH_KMER_MOD;i++) clear_hzmhash(wt->hashs[i]);
	thread_beg_init(midx, ncpu);
	midx->wt = wt;
	midx->mcache = init_midxv(32);
	midx->hcache = init_hzmmv(32);
	midx->beg = beg;
	midx->end = beg;
	midx->task = 0;
	midx->filter_rep_kmer = filter_rep_kmer;
	thread_end_init(midx);
	fprintf(hzm_debug_out, "[%s] - scanning kmers (%d bp, subsampling 1/%d)\n", date(), wt->par->ksize, wt->par->hzmh_kmer_mod / HZMH_KMER_MOD);
	b = e = beg;
	while(b < end){
		thread_beg_iter(midx);
		midx->task = 3;
		midx->beg = b;
		e = num_min(b + batch_size, end);
		midx->end = e;
		b = e;
		thread_wake(midx);
		thread_end_iter(midx);
		thread_wait_all(midx);
		thread_beg_iter(midx);
		midx->task = 5;
		thread_wake(midx);
		thread_end_iter(midx);
		thread_wait_all(midx);
		thread_beg_iter(midx);
		clear_midxv(midx->mcache);
		thread_end_iter(midx);
		if(hzm_debug == 0){ fprintf(hzm_debug_out, "\r%u", e - beg); fflush(hzm_debug_out); }
	}
	if(hzm_debug == 0) fprintf(hzm_debug_out, "\r%u reads\n", end - beg); fflush(hzm_debug_out);
	// estimate kmer_freq_cutoff
	if(wt->par->kmax < 2){
		MAX = 10000;
		kcnts = calloc(MAX, sizeof(u8i));
		ktot = ktyp = 0;
		for(i=0;i<HZMH_KMER_MOD;i++){
			reset_iter_hzmhash(wt->hashs[i]);
			while((u = ref_iter_hzmhash(wt->hashs[i]))){
				ktot += u->cnt;
				kcnts[num_min(MAX - 1, u->cnt)] ++;
				if(u->cnt > 100) ktyp ++;
			}
		}
		ktyp = ktyp * wt->par->ktop;
		for(i=MAX;i>0;i--){
			if(kcnts[i-1] < ktyp){
				ktyp -= kcnts[i-1];
			} else break;
		}
		wt->par->kmax = i;
		free(kcnts);
		fprintf(hzm_debug_out, "[%s] - high frequency kmer depth is set to %d\n", date(), wt->par->kmax);
	}
	ktot = nrem = Nrem = none = nflt = ktyp = 0;
	off = 0;
	thread_apply_all(midx, midx->task = 7);
	thread_beg_iter(midx);
	ktot += midx->ktot;
	nrem += midx->nrem;
	Nrem += midx->Nrem;
	none += midx->none;
	nflt += midx->nflt;
	cnt = midx->offset;
	midx->offset = off;
	off += cnt;
	midx->task = 8;
	thread_wake(midx);
	thread_end_iter(midx);
	thread_wait_all(midx);
	for(i=0;i<HZMH_KMER_MOD;i++){ ktyp += wt->hashs[i]->count; }
	clear_and_encap_hzmv(wt->seeds, off);
	wt->seeds->size = off;
	kavg = ktot / (ktyp + 1);
	fprintf(hzm_debug_out, "[%s] - Total kmers = %llu\n", date(), ktyp);
	fprintf(hzm_debug_out, "[%s] - average kmer depth = %d\n", date(), kavg);
	fprintf(hzm_debug_out, "[%s] - %llu low frequency kmers (<%d)\n", date(), none, wt->par->kmin);
	fprintf(hzm_debug_out, "[%s] - %llu high frequency kmers (>%d)\n", date(), nflt, wt->par->kmax);
	fprintf(hzm_debug_out, "[%s] - indexing %llu kmers, %llu instances\n", date(), nrem, Nrem);
	b = e = beg;
	while(b < end){
		thread_beg_iter(midx);
		midx->task = 4;
		midx->beg = b;
		e = num_min(b + batch_size, end);
		midx->end = e;
		b = e;
		thread_wake(midx);
		thread_end_iter(midx);
		thread_wait_all(midx);
		thread_beg_iter(midx);
		midx->task = 6;
		thread_wake(midx);
		thread_end_iter(midx);
		thread_wait_all(midx);
		thread_beg_iter(midx);
		clear_hzmmv(midx->hcache);
		thread_end_iter(midx);
		if(hzm_debug == 0){ fprintf(hzm_debug_out, "\r%u", e - beg); fflush(hzm_debug_out); }
	}
	if(hzm_debug == 0) fprintf(hzm_debug_out, "\r%u reads\n", end - beg); fflush(hzm_debug_out);
	thread_beg_close(midx);
	free_midxv(midx->mcache);
	free_hzmmv(midx->hcache);
	thread_end_close(midx);
	thread_beg_init(msrt, ncpu);
	msrt->wt = wt;
	thread_end_init(msrt);
	thread_wake_all(msrt);
	thread_wait_all(msrt);
	thread_beg_close(msrt);
	thread_end_close(msrt);
}

/*
 * abundance map
 */

static inline int abundance_map_hzmps(DMOPar *par, hzmpv *rs, int span, wtovlv *hits){
	wt_ovl_t *h;
	hzmp_t *p1, *p2;
	uint32_t x, y, beg, end, key;
	int mat, aln, ret, dir, add;
	ret = 0;
	if(rs->size == 0) return ret;
	sort_array(rs->buffer, rs->size, hzmp_t, num_cmpgtx(a.dir1 ^ a.dir2, b.dir1 ^ b.dir2, a.off2, b.off2));
	key = 0;
	while(key < rs->size){
		if(rs->buffer[key].dir1 ^ rs->buffer[key].dir2) break;
		key ++;
	}
	for(dir=0;dir<2;dir++){
		beg = dir? key : 0;
		end = dir? rs->size : key;
		if(beg >= end) continue;
		x = y = beg;
		mat = rs->buffer[beg].len2;
		h = NULL;
		while(1){
			while(y + 1 < end && rs->buffer[y + 1].off2 + rs->buffer[y + 1].len2 - rs->buffer[x].off2 <= span){
				p1 = ref_hzmpv(rs, y);
				p2 = ref_hzmpv(rs, ++y);
				if(p2->off2 < p1->off2 + p1->len2) mat += (p2->off2 + p2->len2) - (p1->off2 + p1->len2);
				else mat += p2->len2;
			}
			if(mat >= par->min_mat){
				aln = rs->buffer[y].off2 + rs->buffer[y].len2 - rs->buffer[x].off2;
				if(aln >= par->min_aln && aln * par->min_sm <= mat){
					// FOUND
					//TODO TODO TODO: careful about tandom repeats
					//if(h == NULL || num_max(x, h->pb1) + (y - x + 1) / 2 < num_min(y, h->pb2)){
					if(h == NULL || rs->buffer[x].off2 + num_min(aln, h->aln) / 4 > rs->buffer[h->pb2].off2 + rs->buffer[h->pb2].len2){
						h = next_ref_wtovlv(hits);
						ret ++;
						add = 1;
					} else if(h->mat >= mat) add = 0;
					else add = 1;
					if(add){
						h->pb1 = x; // temp use
						h->dir1 = 0;
						h->pb2 = y; // temp use
						h->dir2 = dir;
						h->qb = rs->buffer[x].off2;
						h->qe = rs->buffer[y].off2 + rs->buffer[y].len2;
						h->tb = 0;
						h->te = span;
						//h->score = mat;
						h->score = y + 1 - x;
						h->mat = mat;
						h->aln = aln;
					}
				}
			}
			if(y + 1 >= end) break;
			if(x + 1 >= end) break;
			p1 = ref_hzmpv(rs, x);
			p2 = ref_hzmpv(rs, ++x);
			if(p2->off2 < p1->off2 + p1->len2) mat -= p2->off2 - p1->off2;
			else mat -= p1->len2;
		}
	}
	return ret;
}


/*
 * syntanic map mode
 */

static inline void denoising_hzmps(hzmpv *rs, hzmpv *dst[2], wtseedv *regs[2], int xvar, int yvar, int zmin, diagv *diags, u4v *block, u4v *grps){
	diag_t *d;
	wt_seed_t *seed;
	hzmp_t *p, *p0, P;
	uint32_t i, j, k, doff, dcnt, gid;
	int dir, lst, num;
	int lst_offset, end_offset;
	sort_array(rs->buffer, rs->size, hzmp_t, ((((int64_t)a.off1 - (int64_t)a.off2) << 32) | a.off1) > ((((int64_t)b.off1 - (int64_t)b.off2) << 32) | b.off1));
	memset(&P, 0, sizeof(hzmp_t));
	P.off1 = 0x7FFFFFFFU;
	//denoising
	for(dir=0;dir<2;dir++){
		clear_diagv(diags);
		clear_hzmpv(dst[dir]);
		clear_wtseedv(regs[dir]);
		d = NULL;
		for(i=0;i<rs->size;i++){
			p = ref_hzmpv(rs, i);
			if(p->dir1 ^ p->dir2 ^ dir) continue;
			if(d && d->offset == (int)p->off1 - (int)p->off2){
				d->cnt ++;
			} else {
				d = next_ref_diagv(diags);
				d->offset = (int)p->off1 - (int)p->off2;
				d->off = i;
				d->cnt = 1;
			}
		}
		push_diagv(diags, (diag_t){0x7FFFFFFF, 0, 0});
		doff = 0;
		end_offset = - 0x7FFFFFFF;
		clear_u4v(grps);
		push_u4v(grps, 0);
		while(doff < diags->size){
			// find yvar offsets region
			lst_offset = diags->buffer[doff].offset;
			dcnt = 0;
			while(1){
				if(diags->buffer[dcnt+doff].offset > lst_offset + yvar) break;
				dcnt ++;
			}
			if(dcnt == 0) break;
			if(diags->buffer[doff + dcnt].offset == end_offset){
				doff += dcnt;
				continue;
			}
			end_offset = diags->buffer[doff + dcnt].offset;
			// sort local yvar offsets region by p->off1
			clear_u4v(block);
			for(i=0;i<dcnt;i++){
				d = ref_diagv(diags, i + doff);
				for(j=0;j<d->cnt;j++){
					p = ref_hzmpv(rs, d->off + j);
					if(p->dir1 ^ p->dir2 ^ dir) continue;
					push_u4v(block, d->off + j);
				}
			}
			sort_array(block->buffer, block->size, uint32_t, rs->buffer[a].off1 > rs->buffer[b].off1);
			// scan linear hzmps
			if(block->size){
				p0 = block->size? ref_hzmpv(rs, block->buffer[0]) : NULL;
				num = 1;
			} else {
				num = 0;
				p0 = NULL;
			}
			j = 0;
			for(i=1;i<=block->size;i++){
				p = (i == block->size)? &P : ref_hzmpv(rs, block->buffer[i]);
				if(p->off1 <= p0->off1 + p0->len1 + xvar){
					num ++;
				} else {
					if(num >= zmin){
						gid = 0;
						for(k=j;k<i;k++){
							if(rs->buffer[block->buffer[k]].gid){
								if(gid == 0){
									gid = grps->buffer[rs->buffer[block->buffer[k]].gid];
								} else if(gid > grps->buffer[rs->buffer[block->buffer[k]].gid]){
									gid = grps->buffer[rs->buffer[block->buffer[k]].gid];
								}
							}
						}
						if(gid == 0){
							gid = grps->size;
							push_u4v(grps, gid);
						} else {
							for(k=j;k<i;k++){
								if(rs->buffer[block->buffer[k]].gid){
									grps->buffer[rs->buffer[block->buffer[k]].gid] = gid;
								}
							}
						}
						for(;j<i;j++){
							//fprintf(stderr, "G[%d\t%d]\t%d\t%d\n", gid, len, ref_hzmpv(rs, block->buffer[j])->off1, ref_hzmpv(rs, block->buffer[j])->len1);
							ref_hzmpv(rs, block->buffer[j])->gid = gid;
						}
					}
					j = i;
					num = 1;
				}
				p0 = p;
			}
			// move offset to next
			for(i=doff;i<doff+dcnt;i++){
				if(diags->buffer[i].offset > lst_offset + yvar / 2) break;
			}
			doff = i;
		}
		// tidy gid map
		for(i=1;i<grps->size;i++){
			if(grps->buffer[i] < i) continue;
			for(j=i+1;j<grps->size;j++){
				if(grps->buffer[j] != i) continue;
				for(k=j+1;k<grps->size;k++){
					if(grps->buffer[k] == j){
						grps->buffer[k] = i;
					}
				}
			}
		}
		for(i=0;i<rs->size;i++){
			p = ref_hzmpv(rs, i);
			if(p->dir1 ^ p->dir2 ^ dir) continue;
			if(p->gid == 0) continue;
			// re-map gid
			p->gid = grps->buffer[p->gid];
			push_hzmpv(dst[dir], *p);
		}
		sort_array(dst[dir]->buffer, dst[dir]->size, hzmp_t, num_cmpgtx(a.gid, b.gid, a.off2, b.off2));
		// generate hzmp block: wt_seed_t
		j = 0;
		for(i=1;i<=dst[dir]->size;i++){
			if(i < dst[dir]->size && dst[dir]->buffer[i].gid == dst[dir]->buffer[j].gid) continue;
			seed = next_ref_wtseedv(regs[dir]);
			seed->pb2 = 0;
			seed->closed = 0;
			seed->dir = dir;
			seed->anchors[0] = j;
			seed->anchors[1] = i;
			seed->beg[0] = 0x7FFFFFFF;
			seed->beg[1] = 0x7FFFFFFF;
			seed->end[0] = 0;
			seed->end[1] = 0;
			seed->ovl = 0;
			lst = 0;
			for(k=j;k<i;k++){
				p = ref_hzmpv(dst[dir], k);
				if(p->off1 < seed->beg[0]) seed->beg[0] = p->off1;
				if(p->off1 + p->len1 > seed->end[0]) seed->end[0] = p->off1 + p->len1;
				if(p->off2 < seed->beg[1]) seed->beg[1] = p->off2;
				if(p->off2 + p->len2 > seed->end[1]) seed->end[1] = p->off2 + p->len2;
				//seed->ovl += (p->off1 > lst)? p->len1 : p->off1 + p->len1 - lst;
				//lst = p->off1 + p->len1;
				seed->ovl += (p->off2 > lst)? p->len2 : p->off2 + p->len2 - lst;
				lst = p->off2 + p->len2;
			}
			j = i;
		}
	}
}

static inline void fast_merge_wtseedv(hzmpv *dst, wtseedv *regs, int xvar, int yvar, diagv *diags, u4v *block, u4v *grps){
	diag_t *d;
	wt_seed_t *s, *s0, S;
	uint32_t i, j, k, l, doff, dcnt, gid;
	int lst_offset, end_offset;
	sort_array(regs->buffer, regs->size, wt_seed_t, ((((int64_t)(a.beg[0] - a.beg[1])) << 32) | a.beg[0]) > ((((int64_t)(b.beg[0] - b.beg[1])) << 32) | b.beg[0]));
	memset(&S, 0, sizeof(wt_seed_t));
	S.beg[0] = 0x7FFFFFFFU;
	clear_diagv(diags);
	for(i=0;i<regs->size;i++){
		d = NULL;
		s = ref_wtseedv(regs, i);
		if(d && d->offset == s->beg[0] - s->beg[1]){
			d->cnt ++;
		} else {
			d = next_ref_diagv(diags);
			d->offset = s->beg[0] - s->beg[1];
			d->off = i;
			d->cnt = 1;
		}
	}
	doff = 0;
	end_offset = - 0x7FFFFFFF;
	clear_u4v(grps);
	push_u4v(grps, 0);
	for(i=0;i<regs->size;i++){
		regs->buffer[i].pb2 = i + 1;
		push_u4v(grps, i + 1);
	}
	while(doff < regs->size){
		lst_offset = diags->buffer[doff].offset;
		dcnt = 0;
		while(1){
			if(diags->buffer[dcnt+doff].offset > lst_offset + yvar) break;
			if(dcnt + doff + 1 >= diags->size) break;
			dcnt ++;
		}
		if(dcnt == 0) break;
		if(diags->buffer[doff + dcnt].offset == end_offset){
			doff += dcnt;
			continue;
		}
		end_offset = diags->buffer[doff + dcnt].offset;
		clear_u4v(block);
		for(i=0;i<dcnt;i++){
			d = ref_diagv(diags, i + doff);
			for(j=0;j<d->cnt;j++){
				s = ref_wtseedv(regs, d->off + j);
				push_u4v(block, d->off + j);
			}
		}
		sort_array(block->buffer, block->size, uint32_t, regs->buffer[a].beg[0] > regs->buffer[b].beg[0]);
		if(block->size){
			s0 = block->size? ref_wtseedv(regs, block->buffer[0]) : NULL;
		} else {
			s0 = NULL;
		}
		j = 0;
		for(i=1;i<=block->size;i++){
			s = (i == block->size)? &S : ref_wtseedv(regs, block->buffer[i]);
			if(s->beg[0] <= s0->end[0] + xvar){
			} else {
				gid = 0;
				for(k=j;k<i;k++){
					if(regs->buffer[block->buffer[k]].pb2){
						if(gid == 0){
							gid = grps->buffer[regs->buffer[block->buffer[k]].pb2];
						} else {
							grps->buffer[regs->buffer[block->buffer[k]].pb2] = gid;
						}
					}
				}
				if(gid == 0){
					gid = grps->size;
					push_u4v(grps, gid);
				}
				for(;j<i;j++){ ref_wtseedv(regs, block->buffer[j])->pb2 = gid; }
				j = i;
				s0 = s;
			}
		}
		// move offset to next
		for(i=doff;i<doff+dcnt;i++){
			if(diags->buffer[i].offset > lst_offset + yvar / 2) break;
		}
		doff = i;
	}
	for(i=1;i<grps->size;i++){
		if(grps->buffer[i] < i) continue;
		for(j=i+1;j<grps->size;j++){
			if(grps->buffer[j] != i) continue;
			for(k=j+1;k<grps->size;k++){
				if(grps->buffer[k] == j){
					grps->buffer[k] = i;
				}
			}
		}
	}
	for(i=0;i<regs->size;i++){
		s = ref_wtseedv(regs, i);
		if(s->pb2 == 0) continue;
		// re-map gid
		s->pb2 = grps->buffer[s->pb2];
	}
	sort_array(regs->buffer, regs->size, wt_seed_t, num_cmpgtx(a.pb2, b.pb2, a.beg[0], b.beg[0]));
	for(j=0;j<dst->size;j++) dst->buffer[j].gid = 0;
	for(j=0;j<regs->size;j++) if(regs->buffer[j].pb2) break;
	s0 = NULL;
	for(i=j+1;i<=regs->size;i++){
		if(i < regs->size && regs->buffer[i].pb2 == regs->buffer[j].pb2) continue;
		s0 = ref_wtseedv(regs, j);
		for(l=s0->anchors[0];l<s0->anchors[1];l++) dst->buffer[l].gid = s0->pb2;
		for(k=j+1;k<i;k++){
			s = ref_wtseedv(regs, k);
			s->closed = 1;
			if(s->beg[0] < s0->beg[0]) s0->beg[0] = s->beg[0];
			if(s->end[0] > s0->end[0]) s0->end[0] = s->end[0];
			if(s->beg[1] < s0->beg[1]) s0->beg[1] = s->beg[1];
			if(s->end[1] > s0->end[1]) s0->end[1] = s->end[1];
			s0->ovl += s->ovl;
			for(l=s->anchors[0];l<s->anchors[1];l++) dst->buffer[l].gid = s0->pb2;
		}
		j = i;
	}
	sort_array(regs->buffer, regs->size, wt_seed_t, num_cmpgtx(a.closed, b.closed, a.pb2, b.pb2));
	for(i=0;i<regs->size;i++) if(regs->buffer[i].closed) break;
	regs->size = i;
	sort_array(dst->buffer, dst->size, hzmp_t, num_cmpgtx(a.gid, b.gid, a.off2, b.off2));
	s0 = NULL;
	for(i=j=0;i<regs->size;i++){
		s = ref_wtseedv(regs, i);
		while(j < dst->size && dst->buffer[j].gid < s->pb2) j ++;
		if(s0) s0->anchors[1] = j;
		s0 = s;
		s0->anchors[0] = j;
	}
	if(s0) s0->anchors[1] = dst->size;
}

#define dmo_roundup8x(n) (((n) + 0x7LLU) & 0xFFFFFFFFFFFFFFF8LLU)

static inline int chaining_overhang_wtseedv(int pblen1, int pblen2, int dir, wtseedv *rets, wtseedv *regs, uint32_t beg, uint32_t end, u8list *mem, int tail_margin, int max_overhang, float band_penalty, float gap_penalty){
	typedef struct { int weight:29; uint32_t head:1, tail:1, closed:1; int bt; } node_t;
	node_t *nodes;
	wt_seed_t *r1, *r2;
	uint32_t i, j;
	int mw, bt, band, gap, weight, W, score;
	if(0 && hzm_debug > 2){
		fprintf(hzm_debug_out, "CHAINING\tlen%u\n", pblen2);
	}
	sort_array(regs->buffer + beg, end - beg, wt_seed_t, a.beg[0] > b.beg[0]);
	clear_and_encap_u8list(mem, dmo_roundup8x(sizeof(node_t) * (end - beg)));
	nodes = (node_t*)mem->buffer;
	nodes = nodes - beg;
	for(i=beg;i<end;i++){
		nodes[i].bt = - 1;
		nodes[i].weight = 0;
		nodes[i].head = 0;
		nodes[i].tail = 0;
		nodes[i].closed = regs->buffer[i].closed;
		r1 = ref_wtseedv(regs, i);
		if(r1->beg[0] <= tail_margin || r1->beg[1] <= tail_margin){
			nodes[i].head = 1;
		}
		if(r1->end[0] + tail_margin > pblen1 || r1->end[1] + tail_margin > pblen2){
			nodes[i].tail = 1;
		}
	}
	mw = -1000000;
	bt = -1;
	for(i=beg;i<end;i++){
		r1 = ref_wtseedv(regs, i);
		if(r1->closed) continue;
		//r1->closed = 1;
		nodes[i].weight += r1->ovl + (r1->end[0] - r1->beg[0]) / 4;
		weight = nodes[i].weight * ((nodes[i].head + 3) * (nodes[i].tail + 3)) / 16;
		if(weight > mw){ mw = weight; bt = i; }
		W = nodes[i].weight / gap_penalty;
		//fprintf(hzm_debug_out, "beg[0]=%d\ti=%d\tbt=%d\tweight=%d\tW=%d\n", r1->beg[0], i, nodes[i].bt, nodes[i].weight, W);
		for(j=i+1;j<end;j++){
			r2 = ref_wtseedv(regs, j);
			if(r2->closed) continue;
			if(r2->beg[0] + max_overhang < r1->end[0]) continue;
			if(r2->beg[1] + max_overhang < r1->end[1]) continue;
			if(r2->beg[0] - r1->end[0] > W) break;
			band = num_diff(r2->beg[0] - r1->end[0], r2->beg[1] - r1->end[1]);
			gap  = num_max(r2->beg[0] - r1->end[0], r2->beg[1] - r1->end[1]);
			if(gap < 0) gap = - gap;
			score = band * band_penalty + gap * gap_penalty;
			score = nodes[i].weight - score;
			//fprintf(hzm_debug_out, "beg[0]=%d\ti=%d\tj=%d\t%d\t%d\t%d\t%d\t%d\n", r1->beg[0], i, j, band, gap, score, nodes[i].weight - score, nodes[j].weight);
			if(nodes[j].weight <= score){
				nodes[j].weight = score;
				nodes[j].bt = i;
				nodes[j].head = nodes[i].head;
			}
			//if(score < nodes[i].weight && 2 * band * band_penalty <= regs->buffer[j].ovl) break; // j can replace i in further chaining, heuristically
		}
	}
	mw = 0;
	clear_wtseedv(rets);
	while(bt >= 0){
		r1 = ref_wtseedv(regs, bt);
		//r1->closed = 0;
		//fprintf(hzm_debug_out, "CHAIN[%d->%d]\t%d\t%d\t%d\t%d\n", bt, nodes[bt].bt, r1->beg[0], r1->end[0], r1->beg[1], r1->end[1]);
		//mw += r1->end[0] - r1->beg[0];
		mw += r1->ovl;
		bt = nodes[bt].bt;
		push_wtseedv(rets, *r1);
		r1->closed = 1;
	}
	if(hzm_debug > 2){
		fprintf(hzm_debug_out, "CHAIN-HZMP\tlen%u\tlen%u\t%c\tweight=%d\n", pblen1, pblen2, "+-"[dir], mw);
		j = 0;
		for(i=beg;i<end;i++){
			r1 = ref_wtseedv(regs, i);
			if(r1->closed) continue;
			j ++;
			fprintf(hzm_debug_out, "CHAIN[%d]\t%d\t%d\t%d\t%d", j, r1->beg[0], r1->end[0], r1->beg[1], r1->end[1]);
			fprintf(hzm_debug_out, "\t%d\n", r1->ovl);
		}
	}
	return mw;
}

static inline void dot_plot_hzmps(hzmpv *rs, int dir, int pass, FILE *out){
	hzmp_t *p;
	uint32_t i;
	for(i=0;i<rs->size;i++){
		p = ref_hzmpv(rs, i);
		if(p->dir1 ^ p->dir2 ^ dir) continue;
		if((p->gid > 0) < pass) continue;
		fprintf(out, "%d\t%d\t%d\n", p->off1, p->off2, p->gid);
	}
}

static inline void debug_dot_plot_hzmps(hzmpv *rs){
	{
		FILE *out;
		out = fopen("dot_plot.fwd.src.txt", "w");
		fprintf(out, "%d\t%d\t%d\n", 0, 0, 0);
		dot_plot_hzmps(rs, 0, 0, out);
		fclose(out);
	}
	{
		FILE *out;
		out = fopen("dot_plot.fwd.dst.txt", "w");
		fprintf(out, "%d\t%d\t%d\n", 0, 0, 0);
		dot_plot_hzmps(rs, 0, 1, out);
		fclose(out);
	}
	{
		FILE *out;
		out = fopen("dot_plot.rev.src.txt", "w");
		fprintf(out, "%d\t%d\t%d\n", 0, 0, 0);
		dot_plot_hzmps(rs, 1, 0, out);
		fclose(out);
	}
	{
		FILE *out;
		out = fopen("dot_plot.rev.dst.txt", "w");
		fprintf(out, "%d\t%d\t%d\n", 0, 0, 0);
		dot_plot_hzmps(rs, 1, 1, out);
		fclose(out);
	}
}

static inline wt_ovl_t dot_matrix_align_hzmps(DMOPar *par, uint32_t pb1, int pblen1, uint32_t pb2, int pblen2, kigarv *kigars, hzmpv *rs, hzmpv *dst[2], wtseedv *regs[4], diagv *diags, u4v *block, u4v *grps, u8list *mem){
	wt_seed_t *seed;
	uint32_t i, j, d;
	int weight[2];
	wt_ovl_t ret;
	denoising_hzmps(rs, dst, regs, par->xvar, par->yvar, par->zmin, diags, block, grps);
	if(hzm_debug && pb1 == 0 && pb2 == 1) debug_dot_plot_hzmps(rs);
	fast_merge_wtseedv(dst[0], regs[0], par->xvar, 2 * par->yvar, diags, block, grps);
	fast_merge_wtseedv(dst[1], regs[1], par->xvar, 2 * par->yvar, diags, block, grps);
	weight[0] = chaining_overhang_wtseedv(pblen1, pblen2, 0, regs[2], regs[0], 0, regs[0]->size, mem, par->xvar, par->max_overhang, par->deviation_penalty, par->gap_penalty);
	weight[1] = chaining_overhang_wtseedv(pblen1, pblen2, 1, regs[3], regs[1], 0, regs[1]->size, mem, par->xvar, par->max_overhang, par->deviation_penalty, par->gap_penalty);
	d = (weight[0] < weight[1]);
	ret.pb1 = pb1;
	ret.pb2 = pb2;
	ret.dir1 = 0; ret.dir2 = d;
	ret.score = 0;
	ret.qb = ret.tb = 0x7FFFFFFF;
	ret.qe = ret.te = 0;
	if(kigars) ret.kigar_off = kigars->size;
	else ret.kigar_off = 0;
	ret.kigar_len = 0;
	for(i=0;i<regs[d + 2]->size;i++){
		seed = ref_wtseedv(regs[d + 2], i);
		if(seed->closed == 0){
			ret.score += seed->anchors[1] - seed->anchors[0];
			if(ret.qb > seed->beg[1]) ret.qb = seed->beg[1];
			if(ret.tb > seed->beg[0]) ret.tb = seed->beg[0];
			if(ret.qe < seed->end[1]) ret.qe = seed->end[1];
			if(ret.te < seed->end[0]) ret.te = seed->end[0];
			if(kigars){
				for(j=seed->anchors[0];j<seed->anchors[1];j++){
					push_kigarv(kigars, (kigar_t){dst[d]->buffer[j].off1, dst[d]->buffer[j].off2});
				}
				ret.kigar_len += (seed->anchors[1] - seed->anchors[0]);
			}
		}
	}
	if(kigars){ sort_array(kigars->buffer + ret.kigar_off, ret.kigar_len, kigar_t, num_cmpgt(a.off1, b.off1)); }
	ret.mat = weight[d];
	ret.aln = num_max(ret.te - ret.tb, ret.qe - ret.qb);
	if(hzm_debug > 1){
		//int qb, qe, tb, te, aln;
		for(d=0;d<2;d++){
			fprintf(hzm_debug_out, "\nBLOCKS\t%c\t%d\t%d\n", "+-"[d], (int)regs[d]->size, weight[d]);
			for(i=0;i<regs[d]->size;i++){
				seed = ref_wtseedv(regs[d], i);
				fprintf(hzm_debug_out, "WINDOW[%03d]\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\n", i, "+-"[d], seed->beg[0], seed->end[0], seed->beg[1], seed->end[1], seed->ovl, seed->end[0] - seed->beg[0], seed->end[1] - seed->beg[1], "*x--"[seed->closed], seed->anchors[0], seed->anchors[1]);
			}
			//if(weight[d] == 0) continue;
			//aln = num_max(qe - qb, te - tb);
			//fprintf(hzm_debug_out, "OVL\t%c\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%0.3f\n", '+', qb, qe, "+-"[d], tb, te, weight[d], aln, 1.0 * weight[d] / aln);
		}
	}
	return ret;
}

static inline uint32_t dot_matrix_full_align_hzmps(DMOPar *par, uint32_t pb1, int pblen1, uint32_t pb2, int pblen2, wtovlv *hits, kigarv *kigars, hzmpv *rs, hzmpv *dst[2], wtseedv *regs[4], diagv *diags, u4v *block, u4v *grps, u8list *mem){
	wt_seed_t *seed;
	uint32_t i, j, d, n;
	int weight[2];
	wt_ovl_t ret;
	denoising_hzmps(rs, dst, regs, par->xvar, par->yvar, par->zmin, diags, block, grps);
	if(hzm_debug && pb1 == 0 && pb2 == 1) debug_dot_plot_hzmps(rs);
	fast_merge_wtseedv(dst[0], regs[0], par->xvar, 2 * par->yvar, diags, block, grps);
	fast_merge_wtseedv(dst[1], regs[1], par->xvar, 2 * par->yvar, diags, block, grps);
	n = 0;
	for(d=0;d<2;d++){
		while(1){
			weight[d] = chaining_overhang_wtseedv(pblen1, pblen2, d, regs[d + 2], regs[d], 0, regs[d]->size, mem, par->xvar, par->max_overhang, par->deviation_penalty, par->gap_penalty);
			if(hzm_debug > 1){
				int qb, qe, tb, te, aln;
				fprintf(hzm_debug_out, "\nBLOCKS\t%c\t%d\t%d\n", "+-"[d], (int)regs[d]->size, weight[d]);
				qb = tb = 0x7FFFFFFF;
				qe = te = 0;
				for(i=0;i<regs[d + 2]->size;i++){
					seed = ref_wtseedv(regs[d + 2], i);
					if(seed->closed == 0){
						if(qb > seed->beg[0]) qb = seed->beg[0];
						if(tb > seed->beg[1]) tb = seed->beg[1];
						if(qe < seed->end[0]) qe = seed->end[0];
						if(te < seed->end[1]) te = seed->end[1];
					}
				fprintf(hzm_debug_out, "WINDOW[%03d]\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\n", i, "+-"[d], seed->beg[0], seed->end[0], seed->beg[1], seed->end[1], seed->ovl, seed->end[0] - seed->beg[0], seed->end[1] - seed->beg[1], "*x--"[seed->closed], seed->anchors[0], seed->anchors[1]);
				}
				if(weight[d]){
					aln = num_max(qe - qb, te - tb);
					fprintf(hzm_debug_out, "OVL\t%c\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%0.3f\n", '+', qb, qe, "+-"[d], tb, te, weight[d], aln, 1.0 * weight[d] / aln);
				}
			}
			if(weight[d] == 0 || weight[d] < par->min_mat) break;
			ret.pb1 = pb1;
			ret.pb2 = pb2;
			ret.dir1 = 0; ret.dir2 = d;
			ret.score = 0;
			ret.qb = ret.tb = 0x7FFFFFFF;
			ret.qe = ret.te = 0;
			if(kigars) ret.kigar_off = kigars->size;
			else ret.kigar_off = 0;
			ret.kigar_len = 0;
			for(i=0;i<regs[d+2]->size;i++){
				seed = ref_wtseedv(regs[d+2], i);
				if(seed->closed == 0){
					ret.score += seed->anchors[1] - seed->anchors[0];
					if(ret.qb > seed->beg[1]) ret.qb = seed->beg[1];
					if(ret.tb > seed->beg[0]) ret.tb = seed->beg[0];
					if(ret.qe < seed->end[1]) ret.qe = seed->end[1];
					if(ret.te < seed->end[0]) ret.te = seed->end[0];
					if(kigars){
						for(j=seed->anchors[0];j<seed->anchors[1];j++){
							push_kigarv(kigars, (kigar_t){dst[d]->buffer[j].off1, dst[d]->buffer[j].off2});
						}
						ret.kigar_len += (seed->anchors[1] - seed->anchors[0]);
					}
				}
			}
			if(kigars){ sort_array(kigars->buffer + ret.kigar_off, ret.kigar_len, kigar_t, num_cmpgt(a.off1, b.off1)); }
			ret.mat = weight[d];
			ret.aln = num_max(ret.te - ret.tb, ret.qe - ret.qb);
			if(ret.mat && hzm_debug > 2){
				fprintf(hzm_debug_out, "HIT:rd[%u]\t%c\t%d\t%d\t%d\trd[%u]\t%c\t%d\t%d\t%d\t%d\t%d\t%0.4f\n", pb1, "+-"[ret.dir1], pblen1, ret.tb, ret.te,
					pb2, "+-"[ret.dir2], pblen2, ret.qb, ret.qe, ret.aln, ret.mat, 1.0 * ret.mat / (ret.aln + 1));
			}
			if(num_min(ret.qe - ret.qb, ret.te - ret.tb) < par->min_aln) continue;
			if(ret.aln * par->min_sm > ret.mat) continue;
			if(num_diff(ret.qe - ret.qb, ret.te - ret.tb) > par->aln_var * num_min(ret.qe - ret.qb, ret.te - ret.tb)) continue;
			push_wtovlv(hits, ret);
			n ++;
		}
	}
	sort_array(hits->buffer + hits->size - n, n, wt_ovl_t, b.mat > a.mat);
	return n;
}

static inline DMOAux* init_dmoaux(){
	DMOAux *aux;
	aux = malloc(sizeof(DMOAux));
	aux->refs = init_hzrefv(32);
	aux->heap = init_u32list(32);
	aux->hzoff = init_u32list(32);
	aux->rs = init_hzmpv(32);
	aux->dst[0] = init_hzmpv(32);
	aux->dst[1] = init_hzmpv(32);
	aux->regs[0] = init_wtseedv(32);
	aux->regs[1] = init_wtseedv(32);
	aux->regs[2] = init_wtseedv(32);
	aux->regs[3] = init_wtseedv(32);
	aux->diags = init_diagv(32);
	aux->block = init_u4v(32);
	aux->grps = init_u4v(32);
	aux->mem = init_u8list(32);
	aux->trials = 0;
	return aux;
}


static inline void free_dmoaux(DMOAux *aux){
	free_hzrefv(aux->refs);
	free_u32list(aux->heap);
	free_u32list(aux->hzoff);
	free_hzmpv(aux->rs);
	free_hzmpv(aux->dst[0]);
	free_hzmpv(aux->dst[1]);
	free_wtseedv(aux->regs[0]);
	free_wtseedv(aux->regs[1]);
	free_wtseedv(aux->regs[2]);
	free_wtseedv(aux->regs[3]);
	free_diagv(aux->diags);
	free_u4v(aux->block);
	free_u4v(aux->grps);
	free_u8list(aux->mem);
	free(aux);
}

// store results in aux->refs
static inline void query_index_dmo(DMO *wt, uint32_t pbid, BaseBank *rdseqs, uint64_t pboff, uint32_t pblen, DMOAux *aux){
	hz_ref_t *r;
	hzmh_t *h;
	uint64_t kmer, KMER, krev, kmask, KMASK1, KMASK2, off;
	uint32_t i, j, ksize, kidx, dir, kg;
	uint8_t b, c;
	clear_hzrefv(aux->refs);
	clear_u32list(aux->heap);
	{
		ksize = wt->par->ksize;
		kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - ksize) << 1);
		b = 4;
		kmer = 0;
		off = pboff;
		for(i=j=0;j<pblen;j++){
			c = bits2bit(rdseqs->bits, off); off ++;
			if(wt->par->hk && c == b) continue;
			b = c;
			i ++;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < ksize) continue;
			krev = dna_rev_seq(kmer, ksize);
			if(krev == kmer) continue;
			dir  = krev > kmer? 0 : 1;
			krev = krev > kmer? kmer : krev;
			kidx = hzmh_kmer_smear(krev) % wt->par->hzmh_kmer_mod;
			if(kidx >= HZMH_KMER_MOD) continue;
			h = get_hzmhash(wt->hashs[kidx], (hzmh_t){krev, 0, 0, 0});
			if(h == NULL) continue;
			if(h->flt) continue;
			r = next_ref_hzrefv(aux->refs);
			r->p1.off = j + 1 - ksize;
			r->p1.len = ksize;
			r->p1.dir = dir;
			r->beg = h->off;
			r->end = h->off + h->cnt;
			while(r->beg < r->end){
				r->p2 = ref_hzmv(wt->seeds, r->beg ++);
				if(pbid < wt->n_rd && r->p2->rd_id <= pbid) continue;
				array_heap_push(aux->heap->buffer, aux->heap->size, aux->heap->cap, uint32_t, aux->refs->size - 1, heap_cmp_hz_ref_macro(aux->refs, a, b));
				break;
			}
		}
	}
	for(kg=1;kg<=wt->par->kgap;kg++){
		ksize = wt->par->ksize + kg;
		if(ksize > 32) break;
		kmask  = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - ksize) << 1);
		KMASK1 = 0xFFFFFFFFFFFFFFFFLLU << (((wt->par->ksize >> 1) + kg) << 1);
		KMASK2 = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - (wt->par->ksize >> 1)) << 1);
		b = 4;
		kmer = 0;
		off = pboff;
		for(i=j=0;j<pblen;j++){
			c = bits2bit(rdseqs->bits, off); off ++;
			if(wt->par->hk && c == b) continue;
			b = c;
			i ++;
			kmer = ((kmer << 2) | b) & kmask;
			if(i < ksize) continue;
			KMER = ((kmer & KMASK1) >> (kg << 1)) | (kmer & KMASK2);
			krev = dna_rev_seq(KMER, ksize);
			if(krev == KMER) continue;
			dir  = krev > KMER? 0 : 1;
			krev = krev > KMER? KMER : krev;
			kidx = hzmh_kmer_smear(krev) % wt->par->hzmh_kmer_mod;
			if(kidx >= HZMH_KMER_MOD) continue;
			h = get_hzmhash(wt->hashs[kidx], (hzmh_t){krev, 0, 0, 0});
			if(h == NULL) continue;
			if(h->flt) continue;
			r = next_ref_hzrefv(aux->refs);
			r->p1.off = j - ksize;
			r->p1.len = ksize;
			r->p1.dir = dir;
			r->beg = h->off;
			r->end = h->off + h->cnt;
			while(r->beg < r->end){
				r->p2 = ref_hzmv(wt->seeds, r->beg ++);
				if(pbid < wt->n_rd && r->p2->rd_id <= pbid) continue;
				array_heap_push(aux->heap->buffer, aux->heap->size, aux->heap->cap, uint32_t, aux->refs->size - 1, heap_cmp_hz_ref_macro(aux->refs, a, b));
				break;
			}
		}
	}
}

static inline void query_index_by_kmers_dmo(DMO *wt, uint32_t pbid, u8v *kmers, u1v *kdirs, u4v *koffs, DMOAux *aux){
	hz_ref_t *r;
	hzmh_t *h;
	uint64_t krev;
	uint32_t i, kidx, off;
	int dir;
	clear_hzrefv(aux->refs);
	for(i=0;i<kmers->size;i++){
		krev = kmers->buffer[i];
		dir  = kdirs->buffer[i];
		off  = koffs->buffer[i];
		kidx = hzmh_kmer_smear(krev) % wt->par->hzmh_kmer_mod;
		if(kidx >= HZMH_KMER_MOD) continue;
		h = get_hzmhash(wt->hashs[kidx], (hzmh_t){krev, 0, 0, 0});
		if(h == NULL) continue;
		if(h->flt) continue;
		r = next_ref_hzrefv(aux->refs);
		r->p1.off = off;
		r->p1.len = wt->par->ksize;
		r->p1.dir = dir;
		r->beg = h->off;
		r->end = h->off + h->cnt;
		while(r->beg < r->end){
			r->p2 = ref_hzmv(wt->seeds, r->beg ++);
			if(pbid < wt->n_rd && r->p2->rd_id <= pbid) continue;
			array_heap_push(aux->heap->buffer, aux->heap->size, aux->heap->cap, uint32_t, aux->refs->size - 1, heap_cmp_hz_ref_macro(aux->refs, a, b));
			break;
		}
	}
}

static inline void query_dmo(DMO *wt, uint32_t pbid, BaseBank *rdseqs, uint64_t pboff, uint32_t pblen, int align_mode, wtovlv *hits, kigarv *kigars, DMOAux *aux){
	hz_ref_t *r, R;
	wt_ovl_t hit, *hh;
	hzm_t *p1, P1;
	hzm_t *p2, P2;
	uint32_t i, n, id1, id2, cnt, dir;
	int mat, aln, lst;
	uint8_t has_next;
	if(rdseqs) query_index_dmo(wt, pbid, rdseqs, pboff, pblen, aux);
	// Else you have already prepare aux->refs
	R.p1 = (hzm_t){0xFFFFFFFFU, 0, 0, 0, 0};
	P2   = (hzm_t){0x7FFFFFFFU, 0, 0, 0, 0};
	R.p2 = &P2;
	R.beg = R.end = 0;
	p1 = &P1;
	id2 = 0xFFFFFFFFU;
	mat = aln = cnt = 0; lst = -1;
	dir = 0;
	aux->trials = 0;
	clear_hzmpv(aux->rs);
	while(1){
		if(aux->heap->size){
			id1 = aux->heap->buffer[0];
			r = ref_hzrefv(aux->refs, id1);
			p1 = &r->p1;
			p2 = r->p2;
			has_next = 0;
			while(r->beg < r->end){
				r->p2 = ref_hzmv(wt->seeds, r->beg ++);
				if(pbid < wt->n_rd && r->p2->rd_id <= pbid) continue;
				array_heap_replace(aux->heap->buffer, aux->heap->size, aux->heap->cap, uint32_t, 0, id1, heap_cmp_hz_ref_macro(aux->refs, a, b));
				has_next = 1;
				break;
			}
			if(has_next == 0) array_heap_remove(aux->heap->buffer, aux->heap->size, aux->heap->cap, uint32_t, 0, heap_cmp_hz_ref_macro(aux->refs, a, b));
		} else {
			id1 = 0xFFFFFFFFU;
			p1 = &R.p1;
			p2 = R.p2;
		}
		if(id2 != p2->rd_id){
			if(id2 != 0xFFFFFFFFU){
				if(hzm_debug > 2 && cnt > 2){
					fprintf(hzm_debug_out, "FAXIAN\t%s\t%d\t%d\t%d\n", wt->reads->buffer[id2].tag, aln, mat, cnt);
				}
				if((mat >= wt->par->min_mat && aln >= wt->par->min_aln) || (align_mode == 2 && (int)(cnt * wt->par->ksize) >= wt->par->min_mat)){
					aux->trials ++;
					if(align_mode == 1){
						n = dot_matrix_full_align_hzmps(wt->par, pbid, pblen, id2, wt->reads->buffer[id2].rdlen, hits, kigars, aux->rs, aux->dst, aux->regs, aux->diags, aux->block, aux->grps, aux->mem);
					} else if(align_mode == 0){
						hit = dot_matrix_align_hzmps(wt->par, pbid, pblen, id2, wt->reads->buffer[id2].rdlen, kigars, aux->rs, aux->dst, aux->regs, aux->diags, aux->block, aux->grps, aux->mem);
						hit.pb1 = pbid;
						hit.pb2 = id2;
						if(hit.mat && hzm_debug > 2){
							fprintf(hzm_debug_out, "HIT:%s\t%c\t%d\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%0.4f\n", wt->reads->buffer[hit.pb1].tag, "+-"[hit.dir1], wt->reads->buffer[hit.pb1].rdlen, hit.tb, hit.te,
								wt->reads->buffer[hit.pb2].tag, "+-"[hit.dir2], wt->reads->buffer[hit.pb2].rdlen, hit.qb, hit.qe, hit.aln, hit.mat, 1.0 * hit.mat / (hit.aln + 1));
						}
						if(hit.mat >= wt->par->min_mat
							&& num_min(hit.qe - hit.qb, hit.te - hit.tb) >= wt->par->min_aln
							&& hit.aln * wt->par->min_sm <= hit.mat
							&& num_diff(hit.qe - hit.qb, hit.te - hit.tb) <= wt->par->aln_var * num_min(hit.qe - hit.qb, hit.te - hit.tb))
							push_wtovlv(hits, hit);
					} else {
						n  = abundance_map_hzmps(wt->par, aux->rs, pblen, hits);
						for(i=0;i<n;i++){
							hh = ref_wtovlv(hits, hits->size - 1 - i);
							hh->pb1 = pbid;
							hh->pb2 = id2;
						}
					}
				}
			}
			id2 = p2->rd_id;
			mat = aln = cnt = 0; lst = -1;
			clear_hzmpv(aux->rs);
		}
		if(id1 == 0xFFFFFFFFU) break;
		if(hzm_debug > 3){
			fprintf(hzm_debug_out, "HEAP\t%s\t%c\t%d\t%d\t%c\t%d\t%d\n", wt->reads->buffer[id2].tag, "+-"[p1->dir], p1->off, p1->len, "+-"[p2->dir], p2->off, p2->len);
		}
		if(lst == -1) lst = p1->off;
		if(p1->off >= lst) mat += p1->len;
		else mat += p1->off + p1->len - lst;
		aln += p1->off + p1->len - lst;
		lst = p1->off + p1->len;
		cnt ++;
		if(p1->dir ^ p2->dir){
			push_hzmpv(aux->rs, (hzmp_t){p1->dir, p1->off, p2->dir, wt->reads->buffer[p2->rd_id].rdlen - (p2->off + p2->len), p1->len, p2->len, 0});
		} else {
			push_hzmpv(aux->rs, (hzmp_t){p1->dir, p1->off, p2->dir, p2->off, p1->len, p2->len, 0});
		}
	}
}

static inline void index_reg_dmo(BaseBank *rdseqs, uint64_t pboff, uint32_t pblen, DMOPar *par, hzmhash *hash, hzmrv *seeds){
	hzmh_t U, *u;
	hzmr_t *s;
	uint64_t kmer, krev, kmask, off;
	uint32_t i;
	uint8_t b, c, dir;
	int exists;
	if(par->ksize > 16) par->ksize = 16;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - par->ksize) << 1);
	clear_hzmhash(hash);
	clear_hzmrv(seeds);
	off = pboff;
	b = 4;
	kmer = 0;
	for(i=0;i<pblen;i++){
		c = bits2bit(rdseqs->bits, off); off ++;
		if(par->hk && c == b) continue;
		b = c;
		kmer = ((kmer << 2) | b) & kmask;
		if(i < par->ksize) continue;
		krev = dna_rev_seq(kmer, par->ksize);
		if(krev == kmer) continue;
		dir  = krev > kmer? 0 : 1;
		krev = krev > kmer? kmer : krev;
		push_hzmrv(seeds, (hzmr_t){krev, dir, i - par->ksize, par->ksize, 0});
	}
	sort_array(seeds->buffer, seeds->size, hzmr_t, num_cmpgtx(a.mer, b.mer, a.off, b.off));
	U = (hzmh_t){0xFFFFFFFFFFFFFFFFLLU, 0, 0, 0};
	u = NULL;
	for(i=0;i<seeds->size;i++){
		s = ref_hzmrv(seeds, i);
		if(u == NULL || u->mer != s->mer){
			U.mer = s->mer;
			u = prepare_hzmhash(hash, U, &exists);
			u->mer = s->mer;
			u->off = i;
			u->cnt = 1;
			u->flt = 0;
		} else u->cnt ++;
	}
}

static inline void query_reg_dmo(BaseBank *rdseqs, uint32_t qid, uint64_t qoff, uint32_t qlen, uint32_t tid, uint32_t tlen, hzmhash *hash, hzmrv *seeds, wtovlv *hits, kigarv *kigars, DMOPar *par, DMOAux *aux){
	hzmh_t U, *u;
	hzmr_t *s;
	wt_ovl_t hit, *hh;
	uint64_t kmer, krev, kmask, off;
	uint32_t i, j, n;
	int align_mode = 0;
	uint8_t b, c, dir;
	if(par->ksize > 16) par->ksize = 16;
	U = (hzmh_t){0xFFFFFFFFFFFFFFFFLLU, 0, 0, 0};
	clear_hzmpv(aux->rs);
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - par->ksize) << 1);
	off = qoff;
	b = 4;
	kmer = 0;
	for(i=0;i<qlen;i++){
		c = bits2bit(rdseqs->bits, off); off ++;
		if(par->hk && c == b) continue;
		b = c;
		kmer = ((kmer << 2) | b) & kmask;
		if(i < par->ksize) continue;
		krev = dna_rev_seq(kmer, par->ksize);
		if(krev == kmer) continue;
		dir  = krev > kmer? 0 : 1;
		krev = krev > kmer? kmer : krev;
		U.mer = krev;
		if((u = get_hzmhash(hash, U)) == NULL) continue;
		if(u->flt) continue;
		for(j=0;j<u->cnt;j++){
			s = ref_hzmrv(seeds, u->off + j);
			if(s->dir ^ dir){
				if(par->strand_mask & 2){
					push_hzmpv(aux->rs, (hzmp_t){s->dir, s->off, dir, qlen - i, s->len, par->ksize, 0});
				}
			} else {
				if(par->strand_mask & 1){
					push_hzmpv(aux->rs, (hzmp_t){s->dir, s->off, dir, i - par->ksize, s->len, par->ksize, 0});
				}
			}
		}
	}
	if(align_mode == 1){
		n = dot_matrix_full_align_hzmps(par, tid, tlen, qid, qlen, hits, kigars, aux->rs, aux->dst, aux->regs, aux->diags, aux->block, aux->grps, aux->mem);
	} else if(align_mode == 0){
		hit = dot_matrix_align_hzmps(par, tid, tlen, qid, qlen, kigars, aux->rs, aux->dst, aux->regs, aux->diags, aux->block, aux->grps, aux->mem);
		if(hit.mat && hzm_debug > 2){
			fprintf(hzm_debug_out, "HIT:rd[%u]\t%c\t%d\t%d\t%d\trd[%u]\t%c\t%d\t%d\t%d\t%d\t%d\t%0.4f\n", hit.pb1, "+-"[hit.dir1], tid, hit.tb, hit.te,
				hit.pb2, "+-"[hit.dir2], qlen, hit.qb, hit.qe, hit.aln, hit.mat, 1.0 * hit.mat / (hit.aln + 1));
		}
		if(hit.mat >= par->min_mat
			&& num_min(hit.qe - hit.qb, hit.te - hit.tb) >= par->min_aln
			&& hit.aln * par->min_sm <= hit.mat
			&& num_diff(hit.qe - hit.qb, hit.te - hit.tb) <= par->aln_var * num_min(hit.qe - hit.qb, hit.te - hit.tb))
			push_wtovlv(hits, hit);
	} else {
		n  = abundance_map_hzmps(par, aux->rs, tlen, hits);
		for(i=0;i<n;i++){
			hh = ref_wtovlv(hits, hits->size - 1 - i);
			hh->pb1 = tid;
			hh->pb2 = qid;
		}
	}
}

static inline uint32_t print_hits_dmo(DMO *wt, wtovlv *hits, kigarv *kigars, FILE *out){
	wt_ovl_t *hit;
	uint32_t i, j;
	for(i=0;i<hits->size;i++){
		hit = ref_wtovlv(hits, i);
		if(hit->pb1 >= wt->reads->size)
			fprintf(out, "%s\t%c\t%d\t%d\t%d", "nil", "+-"[hit->dir1], -1, hit->tb, hit->te);
		else
			fprintf(out, "%s\t%c\t%d\t%d\t%d", wt->reads->buffer[hit->pb1].tag, "+-"[hit->dir1], wt->reads->buffer[hit->pb1].rdlen, hit->tb, hit->te);
		fprintf(out, "\t%s\t%c\t%d\t%d\t%d", wt->reads->buffer[hit->pb2].tag, "+-"[hit->dir2], wt->reads->buffer[hit->pb2].rdlen, hit->qb, hit->qe);
		fprintf(out, "\t%d\t%0.3f\t%d\t%d\t%d\t%d\t", hit->score, 1.0 * hit->mat / hit->aln, hit->mat, 0, 0, 0);
		if(kigars && hit->kigar_len){
			for(j=hit->kigar_off;j<hit->kigar_off+hit->kigar_len;j++) fprintf(out, "%d-%d,", kigars->buffer[j].off1, kigars->buffer[j].off2);
		} else fprintf(out, "*");
		fprintf(out, "\n");
	}
	clear_wtovlv(hits);
	return i;
}

#endif
