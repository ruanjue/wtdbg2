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

#ifndef __KMER_BINMAP_RJ_H
#define __KMER_BINMAP_RJ_H

#include "list.h"
#include "hashset.h"
#include "dna.h"
#include "filereader.h"
#include "bitvec.h"
#include "bitsvec.h"
#include "bit2vec.h"
#include "thread.h"
#include "txtplot.h"

//#define __DEBUG__	1
#define TEST_MODE

#define KBM_BIN_BITS	8
#define KBM_BSIZE	256
#define KBM_BIN_SIZE	KBM_BSIZE
#define KBM_MAX_BINCNT	0xFFFFFFFFFFLLU // 40 bits, 1024 G
#define KBM_MAX_RDCNT	0x3FFFFFFF // 30 bits, 1 G
#define KBM_MAX_RDBINCNT	0x7FFFFF // 23 bits
// able to index reference sequences
#define KBM_MAX_RDLEN	0x7FFFFFFFU // 31 bits, 2 G bp
#define KBM_MAX_KSIZE	23
#define KBM_MAX_KCNT	0xFFFF // 16 bits, 65535
#define KBM_N_HASH	4096
#define KBM_LOGF	stderr
#define KBM_LOGFNO	STDERR_FILENO
static int KBM_LOG = 0;
#define KBM_LOG_LOW	1
#define KBM_LOG_MID	2
#define KBM_LOG_HIG	3
#define KBM_LOG_ALL	4
#define KBM_MAX_RDGRP1	0x7FFFFF
#define KBM_MAX_RDGRP2	0xFF

typedef struct {
	u8i seqoff:40, flag:24; // seqoff is counted in BIN(256bp)
	u8i binoff:40, bincnt:23, closed:1;
	char *tag;
} kbm_read_t;
define_list(kbmreadv, kbm_read_t);
const obj_desc_t kbm_read_t_obj_desc = {"kbm_read_t_obj_desc", sizeof(kbm_read_t), 1, {1}, {offsetof(kbm_read_t, tag)}, {&OBJ_DESC_CHAR_ARRAY}, NULL, NULL};
static inline size_t kbmreadv_deep_obj_desc_cnt(void *list, int idx){ if(idx == 0) return ((kbmreadv*)list)->size; else return 0; }
static const obj_desc_t kbmreadv_deep_obj_desc = {.tag = "kbmreadv_deep_obj_desc", .size = sizeof(kbmreadv), .n_child = 1, .mem_type = {1}, .addr = {offsetof(kbmreadv, buffer)}, .desc = {&kbm_read_t_obj_desc}, .cnt = kbmreadv_deep_obj_desc_cnt, .post = NULL};
#define kbm_read_seqoff(rd) ((rd)->seqoff << KBM_BIN_BITS)
#define kbm_read_seqlen(rd) ((rd)->bincnt << KBM_BIN_BITS)

#define KBM_MAX_BIN_DEGREE	0x1FFU
// each BIN takes KBM_BIN_SIZE bp in uncompressed reads
typedef struct {
	u4i ridx:30, off:24, closed:1, degree:9; // off * KBM_BIN_SIZE is the real position
} kbm_bin_t;
define_list(kbmbinv, kbm_bin_t);

typedef struct {
	u4i bidx;
} kbm_bmer_t;
define_list(kbmbmerv, kbm_bmer_t);

typedef struct {
	u1i bidx; // bidx = (kbm_baux_t->bidx << 32 | kbm_bmer_t->bidx)
	u1i dir:1, koff:7; // koff is the real (offset? >> 1), here offset is +0 or +1
} kbm_baux_t;
define_list(kbmbauxv, kbm_baux_t);

#define getval_bidx(kbm, offset) ((((u8i)((kbm)->sauxs->buffer[offset].bidx)) << 32) | (kbm)->seeds->buffer[offset].bidx)

//#define kbm_kmer_smear(K) ((K) ^ ((K) >> 4) ^ ((K) >> 7) ^ ((K) >> 12))
#define kbm_kmer_smear(K) ((K) ^ ((K) >> 4) ^ ((K) >> 7))

typedef struct {
	u8i mer:46, tot:17, flt:1;
} kbm_kmer_t;
define_list(kbmkmerv, kbm_kmer_t);
#define KBM_KMERCODE(E) ((E).mer)
#define KBM_KMEREQUALS(E1, E2) ((E1).mer == (E2).mer)
#define KBM_KEYEQUALS(K, E) ((K) == (E).mer)
define_hashtable(kbmhash, kbm_kmer_t, KBM_KMERCODE, KBM_KMEREQUALS, u8i, ITSELF, KBM_KEYEQUALS, kbm_kmer_t*, ITSELF);

typedef struct {
	u8i off:40, cnt:24;
} kbm_kaux_t;
define_list(kbmkauxv, kbm_kaux_t);

typedef struct {
	kbm_kmer_t *mer;
	kbm_kaux_t *aux;
	u4i kidx;
	u4i off:24, dir:1, pdir:1, fine:1, closed:1, extra_bits1:4;
	u4i qbidx;
	u4i poffs[2];
	u8i bidx, boff, bend;
	//kbm_bmer_t *b, *end;
} kbm_ref_t;
define_list(kbmrefv, kbm_ref_t);
#if 0
#define heap_cmp_kbm_bmer(refs, a, b)	num_cmpx(refs[a].b->bidx, refs[b].b->bidx, refs[a].poffs[refs[a].pdir], refs[b].poffs[refs[b].pdir])
#endif
#define heap_cmp_kbm_bmer(refs, a, b)	num_cmpx(refs[a]->bidx, refs[b]->bidx, refs[a].poffs[refs[a].pdir], refs[b].poffs[refs[b].pdir])

typedef struct {
	u4i koff;
	u4i kcnt:8, kmat:9, boff:15; // offset from the start bin_idx
} kbm_cmer_t;
define_list(kbmcmerv, kbm_cmer_t);
static const kbm_cmer_t KBM_CMER_NULL = {0, 0, 0, 0};

typedef struct {
	u8i beg:46, mat:16, bt:2;
	b2i var;
	u2i gap;
	b4i score;
} kbm_cell_t;
static const kbm_cell_t KBM_CELL_NULL = {0, 0, 0, 0, 0, 0};
define_list(kbmcellv, kbm_cell_t);

typedef struct {
	u8i beg, end;
	u4i mat;
	int score;
} kbm_path_t;
define_list(kbmpathv, kbm_path_t);
#define kbmpath_hashcode(E) E.beg
#define kbmpath_hashequals(E1, E2) (E1).beg == (E2).beg
define_hashset(kbmphash, kbm_path_t, kbmpath_hashcode, kbmpath_hashequals);

typedef struct {
	u4i qidx:31, qdir:1;
	u4i tidx:31, tdir:1;
	u8i cgoff:40, cglen:24;
	int qb, qe, tb, te;
	int mat, cnt, aln, gap; // gap is counted in BINs
} kbm_map_t;
define_list(kbmmapv, kbm_map_t);

typedef struct {
	int rd_len_order; // 0
	int min_bin_degree; // 2
	u4i ksize, psize; // 0, 21
	u4i kmax, kmin, kmer_mod, ksampling; // 1000, 1, 4.0 * KBM_N_HASH, KBM_BSIZE
	float ktop; // 0.05
	// runtime
	u4i strand_mask; // 3. 1: forward; 2: reverse; 3: both
	int self_aln; // 0. 0: map to all; 1: only map to greater read_id; 2: itself but reverse complementary
	int skip_contained; // 1
	u2i max_bgap; // 4
	u2i max_bvar; // 4
	float max_gap; // 0.6
	u8i max_bcnt; // 0xFFFF
	int pgap, pvar; // -7, -21
	u4i max_hit; // 1000
	int min_aln, min_mat; // 2048/256, 200
	float aln_var; // 0.25
	float min_sim; // kmer similarity: 0.05
	int test_mode; // see codes
} KBMPar;

static const obj_desc_t kbmpar_obj_desc = {"kbmpar_obj_desc", sizeof(KBMPar), 0, {}, {}, {}, NULL, NULL};

typedef struct {
	u8i      flags; // 64 bits, 0: mem_load all, 1: mem_load rdseqs+reads; 2: Single Hash Mode;3-63: unused
	KBMPar   *par;
	BaseBank *rdseqs;
	u8i      avail_bases;
	u4i      avail_reads;
	kbmreadv *reads;
	cuhash   *tag2idx;
	kbmbinv  *bins;
	BitVec   *binmarks;
	kbmbmerv *seeds;
	kbmbauxv *sauxs;
	kbmhash  *hashs[KBM_N_HASH];
	kbmkauxv *kauxs[KBM_N_HASH];
} KBM;

static inline size_t kbm_obj_desc_cnt(void *kbm, int idx){
	UNUSED(kbm);
	if(idx == 8 || idx == 9) return KBM_N_HASH;
	else return 1;
}

static inline void rebuild_tag2idx_kbm(void *kbm, size_t aux);

static const obj_desc_t kbm_obj_desc = {.tag = "kbm_obj_desc", .size = sizeof(KBM), .n_child = 10,
		.mem_type = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2},
		.addr = {offsetof(KBM, par), offsetof(KBM, rdseqs), offsetof(KBM, reads), offsetof(KBM, tag2idx), offsetof(KBM, bins), offsetof(KBM, binmarks), offsetof(KBM, seeds), offsetof(KBM, sauxs), offsetof(KBM, hashs), offsetof(KBM, kauxs)},
		.desc = {&kbmpar_obj_desc, &basebank_obj_desc, &kbmreadv_deep_obj_desc, &cuhash_obj_desc, &kbmbinv_obj_desc, &bitvec_obj_desc, &kbmbmerv_obj_desc, &kbmbauxv_obj_desc, &kbmhash_obj_desc, &kbmkauxv_obj_desc},
		kbm_obj_desc_cnt, rebuild_tag2idx_kbm
	};
// Please note that, kbm->tag2idx is not functional after mem_load, because we use cuhash_obj_desc instread of cuhash_deep_obj_desc

typedef struct {
#if 0
	u4i poff, bidx;
	u4i refidx:26, koff:6;
#endif
	u4i poff;
	u4i refidx;
	u8i bidx:40, koff:24;
} kbm_dpe_t;
define_list(kbmdpev, kbm_dpe_t);

typedef struct {
	kbmdpev  *kms; // kmer offset in query and bidx
	u4i      km_len;
	BitVec   *cmask; // bit for kbm_cmer_t
	kbmcmerv *cms;
	u4v      *coffs; // kbm_cmer_t offset for each bin
	BitVec   *rmask[2];
	kbmcellv *cells[2];
	Bit2Vec  *bts; // back trace flag: 0: diagonal, 1: horizontal, 2: vertical
	kbmphash *paths; // storing best unique paths by now
	u8i      boff;
	u8i      last_bidx;
} KBMDP;

typedef struct {
	u8i kmer;
	u4i off;
	u4i kidx:30, dir:1, closed:1;
} kmer_off_t;
define_list(kmeroffv, kmer_off_t);

typedef struct {
	KBM    *kbm;
	KBMPar *par; // can diff from kbm->par
	char     *qtag;
	BaseBank *qseqs;
	u8i      qsoff;
	u4i      qlen, slen, qnbin, qnbit;
	u4i      qidx;
	u8i      bmin, bmax;
	kmeroffv *koffs[2];
	kbmrefv *refs;
	u4v      *rank, **heaps;
	u4i      nheap;
	u4i      hptr;
	u2i      *binmap;
	u4i      bmlen, bmoff, bmcnt;
	kbmdpev  *caches[2];
	KBMDP    *dps[2];
	kbmmapv  *hits;
	BitsVec  *cigars;
	BitVec   *solids;
	String   *str;
} KBMAux;

static inline KBMPar* init_kbmpar(){
	KBMPar *par;
	par = malloc(sizeof(KBMPar));
	par->rd_len_order = 0;
	par->min_bin_degree = 2;
	par->ksize = 0;
	par->psize = 21;
	par->kmax  = 1000;
	par->kmin  = 1;
	par->kmer_mod = 4 * KBM_N_HASH;
	par->ksampling = KBM_BSIZE;
	par->ktop  = 0.05;
	par->strand_mask = 3;
	par->self_aln = 0;
	par->skip_contained = 1;
	par->max_bgap = 4; // 4 * KBM_BIN_SIZE
	par->max_bvar = 4;
	par->max_bcnt = 0xFFFF;
	par->max_gap = 0.6;
	par->pgap = -7;
	par->pvar = -21;
	par->max_hit = 1000;
	par->min_aln = 2048 / KBM_BIN_SIZE;
	par->min_mat = 200;
	par->min_sim = 0.05;
	par->aln_var = 0.25;
	par->test_mode = 0;
	return par;
}

static inline void free_kbmpar(KBMPar *par){ free(par); }

static inline KBM* init_kbm(KBMPar *par){
	KBM *kbm;
	u4i i;
	kbm = malloc(sizeof(KBM));
	kbm->flags = 0;
	kbm->par = par;
	kbm->rdseqs = init_basebank();
	kbm->avail_bases = 0;
	kbm->avail_reads = 0;
	kbm->reads = init_kbmreadv(64);
	kbm->tag2idx = init_cuhash(1023);
	kbm->bins  = init_kbmbinv(64);
	kbm->binmarks = init_bitvec(1024);
	//kbm->kfs = NULL;
	kbm->seeds = init_kbmbmerv(64);
	kbm->sauxs = init_kbmbauxv(64);
	for(i=0;i<KBM_N_HASH;i++) kbm->hashs[i] = init_kbmhash(1023);
	for(i=0;i<KBM_N_HASH;i++) kbm->kauxs[i] = init_kbmkauxv(64);
	return kbm;
}

static inline void free_kbm(KBM *kbm){
	u4i i;
	if(kbm->flags & 0x1) return;
	if(kbm->flags & 0x2){
	} else {
		free_basebank(kbm->rdseqs);
		for(i=0;i<kbm->reads->size;i++){
			if(kbm->reads->buffer[i].tag) free(kbm->reads->buffer[i].tag);
		}
		free_kbmreadv(kbm->reads);
	}
	free_cuhash(kbm->tag2idx);
	free_kbmbinv(kbm->bins);
	free_bitvec(kbm->binmarks);
	//if(kbm->kfs) free(kbm->kfs);
	free_kbmbmerv(kbm->seeds);
	free_kbmbauxv(kbm->sauxs);
	for(i=0;i<KBM_N_HASH;i++) free_kbmhash(kbm->hashs[i]);
	for(i=0;i<KBM_N_HASH;i++) free_kbmkauxv(kbm->kauxs[i]);
	free(kbm);
}

static inline void reset_index_kbm(KBM *kbm){
	u4i i;
	if(kbm->flags & 0x04){
		free_kbmhash(kbm->hashs[0]);
		kbm->hashs[0] = init_kbmhash(1023);
		free_kbmkauxv(kbm->kauxs[0]);
		kbm->kauxs[0] = init_kbmkauxv(64);
	} else {
		for(i=0;i<KBM_N_HASH;i++){
			free_kbmhash(kbm->hashs[i]);
			kbm->hashs[i] = init_kbmhash(1023);
			free_kbmkauxv(kbm->kauxs[i]);
			kbm->kauxs[i] = init_kbmkauxv(64);
		}
	}
	free_kbmbmerv(kbm->seeds);
	kbm->seeds = init_kbmbmerv(64);
	free_kbmbauxv(kbm->sauxs);
	kbm->sauxs = init_kbmbauxv(64);
	//zeros_bitvec(kbm->binmarks);
}

static inline void clear_kbm(KBM *kbm){
	u4i i;
	if(kbm->flags & 0x1) return;
	if(kbm->flags & 0x2){
	} else {
		clear_basebank(kbm->rdseqs);
		for(i=0;i<kbm->reads->size;i++){
			if(kbm->reads->buffer[i].tag) free(kbm->reads->buffer[i].tag);
		}
		clear_kbmreadv(kbm->reads);
		clear_cuhash(kbm->tag2idx);
	}
	clear_kbmbinv(kbm->bins);
	clear_bitvec(kbm->binmarks);
	for(i=0;i<KBM_N_HASH;i++){
		free_kbmhash(kbm->hashs[i]);
		kbm->hashs[i] = init_kbmhash(1023);
		free_kbmkauxv(kbm->kauxs[i]);
		kbm->kauxs[i] = init_kbmkauxv(64);
		if(kbm->flags & 0x04) break;
	}
	clear_kbmbmerv(kbm->seeds);
	clear_kbmbauxv(kbm->sauxs);
}

#define kbm_cvt_length(seqlen) ((seqlen) & (MAX_U8 << KBM_BIN_BITS))

static inline void push_kbm(KBM *kbm, char *tag, int taglen, char *seq, u4i seqlen){
	kbm_read_t *rd;
	char *ptr;
	if(taglen){
		ptr = malloc(taglen + 1);
		memcpy(ptr, tag, taglen);
		ptr[taglen] = 0;
	} else {
		ptr = NULL;
	}
	if(seqlen > KBM_MAX_RDLEN) seqlen = KBM_MAX_RDLEN;
	seqlen = kbm_cvt_length(seqlen);
	rd = next_ref_kbmreadv(kbm->reads);
	rd->seqoff = (kbm->rdseqs->size >> KBM_BIN_BITS);
	rd->flag   = 0;
	rd->binoff = 0;
	rd->bincnt = seqlen / KBM_BIN_SIZE;
	rd->closed = 0;
	rd->tag = ptr;
	seq2basebank(kbm->rdseqs, seq, seqlen);
}

static inline void bitpush_kbm(KBM *kbm, char *tag, int taglen, u8i *seqs, u8i seqoff, u4i seqlen){
	kbm_read_t *rd;
	char *ptr;
	if(taglen){
		ptr = malloc(taglen + 1);
		memcpy(ptr, tag, taglen);
		ptr[taglen] = 0;
	} else {
		ptr = NULL;
	}
	if(seqlen > KBM_MAX_RDLEN) seqlen = KBM_MAX_RDLEN;
	seqlen = kbm_cvt_length(seqlen);
	rd = next_ref_kbmreadv(kbm->reads);
	rd->seqoff = (kbm->rdseqs->size >> KBM_BIN_BITS);
	rd->flag   = 0;
	rd->binoff = 0;
	rd->bincnt = seqlen / KBM_BIN_SIZE;
	rd->closed = 0;
	rd->tag = ptr;
	fast_fwdbits2basebank(kbm->rdseqs, seqs, seqoff, seqlen);
}

// Please call no more than once
static inline u8i filter_reads_kbm(KBM *kbm, u8i retain_size, int strategy){
	u8i m, beg, end, len;
	if(kbm->reads->size == 0) return 0;
	if(retain_size == 0 || retain_size >= kbm->rdseqs->size) return kbm->rdseqs->size;
	if((kbm->flags & 0x2) == 0){
		if(kbm->par->rd_len_order){
			sort_array(kbm->reads->buffer, kbm->reads->size, kbm_read_t, num_cmpgt(b.bincnt, a.bincnt));
			if(strategy == 0){ // longest
				len = 0;
				for(end=0;end<kbm->reads->size;end++){
					len += kbm->reads->buffer[end].bincnt * KBM_BIN_SIZE;
					if(len >= retain_size) break;
				}
				for(m=end;m<kbm->reads->size;m++) free(kbm->reads->buffer[m].tag);
				kbm->reads->size = end;
			} else if(strategy == 1){ // median
				m = kbm->reads->size / 2;
				len = kbm->reads->buffer[m].bincnt * KBM_BIN_SIZE;
				end = m;
				for(beg=0;beg<=m&&len<retain_size;beg++){
					len += kbm->reads->buffer[m - beg].bincnt * KBM_BIN_SIZE;
					len += kbm->reads->buffer[m + beg].bincnt * KBM_BIN_SIZE;
				}
				end = beg * 2;
				beg = m - beg;
				if(beg){
					remove_array_kbmreadv(kbm->reads, 0, beg);
				}
				for(m=end;m<kbm->reads->size;m++) free(kbm->reads->buffer[m].tag);
				kbm->reads->size = end;
			} else {
				return kbm->rdseqs->size;
			}
			return len;
		} else {
			return kbm->rdseqs->size;
		}
	} else {
		return kbm->rdseqs->size;
	}
}

static inline void ready_kbm(KBM *kbm){
	kbm_read_t *rd;
	u4i i, j;
	if((kbm->flags & 0x2) == 0){
		if(kbm->par->rd_len_order){
			sort_array(kbm->reads->buffer, kbm->reads->size, kbm_read_t, num_cmpgt(b.bincnt, a.bincnt));
		}
		encap_basebank(kbm->rdseqs, KBM_BSIZE);
	}
	clear_kbmbinv(kbm->bins);
	for(i=0;i<kbm->reads->size;i++){
		rd = ref_kbmreadv(kbm->reads, i);
		if(rd->tag) put_cuhash(kbm->tag2idx, (cuhash_t){rd->tag, i});
		if((kbm->flags & 0x2) == 0) rd->binoff = kbm->bins->size;
		for(j=0;j<rd->bincnt;j++){
			push_kbmbinv(kbm->bins, (kbm_bin_t){i, j, 0, 0});
		}
		if((kbm->flags & 0x2) == 0) rd->bincnt = j;
	}
	clear_bitvec(kbm->binmarks);
	encap_bitvec(kbm->binmarks, kbm->bins->size);
	kbm->binmarks->n_bit = kbm->bins->size;
	zeros_bitvec(kbm->binmarks);
}

// Share seqs, reads
static inline KBM* clone_seqs_kbm(KBM *src, KBMPar *par){
	KBM *dst;
	dst = init_kbm(par);
	free_basebank(dst->rdseqs);
	free_kbmreadv(dst->reads);
	dst->rdseqs = src->rdseqs;
	dst->reads  = src->reads;
	dst->flags  = 1LLU << 1;
	ready_kbm(dst); // Notice encap_basebank in ready_kbm
	return dst;
}

// rs[0]->n_head MUST >= 1
static inline void split_FIXP_kmers_kbm(BaseBank *rdseqs, u8i offset, u4i length, u1i ksize, u1i psize, u4i kmod, kmeroffv *rs[2]){
	kmer_off_t *kp;
	u8i kmer, krev, hv, npz, kmask, p, pmask;
	u4i ki, npl;
	int i, j;
	u1i c, b, kshift, pshift, kpshf;
	clear_kmeroffv(rs[0]);
	clear_kmeroffv(rs[1]);
	kshift = (ksize - 1) << 1;
	kmask = (1LLU << (ksize << 1)) - 1LLU;
	pshift = (psize - 1) << 1;
	pmask  = (1LLU << (psize << 1)) - 1LLU;
	kpshf  =  psize << 1;
	if(ksize){
		// scan F-part of kmers
		kmer = krev = 0;
		for(i=0;i+1<ksize;i++){
			c = bits2bit(rdseqs->bits, offset + i);
			kmer = (kmer << 2) | c;
			krev = (krev >> 2) | (((u8i)((~c) & 0x3)) << kshift);
		}
		for(;i<(int)length;i++){
			c = bits2bit(rdseqs->bits, offset + i);
			kmer = ((kmer << 2) | c) & kmask;
			krev = (krev >> 2) | (((u8i)((~c) & 0x3)) << kshift);
			if(kmer < krev){
				push_kmeroffv(rs[0], (kmer_off_t){kmer, i - (ksize - 1), 0, 0, 1});
			} else if(krev < kmer){
				push_kmeroffv(rs[1], (kmer_off_t){krev, i - (ksize - 1), 0, 1, 1});
			}
		}
		if(psize){
			// scan P-part of forward kmers
			assert(rs[0]->n_head > 0);
			memset(rs[0]->buffer - 1, 0, sizeof(kmer_off_t)); // rs[0]->n_head > 0
			rs[0]->buffer[-1].off = length;
			npz = 0;
			npl = 0;
			kp = rs[0]->buffer + rs[0]->size - 1;
			b = 4;
			for(i=length-1;i>=0;i--){
				if(kp->off + (ksize - 1) == (u4i)i){
					if(npl >= psize){
						p = npz & pmask;
						kp->closed = 0;
						kp->kmer = (kp->kmer << kpshf) | p;
					}
					kp --;
				}
				c = bits2revbit(rdseqs->bits, offset + i);
				if(c == b){
				} else {
					npz = (npz << 2) | c;
					npl ++;
					b = c;
				}
			}
			// scan P-part of reverse kmers
			encap_kmeroffv(rs[1], 1); memset(rs[1]->buffer + rs[1]->size, 0xFF, sizeof(kmer_off_t));
			npz = 0;
			npl = 0;
			kp = rs[1]->buffer;
			b = 4;
			for(i=0;i<(int)length;i++){
				if(kp->off == (u4i)(i)){
					if(npl >= psize){
						p = npz & pmask;
						kp->closed = 0;
						kp->kmer = (kp->kmer << kpshf) | p;
						kp->off  = kp->off - psize;
					}
					kp ++;
				}
				c = bits2bit(rdseqs->bits, offset + i);
				if(c == b){
				} else {
					npz = (npz << 2) | c;
					npl ++;
					b = c;
				}
			}
		} else {
			for(i=0;(u4i)i<rs[0]->size;i++) rs[0]->buffer[i].closed = 0;
			for(i=0;(u4i)i<rs[1]->size;i++) rs[1]->buffer[i].closed = 0;
		}
	} else if(psize){
		kmer = krev = 0;
		b = 4;
		for(i=j=0;i<(int)length;i++){
			c = bits2bit(rdseqs->bits, offset + i);
			if(b == c) continue;
			b = c;
			kmer = ((kmer << 2) | c) & pmask;
			krev = (krev >> 2) | (((u8i)((~c) & 0x3)) << pshift);
			j ++;
			if(j < psize) continue;
			if(kmer < krev){
				push_kmeroffv(rs[0], (kmer_off_t){kmer, i - (psize - 1), 0, 0, 0});
			} else if(krev < kmer){
				push_kmeroffv(rs[1], (kmer_off_t){krev, i - (psize - 1), 0, 1, 0});
			}
		}
	}
	for(b=0;b<2;b++){
		for(i=0;(u4i)i<rs[b]->size;i++){
			if(rs[b]->buffer[i].closed) continue;
			hv = kbm_kmer_smear(rs[b]->buffer[i].kmer);
			ki = hv % kmod;
			if(ki >= KBM_N_HASH) rs[b]->buffer[i].closed = 1;
			rs[b]->buffer[i].kidx = ki;
		}
	}
}

static inline u8i seed2solid_idx_kbm(KBM *kbm, kbm_dpe_t *p){
	kbm_bin_t *b;
	kbm_read_t *rd;
	u8i seqoff;
	b = kbm->bins->buffer + p->bidx;
	rd = kbm->reads->buffer + b->ridx;
	seqoff = (((rd->seqoff + b->off) * KBM_BSIZE) >> 1) + p->koff;
	return seqoff;
}

static inline u8i rdoff2solid_idx_kbm(KBM *kbm, u4i ridx, u4i roff){
	kbm_read_t *rd;
	u8i seqoff;
	rd = kbm->reads->buffer + ridx;
	seqoff = (rd->seqoff * KBM_BIN_SIZE + roff) >> 1;
	return seqoff;
}

#define binoff2solid_koff_kbm(kbm, bidx, boff) ((boff) >> 1)

typedef struct {
	u8i mer:50, kidx:14;
	u8i bidx;
	u4i cnt:22, koff:8, dir:1, used:1;
} kbm_midx_t;
define_list(kbmmidxv, kbm_midx_t);

typedef struct {
	u8i bidx;
	kbm_baux_t aux;
} kbm_tmp_bmer_t;
define_list(tmpbmerv, kbm_tmp_bmer_t);

thread_beg_def(midx);
KBM *kbm;
u8i beg, end; // (end - beg) * KBM_BSIZE MUST <= KBM_KMEROFF_MAX
u8i ktot, nrem, Nrem, none, nflt, offset;
u8i srem, Srem;
u8i *cnts, n_cnt;
int task;
int cal_degree;
pthread_mutex_t *locks;
thread_end_def(midx);

thread_beg_func(midx);
KBM *kbm;
kbm_bin_t *bin;
kbmmidxv **kidxs;
kbm_midx_t *mx;
kmer_off_t *f;
kbm_kmer_t *u;
kbm_kaux_t *x;
kmeroffv *kmers[2];
tmpbmerv *bms;
u8i bidx, off;
u4i i, j, k, len, ncpu, tidx, kidx;
int exists;
kbm = midx->kbm;
ncpu = midx->n_cpu;
tidx = midx->t_idx;
kmers[0] = adv_init_kmeroffv(64, 0, 1);
kmers[1] = adv_init_kmeroffv(64, 0, 1);
kidxs = malloc(KBM_N_HASH * sizeof(kbmmidxv*));
for(i=0;i<KBM_N_HASH;i++) kidxs[i] = init_kbmmidxv(64);
bms = init_tmpbmerv(KBM_MAX_KCNT);
thread_beg_loop(midx);
if(midx->task == 1){
	// counting kmers
	for(i=0;i<KBM_N_HASH;i++) clear_kbmmidxv(kidxs[i]);
	for(bidx=midx->beg+tidx;bidx<midx->end;bidx+=ncpu){
		if(KBM_LOG == 0 && tidx == 0 && ((bidx - midx->beg) % 100000) == 0){ fprintf(KBM_LOGF, "\r%llu", bidx - midx->beg); fflush(KBM_LOGF); }
		bin = ref_kbmbinv(kbm->bins, bidx);
		if(bin->closed) continue;
		off = (kbm->reads->buffer[bin->ridx].seqoff + bin->off) * KBM_BIN_SIZE;
		len = KBM_BIN_SIZE;
		split_FIXP_kmers_kbm(kbm->rdseqs, off, len, kbm->par->ksize, kbm->par->psize, kbm->par->kmer_mod, kmers);
		for(i=0;i<2;i++){
			for(j=0;j<kmers[i]->size;j++){
				f = ref_kmeroffv(kmers[i], j);
				if(f->closed) continue;
				mx = next_ref_kbmmidxv(kidxs[f->kidx]);
				mx->mer = f->kmer;
				mx->kidx = f->kidx;
				mx->bidx = bidx;
				mx->dir  = i;
				mx->koff = f->off;
				if(kidxs[f->kidx]->size >= 64){
					kidx = f->kidx;
					// lock hashs[kidx]
					pthread_mutex_lock(midx->locks + kidx);
					// hash adding
					for(k=0;k<kidxs[kidx]->size;k++){
						mx = ref_kbmmidxv(kidxs[kidx], k);
						u = prepare_kbmhash(kbm->hashs[kidx], mx->mer, &exists);
						if(exists){
							if(u->tot < KBM_MAX_KCNT) u->tot ++;
						} else {
							u->mer = mx->mer;
							u->tot = 1;
							u->flt = 0;
						}
					}
					// free hashs[f->kidx]
					pthread_mutex_unlock(midx->locks + kidx);
					clear_kbmmidxv(kidxs[kidx]);
				}
			}
		}
	}
	for(kidx=0;kidx<KBM_N_HASH;kidx++){
		if(kidxs[kidx]->size){
			// lock hashs[kidx]
			pthread_mutex_lock(midx->locks + kidx);
			// hash adding
			for(k=0;k<kidxs[kidx]->size;k++){
				mx = ref_kbmmidxv(kidxs[kidx], k);
				u = prepare_kbmhash(kbm->hashs[kidx], mx->mer, &exists);
				if(exists){
					if(u->tot < KBM_MAX_KCNT) u->tot ++;
				} else {
					u->mer = mx->mer;
					u->tot = 1;
					u->flt = 0;
				}
			}
			// free hashs[f->kidx]
			pthread_mutex_unlock(midx->locks + kidx);
			clear_kbmmidxv(kidxs[kidx]);
		}
	}
} else if(midx->task == 2){
	// delete low freq kmers
	for(i=tidx;i<KBM_N_HASH;i+=ncpu){
		reset_iter_kbmhash(kbm->hashs[i]);
		while((u = ref_iter_kbmhash(kbm->hashs[i]))){
			if(u->tot < kbm->par->kmin){
				delete_kbmhash(kbm->hashs[i], u);
			}
		}
	}
} else if(midx->task == 3){
	// stat kmer counts
	memset(midx->cnts, 0, midx->n_cnt * sizeof(u8i));
	for(i=tidx;i<KBM_N_HASH;i+=ncpu){
		reset_iter_kbmhash(kbm->hashs[i]);
		while((u = ref_iter_kbmhash(kbm->hashs[i]))){
			midx->cnts[num_min(midx->n_cnt, u->tot) - 1] ++;
		}
	}
} else if(midx->task == 4){
	// stat counts
	midx->offset = 0;
	for(i=tidx;i<KBM_N_HASH;i+=ncpu){
		reset_iter_kbmhash(kbm->hashs[i]);
		while((u = ref_iter_kbmhash(kbm->hashs[i]))){
			x = ref_kbmkauxv(kbm->kauxs[i], offset_kbmhash(kbm->hashs[i], u));
			x->off = midx->offset;
			x->cnt = 0;
			midx->ktot += u->tot;
			if(u->tot < kbm->par->kmin){
				u->flt = 1;
			} else if(u->tot > kbm->par->kmax){
				u->flt = 1;
			} else {
				midx->offset += u->tot;
			}
		}
	}
} else if(midx->task == 5){
	// revise offset
	for(i=tidx;i<KBM_N_HASH;i+=ncpu){
		for(off=0;off<kbm->kauxs[i]->size;off++) kbm->kauxs[i]->buffer[off].off += midx->offset;
	}
} else if(midx->task == 6){
	// fill seeds
	for(i=0;i<KBM_N_HASH;i++) clear_kbmmidxv(kidxs[i]);
	u4v *chgs;
	chgs = init_u4v(KBM_BSIZE);
	for(bidx=midx->beg+tidx;bidx<midx->end;bidx+=ncpu){
		if(KBM_LOG == 0 && tidx == 0 && ((bidx - midx->beg) % 100000) == 0){ fprintf(KBM_LOGF, "\r%llu", bidx - midx->beg); fflush(KBM_LOGF); }
		bin = ref_kbmbinv(kbm->bins, bidx);
		if(bin->closed) continue;
		bin->degree = 0;
		off = (kbm->reads->buffer[bin->ridx].seqoff + bin->off) * KBM_BIN_SIZE;
		len = KBM_BIN_SIZE;
		split_FIXP_kmers_kbm(kbm->rdseqs, off, len, kbm->par->ksize, kbm->par->psize, kbm->par->kmer_mod, kmers);
		clear_u4v(chgs);
		for(i=0;i<2;i++){
			for(j=0;j<kmers[i]->size;j++){
				f = ref_kmeroffv(kmers[i], j);
				if(f->closed) continue;
				mx = next_ref_kbmmidxv(kidxs[f->kidx]);
				mx->mer = f->kmer;
				mx->kidx = f->kidx;
				mx->bidx = bidx;
				mx->dir  = i;
				mx->koff = f->off;
				if(kidxs[f->kidx]->size == 64){
					push_u4v(chgs, f->kidx);
				}
			}
		}
		for(i=0;i<chgs->size;i++){
			kidx = chgs->buffer[i];
			if(kidxs[kidx]->size == 0) continue;
			pthread_mutex_lock(midx->locks + kidx);
			for(k=0;k<kidxs[kidx]->size;k++){
				mx = ref_kbmmidxv(kidxs[kidx], k);
				u = get_kbmhash(kbm->hashs[kidx], mx->mer);
				if(u && u->flt == 0){
					x = ref_kbmkauxv(kbm->kauxs[kidx], offset_kbmhash(kbm->hashs[kidx], u));
					kbm->bins->buffer[mx->bidx].degree ++;
					if(x->cnt < u->tot){
						if(x->cnt && getval_bidx(kbm, x->off + x->cnt - 1) == mx->bidx && kbm->sauxs->buffer[x->off + x->cnt - 1].dir == mx->dir){
							// repeated kmer within one bin
						} else {
							kbm->seeds->buffer[x->off + x->cnt].bidx = mx->bidx & MAX_U4;
							kbm->sauxs->buffer[x->off + x->cnt].bidx = mx->bidx >> 32;
							kbm->sauxs->buffer[x->off + x->cnt].dir  = mx->dir;
							kbm->sauxs->buffer[x->off + x->cnt].koff = mx->koff >> 1;
							x->cnt ++;
						}
					}
				}
			}
			pthread_mutex_unlock(midx->locks + kidx);
			clear_kbmmidxv(kidxs[kidx]);
		}
	}
	for(kidx=0;kidx<KBM_N_HASH;kidx++){
		if(kidxs[kidx]->size){
			// lock hashs[kidx]
			pthread_mutex_lock(midx->locks + kidx);
			// hash adding
			for(k=0;k<kidxs[kidx]->size;k++){
				mx = ref_kbmmidxv(kidxs[kidx], k);
				u = get_kbmhash(kbm->hashs[kidx], mx->mer);
				if(u && u->flt == 0){
					x = ref_kbmkauxv(kbm->kauxs[kidx], offset_kbmhash(kbm->hashs[kidx], u));
					kbm->bins->buffer[mx->bidx].degree ++;
					if(x->cnt < u->tot){
						if(x->cnt && getval_bidx(kbm, x->off + x->cnt - 1) == mx->bidx && kbm->sauxs->buffer[x->off + x->cnt - 1].dir == mx->dir){
							// repeated kmer within one bin
						} else {
							kbm->seeds->buffer[x->off + x->cnt].bidx = mx->bidx & MAX_U4;
							kbm->sauxs->buffer[x->off + x->cnt].bidx = mx->bidx >> 32;
							kbm->sauxs->buffer[x->off + x->cnt].dir  = mx->dir;
							kbm->sauxs->buffer[x->off + x->cnt].koff = mx->koff >> 1;
							x->cnt ++;
						}
					}
				}
			}
			// free hashs[f->kidx]
			pthread_mutex_unlock(midx->locks + kidx);
			clear_kbmmidxv(kidxs[kidx]);
		}
	}
	free_u4v(chgs);
} else if(midx->task == 7){
	// count added kmers
	midx->srem = midx->Srem = 0;
	for(i=tidx;i<KBM_N_HASH;i+=ncpu){
		reset_iter_kbmhash(kbm->hashs[i]);
		while((u = ref_iter_kbmhash(kbm->hashs[i]))){
			x = ref_kbmkauxv(kbm->kauxs[i], offset_kbmhash(kbm->hashs[i], u));
			if(x->cnt){
				midx->srem ++;
				midx->Srem += x->cnt;
			}
		}
	}
} else if(midx->task == 8){
	// sort seeds within a kmer
	for(i=tidx;i<KBM_N_HASH;i+=ncpu){
		reset_iter_kbmhash(kbm->hashs[i]);
		while((u = ref_iter_kbmhash(kbm->hashs[i]))){
			x = ref_kbmkauxv(kbm->kauxs[i], offset_kbmhash(kbm->hashs[i], u));
			if(x->cnt < 2) continue;
			clear_tmpbmerv(bms);
			for(j=0;j<x->cnt;j++){
				push_tmpbmerv(bms, (kbm_tmp_bmer_t){getval_bidx(kbm, x->off + j), kbm->sauxs->buffer[x->off + j]});
			}
			sort_array(bms->buffer, bms->size, kbm_tmp_bmer_t, num_cmpgt(a.bidx, b.bidx));
			kbm->seeds->buffer[x->off + 0].bidx = bms->buffer[0].bidx & MAX_U4;
			kbm->sauxs->buffer[x->off + 0].bidx = bms->buffer[0].bidx >> 32;
			kbm->sauxs->buffer[x->off + 0]      = bms->buffer[0].aux;
			len = 1;
			for(j=1;j<x->cnt;j++){
				if(bms->buffer[j].bidx < bms->buffer[j - 1].bidx){
					continue;
				}
				kbm->seeds->buffer[x->off + len].bidx = bms->buffer[j].bidx & MAX_U4;
				kbm->sauxs->buffer[x->off + len].bidx = bms->buffer[0].bidx >> 32;
				kbm->sauxs->buffer[x->off + len]      = bms->buffer[j].aux;
				len ++;
			}
			x->cnt = len;
		}
	}
}
thread_end_loop(midx);
free_kmeroffv(kmers[0]);
free_kmeroffv(kmers[1]);
for(i=0;i<KBM_N_HASH;i++) free_kbmmidxv(kidxs[i]);
free(kidxs);
free_tmpbmerv(bms);
thread_end_func(midx);

static inline void index_kbm(KBM *kbm, u8i beg, u8i end, u4i ncpu, FILE *kmstat){
	u8i ktyp, nflt, nrem, Nrem, none, ktot, srem, Srem, off, cnt, *kcnts, MAX;
	u8i i, b, e, n;
	u4i j, kavg, batch_size;
	pthread_mutex_t *hash_locks;
	thread_preprocess(midx);
	batch_size = 10000;
	clear_kbmbmerv(kbm->seeds);
	for(i=0;i<KBM_N_HASH;i++) clear_kbmhash(kbm->hashs[i]);
	kcnts = NULL;
	MAX = KBM_MAX_KCNT;
	hash_locks = calloc(KBM_N_HASH, sizeof(pthread_mutex_t));
	thread_beg_init(midx, ncpu);
	midx->kbm = kbm;
	midx->beg  = beg;
	midx->end  = end;
	midx->cnts = NULL;
	midx->n_cnt = MAX;
	midx->task = 0;
	midx->cal_degree = 0;
	midx->locks = hash_locks;
	thread_end_init(midx);
	fprintf(KBM_LOGF, "[%s] - scanning kmers (K%dP%dS%0.2f) from %llu bins\n", date(), kbm->par->ksize, kbm->par->psize, 1.0 * kbm->par->kmer_mod / KBM_N_HASH, end - beg);
	b = e = beg;
	thread_apply_all(midx, midx->task = 1);
	if(KBM_LOG == 0){ fprintf(KBM_LOGF, "\r%llu bins\n", end - beg); fflush(KBM_LOGF); }
	kcnts = calloc(MAX, sizeof(u8i));
	thread_beg_iter(midx);
	midx->cnts = calloc(MAX, sizeof(u8i));
	midx->task = 3; // counting raw kmers
	thread_wake(midx);
	thread_end_iter(midx);
	thread_beg_iter(midx);
	thread_wait(midx);
	for(i=0;i<MAX;i++) kcnts[i] += midx->cnts[i];
	thread_end_iter(midx);
	if(kmstat){
		fprintf(kmstat, "#Reads: %llu\n", (u8i)kbm->reads->size);
		fprintf(kmstat, "#Bases: %llu bp\n", (u8i)kbm->rdseqs->size);
		fprintf(kmstat, "#K%dP%dS%0.2f\n", kbm->par->ksize, kbm->par->psize, 1.0 * kbm->par->kmer_mod / KBM_N_HASH);
		for(i=0;i+1<MAX;i++){
			fprintf(kmstat, "%llu\t%llu\t%llu\n", i + 1, kcnts[i], (i + 1) * kcnts[i]);
		}
		fflush(kmstat);
	}
	if(kmstat){
		u8i *_kcnts;
		_kcnts = malloc(200 * sizeof(u8i));
		for(i=0;i<200;i++){
			_kcnts[i] = (i + 1) * kcnts[i];
		}
		char *txt = barplot_txt_u8_simple(100, 20, _kcnts, 200, 0);
		fprintf(KBM_LOGF, "********************** Kmer Frequency **********************\n");
		fputs(txt, KBM_LOGF);
		fprintf(KBM_LOGF, "**********************     1 - 201    **********************\n");
		ktyp = 0;
		for(i=0;i<MAX;i++){
			ktyp += (i + 1) * kcnts[i];
		}
		_kcnts[0] = 0.10 * ktyp;
		_kcnts[1] = 0.20 * ktyp;
		_kcnts[2] = 0.30 * ktyp;
		_kcnts[3] = 0.40 * ktyp;
		_kcnts[4] = 0.50 * ktyp;
		_kcnts[5] = 0.60 * ktyp;
		_kcnts[6] = 0.70 * ktyp;
		_kcnts[7] = 0.80 * ktyp;
		_kcnts[8] = 0.90 * ktyp;
		_kcnts[9] = 0.95 * ktyp;
		fprintf(KBM_LOGF, "Quatiles:\n");
		fprintf(KBM_LOGF, "   10%%   20%%   30%%   40%%   50%%   60%%   70%%   80%%   90%%   95%%\n");
		off = 0;
		for(i=j=0;j<10;j++){
			while(off < _kcnts[j] && i < MAX){
				off += kcnts[i] * (i + 1);
				i ++;
			}
			fprintf(KBM_LOGF, "%6llu", i);
		}
		fprintf(KBM_LOGF, "\n");
		//fprintf(KBM_LOGF,
			//"# If the kmer distribution is not good, please kill me and adjust -k, -p, and -K\n"
			//"# Cannot get a good distribution anyway, should adjust -S -s, also -A -e in assembly\n"
		//);
		free(_kcnts);
		free(txt);
	}
	// delete low freq kmer from hash
	thread_apply_all(midx, midx->task = 2);
	// freeze hash to save memory and speed up the query
	for(i=0;i<KBM_N_HASH;i++){
		if(0){
			fprintf(KBM_LOGF, "%12llu ", (u8i)kbm->hashs[i]->count); fflush(KBM_LOGF);
			if((i % 8) == 7){
				fprintf(KBM_LOGF, "\n"); fflush(KBM_LOGF);
			}
		}
		freeze_kbmhash(kbm->hashs[i], 1.0 / 16);
		free_kbmkauxv(kbm->kauxs[i]);
		kbm->kauxs[i] = init_kbmkauxv(kbm->hashs[i]->count);
		kbm->kauxs[i]->size = kbm->hashs[i]->count;
	}
	print_proc_stat_info(0);
	ktot = nrem = Nrem = none = nflt = ktyp = 0;
	for(i=0;i+1<MAX;i++){
		ktot += kcnts[i];
		if(i + 1 < kbm->par->kmin){
			none += kcnts[i];
		} else {
			ktyp += kcnts[i] * (i + 1);
		}
	}
	if(kbm->par->ktop){
		off = 0;
		for(i=MAX;i>kbm->par->kmin;i--){
			off += kcnts[i-1] * i;
			if(off >= (ktyp * kbm->par->ktop)) break;
		}
		if(i > kbm->par->kmax){
			kbm->par->kmax = i;
		}
		fprintf(KBM_LOGF, "[%s] - high frequency kmer depth is set to %d\n", date(), kbm->par->kmax);
	}
	for(i=kbm->par->kmin? kbm->par->kmin - 1 : 0;i<kbm->par->kmax;i++){
		nrem += kcnts[i];
		Nrem += kcnts[i] * (i + 1);
	}
	for(i=kbm->par->kmax;i+1<MAX;i++){
		nflt += kcnts[i];
	}
	off = 0;
	thread_apply_all(midx, midx->task = 4);
	thread_beg_iter(midx);
	cnt = midx->offset;
	midx->offset = off;
	off += cnt;
	midx->task = 5;
	thread_wake(midx);
	thread_end_iter(midx);
	thread_wait_all(midx);
	clear_and_encap_kbmbmerv(kbm->seeds, off + 1);
	kbm->seeds->size = off;
	free_kbmbauxv(kbm->sauxs);
	kbm->sauxs = init_kbmbauxv(off + 1);
	kbm->sauxs->size = off;
	kavg = Nrem / (nrem + 1);
	fprintf(KBM_LOGF, "[%s] - Total kmers = %llu\n", date(), ktot);
	fprintf(KBM_LOGF, "[%s] - average kmer depth = %d\n", date(), kavg);
	fprintf(KBM_LOGF, "[%s] - %llu low frequency kmers (<%d)\n", date(), none, kbm->par->kmin);
	fprintf(KBM_LOGF, "[%s] - %llu high frequency kmers (>%d)\n", date(), nflt, kbm->par->kmax);
	fprintf(KBM_LOGF, "[%s] - indexing %llu kmers, %llu instances (at most)\n", date(), nrem, Nrem);
	thread_apply_all(midx, midx->task = 6);
	if(KBM_LOG == 0){ fprintf(KBM_LOGF, "\r%llu bins\n", end - beg); fflush(KBM_LOGF); }
	thread_apply_all(midx, midx->task = 7);
	srem = Srem = 0;
	thread_beg_iter(midx);
	srem += midx->srem;
	Srem += midx->Srem;
	thread_end_iter(midx);
	fprintf(KBM_LOGF, "[%s] - indexed  %llu kmers, %llu instances\n", date(), srem, Srem);
	{
		n = 0;
		for(i=beg;i<end;i++){
			if(kbm->bins->buffer[i].degree < kbm->par->min_bin_degree){
				kbm->bins->buffer[i].closed = 1;
				one_bitvec(kbm->binmarks, i);
				n ++;
			}
		}
		index_bitvec(kbm->binmarks);
		fprintf(KBM_LOGF, "[%s] - masked %llu bins as closed\n", date(), n);
	}
	fprintf(KBM_LOGF, "[%s] - sorting\n", date());
	thread_apply_all(midx, midx->task = 8);
	if(kcnts) free(kcnts);
	thread_beg_close(midx);
	if(midx->cnts) free(midx->cnts);
	thread_end_close(midx);
	free(hash_locks);
	print_proc_stat_info(0);
}

static inline void simple_index_kbm(KBM *kbm, u8i beg, u8i end){
	kbm_bin_t *bin;
	kbm_kmer_t *u;
	kbm_kaux_t *x;
	kmeroffv *kmers[2];
	tmpbmerv *bms;
	kmer_off_t *f;
	u8i bidx, off;
	u4i i, j, len;
	int exists;
	kbm->flags |= 1LLU << 2;
	kmers[0] = adv_init_kmeroffv(64, 0, 1);
	kmers[1] = adv_init_kmeroffv(64, 0, 1);
	bms = init_tmpbmerv(KBM_MAX_KCNT);
	// count kmers
	{
		for(bidx=beg;bidx<end;bidx++){
			bin = ref_kbmbinv(kbm->bins, bidx);
			if(bin->closed) continue;
			off = (kbm->reads->buffer[bin->ridx].seqoff + bin->off) * KBM_BIN_SIZE;
			split_FIXP_kmers_kbm(kbm->rdseqs, off, KBM_BIN_SIZE, kbm->par->ksize, kbm->par->psize, kbm->par->kmer_mod, kmers);
			for(i=0;i<2;i++){
				for(j=0;j<kmers[i]->size;j++){
					f = ref_kmeroffv(kmers[i], j);
					if(f->closed) continue;
					u = prepare_kbmhash(kbm->hashs[0], f->kmer, &exists);
					if(exists){
						if(u->tot < KBM_MAX_KCNT) u->tot ++;
					} else {
						u->mer = f->kmer;
						u->tot = 1;
						u->flt = 0;
					}
				}
			}
		}
	}
	// delete low freq kmers
	if(kbm->par->kmin){
		reset_iter_kbmhash(kbm->hashs[0]);
		while((u = ref_iter_kbmhash(kbm->hashs[0]))){
			if(u->tot < kbm->par->kmin){
				delete_kbmhash(kbm->hashs[0], u);
			}
		}
	}
	// freeze hash
	{
		freeze_kbmhash(kbm->hashs[0], 1.0 / 16);
		free_kbmkauxv(kbm->kauxs[0]);
		kbm->kauxs[0] = init_kbmkauxv(kbm->hashs[0]->count);
		kbm->kauxs[0]->size = kbm->hashs[0]->count;
	}
	// encap seeds
	{
		off = 0;
		reset_iter_kbmhash(kbm->hashs[0]);
		while((u = ref_iter_kbmhash(kbm->hashs[0]))){
			x = ref_kbmkauxv(kbm->kauxs[0], offset_kbmhash(kbm->hashs[0], u));
			x->off = off;
			x->cnt = 0;
			if(u->tot < kbm->par->kmin){ // Never happens
				u->flt = 1;
			} else if(kbm->par->kmax && u->tot > kbm->par->kmax){
				u->flt = 1;
			} else {
				off += u->tot;
			}
		}
		clear_and_encap_kbmbmerv(kbm->seeds, off + 1);
		kbm->seeds->size = off;
		clear_and_encap_kbmbauxv(kbm->sauxs, off + 1);
		kbm->sauxs->size = off;
	}
	// fill seeds
	{
		for(bidx=beg;bidx<end;bidx++){
			bin = ref_kbmbinv(kbm->bins, bidx);
			if(bin->closed) continue;
			off = (kbm->reads->buffer[bin->ridx].seqoff + bin->off) * KBM_BIN_SIZE;
			split_FIXP_kmers_kbm(kbm->rdseqs, off, KBM_BIN_SIZE, kbm->par->ksize, kbm->par->psize, kbm->par->kmer_mod, kmers);
			for(i=0;i<2;i++){
				for(j=0;j<kmers[i]->size;j++){
					f = ref_kmeroffv(kmers[i], j);
					if(f->closed) continue;
					u = get_kbmhash(kbm->hashs[0], f->kmer);
					if(u == NULL || u->flt) continue;
					x = ref_kbmkauxv(kbm->kauxs[0], offset_kbmhash(kbm->hashs[0], u));
					kbm->bins->buffer[bidx].degree ++;
					if(x->cnt < u->tot){
						if(x->cnt && getval_bidx(kbm, x->off + x->cnt - 1) == bidx && kbm->sauxs->buffer[x->off + x->cnt - 1].dir == i){
							// repeated kmer within one bin
						} else {
							kbm->seeds->buffer[x->off + x->cnt].bidx = bidx & MAX_U4;
							kbm->sauxs->buffer[x->off + x->cnt].bidx = bidx >> 32;
							kbm->sauxs->buffer[x->off + x->cnt].dir  = i;
							kbm->sauxs->buffer[x->off + x->cnt].koff = f->off >> 1;
							x->cnt ++;
						}
					}
				}
			}
		}
	}
	// sort seeds within a kmer
	{
		reset_iter_kbmhash(kbm->hashs[0]);
		while((u = ref_iter_kbmhash(kbm->hashs[0]))){
			x = ref_kbmkauxv(kbm->kauxs[0], offset_kbmhash(kbm->hashs[0], u));
			clear_tmpbmerv(bms);
			for(j=0;j<x->cnt;j++){
				push_tmpbmerv(bms, (kbm_tmp_bmer_t){getval_bidx(kbm, x->off + j), kbm->sauxs->buffer[x->off + j]});
			}
			sort_array(bms->buffer, bms->size, kbm_tmp_bmer_t, num_cmpgt(a.bidx, b.bidx));
			kbm->seeds->buffer[x->off + 0].bidx = bms->buffer[0].bidx & MAX_U4;
			kbm->sauxs->buffer[x->off + 0].bidx = bms->buffer[0].bidx >> 32;
			kbm->sauxs->buffer[x->off + 0]      = bms->buffer[0].aux;
			len = 1;
			for(j=1;j<x->cnt;j++){
				if(bms->buffer[j].bidx < bms->buffer[j - 1].bidx){
					continue;
				}
				kbm->seeds->buffer[x->off + len].bidx = bms->buffer[j].bidx & MAX_U4;
				kbm->sauxs->buffer[x->off + len].bidx = bms->buffer[0].bidx >> 32;
				kbm->sauxs->buffer[x->off + len]      = bms->buffer[j].aux;
				len ++;
			}
			x->cnt = len;
		}
	}
	free_kmeroffv(kmers[0]);
	free_kmeroffv(kmers[1]);
	free_tmpbmerv(bms);
}

static inline KBMDP* init_kbmdp(){
	KBMDP *dp;
	dp = malloc(sizeof(KBMDP));
	dp->kms = init_kbmdpev(1024);
	dp->km_len = 0;
	dp->cmask = init_bitvec(1024);
	dp->cms = init_kbmcmerv(64);
	dp->coffs = init_u4v(32);
	dp->rmask[0] = init_bitvec(256);
	dp->rmask[1] = init_bitvec(256);
	dp->cells[0] = init_kbmcellv(16);
	dp->cells[1] = init_kbmcellv(16);
	dp->bts = init_bit2vec(1024);
	dp->paths = init_kbmphash(13);
	dp->boff = 0;
	dp->last_bidx = 0;
	return dp;
}

static inline void reset_kbmdp(KBMDP *dp, KBMAux *aux, u8i bidx){
	//clear_kbmdpev(dp->kms);
	//dp->km_len = 0;
	clear_bitvec(dp->cmask);
	clear_kbmcmerv(dp->cms);
	clear_u4v(dp->coffs);
	recap_bitvec(dp->rmask[0], aux->qnbit);
	zeros_bitvec(dp->rmask[0]);
	recap_bitvec(dp->rmask[1], aux->qnbit);
	zeros_bitvec(dp->rmask[1]);
	recap_kbmcellv(dp->cells[0], aux->qnbin);
	recap_kbmcellv(dp->cells[1], aux->qnbin);
	clear_bit2vec(dp->bts);
	if(dp->paths->size > 1024 * 10){
		free_kbmphash(dp->paths);
		dp->paths = init_kbmphash(1023);
	} else {
		clear_kbmphash(dp->paths);
	}
	dp->boff = bidx;
	dp->last_bidx = bidx;
}

static inline void clear_kbmdp(KBMDP *dp, KBMAux *aux, u8i bidx){
	reset_kbmdp(dp, aux , bidx);
	clear_kbmdpev(dp->kms);
	dp->km_len = 0;
}

static inline void free_kbmdp(KBMDP *dp){
	free_kbmdpev(dp->kms);
	free_bitvec(dp->cmask);
	free_kbmcmerv(dp->cms);
	free_u4v(dp->coffs);
	free_bitvec(dp->rmask[0]);
	free_bitvec(dp->rmask[1]);
	free_kbmcellv(dp->cells[0]);
	free_kbmcellv(dp->cells[1]);
	free_bit2vec(dp->bts);
	free_kbmphash(dp->paths);
	free(dp);
}

static inline KBMAux* init_kbmaux(KBM *kbm){
	KBMAux *aux;
	aux = malloc(sizeof(KBMAux));
	aux->kbm = kbm;
	aux->par = kbm->par;
	aux->qtag = NULL;
	aux->qseqs = NULL;
	aux->qsoff = 0;
	aux->qlen = 0;
	aux->slen = 0;
	aux->qidx = 0;
	aux->qnbin = 0;
	aux->qnbit = (aux->qnbin + 63) & 0xFFFFFFC0U;
	aux->bmin = 0;
	aux->bmax = MAX_U8;
	aux->koffs[0] = adv_init_kmeroffv(32, 0, 1);
	aux->koffs[1] = adv_init_kmeroffv(32, 0, 1);
	aux->refs = init_kbmrefv(64);
	aux->rank = init_u4v(64);
	aux->nheap = 1024;
	aux->heaps = malloc((aux->nheap + 1) * sizeof(u4v*));
	{
		u4i i;
		for(i=0;i<=aux->nheap;i++){
			aux->heaps[i] = init_u4v(8);
		}
	}
	aux->bmlen = 1;
	aux->bmoff = 0;
	aux->bmcnt = aux->kbm->bins->size;
	aux->binmap = NULL;
	aux->caches[0] = init_kbmdpev(64);
	aux->caches[1] = init_kbmdpev(64);
	aux->dps[0] = init_kbmdp();
	aux->dps[1] = init_kbmdp();
	aux->hits  = init_kbmmapv(16);
	aux->cigars = init_bitsvec(1024, 3);
	aux->solids = NULL;
	aux->str = init_string(1024);
	return aux;
}

static inline void free_kbmaux(KBMAux *aux){
	free_kmeroffv(aux->koffs[0]);
	free_kmeroffv(aux->koffs[1]);
	free_kbmrefv(aux->refs);
	free_u4v(aux->rank);
	if(aux->heaps){
		u4i i;
		for(i=0;i<=aux->nheap;i++){
			if(aux->heaps[i]) free_u4v(aux->heaps[i]);
		}
		free(aux->heaps);
	}
	if(aux->binmap) free(aux->binmap);
	free_kbmdpev(aux->caches[0]);
	free_kbmdpev(aux->caches[1]);
	free_kbmdp(aux->dps[0]);
	free_kbmdp(aux->dps[1]);
	free_kbmmapv(aux->hits);
	free_bitsvec(aux->cigars);
	free_string(aux->str);
	free(aux);
}

static inline void query_index_kbm(KBMAux *aux, char *qtag, u4i qidx, BaseBank *rdseqs, u8i seqoff, u4i seqlen){
	KBM *kbm;
	KBMPar *par;
	kbm_kmer_t *u;
	kbm_kaux_t *x;
	kmer_off_t *f;
	kbm_ref_t *ref;
	u8i sidx, bmin, bmax;
	u4i hidx, next;
	u4i i, j, l, tot, mr, pdir;
	kbm = aux->kbm;
	par = aux->par;
	aux->qtag  = qtag? qtag : kbm->reads->buffer[qidx].tag;
	aux->qseqs = rdseqs;
	aux->qsoff = seqoff;
	aux->qlen = seqlen;
	aux->slen = (seqlen / KBM_BIN_SIZE) * KBM_BIN_SIZE;
	aux->qidx = qidx;
	aux->qnbin = aux->slen / KBM_BIN_SIZE;
	aux->qnbit = (aux->qnbin + 63) & 0xFFFFFFC0U;
	clear_kbmdp(aux->dps[0], aux, 0);
	clear_kbmdp(aux->dps[1], aux, 0);
	clear_kbmrefv(aux->refs);
	clear_u4v(aux->rank);
	clear_kbmdpev(aux->caches[0]);
	clear_kbmdpev(aux->caches[1]);
	clear_kbmmapv(aux->hits);
	clear_bitsvec(aux->cigars);
#ifdef TEST_MODE
	if(par->test_mode >= 7) return;
#endif
	bmin = par->self_aln? kbm->reads->buffer[qidx].binoff + kbm->reads->buffer[qidx].bincnt : 0;
	if(bmin < aux->bmin) bmin = aux->bmin;
	bmax = aux->bmax;
	if(par->self_aln == 2){
		par->strand_mask = 2; // 1: forward; 2: reverse; 3: both
	}
	split_FIXP_kmers_kbm(rdseqs, seqoff, aux->slen, par->ksize, par->psize, par->kmer_mod, aux->koffs);
#ifdef TEST_MODE
	if(par->test_mode >= 6) return;
#endif
	tot = 0;
	for(i=0;i<2;i++){
		next = 0;
		for(j=0;j<aux->koffs[i]->size;j++){
			f = ref_kmeroffv(aux->koffs[i], j);
			if(f->closed) continue;
			if(kbm->flags & (1LLU << 2)) f->kidx = 0;
			u = get_kbmhash(kbm->hashs[f->kidx], f->kmer);
			if(u == NULL || u->flt || u->tot < par->kmin){
				continue;
			}
			x = ref_kbmkauxv(kbm->kauxs[f->kidx], offset_kbmhash(kbm->hashs[f->kidx], u));
			ref = next_ref_kbmrefv(aux->refs);
			ref->mer = u;
			ref->aux = x;
			ref->kidx = f->kidx;
			ref->off = f->off;
			ref->dir = i;
			if(par->self_aln && aux->solids){
				sidx = (kbm->reads->buffer[qidx].seqoff * KBM_BIN_SIZE + ref->off) >> 1;
				ref->fine = get_bitvec(aux->solids, sidx);
			} else {
				ref->fine = 0;
			}
			ref->qbidx = ref->off / KBM_BIN_SIZE;
			ref->poffs[0] = ref->off;
			ref->poffs[1] = aux->slen - (ref->off + (aux->par->ksize + aux->par->psize));
			ref->boff = x->off;
			ref->bend = x->off + x->cnt;
			ref->bidx = getval_bidx(kbm, x->off);
			ref->closed = 0;
			{
				// Refine boundray
				if(par->self_aln == 2){ // reverse complementary only
					while(ref->boff < ref->bend && getval_bidx(kbm, ref->bend - 1) > bmin){
						ref->bend --;
					}
				}
				while(ref->boff < ref->bend){
					if(ref->bidx < bmin){
						ref->boff ++;
						ref->bidx = getval_bidx(kbm, ref->boff);
						continue;
					//} else if(ref->b->bidx >= bmax){
						//break;
					}
					break;
				}
				while(ref->bend > ref->boff){
					if(getval_bidx(kbm, ref->bend - 1) > bmax){
						ref->bend --;
					} else {
						break;
					}
				}
				if(ref->boff >= ref->bend){
					ref->closed = 1;
				}
			}
			if(ref->closed){
				aux->refs->size --;
				continue;
			}
			tot += x->cnt;
		}
	}
	if(0){
		//sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t, num_cmpgt(a.off, b.off));
		for(i=0;i<aux->refs->size;i++){
			ref = ref_kbmrefv(aux->refs, i);
			fprintf(KBM_LOGF, "%s\t%d\t%c\t%d\n", aux->qtag, ref->off, "+-"[ref->dir], (int)ref->aux->cnt);
		}
	}
	if(par->self_aln && aux->solids){
		// Obsolete
		sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t, num_cmpgt(a.off, b.off));
		tot = 0;
		next = 0;
		for(i=0;i<aux->refs->size;i++){
			ref = ref_kbmrefv(aux->refs, i);
			if(ref->closed){
				continue;
			} else if(ref->fine){
				tot += ref->bend - ref->boff;
				next = ref->off + (aux->par->ksize + aux->par->psize) / 2 + 1;
			} else if(ref->off >= next){
				tot += ref->bend - ref->boff;
			} else {
				ref->boff = ref->bend;
				ref->closed = 1;
			}
		}
	} else if(aux->par->ksampling < KBM_BIN_SIZE && aux->refs->size){
		sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t, num_cmpgtx(a.qbidx, b.qbidx, b.bend - b.boff, a.bend - a.boff));
		tot = 0;
		for(i=j=0;i<aux->refs->size;i++){
			if(aux->refs->buffer[i].qbidx != aux->refs->buffer[j].qbidx){
				if(aux->refs->buffer[j].qbidx){ // skip the first and last bin
					if((i - j) > aux->par->ksampling){
						l = j + aux->par->ksampling;
						for(;j<l;j++){
							tot += aux->refs->buffer[j].bend - aux->refs->buffer[j].boff;
						}
						for(;j<i;j++){
							aux->refs->buffer[j].boff = aux->refs->buffer[j].bend;
							aux->refs->buffer[j].bidx = getval_bidx(kbm, aux->refs->buffer[j].boff);
							aux->refs->buffer[j].closed = 1;
						}
					}
				}
				j = i;
			}
		}
		//sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t, num_cmpgt(a.off, b.off));
	}
	sort_array(aux->refs->buffer, aux->refs->size, kbm_ref_t, num_cmpgt(a.off, b.off));
	// estimate binmap
	aux->bmoff = 0;
	if(aux->refs->size){
		mr = aux->par->min_mat / (aux->par->ksize + aux->par->psize);
		if(mr < 512) mr = 512;
		aux->bmlen = tot / mr;
		if(aux->bmlen == 0) aux->bmlen = 1;
		aux->bmcnt = (aux->kbm->bins->size + aux->bmlen - 1) / aux->bmlen;
		if(aux->bmcnt < aux->qnbin * 50){
			aux->bmcnt = aux->qnbin * 50;
			aux->bmlen = (aux->kbm->bins->size + aux->bmcnt - 1) / aux->bmcnt;
		}
	} else {
		aux->bmlen = 1;
		aux->bmcnt = aux->kbm->bins->size;
	}
	if(0 && aux->bmlen > aux->nheap){
		aux->bmlen = aux->nheap;
		aux->bmcnt = (aux->kbm->bins->size + aux->bmlen - 1) / aux->bmlen;
	}
	//fprintf(stderr, " -- %s tot=%d bmlen=%d bmcnt=%d mr=%d in %s -- %s:%d --\n", aux->qtag, tot, aux->bmlen, aux->bmcnt, tot / aux->bmlen, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	for(i=0;i<aux->nheap&&i<aux->bmlen;i++){
		clear_u4v(aux->heaps[i]);
	}
	//fprintf(stderr, " -- %s tot=%d avg=%d bmlen=%d bmcnt=%d mr=%d qnbin=%d in %s -- %s:%d --\n", aux->qtag, tot, tot / aux->bmlen, aux->bmlen, aux->bmcnt, mr, aux->qnbin, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
#ifdef TEST_MODE
		if(par->test_mode >= 5) return;
#endif
	// init heaps
	for(i=0;i<aux->refs->size;i++){
		ref = ref_kbmrefv(aux->refs, i);
		while(ref->boff < ref->bend){
			if(0){
				pdir = (ref->dir ^ kbm->sauxs->buffer[ref->boff].dir);
				if(((aux->par->strand_mask >> pdir) & 0x01) == 0){
					ref->boff ++;
					ref->bidx = getval_bidx(aux->kbm, ref->boff);
					continue;
				}
			}
			hidx = ref->bidx / aux->bmcnt;
			if(hidx - aux->bmoff < aux->nheap){
				push_u4v(aux->heaps[hidx - aux->bmoff], i);
			}
			break;
		}
	}
	aux->hptr = 0;
}

static inline void print_exists_index_kbm(KBM *kbm, char *qtag, BaseBank *rdseqs, u8i seqoff, u4i seqlen, kmeroffv *kmers[2], FILE *out){
	KBMPar *par;
	kbm_kmer_t *u;
	kbm_kaux_t *x;
	kmer_off_t *f;
	u4i i, j;
	par = kbm->par;
	split_FIXP_kmers_kbm(rdseqs, seqoff, seqlen, par->ksize, par->psize, par->kmer_mod, kmers);
	for(i=0;i<2;i++){
		for(j=0;j<kmers[i]->size;j++){
			f = ref_kmeroffv(kmers[i], j);
			if(f->closed) continue;
			if(kbm->flags & 0x04) f->kidx = 0;
			u = get_kbmhash(kbm->hashs[f->kidx], f->kmer);
			if(u == NULL){
				continue;
			}
			x = ref_kbmkauxv(kbm->kauxs[f->kidx], offset_kbmhash(kbm->hashs[f->kidx], u));
			fprintf(out, "%s\t%d\t%c\t0x%llx\t%u\t%u\t%c\n", qtag, f->off, "+-"[f->dir], f->kmer, x->cnt, u->tot, "YN"[u->flt]);
		}
	}
}

static inline int _update_dp_path_kbm(KBMDP *dp, u8i end, kbm_cell_t *c){
	int exists;
	kbm_path_t *p, P;
	P.beg = c->beg;
	p = prepare_kbmphash(dp->paths, P, &exists);
	if(exists){
		if(p->score < c->score){
			p->end = end;
			p->mat = c->mat;
			p->score = c->score;
		} else {
			return 0;
		}
	} else {
		p->beg = c->beg;
		p->end = end;
		p->mat = c->mat;
		p->score = c->score;
	}
	return 1;
}

static inline void _dp_cal_spare_row_kbm(KBMAux *aux, int dir){
	KBMDP *dp;
	kbmcellv *cells[2];
	BitVec *masks[2];
	kbm_cmer_t *m;
	kbm_cell_t D, H, V;
	u4i i, n, ni, rowoff;
	u8i bitoff, celoff;
	int flg, is_gap, score;
	dp = aux->dps[dir];
	flg = (dp->last_bidx - dp->boff) & 0x01;
	cells[0] = dp->cells[flg];
	cells[1] = dp->cells[!flg];
	masks[0] = dp->rmask[flg];
	masks[1] = dp->rmask[!flg];
	rowoff = dp->coffs->buffer[dp->last_bidx - dp->boff];
	bitoff = (dp->last_bidx - dp->boff) * aux->qnbit;
	celoff = (dp->last_bidx - dp->boff) * aux->qnbin;
	// update cells' mask to be calculated by new coming kbm_cmer_t
	for(i=0;i<aux->qnbit;i+=64){
		masks[1]->bits[i >> 6] = masks[0]->bits[i >> 6] | dp->cmask->bits[(i + bitoff) >> 6];
	}
	// dp core
	H = KBM_CELL_NULL;
	i = next_one_bitvec(masks[1], 0);
	D = (i && get_bitvec(masks[0], i - 1))? cells[0]->buffer[i - 1] : KBM_CELL_NULL;
	reg_zeros_bitvec(masks[0], 0, i);
	encap_bit2vec(dp->bts, aux->qnbin);
	n = 0;
	while(i < aux->qnbin){
		if(get_bitvec(dp->cmask, i + bitoff)){
			is_gap = 0;
			m = ref_kbmcmerv(dp->cms, rowoff + n);
			score = m->kmat + m->kcnt;
			n ++;
		} else {
			is_gap = 1;
			m = (kbm_cmer_t*)&KBM_CMER_NULL;
			score = aux->par->pgap;
		}
		// horizontal
		{
			H.score += score;
			H.mat += m->kmat;
			//H.gap = 0;
		}
		if(H.var < 0){
			H.score += - aux->par->pvar; // score increases for abs(var) decreased
		} else {
			H.score += aux->par->pvar;
		}
		H.var ++;
		H.bt = 1; // horizontal backtrace, insertion for query sequence
		// diagonal
		{
			D.score += score;
			D.mat += m->kmat;
			//D.gap = 0;
		}
		D.var = 0; // whether to reset var
		D.bt = 0; // diagonal backtrace
		if(D.score >= H.score){
			H = D;
		}
		// vertical
		V = D = get_bitvec(masks[0], i)? cells[0]->buffer[i] : KBM_CELL_NULL;
		{
			V.score += score;
			V.mat += m->kmat;
			//V.gap = 0;
		}
		if(V.var > 0){
			V.score += - aux->par->pvar; // score increases for abs(var) decreased
		} else {
			V.score += aux->par->pvar;
		}
		V.var --;
		V.bt = 2; // vertical backtrace
		if(V.score > H.score){
			H = V;
		}
		if(is_gap){
			H.gap ++;
		} else {
			H.gap = 0;
		}
		set_bit2vec(dp->bts, dp->bts->size + i, H.bt);
		// init new path ID when there is no progenitor
		if(H.beg == 0){
			H.beg = 1 + i + celoff;
		}
#if __DEBUG__
		if(KBM_LOG >= KBM_LOG_ALL){
			fprintf(KBM_LOGF, "KBMLOG%d [x=%d, y=%llu, beg=%llu, score=%d, gap=%d, var=%d, bt=%d]\n", __LINE__, i, dp->last_bidx, (u8i)H.beg, H.score, H.gap, H.var, H.bt);
		}
#endif
		if(H.score > 0 && H.gap <= aux->par->max_bgap && num_abs(H.var) <= aux->par->max_bvar){
			// set cell, move next
			one_bitvec(masks[1], i);
			cells[1]->buffer[i] = H;
			if(is_gap == 0 && H.mat >= aux->par->min_mat) _update_dp_path_kbm(dp, 1 + i + celoff, &H);
			i ++;
		} else if(D.score > 0){
			// move next, may have diagonal hit
			zero_bitvec(masks[1], i);
			H = KBM_CELL_NULL;
			i ++;
		} else {
			zero_bitvec(masks[1], i);
			ni = next_one_bitvec(masks[1], i + 1);
			if(i + 1 < ni){
				reg_zeros_bitvec(masks[1], i + 1, ni);
				D = H = KBM_CELL_NULL;
			}
			i = ni;
		}
	}
	// prepare next row
	//NEXT_ROW:
	push_u4v(dp->coffs, dp->cms->size);
	zero_bitvec(dp->cmask, dp->cmask->n_bit);
	dp->cmask->n_bit += aux->qnbit;
	one_bitvec(dp->cmask, dp->cmask->n_bit);
	dp->bts->size += aux->qnbin;
	dp->last_bidx ++;
}

static inline int _backtrace_map_kbm(KBMAux *aux, int dir, kbm_path_t *p){
	KBMDP *dp;
	kbm_map_t *hit;
	kbm_cmer_t *c;
	u8i cgoff, sidx;
	u4i i, mat, cnt, gap, cglen;
	int tmp, x, y, bt, lst;
	dp = aux->dps[dir];
	hit = next_ref_kbmmapv(aux->hits);
	hit->qidx = aux->qidx;
	hit->qdir = dir;
	hit->tidx = aux->kbm->bins->buffer[dp->boff + p->beg / aux->qnbin].ridx;
	hit->tdir = 0;
	hit->qb = p->beg % aux->qnbin;
	hit->qe = p->end % aux->qnbin;
	hit->tb = p->beg / aux->qnbin;
	hit->te = p->end / aux->qnbin;
	hit->aln = num_min(hit->qe - hit->qb, hit->te - hit->tb);
	hit->aln ++;
	if(hit->aln < aux->par->min_aln){
		aux->hits->size --;
		return 0;
	}
	cgoff = aux->cigars->size;
	cglen = 0;
	cnt = 0;
	mat = 0;
	gap = 0;
	x = hit->qe;
	y = hit->te;
	lst = 0;
	while(x >= hit->qb && y >= hit->tb){
		bt = get_bit2vec(dp->bts, x + y * aux->qnbin);
		if(get_bitvec(dp->cmask, x + y * aux->qnbit)){
			c = ref_kbmcmerv(dp->cms, rank_bitvec(dp->cmask, x + y * aux->qnbit));
			cnt += c->kcnt;
			mat += c->kmat;
			// try merge 'ID' or 'DI' into 'M'
			if(lst == 0){
				if(bt == 0){
					push_bitsvec(aux->cigars, 0);
				}
				lst = bt;
			} else if(bt == 0){
				push_bitsvec(aux->cigars, lst);
				push_bitsvec(aux->cigars, 0);
				lst = 0;
			} else if(bt == lst){
				push_bitsvec(aux->cigars, lst);
			} else {
				push_bitsvec(aux->cigars, 0);
				lst = 0;
			}
		} else {
			gap ++;
			if(lst){
				push_bitsvec(aux->cigars, lst);
			}
			lst = 0;
			push_bitsvec(aux->cigars, 0x4 | bt);
		}
		switch(bt){
			case 0: x --; y --; break;
			case 1: x --; break;
			default: y --; break;
		}
	}
	if(lst){
		push_bitsvec(aux->cigars, lst);
	}
	cglen = aux->cigars->size - cgoff;
	if(mat < (u4i)aux->par->min_mat || mat < UInt(hit->aln * KBM_BSIZE * aux->par->min_sim)
		|| gap > (u4i)(hit->aln * aux->par->max_gap)
		|| hit->aln < (int)aux->par->min_aln
		|| num_diff(hit->qe - hit->qb, hit->te - hit->tb) > (int)num_max(aux->par->aln_var * hit->aln, 1.0)){
		aux->hits->size --;
		aux->cigars->size = cgoff;
		return 0;
	}
	if(aux->par->self_aln && aux->solids){
		// Obsolete
		x = hit->qe;
		y = hit->te;
		while(x >= hit->qb && y >= hit->tb){
			bt = get_bit2vec(dp->bts, x + y * aux->qnbin);
			if(get_bitvec(dp->cmask, x + y * aux->qnbit)){
				c = ref_kbmcmerv(dp->cms, rank_bitvec(dp->cmask, x + y * aux->qnbit));
				for(i=0;i<c->kcnt;i++){
					sidx = seed2solid_idx_kbm(aux->kbm, dp->kms->buffer + c->koff + i);
					one_bitvec(aux->solids, sidx); // Thread-unsafe, but no hurt
				}
			} else {
			}
			switch(bt){
				case 0: x --; y --; break;
				case 1: x --; break;
				default: y --; break;
			}
		}
	}
	hit->qe = (hit->qe + 1);
	hit->tb = (dp->boff + hit->tb - aux->kbm->reads->buffer[hit->tidx].binoff);
	hit->te = (dp->boff + hit->te - aux->kbm->reads->buffer[hit->tidx].binoff + 1);
	hit->mat = mat;
	hit->cnt = cnt;
	hit->gap = gap;
	hit->cgoff = cgoff;
	hit->cglen = cglen;
	if(dir){
		tmp = aux->qnbin - hit->qb;
		hit->qb = aux->qnbin - hit->qe;
		hit->qe = tmp;
	}
#if __DEBUG__
	if(KBM_LOG){
		fprintf(KBM_LOGF, "HIT\tQ[%d]\t%c\t%d\t%d", hit->qidx, "+-"[hit->qdir], hit->qb, hit->qe);
		fprintf(KBM_LOGF, "\tT[%d]\t%c\t%d\t%d", hit->tidx, "+-"[hit->tdir], hit->tb, hit->te);
		fprintf(KBM_LOGF, "\t%d\t%d\t%d\t%d\n", hit->mat, hit->aln, hit->cnt, hit->gap);
	}
#endif
	return 1;
}

static inline int check_hit_cigar_kbm(kbm_map_t *hit, BitsVec *cigars){
	u4i i, bt;
	int x, y;
	x = y = -1;
	i = hit->cglen;
	while(i){
		bt = get_bitsvec(cigars, hit->cgoff + i - 1);
		bt = bt & 0x03;
		x += (0b0011 >> bt) & 0x01;
		y += (0b0101 >> bt) & 0x01;
		i --;
	}
	return !(x + 1 + hit->qb == hit->qe && y + 1 + hit->tb == hit->te);
}

static inline void print_hit_kbm(KBM *kbm, char *qtag, u4i qlen, kbm_map_t *hit, BitsVec *cigars, String *_str, FILE *out){
	String *str;
	u8i coff;
	u4i clen, len, bt, _bt;
	if(hit->mat == 0) return;
	fprintf(out, "%s\t%c\t%d\t%d\t%d", qtag, "+-"[hit->qdir], qlen, hit->qb * KBM_BIN_SIZE, hit->qe * KBM_BIN_SIZE);
	fprintf(out, "\t%s\t%c\t%d\t%d\t%d", kbm->reads->buffer[hit->tidx].tag, "+-"[hit->tdir], kbm->reads->buffer[hit->tidx].bincnt * KBM_BIN_SIZE, hit->tb * KBM_BIN_SIZE, hit->te * KBM_BIN_SIZE);
	fprintf(out, "\t%d\t%d\t%d\t%d\t", hit->mat, hit->aln * KBM_BIN_SIZE, hit->cnt, hit->gap);
	if(cigars){
		str = _str? _str : init_string(64);
		bt = len = 0;
		coff = hit->cgoff;
		clen = hit->cglen;
		clear_string(str);
		while(clen){
			_bt = get_bitsvec(cigars, coff + clen - 1);
			if(_bt == bt){
				len ++;
			} else {
				if(len > 1){
					add_int_string(str, len);
					add_char_string(str, "MID?mid?"[bt]);
				} else if(len == 1){
					add_char_string(str, "MID?mid?"[bt]);
				}
				bt = _bt;
				len = 1;
			}
			clen --;
		}
		if(len > 1){
			add_int_string(str, len);
			add_char_string(str, "MID?mid?"[bt]);
		} else if(len == 1){
			add_char_string(str, "MID?mid?"[bt]);
		}
		fputs(str->string, out);
		if(_str == NULL) free_string(str);
	} else {
		fputc('*', out);
	}
	fprintf(out, "\n");
}

static inline void fprint_hit_kbm(KBMAux *aux, u4i hidx, FILE *out){
	kbm_map_t *hit;
	hit = ref_kbmmapv(aux->hits, hidx);
	print_hit_kbm(aux->kbm, aux->qtag, aux->slen, hit, aux->cigars, aux->str, out);
}

static inline void flip_hit_kbmaux(KBMAux *dst, KBMAux *src, u4i hidx){
	kbm_map_t *h1, *h2;
	u4i t, i;
	h2 = next_ref_kbmmapv(dst->hits);
	h1 = ref_kbmmapv(src->hits, hidx);
	//h2->qidx = h1->tidx;
	h2->qidx = dst->qidx;
	h2->qdir = h1->qdir;
	h2->tidx = h1->qidx;
	h2->tdir = h1->tdir;
	h2->cgoff = dst->cigars->size;
	h2->cglen = h1->cglen;
	h2->qb = h1->tb;
	h2->qe = h1->te;
	h2->tb = h1->qb;
	h2->te = h1->qe;
	h2->mat = h1->mat;
	h2->cnt = h1->cnt;
	h2->aln = h1->aln;
	h2->gap = h1->gap;
	if(h1->qdir){
		for(i=0;i<h1->cglen;i++){
			t = get_bitsvec(src->cigars, h1->cgoff + h1->cglen - i - 1);
			if((t & 0x03)){
				t = ((~(t & 0x03)) & 0x03) | (t & 0x4);
			}
			push_bitsvec(dst->cigars, t);
		}
	} else {
		for(i=0;i<h1->cglen;i++){
			t = get_bitsvec(src->cigars, h1->cgoff + i);
			if((t & 0x03)){
				t = ((~(t & 0x03)) & 0x03) | (t & 0x4);
			}
			push_bitsvec(dst->cigars, t);
		}
	}
}

static inline int _dp_path2map_kbm(KBMAux *aux, int dir){
	KBMDP *dp;
	kbm_path_t *p;
	u4i ret;
	dp = aux->dps[dir];
	ret = 0;
	index_bitvec_core(dp->cmask, roundup_times(dp->cmask->n_bit, 64 * 8));
	reset_iter_kbmphash(dp->paths);
	while((p = ref_iter_kbmphash(dp->paths))){
		p->beg --; // 1-based to 0-based
		p->end --;
#if __DEBUG__
		if(KBM_LOG >= KBM_LOG_HIG){
			fprintf(KBM_LOGF, "KBMLOG%d\t%d\t%c\tkbm_path_t[%llu(%d:%d),%llu(%d:%d),%d]\n", __LINE__, aux->qidx, "+-"[dir], (u8i)p->beg, (u4i)(p->beg % aux->qnbin), (u4i)((p->beg / aux->qnbin) + dp->boff), (u8i)p->end, (u4i)(p->end % aux->qnbin), (u4i)((p->end / aux->qnbin) + dp->boff), p->score);
		}
#endif
		if(_backtrace_map_kbm(aux, dir, p)){
			ret ++;
		}
	}
	//clear_kbmphash(dp->paths);
	return ret;
}

static inline void push_kmer_match_kbm(KBMAux *aux, int dir, kbm_dpe_t *p){
	KBM *kbm;
	KBMDP *dp;
	kbm_cmer_t *c;
	kbm_dpe_t *e, E;
	u8i bb;
	u4i i, qb, kmat, kcnt, blen;
	kbm = aux->kbm;
	dp = aux->dps[dir];
	if(p == NULL){
		if(dp->kms->size == 0){
			dp->km_len = 0;
			return;
		}
		e = ref_kbmdpev(dp->kms, dp->kms->size - 1);
		if((int)dp->km_len < aux->par->min_mat || e->bidx + 1 < dp->kms->buffer[0].bidx + aux->par->min_aln){
			clear_kbmdpev(dp->kms);
			dp->km_len = aux->par->ksize + aux->par->psize;
			return;
		}
	} else {
		if(dp->kms->size == 0){
			dp->km_len = (aux->par->ksize + aux->par->psize);
			push_kbmdpev(dp->kms, *p);
			return;
		}
		e = ref_kbmdpev(dp->kms, dp->kms->size - 1);
		if(e->bidx == p->bidx){
			if(p->poff <= e->poff + aux->par->ksize + aux->par->psize){
				dp->km_len += p->poff - e->poff;
			} else {
				dp->km_len += aux->par->ksize + aux->par->psize;
			}
			push_kbmdpev(dp->kms, *p);
			return;
#if __DEBUG__
		} else if(p->bidx < e->bidx){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
#endif
		}
		if(e->bidx + aux->par->max_bgap + 1 < p->bidx || dp->kms->buffer[0].bidx + aux->par->max_bcnt < p->bidx || aux->kbm->bins->buffer[e->bidx].ridx != aux->kbm->bins->buffer[p->bidx].ridx){
			if(Int(dp->km_len) < aux->par->min_mat || e->bidx + 1 < dp->kms->buffer[0].bidx + aux->par->min_aln){
				clear_kbmdpev(dp->kms);
				push_kbmdpev(dp->kms, *p);
				dp->km_len = aux->par->ksize + aux->par->psize;
				return;
			}
		} else {
			dp->km_len += aux->par->ksize + aux->par->psize;
			push_kbmdpev(dp->kms, *p);
			return;
		}
	}
	reset_kbmdp(dp, aux, dp->kms->buffer[0].bidx);
#ifdef TEST_MODE
	if(aux->par->test_mode >= 1){
		clear_kbmdpev(dp->kms);
		if(p){
			dp->km_len = aux->par->ksize + aux->par->psize;
			push_kbmdpev(dp->kms, *p);
		}
		return;
	}
#endif
	blen = dp->kms->buffer[dp->kms->size - 1].bidx + 1 - dp->kms->buffer[0].bidx + 2;
	{
		push_u4v(dp->coffs, 0);
		clear_bitvec(dp->cmask);
		encap_bitvec(dp->cmask, aux->qnbit * blen);
		reg_zeros_bitvec(dp->cmask, 0, aux->qnbit * blen);
		dp->cmask->n_bit = aux->qnbit;
		one_bitvec(dp->cmask, dp->cmask->n_bit);
	}
	qb = dp->kms->buffer[0].poff;
	bb = dp->kms->buffer[0].bidx;
	kcnt = 0;
	kmat = aux->par->ksize + aux->par->psize;
	E.bidx = dp->kms->buffer[dp->kms->size - 1].bidx + 1;
	E.poff = 0;
	for(i=0;i<=dp->kms->size;i++){
		e = (i < dp->kms->size)? ref_kbmdpev(dp->kms, i) : &E;
		if(e->bidx == bb && e->poff / KBM_BIN_SIZE == qb / KBM_BIN_SIZE){
			kcnt ++;
			if(qb + aux->par->ksize + aux->par->psize >= e->poff){
				kmat += e->poff - qb;
			} else {
				kmat += aux->par->ksize + aux->par->psize;
			}
		} else {
			one_bitvec(dp->cmask, (bb - dp->boff) * aux->qnbit + qb / KBM_BIN_SIZE);
			c = next_ref_kbmcmerv(dp->cms);
			c->koff = i - kcnt;
			c->kcnt = kcnt;
			c->kmat = kmat;
			c->boff = bb - dp->boff;
#if __DEBUG__
			if(KBM_LOG >= KBM_LOG_ALL){
				fprintf(KBM_LOGF, "KBMLOG%d [x=%d, y=%d(%s:%d), off=%d cnt=%d, mat=%d]\n", __LINE__, qb / KBM_BIN_SIZE, bb,
					aux->kbm->reads->buffer[aux->kbm->bins->buffer[bb].ridx].tag, bb - aux->kbm->reads->buffer[aux->kbm->bins->buffer[bb].ridx].binoff, c->koff, c->kcnt, c->kmat);
			}
#endif
			while(bb < e->bidx){
				_dp_cal_spare_row_kbm(aux, dir);
				bb ++;
			}
			kcnt = 1;
			kmat = aux->par->ksize + aux->par->psize;
		}
		qb = e->poff;
	}
	// flush last row
	_dp_cal_spare_row_kbm(aux, dir);
	//collecting maps
	_dp_path2map_kbm(aux, dir);
	clear_kbmdpev(dp->kms);
	if(p) push_kbmdpev(dp->kms, *p);
	dp->km_len = aux->par->ksize + aux->par->psize;
}

static inline void map_kbm(KBMAux *aux){
	KBM *kbm;
	kbm_ref_t *ref;
	kbm_baux_t *saux;
	u4v *heap;
	u4i idx, hidx, i, j, pdir;
	kbm = aux->kbm;
#ifdef TEST_MODE
	if(aux->par->test_mode >= 4) return;
#endif
	while(aux->hptr < aux->bmlen){
		if(aux->hptr - aux->bmoff >= aux->nheap){
			aux->bmoff += aux->nheap;
			for(i=0;i<aux->nheap;i++){
				clear_u4v(aux->heaps[i]);
			}
			for(i=0;i<aux->refs->size;i++){
				ref = ref_kbmrefv(aux->refs, i);
				while(ref->boff < ref->bend){
					hidx = ref->bidx / aux->bmcnt;
					if(hidx - aux->bmoff < aux->nheap){
						push_u4v(aux->heaps[hidx - aux->bmoff], i);
					}
					break;
				}
			}
		}
		heap = aux->heaps[aux->hptr - aux->bmoff];
		if(heap->size){
			clear_kbmdpev(aux->caches[0]);
			clear_kbmdpev(aux->caches[1]);
			for(i=0;i<heap->size;i++){
				idx = heap->buffer[i];
				ref = ref_kbmrefv(aux->refs, idx);
				while(1){
					saux = ref_kbmbauxv(kbm->sauxs, ref->boff);
					pdir = (ref->dir ^ saux->dir);
					if(((aux->par->strand_mask >> pdir) & 0x01)){
						push_kbmdpev(aux->caches[pdir], (kbm_dpe_t){ref->poffs[pdir], idx, ref->bidx, saux->koff});
					}
					ref->boff ++;
					ref->bidx = getval_bidx(aux->kbm, ref->boff);
					if(ref->boff >= ref->bend) break;
#if __DEBUG__
						if(ref->bidx < getval_bidx(aux->kbm, ref->boff - 1)){
							fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
							abort();
						}
#endif
					hidx = ref->bidx / aux->bmcnt;
					if(hidx > aux->hptr){
						if(hidx - aux->bmoff < aux->nheap){
							push_u4v(aux->heaps[hidx - aux->bmoff], idx);
						}
						break;
					}
				}
			}
			{
#ifdef TEST_MODE
				if(aux->par->test_mode <= 2){
#endif
				if(aux->caches[0]->size * (aux->par->ksize + aux->par->psize) < UInt(aux->par->min_mat)){
					aux->caches[0]->size = 0;
				} else {
					sort_array(aux->caches[0]->buffer, aux->caches[0]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
				}
				if(aux->caches[1]->size * (aux->par->ksize + aux->par->psize) < UInt(aux->par->min_mat)){
					aux->caches[1]->size = 0;
				} else {
					sort_array(aux->caches[1]->buffer, aux->caches[1]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
				}
					//sort_array(aux->caches[0]->buffer, aux->caches[0]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
					//sort_array(aux->caches[1]->buffer, aux->caches[1]->size, kbm_dpe_t, num_cmpgtx(a.bidx, b.bidx, a.poff, b.poff));
					// TODO: sort by bidx+koff is more reasonable, need to modify push_kmer_match_kbm too
#ifdef TEST_MODE
				}
#endif
#ifdef TEST_MODE
				if(aux->par->test_mode <= 1){
#endif
				for(i=0;i<2;i++){
					for(j=0;j<aux->caches[i]->size;j++){
#if __DEBUG__
						if(KBM_LOG >= KBM_LOG_ALL){
							//fprintf(KBM_LOGF, "KBMLOG%d\t%d\t%d\t%c\t%d\t%llu[%d,%d]\t%llu[%d,%d]\n", __LINE__, aux->qidx, ref->poffs[ref->pdir], "+-"[ref->pdir], aux->hptr, ref->bidx, aux->kbm->bins->buffer[ref->bidx].ridx, aux->kbm->bins->buffer[ref->bidx].off * KBM_BIN_SIZE, (u8i)e->bidx, aux->kbm->bins->buffer[e->bidx].ridx, e->poff);
						}
#endif
						push_kmer_match_kbm(aux, i, aux->caches[i]->buffer + j);
					}
				}
				if(aux->hits->size >= aux->par->max_hit) return;
#ifdef TEST_MODE
				}
#endif
			}
		}
		aux->hptr ++;
	}
	if(aux->par->strand_mask & 0x01) push_kmer_match_kbm(aux, 0, NULL);
	if(aux->par->strand_mask & 0x02) push_kmer_match_kbm(aux, 1, NULL);
}

// KBM's tag2idx is wrongly loaded, need to be corrected
static inline void rebuild_tag2idx_kbm(void *_kbm, size_t aux){
	KBM *kbm;
	kbm_read_t *rd;
	u4i i;
	UNUSED(aux);
	kbm = (KBM*)_kbm;
	clear_cuhash(kbm->tag2idx); // hash size is not changed, thus there won't have hash re-size
	for(i=0;i<kbm->reads->size;i++){
		rd = ref_kbmreadv(kbm->reads, i);
		if(rd->tag) put_cuhash(kbm->tag2idx, (cuhash_t){rd->tag, i});
	}
	kbm->flags |= 1LLU << 0;
}

static inline int simple_chain_all_maps_kbm(kbm_map_t *srcs, u4i size, BitsVec *src_cigars, kbm_map_t *dst, BitsVec *dst_cigars, float max_aln_var){
	kbm_map_t *hit;
	u4i i, x, y, z, f;
	if(size < 2) return 0;
	sort_array(srcs, size, kbm_map_t, num_cmpgt(a.tb, b.tb));
	*dst = srcs[size - 1];
	dst->cgoff = dst_cigars->size;
	append_bitsvec(dst_cigars, src_cigars, srcs[size - 1].cgoff, srcs[size - 1].cglen);
	for(i=size-1;i>0;i--){
		hit = srcs + i - 1;
		if(dst->tb < hit->te){
			goto FAILED;
		} else {
			y = dst->tb - hit->te;
			dst->tb = hit->tb;
		}
		if(dst->qdir){
			if(hit->qe > dst->qb){
				goto FAILED;
			} else {
				x = dst->qb - hit->qe;
				dst->qb = hit->qb;
			}
		} else {
			if(dst->qb < hit->qe){
				goto FAILED;
			} else {
				x = dst->qb - hit->qe;
				dst->qb = hit->qb;
			}
		}
		dst->mat += hit->mat;
		dst->cnt += hit->cnt;
		dst->gap += hit->gap;
		dst->gap += num_max(x, y);
		z = num_min(x, y);
		f = 0x4 | 0; // diagonal GAP
		pushs_bitsvec(dst_cigars, f, z);
		x -= z;
		y -= z;
		if(x > y){
			z = x;
			f = 0x4 | 1;
		} else {
			z = y;
			f = 0x4 | 2;
		}
		pushs_bitsvec(dst_cigars, f, z);
		append_bitsvec(dst_cigars, src_cigars, hit->cgoff, hit->cglen);
	}
	dst->aln = num_min(dst->qe - dst->qb, dst->te - dst->tb);
	if(dst->aln * max_aln_var < num_diff(dst->qe - dst->qb, dst->te - dst->tb)){
		goto FAILED;
	}
	dst->cglen = dst_cigars->size - dst->cgoff;
	return 1;
	FAILED:
	dst_cigars->size = dst->cgoff;
	return 0;
}

#endif
