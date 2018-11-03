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
 
#ifndef __BIT_VEC_RJ_H
#define __BIT_VEC_RJ_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "mem_share.h"

static const u1i byte_ones_table[256] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

static inline unsigned int _bitvec_roundup_power2(unsigned int v){
	if(v == 0) return 0;
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v + 1;
}

typedef struct {
	u8i *bits;
	u8i n_bit;
	u8i n_cap;

	u8i *sums;
	u8i sum_size;
	u8i n_ones;
	u8i *hash;
	u8i hash_size;
	u8i hash_mod;
	int64_t iter_idx;
} BitVec;

#if 0

static inline u4i count_ones_bit32(u4i v){
	v = v - ((v >> 1) & 0x55555555U);                        // reuse input as temporary
	v = (v & 0x33333333U) + ((v >> 2) & 0x33333333U);        // temp
	return (((v + (v >> 4)) & 0xF0F0F0FU) * 0x1010101U) >> 24; // count
}

#define ONES_STEP_4 0x1111111111111111ULL
#define ONES_STEP_8 0x0101010101010101ULL

static inline int count_ones_bit64(const u8i x){
	register u8i byte_sums = x - ((x & 0xa * ONES_STEP_4) >> 1);
	byte_sums = (byte_sums & 3 * ONES_STEP_4) + ((byte_sums >> 2) & 3 * ONES_STEP_4);
	byte_sums = (byte_sums + (byte_sums >> 4)) & 0x0f * ONES_STEP_8;
	return byte_sums * ONES_STEP_8 >> 56;
}

#else

#define count_ones_bit32(v) __builtin_popcount(v)
#define count_ones_bit64(v) __builtin_popcountll(v)

#endif

#define reverse_u1i(v) (((((u1i)v) * 0x0202020202ULL) & 0x010884422010ULL) % 1023)

static inline size_t bitvec_obj_desc_cnt(void *bitv, int idx){
	switch(idx){
		case 0: return ((BitVec*)bitv)->n_cap / 64 * 8;
		case 1: return ((BitVec*)bitv)->sums? (((BitVec*)bitv)->sum_size * 2 + 1) * 8 : 0;
		case 2: return ((BitVec*)bitv)->hash? (((BitVec*)bitv)->hash_size) * 8 : 0;
		default: return 0;
	}
}

static const obj_desc_t bitvec_obj_desc = {"bitvec_obj_desc", sizeof(BitVec), 3, {1, 1, 1}, {offsetof(BitVec, bits), offsetof(BitVec, sums), offsetof(BitVec, hash)}, {(obj_desc_t*)&OBJ_DESC_DATA, (obj_desc_t*)&OBJ_DESC_DATA, (obj_desc_t*)&OBJ_DESC_DATA}, bitvec_obj_desc_cnt, NULL};

static inline BitVec* init_bitvec(u8i n_bit){
	BitVec *bitv;
	if(n_bit == 0) n_bit = 64 * 8;
	bitv = (BitVec*)malloc(sizeof(BitVec));
	bitv->n_bit = 0;
	bitv->n_cap = (((n_bit + 63) / 64) + 7) / 8 * 64 * 8;
	bitv->bits  = (u8i*)calloc((bitv->n_cap / 64) + 1, 8);
	bitv->bits[bitv->n_cap / 64] = 0x0000000000000001LLU;
	//memset(bitv->bits, 0, bitv->n_cap / 8);
	bitv->sums = NULL;
	bitv->hash = NULL;
	bitv->sum_size = 0;
	bitv->n_ones = 0;
	bitv->hash_size = 0;
	bitv->hash_mod = 0;
	bitv->iter_idx = 0;
	return bitv;
}

static inline size_t dump_bitvec(BitVec *bitv, FILE *out){
	fwrite(&bitv->n_bit, sizeof(u8i), 1, out);
	fwrite(&bitv->n_cap, sizeof(u8i), 1, out);
	fwrite(bitv->bits, sizeof(u8i), bitv->n_cap / 64, out);
	return sizeof(u8i) * (2 + bitv->n_cap / 64);
}

static inline BitVec* load_bitvec(FILE *inp){
	BitVec *bitv;
	size_t n;
	bitv = (BitVec*)malloc(sizeof(BitVec));
	if((n = fread(&bitv->n_bit, sizeof(u8i), 1, inp)) != 1){
		free(bitv); return NULL;
	}
	if((n = fread(&bitv->n_cap, sizeof(u8i), 1, inp)) != 1){
		free(bitv); return NULL;
	}
	bitv->bits = (u8i*)malloc(bitv->n_cap / 8);
	if(bitv->bits == NULL){
		fprintf(stderr, " Out of memeory in load_bitvec\n "); fflush(stderr); exit(1);
	}
	if((n = fread(bitv->bits, sizeof(u8i), bitv->n_cap / 64, inp)) != bitv->n_cap / 64){
		free(bitv); free(bitv->bits); return NULL;
	}
	bitv->sums = NULL;
	bitv->hash = NULL;
	bitv->hash_size = 0;
	return bitv;
}

#if 0
static inline BitVec* mem_load_bitvec(void *mem, FILE *inp){
	BitVec *bitv;
	size_t off, n;
	bitv = mem;
	off = ((sizeof(BitVec) + 7) / 8) * 8;
	if((n = fread(&bitv->n_bit, sizeof(u8i), 1, inp)) != 1) return NULL;
	if((n = fread(&bitv->n_cap, sizeof(u8i), 1, inp)) != 1) return NULL;
	bitv->sums = NULL;
	bitv->hash = NULL;
	bitv->hash_size = 0;
	bitv->bits = mem + off;
	off += (bitv->n_cap / 64) * 8;
	if((n = fread(bitv->bits, sizeof(u8i), bitv->n_cap / 64, inp)) != bitv->n_cap / 64) return NULL;
	return bitv;
}
#endif

static inline void clear_bitvec(BitVec *bitv){ bitv->n_bit = 0; }

static inline void zeros_bitvec(BitVec *bitv){ memset(bitv->bits, 0, bitv->n_cap / 8); }

// exclusive end
static inline void reg_zeros_bitvec(BitVec *bitv, u8i beg, u8i end){
	u8i b, e;
	if(beg >= end) return;
	b = beg >> 6;
	e = end >> 6;
	if(b == e){
		bitv->bits[b] &= (MAX_U8 << (beg & 0x3FU)) ^ (MAX_U8 >> (64 - (end & 0x3FU)));
	} else {
		bitv->bits[b] &= ~(MAX_U8 << (beg & 0x3FU));
		while(++b < e){ bitv->bits[b] = 0; }
		bitv->bits[b] &= MAX_U8 << (end & 0x3FU);
	}
}

static inline void ones_bitvec(BitVec *bitv){ memset(bitv->bits, 0xFFU, bitv->n_cap / 8); }

// exclusive end
static inline void reg_ones_bitvec(BitVec *bitv, u8i beg, u8i end){
	u8i b, e;
	if(beg >= end) return;
	b = beg >> 6;
	e = end >> 6;
	if(b == e){
		bitv->bits[b] |= (MAX_U8 << (beg & 0x3FU)) & (MAX_U8 >> (64 - (end & 0x3FU)));
	} else {
		bitv->bits[b] |= MAX_U8 << (beg & 0x3FU);
		while(++b < e){ bitv->bits[b] = MAX_U8; }
		bitv->bits[b] |= ~(MAX_U8 << (end & 0x3FU));
	}
}

static inline void flip_bitvec(BitVec *bitv, u8i idx){ bitv->bits[idx>>6] ^= 1LLU << (idx&0x3FU); }

static inline void one_bitvec(BitVec *bitv, u8i idx){ bitv->bits[idx>>6] |= 1LLU << (idx&0x3FU); }

static inline void zero_bitvec(BitVec *bitv, u8i idx){ bitv->bits[idx>>6] &= ~(1LLU << (idx&0x3FU)); }

static inline void set_bitvec(BitVec *bitv, u8i idx, int v){
	if(v){
		one_bitvec(bitv, idx);
	} else {
		zero_bitvec(bitv, idx);
	}
}

static inline u8i get_bitvec(BitVec *bitv, u8i idx){ return (bitv->bits[idx>>6] >> (idx&0x3FU)) & 0x01LLU; }

static inline u8i get64_bitvec(BitVec *bitv, u8i off){
	u8i m, n;
	m = off >> 6;
	n = off & 0x3F;
	if(n){
		return (bitv->bits[m] >> (64 - n)) | (bitv->bits[m + 1] << n);
	} else {
		return bitv->bits[m];
	}
}

static inline void set64_bitvec(BitVec *bitv, u8i off, u8i val){
	u8i m, n;
	m = off >> 6;
	n = off & 0x3F;
	if(n){
		bitv->bits[m] = ((bitv->bits[m] << (64 - n)) >> (64 - n)) | (val << (64 - n));
		m ++;
		bitv->bits[m] = ((bitv->bits[m] >> n) << n) | (val >> (64 - n));
	} else {
		bitv->bits[m] = val;
	}
}

static inline void encap_bitvec(BitVec *bitv, u8i num){
	u8i cap;
	if(bitv->n_bit + num < bitv->n_cap) return;
	cap = bitv->n_cap;
	while(bitv->n_bit + num >= bitv->n_cap){
		if(bitv->n_cap < 1024 * 1024 * 8){
			bitv->n_cap <<= 1;
		} else bitv->n_cap += 1024 * 1024 * 8;
	}
	bitv->bits = (u8i*)realloc(bitv->bits, bitv->n_cap / 8 + 8);
	memset(((void*)bitv->bits) + cap / 8, 0, (bitv->n_cap - cap) / 8 + 8);
	bitv->bits[cap / 64] = 0x0000000000000001LLU;
}

static inline void recap_bitvec(BitVec *bitv, u8i new_cap){
	if(new_cap & 0x3FU) new_cap = (new_cap & 0xFFFFFFFFFFFFFFC0LLU) + 0x40U;
	if(bitv->n_cap == new_cap) return;
	bitv->bits = (u8i*)realloc(bitv->bits, new_cap / 8 + 8);
	if(new_cap > bitv->n_cap){
		memset(((void*)bitv->bits) + bitv->n_cap / 8, 0, (new_cap - bitv->n_cap) / 8 + 8);
	}
	bitv->bits[new_cap / 64] = 0x0000000000000001LLU;
	bitv->n_cap = new_cap;
}

static inline void one2bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); one_bitvec(bitv, bitv->n_bit); bitv->n_bit ++; }

static inline void zero2bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); zero_bitvec(bitv, bitv->n_bit); bitv->n_bit ++; }

static inline u8i get_2bitvec(BitVec *bitv, u8i idx){ return (bitv->bits[idx>>5] >> ((idx&0x1FU) << 1)) & 0x03LLU; }

static inline void set_2bitvec(BitVec *bitv, u8i idx, u8i v){
	bitv->bits[idx>>5] = (bitv->bits[idx>>5] & (~(0x03LLU << ((idx&0x1FU) << 1)))) | ((v&0x03LLU) << ((idx&0x1FU) << 1));
}

static inline void push_2bitvec(BitVec *bitv, u8i v){
	encap_bitvec(bitv, 2);
	set_2bitvec(bitv, bitv->n_bit >> 1, v);
	bitv->n_bit = ((bitv->n_bit >> 1) + 1) << 1;
}

static inline void end_bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); one_bitvec(bitv, bitv->n_bit); }

static inline u8i next_one_bitvec(BitVec *bitv, u8i idx){
	register u8i p, v;
	register u4i s;
	p = idx >> 6;
	s = idx & 0x3F;
	while(!(bitv->bits[p] >> s)){ p ++; s = 0; }
	v = bitv->bits[p] >> s;
	s += __builtin_ctzll(v);
	return (p << 6) + s;
}

static inline u8i reg_count_bitvec(BitVec *bitv, u8i beg, u8i end){
	u8i cnt, b, e, t;
	if(beg >= end) return 0;
	b = beg >> 6;
	e = end >> 6;
	if(b == e){
		t = (bitv->bits[b] & (MAX_U8 >> (64 - (end & 0x3F)))) >> (beg & 0x3F);
		cnt = count_ones_bit64(t);
	} else {
		cnt = count_ones_bit64(bitv->bits[b] >> (beg & 0x3F));
		while(++b < e){
			cnt += count_ones_bit64(bitv->bits[b]);
		}
		if(end & 0x3F){
			cnt += count_ones_bit64(bitv->bits[b] & (MAX_U8 >> (64 - (end & 0x3F))));
		}
	}
	return cnt;
}

static const int Mod37BitPosition[] = // map a bit value mod 37 to its position
{
 32,  0,  1, 26,  2, 23, 27,  0,  3, 16,
 24, 30, 28, 11,  0, 13,  4,  7, 17,  0,
 25, 22, 31, 15, 29, 10, 12,  6,  0, 21,
 14,  9,  5, 20,  8, 19, 18
};

static inline u8i next_one_bitvec2(BitVec *bitv, u8i idx){
	register u8i p;
	register u4i s, v;
	p = idx >> 6;
	s = idx & 0x3F;
	while(!(bitv->bits[p] >> s)){ p ++; s = 0; }
	if(!((bitv->bits[p] >> s) & 0xFFFFFFFFU)) s += 32;
	v = bitv->bits[p] >> s;
	s += Mod37BitPosition[(-v & v) % 37];
	return (p << 6) + s;
}

static inline u8i next_one_bitvec3(BitVec *bitv, u8i idx){
	register u8i p;
	register u4i s;
	p = idx >> 6;
	s = idx & 0x3F;
	while(!(bitv->bits[p] >> s)){ p ++; s = 0; }
	while(!((bitv->bits[p] >> s) & 0xFFU)) s += 8;
	while(!((bitv->bits[p] >> s) & 0x01U)) s ++;
	return (p << 6) + s;
}

//n_cap MUST be times of 64 * 8
static inline void index_bitvec_core(BitVec *bitv, size_t n_cap){
	u8i i, k, s, t, m;
	m = ((n_cap + 63) / 64 + 7) / 8;
	if(bitv->sums) free(bitv->sums);
	bitv->sums = (u8i*)calloc((m * 2 + 1), 8);
	t = 0;
	for(i=0;i<n_cap;i+=64*8){
		k = ((i>>6) >> 3) << 1;
		bitv->sums[k] = t;
		s = 0;
		s += count_ones_bit64(bitv->bits[(i>>6)+0]);
		bitv->sums[k+1] |= s << 0;
		s += count_ones_bit64(bitv->bits[(i>>6)+1]);
		bitv->sums[k+1] |= s << 9;
		s += count_ones_bit64(bitv->bits[(i>>6)+2]);
		bitv->sums[k+1] |= s << 18;
		s += count_ones_bit64(bitv->bits[(i>>6)+3]);
		bitv->sums[k+1] |= s << 27;
		s += count_ones_bit64(bitv->bits[(i>>6)+4]);
		bitv->sums[k+1] |= s << 36;
		s += count_ones_bit64(bitv->bits[(i>>6)+5]);
		bitv->sums[k+1] |= s << 45;
		s += count_ones_bit64(bitv->bits[(i>>6)+6]);
		bitv->sums[k+1] |= s << 54;
		s += count_ones_bit64(bitv->bits[(i>>6)+7]);
		t += s;
	}
	bitv->sums[((i>>6) >> 3) << 1] = t;
	bitv->n_ones = t;
	bitv->sum_size = m;
	bitv->hash_size = (n_cap / 64 / 8) / 2;
	if(bitv->hash_size == 0) bitv->hash_size = 1;
	bitv->hash_mod = (t + bitv->hash_size) / bitv->hash_size;
	if(bitv->hash_mod == 0) bitv->hash_mod = 1;
	if(bitv->hash) free(bitv->hash);
	bitv->hash = (u8i*)malloc(sizeof(u8i) * bitv->hash_size);
	s = 0;
	t = 0;
	for(i=0;i<=m;i++){
		k = bitv->sums[i*2] / bitv->hash_mod;
		if(s < k){
			while(s < k){ bitv->hash[s] = t; s ++; }
			t = i? i - 1 : 0;
		}
	}
	bitv->hash[bitv->sums[m*2] / bitv->hash_mod] = t;
}

static inline void index_bitvec(BitVec *bitv){
	index_bitvec_core(bitv, bitv->n_cap);
}

static inline u8i rank_bitvec(BitVec *bitv, u8i idx){
	u8i p, s, sum;
	p = (idx>>6)>>3;
	s = (idx >> 6) & 0x07U;
	sum = bitv->sums[p<<1];
	if(s) sum += (bitv->sums[(p<<1)+1] >> (9 * (s - 1))) & 0x1FFU;
	if(idx & 0x3FU) sum += count_ones_bit64(bitv->bits[idx>>6]<<(64-(idx&0x3FU)));
	return sum;
}

static inline u1i select_8bytes(u8i word, u1i n_one){
	u1i idx, n, m;
	n = count_ones_bit32((u4i)word);
	if(n >= n_one){
		n = 0;
		idx = 0;
		word = word & 0xFFFFFFFFU;
	} else {
		idx = 32;
		word = word >> 32;
	}
	while(1){
		m = byte_ones_table[(u1i)word];
		if(n + m >= n_one) break;
		n += m;
		idx += 8;
		word >>= 8;
	}
	m = byte_ones_table[(u1i)(word & 0xF)];
	if(n + m < n_one){
		idx += 4;
		word >>= 4;
		n += m;
	}
	while(word){
		idx ++;
		if(word & 0x01){
			n ++;
			if(n == n_one) break;
		}
		word >>= 1;
	}
	return idx;
}

/*
 * To select the 1'st one, use select_bitvec(bitv, 1) - 1
 * */
static inline u8i select_bitvec(BitVec *bitv, u8i idx){
	u8i i, p, s, sum, t;
	p = bitv->hash[idx / bitv->hash_mod];
	while(p + 1 < bitv->sum_size && bitv->sums[(p + 1) << 1] < idx) p ++;
	sum = bitv->sums[p << 1];
	i = 0;
	t = sum;
	while(i < 7){
		s = (bitv->sums[(p << 1) + 1] >> (9 * i)) & 0x1FFU;
		if(sum + s >= idx) break;
		t = sum + s;
		i ++;
	}
	p = p * 8 + i;
	s = idx - t;
	return p * 64 + select_8bytes(bitv->bits[p], s);
}

static inline void begin_iter_bitvec(BitVec *bitv){ bitv->iter_idx = -1; }

static inline u8i iter_bitvec(BitVec *bitv){
	if((u8i)(bitv->iter_idx + 1) > bitv->n_cap) return 0xFFFFFFFFFFFFFFFFLLU;
	bitv->iter_idx = next_one_bitvec(bitv, bitv->iter_idx + 1);
	return (u8i)bitv->iter_idx;
}

static inline void free_bitvec(BitVec *bitv){
	free(bitv->bits);
	if(bitv->sums) free(bitv->sums);
	if(bitv->hash) free(bitv->hash);
	free(bitv);
}

#if 0

static inline size_t mem_size_bitvec(BitVec *bitv){
	size_t m;
	m = (sizeof(BitVec) + 7) / 8 * 8 + ((bitv->n_cap / 64) * 8);
	if(bitv->sums){
		m += (bitv->sum_size * 2 + 1) * 8;
	}
	if(bitv->hash){
		m += bitv->hash_size * 8;
	}
	return m;
}

static inline size_t mem_dump_bitvec(BitVec *bitv, void *mem){
	BitVec *clone;
	size_t off;
	clone = mem;
	memcpy(clone, bitv, sizeof(BitVec));
	off = ((sizeof(BitVec) + 7) / 8) * 8;
	clone->bits = mem + off;
	memcpy(clone->bits, bitv->bits, (bitv->n_cap / 64) * 8);
	off += (bitv->n_cap / 64) * 8;
	if(bitv->sums){
		clone->sums = mem + off;
		memcpy(clone->sums, bitv->sums, (bitv->sum_size * 2 + 1) * 8);
		off += (bitv->sum_size * 2 + 1) * 8;
	}
	if(bitv->hash){
		clone->hash = mem + off;
		memcpy(clone->hash, bitv->hash, bitv->hash_size * 8);
		off += bitv->hash_size * 8;
	}
	return off;
}
#endif

#endif
