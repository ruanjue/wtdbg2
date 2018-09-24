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
 
#ifndef __BITS_VEC_RJ_H
#define __BITS_VEC_RJ_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "mem_share.h"

/* Useful functions when n_bit > 8 */

static inline void u8byte2bits(uint64_t val, uint8_t *dat, uint64_t offset, uint8_t size){
	uint8_t i, *v;
	v = (uint8_t*)&val;
#if __BYTE_ORDER == 1234
	for(i=0;i<size;i++){
		if(v[i >> 3] & (1U << (i & 0x7U))) dat[offset >> 3] |= 1U << (offset & 0x7U);
		else dat[offset >> 3] &= ~(1U << (offset & 0x7U));
		offset ++;
	}
#else
	for(i=0;i<size;i++){
		if(v[7 - (i >> 3)] & (1U << (i & 0x7U))) dat[offset >> 3] |= 1U << (offset & 0x7U);
		else dat[offset >> 3] &= ~(1U << (offset & 0x7U));
		offset ++;
	}
#endif
}

static inline uint64_t bits2u8byte(uint8_t *dat, uint64_t offset, uint8_t size){
	uint64_t ret;
	uint8_t i, *v;
	ret = 0;
	v = (uint8_t*)&ret;
#if __BYTE_ORDER == 1234
	for(i=0;i<size;i++){
		if(dat[offset >> 3] & (1U << (offset & 0x7U))) v[i >> 3] |= 1U << (i & 0x7U);
		offset ++;
	}
#else
	for(i=0;i<size;i++){
		if(dat[offset >> 3] & (1U << (offset & 0x7U))) v[7 - (i >> 3)] |= 1U << (i & 0x7U);
		offset ++;
	}
#endif
	return ret;
}

typedef struct {
	uint8_t *bits;
	uint64_t size;
	uint64_t cap;
	uint8_t n_bit;
	uint32_t mask;
} BitsVec;

static inline BitsVec* init_bitsvec(uint64_t size, uint32_t n_bit){
	BitsVec *vec;
	if(n_bit == 0) n_bit = 1;
	else if(n_bit > 8) n_bit = 8;
	if(size < 8) size = 8;
	vec = calloc(1, sizeof(BitsVec));
	vec->n_bit = n_bit;
	vec->mask  = (1U << n_bit) - 1U;
	vec->size  = 0;
	vec->cap   = size;
	vec->bits  = calloc((size * vec->n_bit + 15) / 8, 1);
	return vec;
}

static inline size_t bitsvec_obj_desc_cnt(void *obj, int idx){
	return (((BitsVec*)obj)->cap * ((BitsVec*)obj)->n_bit + 15) / 8;
	idx = idx;
}

static const obj_desc_t bitsvec_obj_desc = {"bitsvec_obj_desc", sizeof(BitsVec), 1, {1}, {offsetof(BitsVec, bits)}, {(obj_desc_t*)&OBJ_DESC_DATA}, bitsvec_obj_desc_cnt, NULL};

static inline void clear_bitsvec(BitsVec *vec){
	vec->size = 0;
}

static inline void free_bitsvec(BitsVec *vec){
	free(vec->bits);
	free(vec);
}

static inline int encap_bitsvec(BitsVec *vec, u8i n){
	if(vec->size + n <= vec->cap) return 0;
	if(vec->size + n < 0x3FFFFFFFLLU){
		vec->cap = roundup_power2(vec->size + n);
	} else {
		vec->cap = (vec->size + n + 0x3FFFFFFFLLU) & (MAX_U8 << 30);
	}
	vec->bits = realloc(vec->bits, (vec->cap * vec->n_bit + 15) / 8);
	return 1;
}

static inline void set_bitsvec(BitsVec *vec, u8i idx, u1i dat){
	register u8i off;
	register u2i x, d;
	off = (idx * vec->n_bit);
	d = off & 0x07;
	off >>= 3;
	x = (((u2i)vec->bits[off + 1]) << 8) | vec->bits[off + 0];
	x = (x & (~(vec->mask << d))) | ((UInt(dat) & vec->mask) << d);
	vec->bits[off] = x;
	vec->bits[off + 1] = x >> 8;
}

static inline void push_bitsvec(BitsVec *vec, u1i dat){
	encap_bitsvec(vec, 1);
	set_bitsvec(vec, vec->size, dat);
	vec->size ++;
}

// TODO: need to be optimized
static inline void pushs_bitsvec(BitsVec *vec, u1i dat, u8i len){
	u8i i;
	encap_bitsvec(vec, len);
	for(i=0;i<len;i++){
		set_bitsvec(vec, vec->size + i, dat);
	}
	vec->size += len;
}

static inline u2i gets_bitsvec(BitsVec *vec, u8i idx){
	u8i off;
	off = (idx * vec->n_bit);
	return ((((u2i)vec->bits[(off >> 3) + 1]) << 8) | vec->bits[(off >> 3) + 0]) >> (off & 0x07);
}

static inline uint8_t get_bitsvec(BitsVec *vec, u8i idx){
	return gets_bitsvec(vec, idx) & vec->mask;
}

static inline void append_bitsvec(BitsVec *dst, BitsVec *src, u8i off, u8i len){
	u8i i, di, si, se;
	u2i x;
	u1i p, n, sd;
	if(0){ // Assume dst->n_bit == src->n_bit
		if(dst->n_bit != src->n_bit){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
	}
	encap_bitsvec(dst, len);
	p = (dst->n_bit & 0x01)? 8 : (8 / dst->n_bit);
	n = dst->size % p;
	if(n) n = p - n;
	if(len <= n){
		for(i=0;i<len;i++){
			set_bitsvec(dst, dst->size + i, get_bitsvec(src, off + i));
		}
		dst->size += len;
		return;
	} else {
		for(i=0;i<n;i++){
			set_bitsvec(dst, dst->size + i, get_bitsvec(src, off + i));
		}
		dst->size += n;
	}
	di = (dst->size * dst->n_bit) >> 3;
	si = ((off + i) * src->n_bit);
	sd = si & 0x07;
	si >>= 3;
	se = ((off + len) * src->n_bit + 7) >> 3;
	while(si < se){
		x = ((src->bits[si + 1] << 8) | src->bits[si]) >> sd;
		dst->bits[di++] = x;
		si ++;
	}
	dst->size += len - i;
}

static inline int pop_bitsvec(BitsVec *vec){
	if(vec->size == 0) return -1;
	vec->size --;
	return get_bitsvec(vec, vec->size);
}

#endif
