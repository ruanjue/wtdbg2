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
 
#ifndef __BIT2_VEC_RJ_H
#define __BIT2_VEC_RJ_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "mem_share.h"

typedef struct {
	uint8_t *bits;
	uint64_t size;
	uint64_t cap;
} Bit2Vec;

static inline Bit2Vec* init_bit2vec(uint64_t size){
	Bit2Vec *vec;
	vec = malloc(sizeof(Bit2Vec));
	vec->size  = 0;
	vec->cap   = size;
	vec->bits  = calloc(1, (size * 2 + 7) / 8);
	return vec;
}

static inline size_t bit2vec_obj_desc_cnt(void *obj, int idx){
	return (((Bit2Vec*)obj)->cap * 2 + 7) / 8;
	idx = idx;
}

static const obj_desc_t bit2vec_obj_desc = {"bit2vec_obj_desc", sizeof(Bit2Vec), 1, {1}, {offsetof(Bit2Vec, bits)}, {(obj_desc_t*)&OBJ_DESC_DATA}, bit2vec_obj_desc_cnt, NULL};

static inline void clear_bit2vec(Bit2Vec *vec){
	memset(vec->bits, 0, (vec->cap * 2 + 7) / 8);
	vec->size = 0;
}

static inline void free_bit2vec(Bit2Vec *vec){
	free(vec->bits);
	free(vec);
}

static inline int encap_bit2vec(Bit2Vec *vec, uint32_t n){
	uint64_t cap;
	if(vec->size + n <= vec->cap) return 0;
	cap = vec->cap;
	while(vec->size + n > vec->cap){
		if(vec->cap < 1024 * 1024){
			vec->cap <<= 1;
		} else {
			vec->cap += 1024 * 1024;
		}
	}
	vec->bits = realloc(vec->bits, (vec->cap * 2 + 7) / 8);
	memset(vec->bits + (cap * 2 + 7) / 8, 0, (vec->cap * 2 + 7) / 8 - (cap * 2 + 7) / 8);
	return 1;
}

static inline void set_bit2vec(Bit2Vec *vec, uint64_t idx, uint8_t dat){
	vec->bits[idx >> 2] = (vec->bits[idx >> 2] & (~(3U << ((idx & 0x03U) << 1)))) | ((dat & 0x03) << ((idx & 0x03U) << 1));
}

static inline void push_bit2vec(Bit2Vec *vec, uint8_t dat){
	encap_bit2vec(vec, 1);
	set_bit2vec(vec, vec->size, dat);
	vec->size ++;
}

static inline uint8_t get_bit2vec(Bit2Vec *vec, uint64_t idx){
	return (vec->bits[idx >> 2] >> ((idx & 0x03U) << 1)) & 0x03;
}

static inline int pop_bit2vec(Bit2Vec *vec){
	if(vec->size == 0) return -1;
	vec->size --;
	return get_bit2vec(vec, vec->size);
}

#endif
