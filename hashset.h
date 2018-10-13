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
 
#ifndef __HASH_SET_RJ
#define __HASH_SET_RJ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "mem_share.h"
#include "bitvec.h"

static const uint64_t sys_prime_list[61] = {
	0x7LLU, 0xfLLU, 0x1fLLU, 0x43LLU, 0x89LLU,
	0x115LLU, 0x22dLLU, 0x45dLLU, 0x8bdLLU, 0x1181LLU,
	0x2303LLU, 0x4609LLU, 0x8c17LLU, 0x1183dLLU, 0x2307bLLU,
	0x460fdLLU, 0x8c201LLU, 0x118411LLU, 0x230833LLU, 0x461069LLU,
	0x8c20e1LLU, 0x11841cbLLU, 0x2308397LLU, 0x461075bLLU, 0x8c20ecbLLU,
	0x11841da5LLU, 0x23083b61LLU, 0x461076c7LLU, 0x8c20ed91LLU, 0x11841db31LLU,
	0x23083b673LLU, 0x461076d1bLLU, 0x8c20eda41LLU, 0x11841db48dLLU, 0x23083b6937LLU,
	0x461076d27fLLU, 0x8c20eda50dLLU, 0x11841db4a59LLU, 0x23083b694ebLLU, 0x461076d29f1LLU,
	0x8c20eda5441LLU, 0x11841db4a887LLU, 0x23083b69511fLLU, 0x461076d2a2c1LLU, 0x8c20eda54591LLU,
	0x11841db4a8b55LLU, 0x23083b69516c1LLU, 0x461076d2a2da5LLU, 0x8c20eda545b55LLU, 0x11841db4a8b6b5LLU,
	0x23083b69516d91LLU, 0x461076d2a2db3bLLU, 0x8c20eda545b69dLLU, 0x11841db4a8b6d5dLLU, 0x23083b69516daf5LLU,
	0x461076d2a2db5edLLU, 0x8c20eda545b6c5fLLU, 0x11841db4a8b6d8ebLLU, 0x23083b69516db1ffLLU, 0x461076d2a2db643fLLU,
	0x8c20eda545b6c8f3LLU
};

static inline uint64_t _rj_hashset_find_prime(uint64_t n){
	uint32_t i;
	i = 0;
	while(i < 60 && n > sys_prime_list[i]) i ++;
	return sys_prime_list[i];
}

#define init_hashset_macro(hash_type, hash_ele_type) \
typedef struct { hash_ele_type *array;  BitVec *ones, *dels; size_t e_size; size_t ocp; size_t size; size_t count; size_t max; float load_factor; size_t iter_ptr; void *userdata; } hash_type; \
static inline size_t hash_type##_obj_desc_cnt(void *obj, int idx){	\
	hash_type *set;	\
	set = (hash_type*)obj;	\
	if(set->dels){	\
		switch(idx){	\
			case 0: return ((hash_type*)obj)->size * sizeof(hash_ele_type);	\
			default: return 1;	\
		}	\
	} else {	\
		switch(idx){	\
			case 0: return ((hash_type*)obj)->count * sizeof(hash_ele_type);	\
			case 1: return 1;	\
			default: return 0;	\
		}	\
	}	\
}	\
static const obj_desc_t hash_type##_obj_desc = {TOSTR(_hashset_##hash_type), sizeof(hash_type), 3, {1, 1, 1}, {offsetof(hash_type, array), offsetof(hash_type, ones), offsetof(hash_type, dels)}, {(obj_desc_t*)&OBJ_DESC_DATA, (obj_desc_t*)&bitvec_obj_desc, (obj_desc_t*)&bitvec_obj_desc}, hash_type##_obj_desc_cnt, NULL};	\
static inline int hash_type##_is_prime(uint64_t num){                          \
	uint64_t i, max;                                                           \
	if(num < 4) return 1;                                                      \
	if(num % 2 == 0) return 0;                                                 \
	max = (uint64_t)sqrt((double)num);                                         \
	for(i=3;i<max;i+=2){ if(num % i == 0) return 0; }                          \
	return 1;                                                                  \
}                                                                              \
static inline uint64_t hash_type##_find_next_prime(uint64_t num){              \
	if(num % 2 == 0) num ++;                                                   \
	while(1){ if(hash_type##_is_prime(num)) return num; num += 2; }            \
}                                                                              \
static inline hash_type* init2_##hash_type(uint32_t size, float factor){       \
	hash_type *set;                                                            \
	set = (hash_type*)calloc(1, sizeof(hash_type));                            \
	set->e_size = sizeof(hash_ele_type);                                       \
	set->size   = _rj_hashset_find_prime(size);                                \
	set->count  = 0;                                                           \
	set->ocp    = 0;                                                           \
	set->load_factor = factor;                                                 \
	set->max    = set->size * set->load_factor;                                \
	set->iter_ptr    = 0;                                                      \
	set->array       = calloc(set->size, set->e_size);                         \
	set->ones = init_bitvec(set->size);                                        \
	set->dels = init_bitvec(set->size);                                        \
	set->userdata = NULL;                                                      \
	return set;                                                                \
}                                                                              \
static inline void set_userdata_##hash_type(hash_type *set, void *userdata){ set->userdata = userdata; } \
static inline hash_type* init_##hash_type(uint32_t size){ return init2_##hash_type(size, 0.67f); }

#define get_hashset_macro(hash_type, hash_ele_type, hash_key_type, hash_key_code, hash_key_equal, hash_val_type, hash_ele2val) \
static inline hash_ele_type* get_##hash_type(hash_type *set, hash_key_type key){\
	hash_ele_type *e;                                                          \
	size_t hc, hi;                                                             \
	hc = hash_key_code(key) % set->size;                                       \
	if(set->dels){	\
		while(1){                                                                  \
			if(get_bitvec(set->ones, hc) == 0){	\
				return NULL;	\
			} else if(get_bitvec(set->dels, hc)){	\
			} else {	\
				e = ((hash_ele_type*)set->array) + hc;                             \
				if(hash_key_equal((key), (*e))) return e;                          \
			}	\
			hc = (hc + 1) % set->size;                                             \
		}                                                                          \
	} else {	\
		hi = MAX_U8;	\
		while(1){	\
			if(get_bitvec(set->ones, hc)){	\
				if(hi == MAX_U8){	\
					hi = rank_bitvec(set->ones, hc);	\
				}	\
				e = ((hash_ele_type*)set->array) + hi;	\
				if(hash_key_equal((key), (*e))) return e;	\
			} else {	\
				return NULL;	\
			}	\
			hc ++; \
			hi ++;	\
		}	\
	}	\
	return NULL;                                                               \
}                                                                              \
static inline size_t offset_##hash_type(hash_type *set, hash_ele_type *ptr){   \
	return ptr - set->array;                                                   \
}	\
static inline hash_ele_type* ref_##hash_type(hash_type *set, size_t off){ return set->array + off; }	\
static inline hash_val_type getval_##hash_type(hash_type *set, hash_key_type key){	\
	hash_ele_type *e;	\
	e = get_##hash_type(set, key);	\
	return hash_ele2val(e);	\
}

#define prepare_hashset_macro(hash_type, hash_ele_type, hash_key_type, hash_key_code, hash_key_equal) \
static inline void encap_##hash_type(hash_type *set, size_t num);              \
static inline hash_ele_type* prepare_##hash_type(hash_type *set, hash_key_type key, int *exists){\
	hash_ele_type *e;                                                          \
	size_t hc, d;                                                              \
	if(set->dels == NULL){ *exists = 0;  return NULL; }	\
	encap_##hash_type(set, 1);                                                 \
	hc = hash_key_code((key)) % set->size;                                     \
	d = set->size;                                                             \
	while(1){	\
		if(get_bitvec(set->ones, hc) == 0){	\
			if(d == set->size){	\
				one_bitvec(set->ones, hc);	\
				set->ocp ++;	\
			} else {	\
				hc = d;	\
				zero_bitvec(set->dels, hc);	\
			}	\
			if(exists) *exists = 0;	\
			set->count ++;	\
			e = ((hash_ele_type*)set->array) + hc;	\
			return e;	\
		} else if(get_bitvec(set->dels, hc)){	\
			if(d == set->size) d = hc;	\
		} else {	\
			e = ((hash_ele_type*)set->array) + hc;	\
			if(hash_key_equal((key), (*e))){	\
				if(exists) *exists = 1;	\
				return e;	\
			}	\
		}	\
		hc = (hc + 1) % set->size;	\
	}	\
	return NULL;                                                               \
}

#define exists_hashset_macro(hash_type, hash_ele_type, hash_key_type, hash_key_code, hash_key_equal) \
static inline int exists_##hash_type(hash_type *set, hash_key_type key){       \
	return get_##hash_type(set, key) != NULL;	\
}

#define add_hashset_macro(hash_type, hash_ele_type, hash_code_macro, hash_equal_macro) \
static inline hash_ele_type* add_##hash_type(hash_type *set, hash_ele_type ele){       \
	hash_ele_type *e;                                                          \
	size_t d, hc;                                                              \
	if(set->dels == NULL) return NULL;	\
	hc = hash_code_macro(ele) % set->size;                                     \
	d  = set->size;                                                            \
	do {	\
		if(get_bitvec(set->ones, hc) == 0){	\
			if(d == set->size){	\
				one_bitvec(set->ones, hc);	\
				set->ocp ++;	\
			} else {	\
				hc = d;	\
				zero_bitvec(set->dels, hc);	\
			}	\
			set->count ++;	\
			e = ((hash_ele_type*)set->array) + hc;	\
			*e = ele;	\
			return e;	\
		} else if(get_bitvec(set->dels, hc)){	\
			if(d == set->size) d = hc;	\
		} else {	\
			e = ((hash_ele_type*)set->array) + hc;	\
			if(hash_equal_macro((ele), (*e))){	\
				*e = ele;	\
				return e;	\
			}	\
		}	\
		hc = (hc + 1) % set->size;	\
	} while(1);	\
	return NULL;                                                                  \
}

#define put_hashset_macro(hash_type, hash_ele_type) \
static inline hash_ele_type* put_##hash_type(hash_type *set, hash_ele_type ele){         \
	encap_##hash_type(set, 1);                                                 \
	return add_##hash_type(set, ele);                                          \
}

#define remove_hashset_macro(hash_type, hash_ele_type, hash_key_type, hash_key_code, hash_key_equal) \
static inline int delete_##hash_type(hash_type *set, hash_ele_type *ele){	\
	size_t hc;	\
	if(set->dels == NULL) return 0;	\
	hc = offset_##hash_type(set, ele);	\
	if(get_bitvec(set->ones, (hc + 1) % set->size) == 0){	\
		zero_bitvec(set->ones, hc);	\
		set->ocp --;	\
	} else {	\
		one_bitvec(set->dels, hc);	\
	}	\
	set->count --;	\
	return 1;	\
}	\
	\
static inline int remove_##hash_type(hash_type *set, hash_key_type key){       \
	hash_ele_type *e;                                                          \
	size_t hc;                                                                 \
	if(set->dels == NULL) return 0;                                           \
	hc = hash_key_code(key) % set->size;                                       \
	while(1){	\
		if(get_bitvec(set->ones, hc) == 0){	\
			return 0;	\
		} else if(get_bitvec(set->dels, hc)){	\
		} else {	\
			e = ((hash_ele_type*)set->array) + hc;	\
			if(hash_key_equal((key), (*e))){	\
				if(get_bitvec(set->ones, (hc + 1) % set->size) == 0){	\
					zero_bitvec(set->ones, hc);	\
					set->ocp --;	\
				} else {	\
					one_bitvec(set->dels, hc);	\
				}	\
				set->count --;	\
				return 1;	\
			}	\
		}	\
		hc = (hc + 1) % set->size;	\
	}	\
	return 0;                                                                  \
}

#define reset_iter_hashset_macro(hash_type) static inline void reset_iter_##hash_type(hash_type *set){ set->iter_ptr = 0; }

#define ref_iter_hashset_macro(hash_type, hash_ele_type) \
static inline hash_ele_type* ref_iter2_##hash_type(hash_type *set, size_t *iter_ptr){             \
	if(set->dels){	\
		while(((*iter_ptr) = next_one_bitvec(set->ones, (*iter_ptr))) < set->size){	\
			if(get_bitvec(set->dels, (*iter_ptr))){	\
				(*iter_ptr) ++;	\
			} else {	\
				return (((hash_ele_type*)set->array) + (*iter_ptr)++);	\
			}	\
		}	\
	} else {	\
		while((*iter_ptr) < set->count){                                          \
				return (((hash_ele_type*)set->array) + (*iter_ptr)++);           \
		}                                                                          \
	}	\
	return NULL;                                                               \
}	\
static inline hash_ele_type* ref_iter_##hash_type(hash_type *set){             \
	return ref_iter2_##hash_type(set, &(set->iter_ptr));	\
}

#define count_hashset_macro(hash_type) static inline int64_t count_##hash_type(hash_type *set){ return set->count; }

#define freeze_hashset_macro(hash_type, hash_ele_type, hash_code_macro)	\
static inline int freeze_##hash_type(hash_type *set, float load_factor){	\
	size_t *hvs, i, j, sz;	\
	if(set->dels == NULL) return 0;	\
	if(load_factor == 0) load_factor = set->load_factor;	\
	sz = set->count / load_factor;	\
	sz = _rj_hashset_find_prime(sz);	\
	for(i=j=0;(i=next_one_bitvec(set->ones, i))<set->size;i++){	\
		if(get_bitvec(set->dels, i)) continue;	\
		if(j < i){	\
			set->array[j] = set->array[i];	\
		}	\
		j ++;	\
	}	\
	free_bitvec(set->ones);	\
	set->ones = NULL;	\
	free_bitvec(set->dels);	\
	set->dels = NULL;	\
	set->size = sz;	\
	set->load_factor = load_factor;	\
	set->ocp = set->count;	\
	set->array = realloc(set->array, (set->count + 1) * sizeof(hash_ele_type));	\
	memset(set->array + set->count, 0, sizeof(hash_ele_type));	\
	hvs = malloc(set->count * sizeof(size_t));	\
	for(i=0;i<set->count;i++){	\
		hvs[i] = hash_code_macro(set->array[i]) % sz;	\
	}	\
	sort_array_adv(set->count, hvs[a] > hvs[b], swap_var(hvs[a], hvs[b]); swap_var(set->array[a], set->array[b]));	\
	for(i=j=0;i<set->count;i++){	\
		if(j < hvs[i]) j = hvs[i];	\
		j ++;	\
	}	\
	if(j < sz) j = sz;	\
	set->ones = init_bitvec(j + 1);	\
	for(i=j=0;i<set->count;i++){	\
		if(j < hvs[i]) j = hvs[i];	\
		one_bitvec(set->ones, j);	\
		j ++;	\
	}	\
	free(hvs);	\
	index_bitvec(set->ones);	\
	return 1;	\
}

#define clear_hashset_macro(hash_type) \
static inline void clear_##hash_type(hash_type *set){                          \
	if(set->dels == NULL){	\
		return;	\
	}	\
	zeros_bitvec(set->ones);	\
	zeros_bitvec(set->dels);	\
	set->count = 0;                                                            \
	set->ocp   = 0;                                                            \
	set->iter_ptr = 0;                                                         \
}

#define free_hashset_macro(hash_type) \
static inline void free_##hash_type(hash_type *set){                           \
	free(set->array);                                                          \
	if(set->ones) free_bitvec(set->ones);                                      \
	if(set->dels) free_bitvec(set->dels);                                      \
	free(set);                                                                 \
}

#define encap_hashset_macro(hash_type, hash_ele_type, hash_code_macro) \
static inline void encap_##hash_type(hash_type *set, size_t num){             \
	BitVec *ones, *dels;	\
	size_t i, n, hc;                                                  \
	hash_ele_type key;                                                        \
	if(set->dels == NULL) return;	\
	if(set->ocp + num <= set->max) return;                                    \
	n = set->size;                                                            \
	do{ n = _rj_hashset_find_prime(n * 2); } while(n * set->load_factor < set->count + num);    \
	set->array = realloc(set->array, n * set->e_size);                        \
	if(set->array == NULL){                                                   \
		fprintf(stderr, "-- Out of memory --\n");                             \
		print_backtrace(stderr, 20);                                          \
		exit(1);                                                              \
	}                                                                         \
	ones = init_bitvec(n);	\
	dels = init_bitvec(n);	\
	set->ocp  = set->count;                                                   \
	set->max = n * set->load_factor;                                          \
	for(i=0;(i=next_one_bitvec(set->ones, i))<set->size;i++){	\
		if(get_bitvec(set->dels, i)) continue;	\
		key = ((hash_ele_type*)set->array)[i];	\
		one_bitvec(set->dels, i);	\
		while(1){	\
			hc = hash_code_macro(key) % n;	\
			while(get_bitvec(ones, hc)){	\
				hc = (hc + 1) % n;	\
			}	\
			one_bitvec(ones, hc);	\
			if(hc < set->size && get_bitvec(set->ones, hc) && get_bitvec(set->dels, hc) == 0){	\
				swap_var(key, ((hash_ele_type*)set->array)[hc]);	\
				one_bitvec(set->dels, hc);	\
			} else {	\
				((hash_ele_type*)set->array)[hc] = key;	\
				break;	\
			}	\
		}	\
	}	\
	swap_var(ones, set->ones);	\
	swap_var(dels, set->dels);	\
	set->size = n;	\
	free_bitvec(ones);	\
	free_bitvec(dels);	\
}                                                                             \
static inline size_t offsetof_##hash_type(hash_type *set, hash_ele_type *ptr){ return ptr - set->array; }	\


#define ITSELF(E) (E)
#define NUM_EQUALS(E1, E2) ((E1) == (E2))

#define define_hashtable(hash_type, hash_ele_type, hash_code_macro, hash_equal_macro, hash_key_type, hash_key_code, hash_key_equal, hash_val_type, hash_ele2val)    \
	init_hashset_macro(hash_type, hash_ele_type);                              \
	get_hashset_macro(hash_type, hash_ele_type, hash_key_type, hash_key_code, hash_key_equal, hash_val_type, hash_ele2val);    \
	prepare_hashset_macro(hash_type, hash_ele_type, hash_key_type, hash_key_code, hash_key_equal);    \
	exists_hashset_macro(hash_type, hash_ele_type, hash_key_type, hash_key_code, hash_key_equal);    \
	add_hashset_macro(hash_type, hash_ele_type, hash_code_macro, hash_equal_macro);    \
	put_hashset_macro(hash_type, hash_ele_type);                               \
	remove_hashset_macro(hash_type, hash_ele_type, hash_key_type, hash_key_code, hash_key_equal);    \
	ref_iter_hashset_macro(hash_type, hash_ele_type);                          \
	reset_iter_hashset_macro(hash_type);                                       \
	count_hashset_macro(hash_type);                                            \
	clear_hashset_macro(hash_type);                                            \
	freeze_hashset_macro(hash_type, hash_ele_type, hash_code_macro);           \
	free_hashset_macro(hash_type);                                             \
	encap_hashset_macro(hash_type, hash_ele_type, hash_code_macro);

#define define_hashset(hash_type, hash_ele_type, hash_code_macro, hash_equal_macro)  define_hashtable(hash_type, hash_ele_type, hash_code_macro, hash_equal_macro, hash_ele_type, hash_code_macro, hash_equal_macro, hash_ele_type*, ITSELF)

/* ------------------ Useful functions ------------------------------------- */

static inline uint32_t __lh3_Jenkins_hash_int(uint32_t key){
	key += (key << 12);
	key ^= (key >> 22);
	key += (key << 4);
	key ^= (key >> 9);
	key += (key << 10);
	key ^= (key >> 2);
	key += (key << 7);
	key ^= (key >> 12);
	return key;
}

static inline uint64_t __lh3_Jenkins_hash_64(uint64_t key){
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

static inline uint32_t jenkins_one_at_a_time_hash(char *key, size_t len){
	uint32_t hash, i;
	for(hash = i = 0; i < len; ++i){
		hash += key[i];
		hash += (hash << 10);
		hash ^= (hash >> 6);
	}
	hash += (hash << 3);
	hash ^= (hash >> 11);
	hash += (hash << 15);
	return hash;
}

static inline u8i invertible_hashcode(u8i x, int p){
	u8i m;
	m = 0xFFFFFFFFFFFFFFFFLLU >> (64 - p);
	x = ((~x) + (x << 21)) & m;
	x = x ^ (x >> 24);
	x = (x + (x << 3) + (x << 8)) & m;
	x = x ^ (x >> 14);
	x = (x + (x << 2) + (x << 4)) & m;
	x = x ^ (x >> 28);
	x = (x + (x << 31)) & m;
	return x;
}

static inline uint64_t hash64shift(uint64_t key){
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}


static inline uint64_t MurmurHash64A(const void * key, int len, uint32_t seed){
	const uint64_t m = 0xc6a4a7935bd1e995LLU;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end){
		uint64_t k = *data++;

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7){
	case 7: h ^= ((uint64_t)data2[6]) << 48;
	case 6: h ^= ((uint64_t)data2[5]) << 40;
	case 5: h ^= ((uint64_t)data2[4]) << 32;
	case 4: h ^= ((uint64_t)data2[3]) << 24;
	case 3: h ^= ((uint64_t)data2[2]) << 16;
	case 2: h ^= ((uint64_t)data2[1]) << 8;
	case 1: h ^= ((uint64_t)data2[0]);
	        h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}

#define u32hashcode(key) __lh3_Jenkins_hash_int(key)
#define u64hashcode(key) __lh3_Jenkins_hash_64(key)

static inline uint32_t __string_hashcode(const char *s){
	uint32_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}

#define u32hash_code(e) u32hashcode(e)
#define u64hash_code(e) u64hashcode(e)
#define uxxhash_equals(e1, e2) ((e1) == (e2))
define_hashset(u32hash, uint32_t, u32hash_code, uxxhash_equals);
define_hashset(u64hash, uint64_t, u64hash_code, uxxhash_equals);

#define i32hash_code(e) u32hashcode((uint32_t)(e))
#define i32hash_equals(e1, e2) ((e1) == (e2))
define_hashset(i32hash, int, i32hash_code, i32hash_equals);

#define chash_code(e) __string_hashcode(e)
#define chash_equals(e1, e2) (strcmp(e1, e2) == 0)
define_hashset(chash, char*, chash_code, chash_equals);

#define KV_HASH_GET_VAL(e) (e)? (e)->val : ((typeof(e->val))MAX_U8)

typedef struct { u4i key, val; } uuhash_t;
#define uuhash_code(e) u32hashcode((e).key)
#define uuhash_equals(e1, e2) ((e1).key == (e2).key)
#define uuhash_key_equals(e1, e2) ((e1) == (e2).key)
define_hashtable(uuhash, uuhash_t, uuhash_code, uuhash_equals, u4i, u32hashcode, uuhash_key_equals, u4i, KV_HASH_GET_VAL);

typedef struct { u4i key; int val; } uihash_t;
#define uihashcode(E) u32hashcode((E).key)
#define uihashequals(E1, E2) (E1).key == (E2).key
#define uihashkeyequals(E1, E2) (E1) == (E2).key
define_hashtable(uihash, uihash_t, uihashcode, uihashequals, u4i, u32hashcode, uihashkeyequals, b4i, KV_HASH_GET_VAL);

typedef struct { u8i key, val; } UUhash_t;
#define UUhashcode(E) u64hashcode((E).key)
#define UUhashequals(E1, E2) (E1).key == (E2).key
#define UUhashkeyequals(E1, E2) (E1) == (E2).key
define_hashtable(UUhash, UUhash_t, UUhashcode, UUhashequals, u8i, u64hashcode, UUhashkeyequals, u8i, KV_HASH_GET_VAL);

typedef struct { char *key; u4i val; } cuhash_t;
#define cuhash_code(e) __string_hashcode((e).key)
#define cuhash_equals(e1, e2) (strcmp((e1).key, (e2).key) == 0)
#define cuhash_key_equals(e1, e2) (strcmp((char*)(e1), (e2).key) == 0)
define_hashtable(cuhash, cuhash_t, cuhash_code, cuhash_equals, char*, __string_hashcode, cuhash_key_equals, u4i, KV_HASH_GET_VAL);
static const obj_desc_t cuhash_struct_deep_obj_desc = {"cuhash_struct_deep_obj_desc", sizeof(cuhash_t), 1, {1}, {offsetof(cuhash_t, key)}, {(obj_desc_t*)&OBJ_DESC_CHAR_ARRAY}, NULL, NULL};
static const obj_desc_t cuhash_deep_obj_desc = {"cuhash_deep_obj_desc", sizeof(cuhash), 3, {1, 1, 1}, {offsetof(cuhash, array), offsetof(cuhash, ones), offsetof(cuhash, dels)}, {(obj_desc_t*)&cuhash_struct_deep_obj_desc, (obj_desc_t*)&bitvec_obj_desc, &bitvec_obj_desc}, cuhash_obj_desc_cnt, NULL};

typedef struct { char *key; int val; } cihash_t;
#define cihash_code(e) __string_hashcode((e).key)
#define cihash_equals(e1, e2) (strcmp((e1).key, (e2).key) == 0)
#define cihash_key_equals(e1, e2) (strcmp((char*)(e1), (e2).key) == 0)
define_hashtable(cihash, cihash_t, cihash_code, cihash_equals, char*, __string_hashcode, cihash_key_equals, b4i, KV_HASH_GET_VAL);

typedef struct { char *key; unsigned long long val; } clhash_t;
#define clhash_code(e) __string_hashcode((e).key)
#define clhash_equals(e1, e2) (strcmp((e1).key, (e2).key) == 0)
#define clhash_key_equals(e1, e2) (strcmp((char*)(e1), (e2).key) == 0)
define_hashtable(clhash, clhash_t, clhash_code, clhash_equals, char*, __string_hashcode, clhash_key_equals, u8i, KV_HASH_GET_VAL);

typedef struct { char *key; char *val; } cchash_t;
#define cchash_code(e) __string_hashcode((e).key)
#define cchash_equals(e1, e2) (strcmp((e1).key, (e2).key) == 0)
#define cchash_key_equals(e1, e2) (strcmp((char*)(e1), (e2).key) == 0)
#define KV_CCHASH_GET_VAL(e) ((e)? (e)->val : NULL)
define_hashtable(cchash, cchash_t, cchash_code, cchash_equals, char*, __string_hashcode, cchash_key_equals, char*, KV_CCHASH_GET_VAL);

/**
* Example of using userdata in thread-safe mode
* char **strs;
* ... codes init strs
* #define test_hc(E) __string_hashcode(((char**)set->userdata)[E])
* #define test_he(E1, E2) (strcmp(((char**)set->userdata)[E1], ((char**)set->userdata)[E2]) == 0)
* define_hashset(testhash, uint32_t, test_hc, test_he);
* testhash *hash = init_testhash(13);
* set_userdata_testhash(hash, strs);
* ... now, the key of testhash is uint32_t, but refer to strs
*/

#endif
