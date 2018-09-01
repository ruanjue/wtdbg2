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
 
#ifndef __SORT_RJ_H
#define __SORT_RJ_H

#include <stdio.h>
#include <stdlib.h>

#define cmp_2nums_proc(a, b) if((a) < (b)) return -1; else if((a) > (b)) return 1;
#define num_cmp_script(e1, e2, obj, val_macro) ((val_macro(e1, obj) == val_macro(e2, obj))? 0 : ((val_macro(e1, obj) < val_macro(e2, obj))? -1 : 1))
#define cmpgt_2nums_proc(a, b) if((a) < (b)) return 0; else if((a) > (b)) return 1;

#define define_bubble_sort(name, e_type, is_greater_func)	\
static inline void name(e_type* list, size_t size, void *ref){	\
	size_t i, j, n;	\
	e_type t;	\
	i = 0;	\
	while(i < size){	\
		n = 0;	\
		for(j=size-1;j>i;j--){	\
			if(is_greater_func(list[j-1], list[j], ref) > 0){	\
				t = list[j-1]; list[j-1] = list[j]; list[j] = t;	\
				n = 1;	\
			}	\
		}	\
		if(n == 0) break;	\
		i ++;	\
	}	\
	if(ref == ref) return;	\
}

#define bubble_sort_array(rs, size, e_type, is_a_greater_than_b)	\
do {	\
	size_t bubble_i, bubble_j, bubble_n, bubble_size;	\
	e_type a, b;	\
	bubble_size = size;	\
	for(bubble_i=0;bubble_i<bubble_size;bubble_i++){	\
		bubble_n = 0;	\
		for(bubble_j=bubble_size-1;bubble_j>bubble_i;bubble_j--){	\
			a = (rs)[bubble_j - 1];	\
			b = (rs)[bubble_j];	\
			if((int)(is_a_greater_than_b) > 0){	\
				(rs)[bubble_j] = a; (rs)[bubble_j - 1] = b;	\
				bubble_n = 1;	\
			}	\
		}	\
		if(bubble_n == 0) break;	\
	}	\
} while(0)

#define divide_array(rs_ary, rs_size, e_type, is_a_greater_than_b, ret_val)	\
do {	\
	e_type *_rs;	\
	_rs = (e_type*)(rs_ary);	\
	size_t s, e, i, j, m;	\
	e_type p, t, a, b;	\
	if((rs_size) < 2){ (ret_val) = 0; break; }	\
	{	\
		s = 0;	\
		e = (rs_size) - 1;	\
		m = s + (e - s) / 2;	\
		a = _rs[s]; b = _rs[m];	\
		if(is_a_greater_than_b){ t = _rs[s]; _rs[s] = _rs[m]; _rs[m] = t; }	\
		a = _rs[m]; b = _rs[e];	\
		if(is_a_greater_than_b){	\
			t = _rs[e]; _rs[e] = _rs[m]; _rs[m] = t;	\
			a = _rs[s]; b = _rs[m];	\
			if(is_a_greater_than_b){ t = _rs[s]; _rs[s] = _rs[m]; _rs[m] = t; }	\
		}	\
		p = _rs[m];	\
		i = s + 1; j = e - 1;	\
		while(1){	\
			a = p;	\
			while(b = _rs[i], (is_a_greater_than_b)) i ++;	\
			b = p;	\
			while(a = _rs[j], (is_a_greater_than_b)) j --;	\
			if(i < j){	\
				t = _rs[i]; _rs[i] = _rs[j]; _rs[j] = t;	\
				i ++; j --;	\
			} else break;	\
		}	\
		if(i == j){ i ++; j --; }	\
		(ret_val) = i;	\
	}	\
} while(0)

#define sort_array(rs_ary, rs_size, e_type, is_a_greater_than_b)	\
do {	\
	e_type *_rs;	\
	_rs = (e_type*)(rs_ary);	\
	size_t _qsort_n;	\
	_qsort_n = rs_size;	\
	size_t s, e, i, j, m, stack[64][2], x;	\
	e_type p, t, a, b;	\
	if(_qsort_n < 2) break;	\
	x = 0;	\
	stack[x][0] = 0; stack[x][1] = _qsort_n - 1; x ++;	\
	while(x){	\
		x --; s = stack[x][0]; e = stack[x][1];	\
		m = s + (e - s) / 2;	\
		a = _rs[s]; b = _rs[m];	\
		if(is_a_greater_than_b){ t = _rs[s]; _rs[s] = _rs[m]; _rs[m] = t; }	\
		a = _rs[m]; b = _rs[e];	\
		if(is_a_greater_than_b){	\
			t = _rs[e]; _rs[e] = _rs[m]; _rs[m] = t;	\
			a = _rs[s]; b = _rs[m];	\
			if(is_a_greater_than_b){ t = _rs[s]; _rs[s] = _rs[m]; _rs[m] = t; }	\
		}	\
		p = _rs[m];	\
		i = s + 1; j = e - 1;	\
		while(1){	\
			a = p;	\
			while(b = _rs[i], (is_a_greater_than_b)) i ++;	\
			b = p;	\
			while(a = _rs[j], (is_a_greater_than_b)) j --;	\
			if(i < j){	\
				t = _rs[i]; _rs[i] = _rs[j]; _rs[j] = t;	\
				i ++; j --;	\
			} else break;	\
		}	\
		if(i == j){ i ++; j --; }	\
		if(j - s > e - i){	\
			if(s + 4 < j){ stack[x][0] = s; stack[x][1] = j; x ++; }	\
			if(i + 4 < e){ stack[x][0] = i; stack[x][1] = e; x ++; }	\
		} else {	\
			if(i + 4 < e){ stack[x][0] = i; stack[x][1] = e; x ++; }	\
			if(s + 4 < j){ stack[x][0] = s; stack[x][1] = j; x ++; }	\
		}	\
	}	\
	for(i=0;i<_qsort_n;i++){	\
		x = 0;	\
		for(j=_qsort_n-1;j>i;j--){	\
			a = _rs[j - 1]; b = _rs[j];	\
			if(is_a_greater_than_b){ t = _rs[j - 1]; _rs[j - 1] = _rs[j]; _rs[j] = t; x = 1; }	\
		}	\
		if(x == 0) break;	\
	}	\
} while(0)

#define sort_array_adv(rs_size, is_a_greater_than_b, swap_expr)	\
do {	\
	size_t _qsort_n;	\
	_qsort_n = rs_size;	\
	size_t s, e, i, j, m, stack[64][2], x, a, b;	\
	if(_qsort_n < 2) break;	\
	x = 0;	\
	stack[x][0] = 0; stack[x][1] = _qsort_n - 1; x ++;	\
	while(x){	\
		x --; s = stack[x][0]; e = stack[x][1];	\
		m = s + (e - s) / 2;	\
		a = s; b = m;	\
		if(is_a_greater_than_b){ swap_expr; }	\
		a = m; b = e;	\
		if(is_a_greater_than_b){	\
			swap_expr;	\
			a = s; b = m;	\
			if(is_a_greater_than_b){ swap_expr; }	\
		}	\
		i = s + 1; j = e - 1;	\
		while(1){	\
			a = m;	\
			while(b = i, (is_a_greater_than_b)) i ++;	\
			b = m;	\
			while(a = j, (is_a_greater_than_b)) j --;	\
			if(i < j){	\
				if(m == i) m = j;	\
				else if(m == j) m = i;	\
				a = i; b = j; \
				swap_expr; \
				i ++; j --;	\
			} else break;	\
		}	\
		if(i == j){ i ++; j --; }	\
		if(s + 4 < j){ stack[x][0] = s; stack[x][1] = j; x ++; }	\
		if(i + 4 < e){ stack[x][0] = i; stack[x][1] = e; x ++; }	\
	}	\
	for(i=0;i<_qsort_n;i++){	\
		x = 0;	\
		for(j=_qsort_n-1;j>i;j--){	\
			a = j - 1; b = j;	\
			if(is_a_greater_than_b){ swap_expr; x = 1; }	\
		}	\
		if(x == 0) break;	\
	}	\
} while(0)

// Must #include "thread.h" and "list.h"
#define psort_array(rs_ary, rs_size, e_type, ncpu, is_a_greater_than_b)	\
do {	\
	thread_beg_def(psrt);	\
	e_type *rs;	\
	size_t beg, end, div;	\
	int divide;	\
	thread_end_def(psrt);	\
	thread_beg_func_inline(psrt);	\
	thread_beg_loop(psrt);	\
	if(psrt->divide){	\
		divide_array(psrt->rs + psrt->beg, psrt->end - psrt->beg, e_type, is_a_greater_than_b, psrt->div);	\
		psrt->div += psrt->beg;	\
	} else {	\
		sort_array(psrt->rs + psrt->beg, psrt->end - psrt->beg, e_type, is_a_greater_than_b);	\
	}	\
	thread_end_loop(psrt);	\
	thread_end_func(psrt);	\
	thread_preprocess(psrt);	\
	e_type *_rs;	\
	_rs = (e_type*)(rs_ary);	\
	size_t _qsort_n, _min_blk_size;	\
	_qsort_n = rs_size;	\
	_min_blk_size = _qsort_n / ncpu / 8;	\
	if(_min_blk_size < 1024) _min_blk_size = 1024;	\
	if(_qsort_n < ((uint32_t)(ncpu)) * 32){	\
		sort_array(rs_ary, rs_size, e_type, is_a_greater_than_b);	\
		break;	\
	}	\
	thread_beg_init(psrt, (int)(ncpu));	\
	psrt->rs = _rs; psrt->beg = psrt->end = psrt->div = 0; psrt->divide = 0;	\
	thread_end_init(psrt);	\
	u64list *stacks[2];	\
	int x;	\
	stacks[0] = init_u64list(32);	\
	stacks[1] = init_u64list(32);	\
	push_u64list(stacks[0], 0);	\
	push_u64list(stacks[1], _qsort_n);	\
	x = 0;	\
	while(stacks[0]->size || x > 0){	\
		thread_waitfor_one_idle(psrt);	\
		if(psrt->divide){	\
			if(psrt->div - psrt->beg <= psrt->end - psrt->div){	\
				push_u64list(stacks[0], psrt->beg);	\
				push_u64list(stacks[1], psrt->div);	\
				push_u64list(stacks[0], psrt->div);	\
				push_u64list(stacks[1], psrt->end);	\
			} else {	\
				push_u64list(stacks[0], psrt->div);	\
				push_u64list(stacks[1], psrt->end);	\
				push_u64list(stacks[0], psrt->beg);	\
				push_u64list(stacks[1], psrt->div);	\
			}	\
			x --;	\
			psrt->divide = 0;	\
		} else if(stacks[0]->size){	\
			psrt->beg = stacks[0]->buffer[--stacks[0]->size];	\
			psrt->end = stacks[1]->buffer[--stacks[1]->size];	\
			psrt->divide = (psrt->end - psrt->beg > _min_blk_size);	\
			if(psrt->divide) x ++;	\
			thread_wake(psrt);	\
		}	\
	}	\
	thread_waitfor_all_idle(psrt);	\
	thread_beg_close(psrt);	\
	thread_end_close(psrt);	\
	free_u64list(stacks[0]);	\
	free_u64list(stacks[1]);	\
} while(0)

#define quick_median_array(_rs, _rs_size, e_type, expr)	\
({	\
	e_type key;	\
	do {	\
		e_type *rs;	\
		rs = (e_type*)(_rs);	\
		size_t size;	\
		size = (size_t)(_rs_size);	\
		size_t i, j, beg, mid, end;	\
		if(size == 0){	\
			memset(&key, 0, sizeof(e_type));	\
			break;	\
		}	\
		beg = 0;	\
		end = size - 1;	\
		e_type tmp, a, b;	\
		while(beg < end){	\
			mid = beg + (end - beg) / 2;	\
			a = rs[beg]; b = rs[mid];	\
			if(expr){ tmp = rs[beg]; rs[beg] = rs[mid]; rs[mid] = tmp; }	\
			a = rs[mid]; b = rs[end];	\
			if(expr){	\
				tmp = rs[end]; rs[end] = rs[mid]; rs[mid] = tmp;	\
				a = rs[beg]; b = rs[mid];	\
				if(expr){ tmp = rs[beg]; rs[beg] = rs[mid]; rs[mid] = tmp; }	\
			}	\
			key = rs[mid];	\
			i = beg + 1; j = end - 1;	\
			while(1){	\
				a = key;	\
				while(b = rs[i], (expr)) i ++;	\
				b = key;	\
				while(a = rs[j], (expr)) j --;	\
				if(i < j){	\
					tmp = rs[i]; rs[i] = rs[j]; rs[j] = tmp;	\
					i ++; j --;	\
				} else break;	\
			}	\
			if(i == j){ i ++; j --; }	\
			if(i <= size / 2) beg = i;	\
			else end = j;	\
		}	\
		key = rs[size/2];	\
	} while(0);	\
	key;	\
})


#define apply_array(rs, rs_size, e_type, expression)	\
do {	\
	size_t _i, _rs_size;	\
	e_type a;	\
	_rs_size = rs_size;	\
	for(_i=0;_i<_rs_size;_i++){	\
		a = (rs)[_i];	\
		expression;	\
	}	\
} while(0)

#define ref_apply_array(rs, rs_size, e_type, expression)	\
do {	\
	size_t _i, _rs_size;	\
	e_type *a;	\
	_rs_size = rs_size;	\
	for(_i=0;_i<_rs_size;_i++){	\
		a = (rs) + _i;	\
		(expression);	\
	}	\
} while(0)

#define locate_array(rs, rs_size, e_type, expr)		\
({	\
	 size_t _i, _size;	\
	 e_type a;	\
	 _size = rs_size;	\
	 for(_i=0;_i<_size;_i++){	\
		 a = rs[_i];	\
		 if(expr) break;	\
	 };	\
	 _i;	\
 })

// sort the array according to bool value (true then flase), and return the size of trues
#define apply_xchg_array(rs, rs_size, e_type, expr)	\
({	\
	size_t _i, _j, _size;	\
	e_type a;	\
	_size = rs_size;	\
	for(_i=_j=0;_i<_size;_i++){	\
		a = (rs)[_i];	\
		if(!(expr)) continue;	\
		if(_j < _i){	\
			a = (rs)[_j];	\
			(rs)[_j] = (rs)[_i];	\
			(rs)[_i] = a;	\
		}	\
		_j ++;	\
	}	\
	_j;	\
})

#define define_quick_sort(name, e_type, is_greater_func)	\
static inline void name(e_type *rs, size_t n, void *obj){	\
	size_t s, e, i, j, m, stack[64][2], x;	\
	e_type p, t;	\
	if(n < 2) return;	\
	x = 0;	\
	stack[x][0] = 0; stack[x][1] = n - 1; x ++;	\
	while(x){	\
		x --; s = stack[x][0]; e = stack[x][1];	\
		m = s + (e - s) / 2;	\
		if(is_greater_func(rs[s], rs[m], obj) > 0){ t = rs[s]; rs[s] = rs[m]; rs[m] = t; }	\
		if(is_greater_func(rs[m], rs[e], obj) > 0){	\
			t = rs[e]; rs[e] = rs[m]; rs[m] = t;	\
			if(is_greater_func(rs[s], rs[m], obj) > 0){ t = rs[s]; rs[s] = rs[m]; rs[m] = t; }	\
		}	\
		p = rs[m];	\
		i = s + 1; j = e - 1;	\
		while(1){	\
			while(is_greater_func(p, rs[i], obj) > 0) i ++;	\
			while(is_greater_func(rs[j], p, obj) > 0) j --;	\
			if(i < j){	\
				t = rs[i]; rs[i] = rs[j]; rs[j] = t;	\
				i ++; j --;	\
			} else break;	\
		}	\
		if(i == j){ i ++; j --; }	\
		if(j - s > e - i){	\
			if(s + 4 < j){ stack[x][0] = s; stack[x][1] = j; x ++; }	\
			if(i + 4 < e){ stack[x][0] = i; stack[x][1] = e; x ++; }	\
		} else {	\
			if(i + 4 < e){ stack[x][0] = i; stack[x][1] = e; x ++; }	\
			if(s + 4 < j){ stack[x][0] = s; stack[x][1] = j; x ++; }	\
		}	\
	}	\
	for(i=0;i<n;i++){	\
		x = 0;	\
		for(j=n-1;j>i;j--){	\
			if(is_greater_func(rs[j - 1], rs[j], obj) > 0){ t = rs[j - 1]; rs[j - 1] = rs[j]; rs[j] = t; x = 1; }	\
		}	\
		if(x == 0) break;	\
	}	\
	if(obj == obj) return;	\
}

#define define_merge(name, e_type, cmp_func, output_func)	\
static inline void name(e_type *list1, size_t size1, e_type *list2, size_t size2, void *ref){	\
	size_t i, j;	\
	i = j = 0;	\
	while(i < size1 && j < size2){	\
		if(cmp_func(list1[i], list2[j], ref) <= 0){	\
			output_func(list1[i], ref);	\
			i ++;	\
		} else {	\
			output_func(list2[j], ref);	\
			j ++;	\
		}	\
	}	\
	while(i < size1){ output_func(list1[i++], ref); }	\
	while(j < size2){ output_func(list2[j++], ref); }	\
}	\
	\
static inline size_t name##_files(FILE **files, int n, void *ref){	\
	e_type *es;	\
	int *flags, i, min;	\
	size_t ret;	\
	ret = 0;	\
	es = malloc(sizeof(e_type) * n);	\
	flags = malloc(sizeof(int) * n);	\
	for(i=0;i<n;i++) flags[i] = 0;	\
	while(1){	\
		min = -1;	\
		for(i=0;i<n;i++){	\
			if(flags[i] == 0){	\
				flags[i] = (fread(es + i, sizeof(e_type), 1, files[i]) == 1)? 1 : 2;	\
			}	\
			if(flags[i] == 1){	\
				if(min == -1){	\
					min = i;	\
				} else if(cmp_func(es[i], es[min], ref) <= 0){	\
					min = i;	\
				}	\
			}	\
		}	\
		if(min == -1) break;	\
		output_func(es[min], ref);	\
		flags[min] = 0;	\
		ret ++;	\
	}	\
	free(es);	\
	free(flags);	\
	return ret;	\
}

#define define_unique_merge(name, e_type, cmp_func, output_func)	\
static inline void name(e_type *list1, size_t size1, e_type *list2, size_t size2, void *ref){	\
	size_t i, j;	\
	i = j = 0;	\
	while(i < size1 && j < size2){	\
		switch(cmp_func(list1[i], list2[j])){	\
			case 0:  output_func(list1[i++], ref); j ++; break;	\
			case -1: output_func(list1[i++], ref); break;	\
			default: output_func(list2[j++], ref);	\
		}	\
	}	\
	while(i < size1){ output_func(list1[i++], ref); }	\
	while(j < size2){ output_func(list2[j++], ref); }	\
}

#define reverse_array(_rs, _rs_size, e_type)	\
do {	\
	size_t i, n;	\
	n = (_rs_size);	\
	e_type* __rs = (_rs);	\
	e_type e;	\
	for(i=0;i<n>>1;i++){	\
		e = __rs[i]; __rs[i] = __rs[n-1-i]; __rs[n-1-i] = e;	\
	}	\
} while(0)

#define define_reverse_array(name, e_type)	\
static inline void name(e_type *list, size_t size){	\
	size_t i, j;	\
	e_type t;	\
	if(size == 0) return;	\
	i = 0;	\
	j = size - 1;	\
	while(i < j){	\
		t = list[i]; list[i] = list[j]; list[j] = t;	\
		i ++; j --;	\
	}	\
}

#define define_apply_array(name, e_type, apply_func)	\
static inline size_t name(e_type *list, size_t size, void *ref){	\
	size_t i, ret;	\
	ret = 0;	\
	for(i=0;i<size;i++){	\
		ret += apply_func(list[i], ref);	\
	}	\
	return ret;	\
	ref = NULL;	\
}

#define define_search_array(name, e_type, cmp_func)	\
static inline long long name(e_type *array, long long size, e_type key, void *ref){	\
	long long i, j, m;	\
	i = 0;	\
	j = size;	\
	while(i < j){	\
		m = i + (j - i) / 2;	\
		if(cmp_func(array[m], key, ref) < 0){	\
			i = m + 1;	\
		} else {	\
			j = m;	\
		}	\
	}	\
	if(i < size && cmp_func(array[i], key, ref) == 0) return i;	\
	else return - (i + 1);	\
	if(ref) return 0;	\
}

#define bsearch_array(_rs_array, _rs_size, e_type, _bs_ret, is_a_less_than_your)	\
	\
do {	\
	e_type*_bs_rs;	\
	e_type a;	\
	_bs_rs = (e_type*)(_rs_array);	\
	size_t _bs_size;	\
	_bs_size = _rs_size;	\
	size_t _bs_i, _bs_j, _bs_m;	\
	_bs_i = 0;	\
	_bs_j = _bs_size;	\
	while(_bs_i < _bs_j){	\
		_bs_m = _bs_i + (_bs_j - _bs_i) / 2;	\
		a = _bs_rs[_bs_m];	\
		if(is_a_less_than_your){	\
			_bs_i = _bs_m + 1;	\
		} else {	\
			_bs_j = _bs_m;	\
		}	\
	}	\
	_bs_ret = _bs_i;	\
} while(0)

#endif
