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
 
#ifndef __DNA_RJ_H
#define __DNA_RJ_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "list.h"
#include "bitvec.h"
#include "hashset.h"
#include "thread.h"

static const u1i base_bit_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static const u1i base_bit4_table[256] = {
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,

	15,  1, 14,  2,  13, 15, 15,  4,  11, 15, 15, 12,  15,  3, 15, 15,
	15, 15,  5,  6,   8, 15,  7,  9,  15, 10, 15, 15,  15, 15, 15, 15,
	15,  1, 14,  2,  13, 15, 15,  4,  11, 15, 15, 12,  15,  3, 15, 15,
	15, 15,  5,  6,   8, 15,  7,  9,  15, 10, 15, 15,  15, 15, 15, 15,

	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,

	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15
};

static const u1i bit4_bit_table[16] = { 4, 0, 1, 4,  2, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4 };

static const char bit_base_table[12] = "ACGTN-acgtn*";
static const char bit4_base_table[16] = "-ACMGRSVTWYHKDBN";

// u8i = 0|1|2|3|4|5|6|...
#define bits2bit(bits, off) (((bits)[(off) >> 5] >> (((~(off)) & 0x1FU) << 1)) & 0x03U)
#define bits2revbit(bits, off) ((~((bits)[(off) >> 5] >> (((~(off)) & 0x1FU) << 1))) & 0x03U)

static inline u8i dna_xor2ones(u8i seq){
	return ((seq & 0xAAAAAAAAAAAAAAAALLU) >> 1) | (seq & 0x5555555555555555LLU);
}

static inline u8i dna_rev_seq32(u8i seq){
	seq = ~seq;
	seq = ((seq & 0x3333333333333333LLU)<< 2) | ((seq & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
#if 0
	seq = ((seq & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq & 0xFF00FF00FF00FF00LLU)>> 8);
	seq = ((seq & 0x0000FFFF0000FFFFLLU)<<16) | ((seq & 0xFFFF0000FFFF0000LLU)>>16);
	seq = ((seq & 0x00000000FFFFFFFFLLU)<<32) | ((seq & 0xFFFFFFFF00000000LLU)>>32);
#else
	seq = __builtin_bswap64(seq);
#endif
	return seq;
}

static inline u8i dna_rev_seq(u8i seq, u1i seq_size){
	return dna_rev_seq32(seq) >> (64 - (seq_size<<1));
}

// order of 2-bit in this->seqs is different with that in dna_rev_seq->seq
static inline void dna_rev_seqs(u8i *seqs, u8i seq_size){
	register u8i t;
	int i, j;
	register u1i d, e;
	j = (seq_size + 31) >> 5;
	// Swap within 64bit
	for(i=0;i<j;i++){
		seqs[i] = ~seqs[i];
		seqs[i] = ((seqs[i] & 0x3333333333333333LLU)<< 2) | ((seqs[i] & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
		seqs[i] = ((seqs[i] & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seqs[i] & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
		seqs[i] = __builtin_bswap64(seqs[i]);
	}
	// Swap 64bit blocks
	for(i=0;i<j>>1;i++){
		t = seqs[i]; seqs[i] = seqs[j - i - 1]; seqs[j - i - 1] = t;
	}
	// left-align seqs
	if((d = ((j << 5) - seq_size) << 1)){
		e = 64 - d;
		for(i=0;i<j-1;i++){
			seqs[i] = (seqs[i] << d) | (seqs[i+1] >> e);
		}
		seqs[i] = (seqs[i] << d) | 0;
	}
}

//shift one base, and append one base, useful to build big-kmer
static inline void dna_shl_seqs(u8i *seqs, u8i seq_size, u1i base_val){
	const u1i d = 2;
	const u1i e = 62;
	int i, j;
	j = (seq_size + 31) >> 5;
	for(i=0;i<j-1;i++){
		seqs[i] = (seqs[i] << d) | (seqs[i+1] >> e);
	}
	seqs[i] = (seqs[i] << d) | (((u8i)(base_val & 0x03U)) << ((32 - (seq_size & 0x1FU)) << 1));
}

static inline int dna_cmp_seqs(u8i *seqs1, u8i *seqs2, u8i seq_size){
	int i, j;
	j = (seq_size + 31) >> 5;
	for(i=0;i<j;i++){
		if(seqs1[i] < seqs2[i]) return -1;
		if(seqs1[i] > seqs2[i]) return 1;
	}
	return 0;
}

static inline int dna_cmpx_seqs(u8i *seqs, u8i seq_size){
	register int i, j;
	register u1i a, b;
	j = (seq_size + 1) >> 1;
	for(i=0;i<j;i++){
		a = bits2bit(seqs, i);
		b = (~bits2bit(seqs, seq_size - 1 - i)) & 0x03U;
		if(a < b) return -1;
		if(a > b) return 1;
	}
	return 0;
}

static inline u8i seq2kmer(char *seq, u4i ksize){
	u8i kmer;
	u4i i;
	kmer = 0;
	for(i=0;i<ksize;i++) kmer = (kmer << 2) | base_bit_table[(int)seq[i]];
	return kmer;
}

static inline u8i seq2revkmer(char *seq, u4i ksize){
	u8i kmer;
	u4i i;
	kmer = 0;
	for(i=0;i<ksize;i++) kmer = (kmer << 2) | ((~base_bit_table[(int)seq[ksize - 1 - i]]) & 0x03);
	return kmer;
}

static inline void kmer2seq(char *seq, u8i kmer, u4i ksize){
	u4i i;
	for(i=0;i<ksize;i++){
		seq[i] = bit_base_table[(kmer >> ((ksize - 1 - i) << 1)) & 0x03];
	}
	seq[i] = 0;
}

static inline void kmer2revseq(char *seq, u8i kmer, u4i ksize){
	u4i i;
	kmer = ~kmer;
	for(i=0;i<ksize;i++){
		seq[i] = bit_base_table[(kmer >> (i << 1)) & 0x03];
	}
	seq[i] = 0;
}

static inline void print_kmer_seq(u8i kmer, u4i ksize, FILE *out){
	char seq[33];
	kmer2seq(seq, kmer, ksize);
	fputs(seq, out);
}

static inline void print_kmer_revseq(u8i kmer, u4i ksize, FILE *out){
	char seq[33];
	kmer2revseq(seq, kmer, ksize);
	fputs(seq, out);
}

#define kmer_mask(ksize) (0xFFFFFFFFFFFFFFFFLLU >> ((32 - (ksize)) * 2))

#define beg_seq2kmers(seq, seqlen, ksize, kmask, kmer, idx) {	\
u1i beg_seq2kmers_v;	\
kmer = 0;	\
for(idx=0;(int)idx+1<(int)ksize;idx++){	\
	beg_seq2kmers_v = base_bit_table[(int)(seq)[idx]];	\
	if(beg_seq2kmers_v == 4) beg_seq2kmers_v = lrand48() & 0x03;	\
	kmer = (((kmer) << 2) | beg_seq2kmers_v);	\
}	\
for(idx=0;(int)idx<=(int)(seqlen-ksize);idx++){	\
	beg_seq2kmers_v = base_bit_table[(int)(seq)[idx + (ksize) - 1]];	\
	if(beg_seq2kmers_v == 4) beg_seq2kmers_v = lrand48() & 0x03;	\
	kmer = ((kmer << 2) | beg_seq2kmers_v) & kmask;
#define end_seq2kmers } }

#define beg_seq2revkmers(seq, seqlen, ksize, kmask, kmer, idx) {	\
u1i beg_seq2revkmers_v;	\
kmer = 0;	\
for(idx=0;(int)idx+1<(int)ksize;idx++){	\
	beg_seq2revkmers_v = base_bit_table[(int)(seq)[seqlen - 1 - idx]];	\
	if(beg_seq2revkmers_v == 4) beg_seq2revkmers_v = lrand48() & 0x03;	\
	kmer = (((kmer) << 2) | beg_seq2revkmers_v);	\
}	\
for(idx=0;(int)idx<=(int)seqlen-ksize;idx++){	\
	beg_seq2revkmers_v = base_bit_table[(int)(seq)[seqlen - idx - (ksize)]];	\
	if(beg_seq2revkmers_v == 4) beg_seq2revkmers_v = lrand48() & 0x03;	\
	kmer = ((kmer << 2) | beg_seq2revkmers_v) & kmask;
#define end_seq2revkmers } }

static inline char reverse_dna_base(char b){
	switch(b){
		case 'a': return 't';
		case 'A': return 'T';
		case 'c': return 'g';
		case 'C': return 'G';
		case 'g': return 'c';
		case 'G': return 'C';
		case 't': return 'a';
		case 'T': return 'A';
		default: return 'N';
	}
}

static inline void reverse_dna(char *seq, int len){
	int i, j;
	char c;
	i = 0;
	j = len - 1;
	while(i < j){
		c = seq[i]; seq[i] = seq[j]; seq[j] = c;
		i ++; j --;
	}
	for(i=0;i<len;i++){
		switch(seq[i]){
			case 'a': seq[i] = 't'; break;
			case 'A': seq[i] = 'T'; break;
			case 'c': seq[i] = 'g'; break;
			case 'C': seq[i] = 'G'; break;
			case 'g': seq[i] = 'c'; break;
			case 'G': seq[i] = 'C'; break;
			case 't': seq[i] = 'a'; break;
			case 'T': seq[i] = 'A'; break;
		}
	}
}

#define reverse_dna_coord(x, y, tot_len) { x = x ^ y; y = x ^ y; x = x ^ y; x = tot_len - x; y = tot_len - y; }

#define old_bit2bits(bits, off, bit) { if(((off) & 0x1FU) == 0) (bits)[(off) >> 5] = 0; (bits)[(off) >> 5] |= ((u8i)(bit)) << (((~(off)) & 0x1FU) << 1); }

#define bit2bits(bits, off, bit) {	\
	u8i __off1;	\
	u4i __off2;	\
	__off1 = (off) >> 5;	\
	__off2 = (((~(off)) & 0x1FU) << 1);	\
	(bits)[__off1] = ((bits)[__off1] & (~(0x3LLU << __off2))) | (((u8i)(bit)) << __off2);	\
}

static inline void seq2bits(u8i *bits, u8i bitoff, char *seq, u4i seqlen){
	u8i i, c;
	for(i=0;i<seqlen;i++){
		c = base_bit_table[(int)seq[i]];
		if(c == 4) c = lrand48() & 0x03;
		bit2bits(bits, bitoff + i, c);
	}
}

static inline void revseq2bits(u8i *bits, u8i bitoff, char *seq, u4i seqlen){
	u8i i, c;
	for(i=0;i<seqlen;i++){
		c = base_bit_table[(int)seq[seqlen - i - 1]];
		if(c == 4) c = lrand48();
		c = (~c) & 0x03;
		bit2bits(bits, bitoff + i, c);
	}
}

static inline void bits2seq(char *seq, u8i *bits, u8i off, u4i len){
	u4i i, c;
	for(i=0;i<len;i++){
		c = bits2bit(bits, off + i);
		seq[i] = bit_base_table[c];
	}
	seq[i] = 0;
}

static inline void bits2revseq(char *seq, u8i *bits, u8i off, u4i len){
	u4i i, c;
	for(i=0;i<len;i++){
		c = (bits[(off + i)>>5] >> (((~(off + i)) & 0x1FU) << 1)) & 0x03;
		seq[len - i - 1] = bit_base_table[(~c)&0x03];
	}
	seq[i] = 0;
}

static inline u8i sub32seqbits(u8i *src, u8i off){
	u8i m;
	u4i n;
	m = off >> 5;
	n = (off & 0x1F) << 1;
	return (src[m] << n) | (((src[m + 1] >> (62 - n)) >> 2));
	//n = off & 0x1F;
	//if(n){
		//return (src[m] << (n << 1)) | (src[m + 1] >> ((32 - n) << 1));
	//} else {
		//return src[m];
	//}
}

static inline u8i sub8seqbits(u8i *src, u8i off){
	u8i off1;
	u4i off2;
	off1 = off >> 5;
	off2 = (off & 0x1FU) << 1;
	return ((src[off1] << off2) | (((src[off1 + 1] >> (62 - off2)) >> 2))) >> 48;
}

static inline u8i sub4seqbits(u8i *src, u8i off){
	u8i off1;
	u4i off2;
	off1 = off >> 5;
	off2 = (off & 0x1FU) << 1;
	return ((src[off1] << off2) | (((src[off1 + 1] >> (62 - off2)) >> 2))) >> 56;
}

static inline u8i sub2seqbits(u8i *src, u8i off){
	u8i off1;
	u4i off2;
	off1 = off >> 5;
	off2 = (off & 0x1FU) << 1;
	return ((src[off1] << off2) | (((src[off1 + 1] >> (62 - off2)) >> 2))) >> 60;
}

static inline u8i sub_seqbits(u8i *src, u8i off, u1i len){
	u8i off1;
	u4i off2;
	off1 = off >> 5;
	off2 = (off & 0x1FU) << 1;
	return ((src[off1] << off2) | (((src[off1 + 1] >> (62 - off2)) >> 2))) >> ((32 - len) << 1);
}
#define subseqbits(src, off, len) sub_seqbits(src, off, len)

static inline int cmpgt_seqbits(u8i *bits, u8i off1, u8i off2, u4i _len){
	u8i idxs[2], v[2];
	u4i offs[2], i, len;
	idxs[0] = off1 >> 5;
	idxs[1] = off2 >> 5;
	offs[0] = (off1 & 0x1FU) << 1;
	offs[1] = (off2 & 0x1FU) << 1;
	len = roundup_times(_len, 32);
	for(i=0;i<len;i+=32){
		v[0] = (bits[idxs[0]] << offs[0]) | (((bits[idxs[0] + 1] >> (62 - offs[0])) >> 2));
		v[1] = (bits[idxs[1]] << offs[1]) | (((bits[idxs[1] + 1] >> (62 - offs[1])) >> 2));
		if(v[0] > v[1]){
			return 1;
		} else if(v[0] < v[1]){
			return 0;
		}
		idxs[0] ++;
		idxs[1] ++;
	}
	return 0;
}

#if __BYTE_ORDER == 1234
static const u4i spare_2bits_table[256] = {
         0,  16777216,  33554432,  50331648,     65536,  16842752,  33619968,  50397184,
    131072,  16908288,  33685504,  50462720,    196608,  16973824,  33751040,  50528256,
       256,  16777472,  33554688,  50331904,     65792,  16843008,  33620224,  50397440,
    131328,  16908544,  33685760,  50462976,    196864,  16974080,  33751296,  50528512,
       512,  16777728,  33554944,  50332160,     66048,  16843264,  33620480,  50397696,
    131584,  16908800,  33686016,  50463232,    197120,  16974336,  33751552,  50528768,
       768,  16777984,  33555200,  50332416,     66304,  16843520,  33620736,  50397952,
    131840,  16909056,  33686272,  50463488,    197376,  16974592,  33751808,  50529024,
         1,  16777217,  33554433,  50331649,     65537,  16842753,  33619969,  50397185,
    131073,  16908289,  33685505,  50462721,    196609,  16973825,  33751041,  50528257,
       257,  16777473,  33554689,  50331905,     65793,  16843009,  33620225,  50397441,
    131329,  16908545,  33685761,  50462977,    196865,  16974081,  33751297,  50528513,
       513,  16777729,  33554945,  50332161,     66049,  16843265,  33620481,  50397697,
    131585,  16908801,  33686017,  50463233,    197121,  16974337,  33751553,  50528769,
       769,  16777985,  33555201,  50332417,     66305,  16843521,  33620737,  50397953,
    131841,  16909057,  33686273,  50463489,    197377,  16974593,  33751809,  50529025,
         2,  16777218,  33554434,  50331650,     65538,  16842754,  33619970,  50397186,
    131074,  16908290,  33685506,  50462722,    196610,  16973826,  33751042,  50528258,
       258,  16777474,  33554690,  50331906,     65794,  16843010,  33620226,  50397442,
    131330,  16908546,  33685762,  50462978,    196866,  16974082,  33751298,  50528514,
       514,  16777730,  33554946,  50332162,     66050,  16843266,  33620482,  50397698,
    131586,  16908802,  33686018,  50463234,    197122,  16974338,  33751554,  50528770,
       770,  16777986,  33555202,  50332418,     66306,  16843522,  33620738,  50397954,
    131842,  16909058,  33686274,  50463490,    197378,  16974594,  33751810,  50529026,
         3,  16777219,  33554435,  50331651,     65539,  16842755,  33619971,  50397187,
    131075,  16908291,  33685507,  50462723,    196611,  16973827,  33751043,  50528259,
       259,  16777475,  33554691,  50331907,     65795,  16843011,  33620227,  50397443,
    131331,  16908547,  33685763,  50462979,    196867,  16974083,  33751299,  50528515,
       515,  16777731,  33554947,  50332163,     66051,  16843267,  33620483,  50397699,
    131587,  16908803,  33686019,  50463235,    197123,  16974339,  33751555,  50528771,
       771,  16777987,  33555203,  50332419,     66307,  16843523,  33620739,  50397955,
    131843,  16909059,  33686275,  50463491,    197379,  16974595,  33751811,  50529027
};
#else
static const u4i spare_2bits_table[256] = {
         0,         1,         2,         3,       256,       257,       258,       259,
       512,       513,       514,       515,       768,       769,       770,       771,
     65536,     65537,     65538,     65539,     65792,     65793,     65794,     65795,
     66048,     66049,     66050,     66051,     66304,     66305,     66306,     66307,
    131072,    131073,    131074,    131075,    131328,    131329,    131330,    131331,
    131584,    131585,    131586,    131587,    131840,    131841,    131842,    131843,
    196608,    196609,    196610,    196611,    196864,    196865,    196866,    196867,
    197120,    197121,    197122,    197123,    197376,    197377,    197378,    197379,
  16777216,  16777217,  16777218,  16777219,  16777472,  16777473,  16777474,  16777475,
  16777728,  16777729,  16777730,  16777731,  16777984,  16777985,  16777986,  16777987,
  16842752,  16842753,  16842754,  16842755,  16843008,  16843009,  16843010,  16843011,
  16843264,  16843265,  16843266,  16843267,  16843520,  16843521,  16843522,  16843523,
  16908288,  16908289,  16908290,  16908291,  16908544,  16908545,  16908546,  16908547,
  16908800,  16908801,  16908802,  16908803,  16909056,  16909057,  16909058,  16909059,
  16973824,  16973825,  16973826,  16973827,  16974080,  16974081,  16974082,  16974083,
  16974336,  16974337,  16974338,  16974339,  16974592,  16974593,  16974594,  16974595,
  33554432,  33554433,  33554434,  33554435,  33554688,  33554689,  33554690,  33554691,
  33554944,  33554945,  33554946,  33554947,  33555200,  33555201,  33555202,  33555203,
  33619968,  33619969,  33619970,  33619971,  33620224,  33620225,  33620226,  33620227,
  33620480,  33620481,  33620482,  33620483,  33620736,  33620737,  33620738,  33620739,
  33685504,  33685505,  33685506,  33685507,  33685760,  33685761,  33685762,  33685763,
  33686016,  33686017,  33686018,  33686019,  33686272,  33686273,  33686274,  33686275,
  33751040,  33751041,  33751042,  33751043,  33751296,  33751297,  33751298,  33751299,
  33751552,  33751553,  33751554,  33751555,  33751808,  33751809,  33751810,  33751811,
  50331648,  50331649,  50331650,  50331651,  50331904,  50331905,  50331906,  50331907,
  50332160,  50332161,  50332162,  50332163,  50332416,  50332417,  50332418,  50332419,
  50397184,  50397185,  50397186,  50397187,  50397440,  50397441,  50397442,  50397443,
  50397696,  50397697,  50397698,  50397699,  50397952,  50397953,  50397954,  50397955,
  50462720,  50462721,  50462722,  50462723,  50462976,  50462977,  50462978,  50462979,
  50463232,  50463233,  50463234,  50463235,  50463488,  50463489,  50463490,  50463491,
  50528256,  50528257,  50528258,  50528259,  50528512,  50528513,  50528514,  50528515,
  50528768,  50528769,  50528770,  50528771,  50529024,  50529025,  50529026,  50529027
};
#endif

static inline void spare_2bits(u1i bs[32], u8i v){
	((u4i*)bs)[0] = spare_2bits_table[((v >> 56) & 0xFF)];
	((u4i*)bs)[1] = spare_2bits_table[((v >> 48) & 0xFF)];
	((u4i*)bs)[2] = spare_2bits_table[((v >> 40) & 0xFF)];
	((u4i*)bs)[3] = spare_2bits_table[((v >> 32) & 0xFF)];
	((u4i*)bs)[4] = spare_2bits_table[((v >> 24) & 0xFF)];
	((u4i*)bs)[5] = spare_2bits_table[((v >> 16) & 0xFF)];
	((u4i*)bs)[6] = spare_2bits_table[((v >>  8) & 0xFF)];
	((u4i*)bs)[7] = spare_2bits_table[((v >>  0) & 0xFF)];
}

typedef struct {
	u8i *bits;
	u8i size;
	u8i cap;
} BaseBank;

static inline size_t basebank_obj_desc_cnt(void *obj, int idx){ return ((((BaseBank*)obj)->size + 31) / 32 + 1) * 8; idx = idx; }

static inline void basebank_obj_desc_post_load(void *obj, size_t aux_data){
	BaseBank *bnk;
	UNUSED(aux_data);
	bnk = (BaseBank*)obj;
	bnk->cap = ((bnk->size + 31) / 32) * 32;
}

static const obj_desc_t basebank_obj_desc = {"BaseBank", sizeof(BaseBank), 1, {1}, {offsetof(BaseBank, bits)}, {(obj_desc_t*)&OBJ_DESC_DATA}, basebank_obj_desc_cnt, basebank_obj_desc_post_load};

static inline BaseBank* init_basebank(){
	BaseBank *bnk;
	bnk = malloc(sizeof(BaseBank));
	bnk->size = 0;
	bnk->cap  = 256;
	bnk->bits = calloc(bnk->cap / 32 + 1, 8);
	return bnk;
}

static inline void free_basebank(BaseBank *bnk){
	free(bnk->bits);
	free(bnk);
}

static inline void encap_basebank(BaseBank *bnk, u8i inc){
	u8i old;
	u8i *bits;
	if(bnk->cap - bnk->size >= inc) return;
	old = bnk->cap;
	if(MAX_U8 - inc <= bnk->size){
		fprintf(stderr, " -- Overflow(64bits) %llu + %llu, in %s -- %s:%d --\n", (u8i)bnk->size, (u8i)inc, __FUNCTION__, __FILE__, __LINE__);
		print_backtrace(stderr, 20);
		abort();
	}
	if(MAX_U8 - inc < 0x3FFFFFFFLLU){
		fprintf(stderr, " -- Overflow(64bits) %llu + %llu, in %s -- %s:%d --\n", (u8i)bnk->size, (u8i)inc, __FUNCTION__, __FILE__, __LINE__);
		print_backtrace(stderr, 20);
		abort();
	}
	if(bnk->size + inc <= 0x3FFFFFFFLLU){
		bnk->cap = roundup_times(2 * (bnk->size + inc), 32);
	} else {
		//bnk->cap = ((bnk->size + inc + 0xFFFFFFFLLU - 1LLU) / 0xFFFFFFFLLU) * 0xFFFFFFFLLU;
		bnk->cap = (bnk->size + inc + 0x3FFFFFFFLLU) & (MAX_U8 << 30);
	}
	if(bnk->cap < 32) bnk->cap = 32;
	bits = realloc(bnk->bits, ((bnk->cap >> 5) + 1) << 3);
	if(bits == NULL){
		fprintf(stderr, " -- Out of memory, try to allocate %llu bytes, old size %llu, in %s -- %s:%d --\n", (u8i)bnk->cap >> 2, old >> 2, __FUNCTION__, __FILE__, __LINE__);
		print_backtrace(stderr, 20);
		abort();
	}
	bnk->bits = bits;
	memset(bnk->bits + (old / 32), 0, (bnk->cap + 32 - old) / 4);
}

static inline void clear_basebank(BaseBank *bnk){
	//memset(bnk->bits, 0, ((bnk->size + 31) / 32) * 8);
	bnk->size = 0;
}

static inline void normalize_basebank(BaseBank *bnk){
	if(bnk->size < bnk->cap){
		if(bnk->size & 0x1FU){
			bnk->bits[bnk->size>>5] = bnk->bits[bnk->size>>5] & (MAX_U8 << (64 - ((bnk->size & 0x1FU) << 1)));
		}
	}
}

static inline void pack_basebank(BaseBank *bnk){
	u8i size;
	size = (bnk->size + 31) & (~0x1FLLU);
	if(size == 0) size = 32;
	if(size >= bnk->cap) return;
	bnk->cap = ((size + 31) / 32) * 32;
	bnk->bits = realloc(bnk->bits, ((bnk->cap >> 5) + 1) << 3);
	memset(bnk->bits + (bnk->cap >> 5), 0, 8);
}

static inline void bit2basebank(BaseBank *bnk, u1i v){
	encap_basebank(bnk, 1);
	bit2bits(bnk->bits, bnk->size, (v & 0x03));
	bnk->size ++;
}

static inline void bits2basebank(BaseBank *bnk, u8i *bits, u8i off, u8i len){
	u8i offset;
	encap_basebank(bnk, len);
	for(offset=off;offset<off+len;offset++){
		bit2bits(bnk->bits, bnk->size, bits2bit(bits, offset));
		bnk->size ++;
	}
}

#define fwdbits2basebank(bnk, bits, off, len) bits2basebank(bnk, bits, off, len)

static inline void fast_bits2basebank(BaseBank *bnk, u8i *bits, u8i off, u8i len){
	u8i end, dat;
	u4i gap;
	encap_basebank(bnk, len);
	if(len == 0) return;
	if(bnk->size & 0x1FU){
		gap = 32 - (bnk->size & 0x1FU);
		if(len <= gap){
			dat = subseqbits(bits, off, len);
			bnk->bits[bnk->size >> 5] |= dat << ((gap - len) << 1);
			bnk->size += len;
			return;
		} else {
			dat = subseqbits(bits, off, gap);
			bnk->bits[bnk->size >> 5] |= dat;
			bnk->size += gap;
			off += gap;
			len -= gap;
		}
	}
	end = off + len;
	for(;off+32<=end;off+=32){
		dat = sub32seqbits(bits, off);
		bnk->bits[bnk->size >> 5] = dat;
		bnk->size += 32;
	}
	if(off < end){
		dat = sub32seqbits(bits, off);
		bnk->bits[bnk->size >> 5] = dat & (MAX_U8 << ((32 - (end - off)) << 1));
		bnk->size += end - off;
	} else {
		bnk->bits[bnk->size >> 5] = 0;
	}
}

#define fast_fwdbits2basebank(bnk, bits, off, len) fast_bits2basebank(bnk, bits, off, len)

static inline void revbits2basebank(BaseBank *bnk, u8i *bits, u8i off, u8i len){
	u8i i;
	encap_basebank(bnk, len);
	for(i=1;i<=len;i++){
		bit2bits(bnk->bits, bnk->size, bits2revbit(bits, (off + len - i)));
		bnk->size ++;
	}
}

static inline void fast_revbits2basebank(BaseBank *bnk, u8i *bits, u8i off, u8i len){
	u8i end, dat;
	u4i gap;
	if(len == 0) return;
	encap_basebank(bnk, len);
	if(bnk->size & 0x1FU){
		gap = 32 - (bnk->size & 0x1FU);
		if(len <= gap){
			dat = subseqbits(bits, off, len);
			dat = dna_rev_seq(dat, len);
			bnk->bits[bnk->size >> 5] |= dat << ((gap - len) << 1);
			bnk->size += len;
			return;
		} else {
			dat = subseqbits(bits, off + len - gap, gap);
			dat = dna_rev_seq(dat, gap);
			bnk->bits[bnk->size >> 5] |= dat;
			bnk->size += gap;
			//off += gap;
			len -= gap;
		}
	}
	end = off + len;
	for(;off+32<=end;){
		end -= 32;
		dat = sub32seqbits(bits, end);
		dat = dna_rev_seq32(dat);
		bnk->bits[bnk->size >> 5] = dat;
		bnk->size += 32;
	}
	if(off < end){
		dat = sub32seqbits(bits, off);
		dat = dna_rev_seq32(dat);
		//bnk->bits[bnk->size >> 5] = dat & (MAX_U8 << ((32 - (end - off)) << 1));
		bnk->bits[bnk->size >> 5] = dat << ((32 - (end - off)) << 1);
		bnk->size += end - off;
	} else {
		bnk->bits[bnk->size >> 5] = 0;
	}
}


static inline void seq2basebank(BaseBank *bnk, char *seq, u8i len){
	u8i idx1, i, c;
	u1i idx2;
	encap_basebank(bnk, len);
	idx1 = bnk->size >> 5;
	idx2 = ((bnk->size) & 0x1FU) << 1;
	bnk->size += len;
	if(idx2 == 0) bnk->bits[idx1] = 0;
	for(i=0;i<len;i++){
		c = base_bit_table[(int)seq[i]] & 0x03;
		bnk->bits[idx1] |= c << (62 - idx2);
		idx2 = (idx2 + 2) & 0x3F;
		if(idx2 == 0){
			bnk->bits[++idx1] = 0;
		}
	}
}

#define fwdseq2basebank(bnk, seq, len) seq2basebank(bnk, seq, len)

static inline void revseq2basebank(BaseBank *bnk, char *seq, u8i len){
	char *p;
	u1i c;
	p = seq + len;
	encap_basebank(bnk, len);
	while(p > seq){
		p --;
		c = base_bit_table[(int)*p];
		c = (~c) & 0x03;
		bit2bits(bnk->bits, bnk->size, c);
		bnk->size ++;
	}
}

static inline void seq2basebank2(BaseBank *bnk, char *seq, u8i len){
	char *p;
	u1i c;
	p = seq;
	seq = seq + len;
	encap_basebank(bnk, len);
	while(p < seq){
		c = base_bit_table[(int)*p];
		if(c == 4) c = lrand48() & 0x03;
		bit2bits(bnk->bits, bnk->size, c);
		bnk->size ++;
		p ++;
	}
}

static inline void revseq2basebank2(BaseBank *bnk, char *seq, u8i len){
	char *p;
	u1i c;
	p = seq + len;
	encap_basebank(bnk, len);
	while(p > seq){
		p --;
		c = base_bit_table[(int)*p];
		if(c == 4) c = lrand48() & 0x03;
		c = (~c) & 0x03;
		bit2bits(bnk->bits, bnk->size, c);
		bnk->size ++;
	}
}

static inline u1i get_basebank(BaseBank *bnk, u8i off){ return bits2bit(bnk->bits, off); }

static inline void seq_basebank(BaseBank *bnk, u8i off, u8i len, char *seq){
	u8i i;
	for(i=0;i<len;i++){
		seq[i] = bit_base_table[bits2bit(bnk->bits, off + i)];
	}
	seq[i] = 0;
}

#define fwdseq_basebank(bnk, off, len, seq) seq_basebank(bnk, off, len, seq)

static inline void bitseq_basebank(BaseBank *bnk, u8i off, u8i len, u1i *seq){
	u8i i;
	for(i=0;i<len;i++){
		seq[i] = bits2bit(bnk->bits, off + i);
	}
}

static inline void revseq_basebank(BaseBank *bnk, u8i off, u8i len, char *seq){
	u8i i;
	for(i=0;i<len;i++){
		seq[i] = bit_base_table[(~bits2bit(bnk->bits, off + len - 1 - i)) & 0x03];
	}
	seq[i] = 0;
}

static inline void revbitseq_basebank(BaseBank *bnk, u8i off, u8i len, u1i *seq){
	u8i i;
	for(i=0;i<len;i++){
		seq[i] = (~bits2bit(bnk->bits, off + len - 1 - i)) & 0x03;
	}
}

static inline void reverse_basebank(BaseBank *bnk){
	u8i size, rsize;
	size = bnk->size;
	rsize = (bnk->size + 31) & (~0x1FLLU);
	encap_basebank(bnk, rsize + 32);
	memcpy(bnk->bits + (rsize >> 5), bnk->bits + 0, (rsize >> 5) << 3);
	bnk->size = 0;
	fast_revbits2basebank(bnk, bnk->bits, rsize, size);
}

static inline void print_seq_basebank(BaseBank *bnk, u8i off, u8i len, FILE *out){
	u8i i, b, e;
	char buf[101];
	for(b=off;b<off+len;){
		e = num_min(b + 100, off + len);
		for(i=b;i<e;i++){
			buf[i - b] = bit_base_table[bits2bit(bnk->bits, i)];
		}
		buf[e - b] = '\0';
		fputs(buf, out);
		//fputc('\n', out);
		b = e;
	}
}

static inline void print_lines_basebank(BaseBank *bnk, u8i off, u8i len, FILE *out, int linewidth){
	u8i i, b, e;
	char *buf;
	if(linewidth < 1) linewidth = 100;
	buf = malloc(linewidth + 1);
	for(b=off;b<off+len;){
		e = num_min(b + linewidth, off + len);
		for(i=b;i<e;i++){
			buf[i - b] = bit_base_table[bits2bit(bnk->bits, i)];
		}
		buf[e - b] = '\0';
		fputs(buf, out);
		fputc('\n', out);
		b = e;
	}
	free(buf);
}

#define print_fwdseq_basebank(bnk, off, len, out) print_seq_basebank(bnk, off, len, out)

static inline void println_seq_basebank(BaseBank *bnk, u8i off, u8i len, FILE *out){
	print_seq_basebank(bnk, off, len, out);
	fputc('\n', out);
}

#define println_fwdseq_basebank(bnk, off, len, out) println_seq_basebank(bnk, off, len, out)

static inline void print_revseq_basebank(BaseBank *bnk, u8i off, u8i len, FILE *out){
	u8i i;
	char buf[65];
	buf[64] = '\0';
	for(i=0;i<len;){
		buf[i & 0x3F] = bit_base_table[bits2revbit(bnk->bits, off + len - 1 - i)];
		i ++;
		if((i & 0x3F) == 0){
			fprintf(out, "%s", buf);
		}
	}
	if(i & 0x3F){
		buf[i & 0x3F] = '\0';
		fprintf(out, "%s", buf);
	}
}

static inline u8i sub32_basebank(BaseBank *bnk, u8i off){ return sub32seqbits(bnk->bits, off); }

static inline u8i sub4_basebank(BaseBank *bnk, u8i off){ return sub4seqbits(bnk->bits, off); }

// assert(len > 0 && len <= 32)
static inline u8i subbits_basebank(BaseBank *bnk, u8i off, u1i len){ return sub_seqbits(bnk->bits, off, len); }

static inline void println_revseq_basebank(BaseBank *bnk, u8i off, u8i len, FILE *out){
	print_revseq_basebank(bnk, off, len, out);
	fputc('\n', out);
}

static inline u8i hzsubbits_basebank(BaseBank *bnk, u8i off, u1i len){
	u8i k;
	u1i i, b, c;
	k = 0;
	b = 4;
	for(i=0;i<len;off++){
		c = bits2bit(bnk->bits, off);
		if(c == b) continue;
		i ++;
		b = c;
		k = (k << 2) | b;
	}
	return k;
}

static inline int bitsearch_basebank(BaseBank *bnk, u8i *_off, u8i len, u8i bits, u1i size, int max_occ){
	u8i off, end, k, mask;
	u1i b;
	int ret;
	off = *_off;
	end = off + len;
	mask = MAX_U8 >> ((32 - size) << 1);
	k = subbits_basebank(bnk, off, size - 1);
	off += size - 1;
	ret = 0;
	for(;off<end;off++){
		b = bits2bit(bnk->bits, off);
		k = ((k << 2) | b) & mask;
		if(k == bits){
			_off[ret++] = off - (size - 1);
			if(ret >= max_occ) break;
		}
	}
	return ret;
}

static inline int hzbitsearch_basebank(BaseBank *bnk, u8i *_off, u8i len, u8i bits, u1i size, int max_occ){
	u8i off, h, end, k, mask;
	u1i b, c;
	int ret;
	off = *_off;
	end = off + len;
	mask = MAX_U8 >> ((32 - size) << 1);
	k = 0;
	h = 0;
	b = 4;
	ret = 0;
	for(;off<end;off++){
		c = bits2bit(bnk->bits, off);
		if(c == b) continue;
		b = c;
		h ++;
		k = ((k << 2) | b) & mask;
		if(h >= size && k == bits){
			_off[ret++] = off - (size - 1);
			if(ret >= max_occ) break;
		}
	}
	return ret;
}

static inline u4i mismatch_basebank(BaseBank *bnk, u8i off1, u8i off2, u4i len){
	u8i seq1, seq2;
	u4i mm, i;
	mm = 0;
	for(i=0;i+32<=len;i+=32){
		seq1 = sub32seqbits(bnk->bits, off1 + i);
		seq2 = sub32seqbits(bnk->bits, off2 + i);
		mm += count_ones_bit64(dna_xor2ones(seq1 ^ seq2));
	}
	if(i < len){
		seq1 = sub32seqbits(bnk->bits, off1 + i);
		seq2 = sub32seqbits(bnk->bits, off2 + i);
		mm += count_ones_bit64((dna_xor2ones(seq1 ^ seq2)) >> ((32 - (len - i)) << 1));
	}
	return mm;
}

thread_beg_def(_mradix);
BaseBank *bb;
u4i *counts[2];
u4i *offs;
u1v *lcps;
u4i size, klen;
int task;
FILE *log;
thread_end_def(_mradix);

thread_beg_func(_mradix);
BaseBank *bb;
u4i *offs, *counts[2];
u4i i, j, size, klen, m, n, v, t;
u4i ncpu, tidx;
bb = _mradix->bb;
counts[0] = calloc((MAX_U2 + 1), sizeof(u4i)); // used in twice
counts[1] = calloc((MAX_U2 + 1), sizeof(u4i));
ncpu = _mradix->n_cpu;
tidx = _mradix->t_idx;
thread_beg_loop(_mradix);
if(_mradix->task == 1){
	size = _mradix->size;
	for(i=_mradix->t_idx;i<size;i+=_mradix->n_cpu){
		v = sub8seqbits(bb->bits, i);
		counts[1][v] ++;
	}
	_mradix->counts[0] = counts[0];
	_mradix->counts[1] = counts[1];
} else if(_mradix->task == 11){
	offs = _mradix->offs;
	size = _mradix->size;
	for(i=0;i<size;i++){
		v = sub8seqbits(bb->bits, i);
		if((v % _mradix->n_cpu) == (u4i)_mradix->t_idx){
			offs[_mradix->counts[0][v]++] = i;
		}
		if(_mradix->t_idx == 0 && _mradix->log && (i % 1000000) == 0){
			fprintf(_mradix->log, "\r%u", i); fflush(_mradix->log);
		}
	}
	if(_mradix->t_idx == 0 && _mradix->log){
		fprintf(_mradix->log, "\r%u\n", size);
	}
} else if(_mradix->task == 2) {
	offs = _mradix->offs;
	size = _mradix->size;
	klen = _mradix->klen - 8;
	if(size <= MAX_U1){
		sort_array(offs, size, u4i, cmpgt_seqbits(bb->bits, a + 8, b + 8, klen)); // 8 bp already sorted
	} else {
		memset(counts[1], 0, (MAX_U1 + 1) * sizeof(u4i));
		for(i=0;i<size;i++){
			v = sub4seqbits(bb->bits, offs[i] + 8);
			counts[1][v] ++;
		}
		m = 0;
		for(i=0;i<=MAX_U1;i++){
			counts[0][i] = m;
			m += counts[1][i];
			counts[1][i] = m;
		}
		for(m=0;m<=MAX_U1;m++){
			while(counts[0][m] < counts[1][m]){
				v = offs[counts[0][m]];
				n = sub4seqbits(bb->bits, v + 8);
				while(n > m){
					t = offs[counts[0][n]];
					offs[counts[0][n]] = v;
					counts[0][n] ++;
					v = t;
					n = sub4seqbits(bb->bits, v + 8);
				}
				offs[counts[0][m]++] = v;
			}
		}
		n = 0;
		klen -= 4;
		for(m=0;m<=MAX_U1;m++){
			if(counts[0][m] - n < 2){
				// nothing to do
			} else {
				sort_array(offs + n, counts[0][m] - n, u4i, cmpgt_seqbits(bb->bits, a + 8 + 4, b + 8 + 4, klen));
			}
			n = counts[0][m];
		}
	}
} else if(_mradix->task == 3){
	u1v *lcps;
	u8i mask;
	u4i beg, end;
	u1i lcp;
	offs = _mradix->offs;
	lcps = _mradix->lcps;
	beg = ((_mradix->size + ncpu - 1) / ncpu) * tidx;
	end = ((_mradix->size + ncpu - 1) / ncpu) * (tidx + 1);
	if(end > _mradix->size) end = _mradix->size;
	klen = _mradix->klen;
	if(beg == 0){
		lcps->buffer[0] = 0;
		beg = 1;
	}
	for(i=beg;i<end;i++){
		lcp = 0;
		for(j=0;j<klen;j+=32){
			mask = sub32seqbits(bb->bits, offs[i - 1] + j) ^ sub32seqbits(bb->bits, offs[i] + j);
			if(mask == 0){
				lcp += 32;
			} else {
				lcp += __builtin_clzll(mask) >> 1;
				break;
			}
		}
		lcps->buffer[i] = lcp;
	}
}
thread_end_loop(_mradix);
free(counts[0]);
free(counts[1]);
thread_end_func(_mradix);

// bullet size = 4^8 = (MAX_U2 + 1)
static inline void msd_radix_sort_u4_basebank(BaseBank *bb, u4v *offs, u1v *lcps, u1i _klen, u4i ncpu, FILE *log){
	u4i *counts[3]; // off, end, off
	u4i klen, i, j, size, m, n;
	thread_preprocess(_mradix);
	klen = roundup_times(_klen, 32);
	size = num_max(bb->size, klen) - klen;
	clear_u4v(offs);
	encap_u4v(offs, size);
	offs->size = size;
	counts[0] = calloc(MAX_U2 + 1, sizeof(u4i));
	counts[1] = calloc(MAX_U2 + 1, sizeof(u4i));
	if(log) fprintf(log, "[%s] msd_radix_sort length=%u depth=%u\n", date(), (u4i)bb->size, klen);
	thread_beg_init(_mradix, ncpu);
	_mradix->bb = bb;
	_mradix->counts[0] = NULL;
	_mradix->counts[1] = NULL;
	_mradix->size = size;
	_mradix->offs = NULL;
	_mradix->klen = klen;
	_mradix->task = 0;
	_mradix->log = log;
	thread_end_init(_mradix);
	thread_apply_all(_mradix, _mradix->task = 1);
	thread_beg_iter(_mradix);
	for(j=0;j<=MAX_U2;j++){
		counts[1][j] += _mradix->counts[1][j];
	}
	thread_end_iter(_mradix);
	m = 0;
	for(i=0;i<=MAX_U2;i++){
		counts[0][i] = m;
		m += counts[1][i];
		counts[1][i] = m;
	}
	thread_beg_iter(_mradix);
	_mradix->offs = offs->buffer;
	_mradix->counts[0] = counts[0];
	_mradix->counts[1] = counts[1];
	_mradix->task = 11;
	thread_wake(_mradix);
	thread_end_iter(_mradix);
	thread_wait_all(_mradix);
	if(log) fprintf(log, "[%s] msd_radix_sort sorted by first 8 bp\n", date());
	n = 0;
	for(m=0;m<=MAX_U2;m++){
		if(log && (m % 100) == 0){
			fprintf(log, "\r%u", counts[0][m]); fflush(log);
		}
		if(counts[0][m] - n < 2){
			// nothing to do
		} else {
			thread_wait_one(_mradix);
			_mradix->offs = offs->buffer + n;
			_mradix->size = counts[0][m] - n;
			_mradix->task = 2;
			thread_wake(_mradix);
		}
		n = counts[0][m];
	}
	thread_wait_all(_mradix);
	if(log) fprintf(log, "\r%u\n", n);
	if(log) fprintf(log, "[%s] msd_radix_sort sorted %u bases\n", date(), klen - 8);
	free(counts[0]);
	free(counts[1]);
	if(lcps){
		clear_u1v(lcps);
		encap_u1v(lcps, size);
		lcps->size = size;
		thread_beg_iter(_mradix);
		_mradix->offs = offs->buffer;
		_mradix->size = size;
		_mradix->lcps = lcps;
		_mradix->task = 3;
		thread_wake(_mradix);
		thread_end_iter(_mradix);
		thread_wait_all(_mradix);
		if(log) fprintf(log, "[%s] msd_radix_sort calculated LCP\n", date());
	}
	thread_beg_close(_mradix);
	thread_end_close(_mradix);
}

/*
 * Sequence DB
 */

typedef struct {
	u4i      nseq;
	BaseBank *rdseqs;
	cplist   *rdtags;
	u8v      *rdoffs;
	u4v      *rdlens;
	cuhash   *rdhash;
} SeqBank;
static inline void rebuild_rdhash_seqbank(void *sb, size_t aux);
static const obj_desc_t seqbank_obj_desc = {"SeqBank", sizeof(SeqBank), 5, {1, 1, 1, 1, 1}, 
	{offsetof(SeqBank, rdseqs), offsetof(SeqBank, rdtags), offsetof(SeqBank, rdoffs), offsetof(SeqBank, rdlens), offsetof(SeqBank, rdhash)},
	{&basebank_obj_desc, &cplist_deep_obj_desc, &u8v_obj_desc, &u4v_obj_desc, &cuhash_obj_desc},
	NULL, rebuild_rdhash_seqbank};

static inline SeqBank* init_seqbank(){
	SeqBank *sb;
	sb = malloc(sizeof(SeqBank));
	sb->nseq   = 0;
	sb->rdseqs = init_basebank();
	sb->rdtags = init_cplist(16);
	sb->rdoffs = init_u8v(16);
	sb->rdlens = init_u4v(16);
	sb->rdhash = init_cuhash(1023);
	return sb;
}

static inline void free_seqbank(SeqBank *sb){
	u4i i;
	for(i=0;i<sb->rdtags->size;i++) if(sb->rdtags->buffer[i]) free(sb->rdtags->buffer[i]);
	free_basebank(sb->rdseqs);
	free_cplist(sb->rdtags);
	free_u8v(sb->rdoffs);
	free_u4v(sb->rdlens);
	free_cuhash(sb->rdhash);
	free(sb);
}

static inline void clear_seqbank(SeqBank *sb){
	u4i i;
	for(i=0;i<sb->rdtags->size;i++) if(sb->rdtags->buffer[i]) free(sb->rdtags->buffer[i]);
	clear_basebank(sb->rdseqs);
	clear_cplist(sb->rdtags);
	clear_u8v(sb->rdoffs);
	clear_u4v(sb->rdlens);
	clear_cuhash(sb->rdhash);
	sb->nseq = 0;
}


// SeqBank's rdhash is wrongly loaded, need to be corrected
static inline void rebuild_rdhash_seqbank(void *_sb, size_t aux){
	SeqBank *sb;
	u4i i;
	UNUSED(aux);
	sb = (SeqBank*)_sb;
	clear_cuhash(sb->rdhash); // hash size is not changed, thus there won't have hash re-size
	for(i=0;i<sb->rdtags->size;i++){
		put_cuhash(sb->rdhash, (cuhash_t){get_cplist(sb->rdtags, i), i});
	}
}


static inline void push_seqbank(SeqBank *sb, char *tag, int tag_len, char *seq, int seq_len){
	char *ptr;
	if(tag && tag_len){
		ptr = malloc(tag_len + 1);
		memcpy(ptr, tag, tag_len);
		ptr[tag_len] = 0;
	} else {
		ptr = NULL;
	}
	push_cplist(sb->rdtags, ptr);
	push_u8v(sb->rdoffs, sb->rdseqs->size);
	seq2basebank(sb->rdseqs, seq, seq_len);
	push_u4v(sb->rdlens, seq_len);
	if(ptr) put_cuhash(sb->rdhash, (cuhash_t){ptr, sb->nseq});
	sb->nseq ++;
}

static inline void fwdbitpush_seqbank(SeqBank *sb, char *tag, int tag_len, u8i *bits, u8i off, u4i len){
	char *ptr;
	if(tag && tag_len){
		ptr = malloc(tag_len + 1);
		memcpy(ptr, tag, tag_len);
		ptr[tag_len] = 0;
	} else {
		ptr = NULL;
	}
	push_cplist(sb->rdtags, ptr);
	push_u8v(sb->rdoffs, sb->rdseqs->size);
	fast_fwdbits2basebank(sb->rdseqs, bits, off, len);
	push_u4v(sb->rdlens, len);
	if(ptr) put_cuhash(sb->rdhash, (cuhash_t){ptr, sb->nseq});
	sb->nseq ++;
}

static inline void revbitpush_seqbank(SeqBank *sb, char *tag, int tag_len, u8i *bits, u8i off, u4i len){
	char *ptr;
	if(tag && tag_len){
		ptr = malloc(tag_len + 1);
		memcpy(ptr, tag, tag_len);
		ptr[tag_len] = 0;
	} else {
		ptr = NULL;
	}
	push_cplist(sb->rdtags, ptr);
	push_u8v(sb->rdoffs, sb->rdseqs->size);
	fast_revbits2basebank(sb->rdseqs, bits, off, len);
	push_u4v(sb->rdlens, len);
	if(ptr) put_cuhash(sb->rdhash, (cuhash_t){ptr, sb->nseq});
	sb->nseq ++;
}

static inline u4i find_seqbank(SeqBank *sb, char *tag){ cuhash_t *e; if((e = get_cuhash(sb->rdhash, tag))) return e->val; else return MAX_U4; }

static inline u4i off2idx_seqbank(SeqBank *sb, u8i off){
	u4i ret;
	bsearch_array(sb->rdoffs->buffer, sb->rdoffs->size, u8i, ret, a < off);
	return ret? ret - 1 : 0;
}

static inline u4i num_n50(u4v *lens, FILE *out){
	u8i tot, cum;
	u4i i, max, min, n50, l50, n90, l90, avg;
	if(lens->size == 0) return 0;
	sort_array(lens->buffer, lens->size, u4i, num_cmpgt(b, a));
	tot = 0;
	max = lens->buffer[0];
	min = lens->buffer[lens->size - 1];
	for(i=0;i<lens->size;i++){
		tot += lens->buffer[i];
	}
	avg = (tot + lens->size - 1) / lens->size;
	cum = 0;
	i = 0;
	while(i < lens->size){
		cum += lens->buffer[i];
		if((b8i)cum >= tot * 0.5) break;
		i ++;
	}
	n50 = i < lens->size? lens->buffer[i] : min;
	l50 = i < lens->size? i + 1 : i;
	i ++;
	while(i < lens->size){
		cum += lens->buffer[i];
		if((b8i)cum >= tot * 0.9) break;
		i ++;
	}
	n90 = i < lens->size? lens->buffer[i] : min;
	l90 = i < lens->size? i + 1 : i;
	if(out){
		fprintf(out, "TOT %llu, CNT %u, AVG %u, MAX %u, N50 %u, L50 %u, N90 %u, L90 %u, Min %u", tot, (u4i)lens->size, avg, max, n50, l50, n90, l90, min);
		fflush(out);
	}
	return n50;
}

#endif
