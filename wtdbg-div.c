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

#include "dna.h"
#include "list.h"
#include "bitvec.h"
#include "upgma.h"
#include "file_reader.h"

#define WTDIV_MAX_KSIZE	16

typedef struct {
	uint32_t kmer;
	uint32_t seqid:31, closed:1;
} seq_kmer_t;
define_list(seqkmerv, seq_kmer_t);

typedef struct {
	uint32_t off, len;
} frag_t;
define_list(fragv, frag_t);

#define signsof8bytes(v) (((v) * 567382630219905LLU) >> 56)

typedef struct {
	uint32_t nseq;
	uint64_t nk, nx;
	uint32_t ksize, kmask;
	float    min_kc;
	BitVec   *kbits;
	seqkmerv *sks;
	fragv    *frags;
	u4v      *kmers;
	u2v      *kcnts;
	u2v      *kcoes;
	UPGMA    *upgma;
} WTDIV;

WTDIV* init_wtdiv(uint32_t ksize, float min_kmer_cnt){
	WTDIV *wt;
	wt = malloc(sizeof(WTDIV));
	wt->nseq = 0;
	wt->nk = 0;
	wt->nx = 64;
	wt->ksize = ksize;
	wt->kmask = 0xFFFFFFFFU >> ((16 - ksize) << 1);
	wt->min_kc = min_kmer_cnt;
	wt->kbits = NULL;
	wt->sks = init_seqkmerv(32);
	wt->frags = init_fragv(32);
	wt->kmers = init_u4v(32);
	wt->kcnts = init_u2v(32);
	wt->kcoes = init_u2v(32);
	wt->upgma = init_upgma();
	return wt;
}

void clear_wtdiv(WTDIV *wt){
	wt->nseq = 0;
	wt->nk = 0;
	wt->nx = 64;
	if(wt->kbits) free_bitvec(wt->kbits);
	wt->kbits = NULL;
	clear_seqkmerv(wt->sks);
	clear_fragv(wt->frags);
	clear_u4v(wt->kmers);
	clear_u2v(wt->kcnts);
	clear_u2v(wt->kcoes);
}

void free_wtdiv(WTDIV *wt){
	if(wt->kbits) free_bitvec(wt->kbits);
	free_seqkmerv(wt->sks);
	free_fragv(wt->frags);
	free_u4v(wt->kmers);
	free_u2v(wt->kcnts);
	free_u2v(wt->kcoes);
	free_upgma(wt->upgma);
	free(wt);
}

void push_wtdiv(WTDIV *wt, char *seq, int len){
	uint32_t kmer;
	int i;
	beg_seq2kmers(seq, len, wt->ksize, wt->kmask, kmer, i);
	push_seqkmerv(wt->sks, (seq_kmer_t){kmer, wt->nseq, 0});
	end_seq2kmers;
	wt->nseq ++;
}

void ready_wtdiv(WTDIV *wt){
	frag_t *f;
	uint32_t i, j, kc, nk, discard;
	int kidx;
	if(wt->sks->size == 0) return;
	sort_array(wt->sks->buffer, wt->sks->size, seq_kmer_t, num_cmpgtx(a.kmer, b.kmer, a.seqid, b.seqid));
	encap_seqkmerv(wt->sks, 1);
	wt->sks->buffer[wt->sks->size].kmer = 0xFFFFFFFFU;
	clear_fragv(wt->frags);
	j = 0;
	nk = 0;
	kc = 1;
	discard = 0;
	for(i=1;i<=wt->sks->size;i++){
		if(i == wt->sks->size || wt->sks->buffer[i].kmer != wt->sks->buffer[i-1].kmer){
			if(kc < 2 || kc < wt->min_kc * wt->nseq){
				while(j < i){ wt->sks->buffer[j++].closed = 1; }
			} else {
				nk ++;
				push_fragv(wt->frags, (frag_t){j, i - j});
			}
			j = i;
			kc = 1;
			discard = 0;
		} else if(!discard){
			if(wt->sks->buffer[i].seqid == wt->sks->buffer[i-1].seqid){ discard = 1; kc = 0; }
			else kc ++;
		}
	}
	wt->nk = nk;
	wt->nx = (nk + 63) & 0xFFFFFFC0U;
	wt->kbits = init_bitvec(wt->nx * wt->nseq);
	sort_array(wt->frags->buffer, wt->frags->size, frag_t, num_cmpgt(b.len, a.len));
	clear_u4v(wt->kmers);
	clear_u2v(wt->kcnts);
	kidx = -1;
	for(i=0;i<wt->frags->size;i++){
		f = ref_fragv(wt->frags, i);
		kidx ++;
		push_u4v(wt->kmers, wt->sks->buffer[f->off].kmer);
		push_u2v(wt->kcnts, f->len);
		for(j=0;j<f->len;j++){
			one_bitvec(wt->kbits, wt->nx * wt->sks->buffer[f->off + j].seqid + kidx);
		}
	}
}

void gen_co_existence_kmer_matrix_wtdiv(WTDIV *wt){
	uint32_t i, j, idx, cnt;
	clear_and_encap_u2v(wt->kcoes, wt->nk * wt->nk);
	wt->kcoes->size = wt->nk * wt->nk;
	for(i=0;i<wt->nk;i++){
		wt->kcoes->buffer[i * wt->nk + i] = wt->kcnts->buffer[i];
		for(j=i+1;j<wt->nk;j++){
			cnt = 0;
			for(idx=0;idx<wt->nseq;idx++){
				if(get_bitvec(wt->kbits, idx * wt->nx + i) && get_bitvec(wt->kbits, idx * wt->nx + j)) cnt ++;
			}
			wt->kcoes->buffer[i * wt->nk + j] = cnt;
			wt->kcoes->buffer[j * wt->nk + i] = cnt;
		}
	}
}

void print_seq_kmer_matrix_wtdiv(WTDIV *wt, FILE *out){
	uint32_t i, j, off;
	for(i=0;i<wt->nseq;i++){
		off = i * wt->nx;
		for(j=0;j<wt->nk;j++){
			if(j) fputc(' ', out);
			fputc(get_bitvec(wt->kbits, off + j) + '0', out);
		}
		fputc('\n', out);
	}
}

void print_kmer_kmer_jaccard_index_matrix_wtdiv(WTDIV *wt, FILE *out){
	uint32_t i, j, x, y, z;
	for(i=0;i<wt->nk;i+=50){
		for(j=0;j<wt->nk;j++){
			if(j) fputc('\t', out);
			x = wt->kcnts->buffer[i];
			y = wt->kcnts->buffer[j];
			z = wt->kcoes->buffer[i * wt->nk + j];
			fprintf(out, "%0.4f", 1.0 * (z * z) / ((x + y - z) * wt->nseq));
		}
		fputc('\n', out);
	}
}

static inline uint32_t cal_seq_matched_kmers(WTDIV *wt, uint32_t a, uint32_t b){
	uint64_t *A, *B, C;
	uint32_t i, n, m;
	n = wt->nx >> 6;
	A = wt->kbits->bits + ((wt->nx * a) >> 6);
	B = wt->kbits->bits + ((wt->nx * b) >> 6);
	m = 0;
	for(i=0;i<n;i++){
		C = A[i] & B[i];
		m += count_ones_bit64(C);
	}
	//printf("[%3d, %3d] %d\n", a, b, m);
	return m;
}

void upgma_wtdiv(WTDIV *wt){
	uint32_t i, j;
	reset_upgma(wt->upgma, wt->nseq);
	for(i=0;i<wt->nseq;i++){
		for(j=i+1;j<wt->nseq;j++){
			set_upgma(wt->upgma, i, j, cal_seq_matched_kmers(wt, i, j));
		}
	}
	do_upgma(wt->upgma);
}

/*
void beg_wtdiv(WTDIV *wt){
	u4v *cnts;
	uint32_t iter_idx, idx, cnt, best_idx, best_cnt;
	int flag;
	cnts = init_u4v(wt->upgma->tre->size);
	iter_idx = 0xFFFFFFFFU; flag = 1; cnt = 0;
	while((flag = iter_upgma(wt->upgma, &iter_idx, flag))){
		switch(flag){
			case 2: cnts->buffer[iter_idx] = 1; break;
			case 4:
				idx = wt->upgma->tre->buffer[iter_idx].child;
				while(idx != 0xFFFFFFFFU){
					cnt += cnts->buffer[idx];
					idx = wt->upgma->tre->buffer[idx].sibling;
				}
				cnts->buffer[iter_idx] = cnt;
				break;
		}
	}
	cnt = wt->nseq / 2;
	for(idx=0;idx<wt->upgma->tre->size;idx++){
	}
}
*/

int usage(){
	printf(
	"No Usage\n"
	);
	return 1;
}

int main(int argc, char **argv){
	WTDIV *wt;
	FileReader *fr;
	Sequence *seq;
	int ksize;
	float min_kc;
	ksize = 11;
	min_kc = 0.02;
	fr = stdin_filereader();
	wt = init_wtdiv(ksize, min_kc);
	seq = NULL;
	while(fread_seq(&seq, fr)){
		push_wtdiv(wt, seq->seq.string, seq->seq.size);
	}
	fclose_filereader(fr);
	ready_wtdiv(wt);
	gen_co_existence_kmer_matrix_wtdiv(wt);
	print_kmer_kmer_jaccard_index_matrix_wtdiv(wt, stdout);
	//print_seq_kmer_matrix_wtdiv(wt, stdout);
	//upgma_wtdiv(wt);
	//print_upgma(wt->upgma, stdout);
	free_wtdiv(wt);
	return 0;
}

