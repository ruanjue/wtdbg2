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

#include "dmo.h"
#include "file_reader.h"


wt_ovl_t find_ttr(u1v *seqs, DMOPar *par, hzmrv *kmers, DMOAux *aux){
	u8i kmer, kmask;
	u4i i, j, b, c;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - par->ksize) << 1);
	b = 4;
	kmer = 0;
	clear_hzmrv(kmers);
	for(i=j=0;i<seqs->size;i++){
		if(par->hk && seqs->buffer[i] == b) continue;
		b = seqs->buffer[i];
		kmer = ((kmer << 2) | b) & kmask;
		j ++;
		if(j < par->ksize) continue;
		push_hzmrv(kmers, (hzmr_t){kmer, 0, i, par->ksize, 0});
	}
	sort_array(kmers->buffer, kmers->size, hzmr_t, num_cmpgt(a.mer, b.mer));
	clear_hzmpv(aux->rs);
	for(i=j=0;i<=kmers->size;i++){
		if(i<kmers->size && kmers->buffer[i].mer == kmers->buffer[j].mer){
			continue;
		}
		for(b=j;b+1<i;b++){
			for(c=b+1;c<i;c++){
				push_hzmpv(aux->rs, (hzmp_t){0, kmers->buffer[b].off, 0, kmers->buffer[c].off, par->ksize, par->ksize, 0});
			}
		}
		j = i;
	}
	return dot_matrix_align_hzmps(par, 0, seqs->size, 1, seqs->size, NULL, aux->rs, aux->dst, aux->regs, aux->diags, aux->block, aux->grps, aux->mem);
}

int main(){
	FileReader *fr;
	Sequence *seq;
	u1v *seqs;
	hzmrv *kmers;
	DMOPar *par;
	DMOAux *aux;
	wt_ovl_t hit;
	u4i i;
	par = init_dmopar(15, 0, 1, 1, 256, 96, 3, 600, 1.0, 1.0);
	par->min_mat = 400;
	par->min_aln = 4000;
	aux = init_dmoaux();
	seqs = init_u1v(1024);
	kmers = init_hzmrv(1024);
	fr = stdin_filereader();
	seq = NULL;
	while(fread_seq(&seq, fr)){
		clear_u1v(seqs);
		for(i=0;i<(u4i)seq->seq.size;i++) push_u1v(seqs, base_bit_table[(int)seq->seq.string[i]]);
		hit = find_ttr(seqs, par, kmers, aux);
		if(hit.mat >= par->min_mat && num_min(hit.qe - hit.qb, hit.te - hit.tb) >= par->min_aln && hit.aln * par->min_sm <= hit.mat
			&& num_diff(hit.qe - hit.qb, hit.te - hit.tb) <= par->aln_var * num_min(hit.qe - hit.qb, hit.te - hit.tb)){
			fprintf(stdout, "%s\t+\t%d\t%d\t%d\t", seq->tag.string, seq->seq.size, hit.tb, hit.te);
			fprintf(stdout, "%s\t+\t%d\t%d\t%d\t", seq->tag.string, seq->seq.size, hit.qb, hit.qe);
			fprintf(stdout, "%d\t%0.3f\t%d\t%d\t%d\t%d\t*\n", hit.score, 1.0 * hit.mat / hit.aln, hit.mat, 0, 0, 0);
		}
	}
	fclose_filereader(fr);
	free_u1v(seqs);
	free_dmopar(par);
	free_dmoaux(aux);
	free_hzmrv(kmers);
	return 0;
}
