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
#include "timer.h"
#include <getopt.h>

thread_beg_def(mdmo);
DMO *wt;
uint32_t pbid;
wtovlv *hits;
kigarv *kigars;
DMOPar *par2;
int align_mode;
BaseBank *rdseqs;
pbread_t rd;
thread_end_def(mdmo);

thread_beg_func(mdmo);
DMO *wt;
pbread_t *rd;
wtovlv *hits, *hs;
wt_ovl_t *hit;
hzmhash *hash;
hzmrv *seeds;
DMOAux *aux;
uint32_t i;
wt = mdmo->wt;
hits = mdmo->hits;
hs = init_wtovlv(32);
hash = init_hzmhash(1023);
seeds = init_hzmrv(32);
aux = init_dmoaux();

thread_beg_loop(mdmo);
if(hzm_debug){
	fprintf(hzm_debug_out, "=%s\n", wt->reads->buffer[mdmo->pbid].tag);
}
rd = (pbread_t*)&mdmo->rd;
if(mdmo->par2){
	clear_wtovlv(hs);
	query_dmo(wt, mdmo->pbid, mdmo->rdseqs, rd->rdoff, rd->rdlen, mdmo->align_mode, hs, NULL, aux);
	if(hs->size == 0) continue;
	index_reg_dmo(mdmo->rdseqs, rd->rdoff, rd->rdlen, mdmo->par2, hash, seeds);
	for(i=0;i<hs->size;i++){
		hit = ref_wtovlv(hs, i);
		query_reg_dmo(wt->rdseqs, hit->pb2, wt->reads->buffer[hit->pb2].rdoff, wt->reads->buffer[hit->pb2].rdlen, mdmo->pbid, rd->rdlen, hash, seeds, hits, mdmo->kigars, mdmo->par2, aux);
	}
} else {
	query_dmo(wt, mdmo->pbid, mdmo->rdseqs, rd->rdoff, rd->rdlen, mdmo->align_mode, hits, mdmo->kigars, aux);
}
thread_end_loop(mdmo);

free_dmoaux(aux);
free_hzmhash(hash);
free_hzmrv(seeds);
free_wtovlv(hs);
thread_end_func(mdmo);

uint64_t overlap_dmo(DMO *wt, DMOPar *par2, int online, int align_mode, int has_kigar, int ncpu, uint32_t n_job, uint32_t i_job, FILE *out){
	FileReader *fr;
	Sequence *seq;
	uint64_t ret;
	uint32_t i;
	int n_cpu;
	def_clock(secs);
	thread_preprocess(mdmo);
	if(hzm_debug) n_cpu = 1;
	else n_cpu = ncpu;
	beg_clock(secs);
	thread_beg_init(mdmo, n_cpu);
	mdmo->wt = wt;
	mdmo->rdseqs = wt->rdseqs;
	mdmo->pbid = 0xFFFFFFFFU;
	memset((void*)&mdmo->rd, 0, sizeof(pbread_t));
	mdmo->hits    = init_wtovlv(64);
	mdmo->kigars = has_kigar? init_kigarv(32) : NULL;
	mdmo->par2 = par2;
	mdmo->align_mode = align_mode;
	thread_end_init(mdmo);
	ret = 0;
	if(online){
		fr = stdin_filereader();
		seq = NULL;
		thread_beg_iter(mdmo);
		mdmo->rdseqs = init_basebank();
		thread_end_iter(mdmo);
		for(i=0;fread_seq(&seq, fr);i++){
			if((i % n_job) != i_job) continue;
			thread_waitfor_one_idle(mdmo);
			ret += print_hits_dmo(wt, mdmo->hits, mdmo->kigars, out);
			mdmo->pbid = wt->reads->size;
			clear_basebank(mdmo->rdseqs);
			mdmo->rd.tag = strdup(seq->tag.string);
			mdmo->rd.rdoff = 0;
			mdmo->rd.rdlen = seq->seq.size;
			push_pbreadv(wt->reads, mdmo->rd);
			seq2basebank(mdmo->rdseqs, seq->seq.string, seq->seq.size);
			thread_wake(mdmo);
		}
		fclose_filereader(fr);
	} else {
		for(i=0;i<wt->n_rd;i++){
			if(hzm_debug == 0 && (i % 100) == 0){
				fprintf(hzm_debug_out, "\r%012u\t%llu\t%02f%%\t%0.2f secs", i, (unsigned long long)ret, 100.0 * i / wt->n_rd, 1.0 * count_clock(secs) / CLOCKS_PER_SEC); fflush(hzm_debug_out);
			}
			if((i % n_job) != i_job) continue;
			thread_waitfor_one_idle(mdmo);
			ret += print_hits_dmo(wt, mdmo->hits, mdmo->kigars, out);
			mdmo->pbid = i;
			mdmo->rd = wt->reads->buffer[i];
			thread_wake(mdmo);
		}
	}
	thread_waitfor_all_idle(mdmo);
	thread_beg_close(mdmo);
	ret += print_hits_dmo(wt, mdmo->hits, mdmo->kigars, out);
	if(mdmo->rdseqs != wt->rdseqs) free_basebank(mdmo->rdseqs);
	free_wtovlv(mdmo->hits);
	if(mdmo->kigars) free_kigarv(mdmo->kigars);
	thread_end_close(mdmo);
	fprintf(hzm_debug_out, "\nprogress: %u %llu 100.00%%, %0.2f CPU seconds\n", (int)i, (unsigned long long)ret, 1.0 * end_clock(secs) / CLOCKS_PER_SEC); fflush(hzm_debug_out);
	return ret;
}

int usage(){
	printf(
	"DMO: Overlaper of long reads using dot matrix alignment\n"
	"SMARTdenovo: Ultra-fast de novo assembler for high noisy long reads\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: dmo [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -P <int>    Total parallel jobs, [1]\n"
	" -p <int>    Index of current job (0-based), [0]\n"
	"             Suppose to run dmo parallelly in 60 nodes. For node1, -P 60 -p 0; node2, -P 60 -p 1, ...\n"
	" -i <string> Long reads sequences file, + *\n"
	" -I          Read query sequences from stdin\n"
	" -R          Reads in original order, default in length DSC\n"
	" -J <int>    Jack knife of original read length, [0]\n"
	" -o <string> Output file of alignments, *\n"
	" -Q          Append kigar to alignment\n"
	" -f          Force overwrite\n"
	"----------------------Align Parameter------------------------------------\n"
	" -H          Trun off homopolymer compression\n"
	" -k <int>    Kmer size, 5 <= value <= 27, [16]\n"
	" -W <int>    Max size of insertion in the middle of kmer when querying, [0]\n"
	"             PART1|ins|PART2, PART1 + PART2 = ksize, PART2 = ksize/2, ins <= max_kgap, max_kgap + ksize <= 32\n"
	" -K <float>  Filter high frequency kmers, maybe repetitive, [0.05]\n"
	"             if K >= 1, take the integer value as cutoff\n"
	"             else, mask the top fraction part high frequency kmers\n"
	" -E <int>    Filter low frequency kmers less than it, [2]\n"
	" -S <int>    Subsampling kmers, 1/<-S> kmers are indexed, [2]\n"
	" -x <int>    Intra-block: Max gap  [256]\n"
	" -y <int>    Intra-block: Max deviation [96]\n"
	" -z <int>    Intra-block: Min kmer [3]\n"
	" -w <int>    Inter-block: deviation penalty [1.0]\n"
	" -u <int>    Inter-block: gap penalty [0.1]\n"
	" -l <int>    Min length of alignment, [1000]\n"
	" -m <int>    Min matched of alignment, [100]\n"
	" -s <float>  Max length variation of two aligned fragments, [0.4]\n"
	"-------------------------------------------------------------------------\n"
	" -Z <int>    Perform alignment on matched pairs based on different kmer-size and other paramerter\n"
	"             5 <= value <= 27, like -k\n"
	"             All align parameters after -Z is set to realignment parameters\n"
	"             default: -K 0 -S 1 -x 256 -y 64 -z 3 -w 1.0 -u 0.1 -l 1000 -m 100 -s 0.4\n"
	"             Align mode for realignment is fixed at 0\n"
	" -F <int>    Align mode: [0]\n"
	"                0, using kmer synteny, one hit per pair of reads;\n"
	"                1, using kmer synteny, may more than one hits;\n"
	"                2, using kmer abundance, more than one hits; [0]\n"
	" -v          Verbose, BE careful, HUGEEEEEEEE output on STDOUT\n"
	"\n"
	);
	return 1;
}

int main(int argc, char **argv){
	DMOPar pars[2];
	DMO *wt;
	cplist *pbs;
	FileReader *fr;
	Sequence *seq;
	char *output;
	FILE *out;
	unsigned long long tot_bp;
	int c, ncpu, min_rdlen, overwrite, n_job, i_job, align_mode, read_order, has_kigar, zalign, online;
	output = NULL;
	overwrite = 0;
	n_job = 1;
	i_job = 0;
	ncpu = 1;
	min_rdlen = 0;
	pars[0].min_mat = 100;
	pars[0].min_aln = 1000;
	pars[0].min_sm  = 0.05;
	pars[0].aln_var = 0.4;
	pars[0].hk = 1;
	pars[0].ksize = 16;
	pars[0].kmax  = 0;
	pars[0].ktop  = 0.05;
	pars[0].kmin  = 2;
	pars[0].kgap  = 0;
	pars[0].max_hit = 0;
	pars[0].hzmh_kmer_mod = HZMH_KMER_MOD * 2;
	pars[0].hzmh_kmer_win = 1;
	pars[0].xvar = 256;
	pars[0].yvar = 96;
	pars[0].zmin = 3;
	pars[0].max_overhang = 600;
	pars[0].deviation_penalty = 1.0;
	pars[0].gap_penalty = 0.1;
	pars[0].ttr_len = 5000;
	pars[0].ttr_cnt = 2;
	pars[0].ttr_win = 100;
	pars[0].ttr_drate = 0.05;
	pars[0].ttr_crate = 0.20;
	pars[0].ttr_xrate = 0.80;
	pars[1].min_mat = 100;
	pars[1].min_aln = 1000;
	pars[1].min_sm  = 0.05;
	pars[1].aln_var = 0.4;
	pars[1].hk = 1;
	pars[1].ksize = 16;
	pars[1].kmax  = 0;
	pars[1].ktop  = 0.05;
	pars[1].kmin  = 2;
	pars[1].kgap  = 0;
	pars[1].max_hit = 0;
	pars[1].hzmh_kmer_mod = HZMH_KMER_MOD;
	pars[1].hzmh_kmer_win = 1;
	pars[1].xvar = 256;
	pars[1].yvar = 64;
	pars[1].zmin = 3;
	pars[1].max_overhang = 600;
	pars[1].deviation_penalty = 1.0;
	pars[1].gap_penalty = 0.1;
	pars[1].ttr_len = 5000;
	pars[1].ttr_cnt = 2;
	pars[1].ttr_win = 100;
	pars[1].ttr_drate = 0.05;
	pars[1].ttr_crate = 0.20;
	pars[1].ttr_xrate = 0.80;
	align_mode = 0;
	read_order = 1;
	has_kigar = 0;
	online = 0;
	pbs = init_cplist(4);
	zalign = 0;
	while((c = getopt(argc, argv, "ht:P:p:i:IRb:J:o:QS:fHk:W:K:E:x:y:z:w:u:l:m:s:Z:F:v")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'P': n_job = atoi(optarg); break;
			case 'p': i_job = atoi(optarg); break;
			case 'i': push_cplist(pbs, optarg); break;
			case 'I': online = 1; break;
			case 'R': read_order = 0; break;
			case 'J': min_rdlen = atoi(optarg); break;
			case 'o': output = optarg; break;
			case 'Q': has_kigar = 1; break;
			case 'f': overwrite = 1; break;
			case 'H': pars[zalign].hk = 0; break;
			case 'k': pars[zalign].ksize = atoi(optarg); break;
			case 'W': pars[zalign].kgap  = atoi(optarg); break;
			case 'K': pars[zalign].kmax = atoi(optarg); break;
			case 'E': pars[zalign].kmin = atoi(optarg); break;
			case 'S': pars[zalign].hzmh_kmer_mod = atoi(optarg) * HZMH_KMER_MOD; break;
			case 'x': pars[zalign].xvar = atoi(optarg); break;
			case 'y': pars[zalign].yvar = atoi(optarg); break;
			case 'z': pars[zalign].zmin = atoi(optarg); break;
			case 'w': pars[zalign].deviation_penalty = atof(optarg); break;
			case 'u': pars[zalign].gap_penalty = atof(optarg); break;
			case 'l': pars[zalign].min_aln = atoi(optarg); break;
			case 'm': pars[zalign].min_mat = atoi(optarg); break;
			case 's': pars[zalign].aln_var = atof(optarg); break;
			case 'Z': zalign = 1; pars[zalign].ksize = atoi(optarg); break;
			case 'F': align_mode = atoi(optarg); break;
			case 'v': hzm_debug ++; break;
			default: return usage();
		}
	}
	if(output == NULL) return usage();
	if(!overwrite && strcmp(output, "-") && file_exists(output)){
		fprintf(hzm_debug_out, "File exists! '%s'\n\n", output);
		return usage();
	}
	if(pbs->size == 0) return usage();
	if(pars[0].ksize > DMO_MAX_KSIZE || pars[0].ksize < 5) return usage();
	if(pars[1].ksize == 0) zalign = 0;
	if(zalign && (pars[1].ksize > DMO_MAX_KSIZE || pars[1].ksize < 5)) return usage();
	pars[0].max_overhang = 2 * pars[0].xvar;
	pars[1].max_overhang = 2 * pars[1].xvar;
	wt = init_dmo(pars + 0);
	wt->len_order = read_order;
	if((fr = fopen_m_filereader(pbs->size, pbs->buffer)) == NULL){
		fprintf(hzm_debug_out, " -- Cannot open %s in %s -- %s:%d --\n", pbs->buffer[0], __FUNCTION__, __FILE__, __LINE__); exit(1);
	}
	fprintf(hzm_debug_out, "[%s] loading long reads\n", date());
	seq = NULL;
	tot_bp = 0;
	while(fread_seq(&seq, fr)){
		if(seq->seq.size < min_rdlen) continue;
		push_long_read_dmo(wt, seq->name.string, seq->name.size, seq->seq.string, seq->seq.size);
		tot_bp += seq->seq.size;
		if(hzm_debug == 0 && (wt->n_rd % 1000) == 0){
			fprintf(hzm_debug_out, "\r%u", wt->n_rd); fflush(hzm_debug_out);
		}
	}
	fclose_filereader(fr);
	if(hzm_debug == 0) fprintf(hzm_debug_out, "\r");
	fprintf(hzm_debug_out, "[%s] Done, %u reads (length >= %d), %llu bp\n", date(), (unsigned)wt->n_rd, min_rdlen, tot_bp);
	ready_dmo(wt);
	index_dmo(wt, 0, wt->n_rd, 0, ncpu);
	out = strcmp(output, "-")? fopen(output, "w") : stdout;
	fprintf(hzm_debug_out, "[%s] calculating overlaps, %d threads\n", date(), ncpu);
	overlap_dmo(wt, zalign? pars + 1 : NULL, online, align_mode, has_kigar, ncpu, n_job, i_job, out);
	if(strcmp(output, "-")) fclose(out);
	fprintf(hzm_debug_out, "[%s] Done\n", date());
	free_dmo(wt);
	free_cplist(pbs);
	return 0;
}
