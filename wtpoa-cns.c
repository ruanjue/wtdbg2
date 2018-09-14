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

#include "tripoa.h"
#include "kswx.h"
#include "file_reader.h"

#define MCNS_TASK_POACNS	1
#define MCNS_TASK_OVERLAP	2

thread_beg_def(mcns);
uint32_t eid;
TriPOG *g;
u1v *seq1, *seq2;
u4i node1, node2;
int reglen;
kswx_t ret;
u32list *cigars;
int task;
thread_end_def(mcns);

thread_beg_func(mcns);
TriPOG *g;
u8list *mem_cache;
u1v    *mem_buffer;
u32list *cigars[2];
kswx_t *xs[2];
int qb, qe, tb, te;
mem_cache = init_u8list(1024);
cigars[0] = mcns->cigars;
cigars[1] = NULL;
xs[0] = malloc(sizeof(kswx_t));
xs[1] = NULL;
g = mcns->g;
mem_buffer = init_u1v(1024);
thread_beg_loop(mcns);
if(mcns->task == MCNS_TASK_POACNS){
	if(g->seqs->nseq){
	end_tripog(g);
	}
} else if(mcns->task == MCNS_TASK_OVERLAP){
	if(mcns->seq2->size == 0){
		mcns->ret = KSWX_NULL;
		continue;
	}
	qb = 0; qe = mcns->seq1->size;
	tb = 0; te = mcns->seq2->size;
	if(qe > mcns->reglen) qb = qe - mcns->reglen;
	if(te > mcns->reglen) te = mcns->reglen;
	kswx_overlap_align_core(xs, cigars, qe - qb, mcns->seq1->buffer + qb, te - tb, mcns->seq2->buffer + tb, 1, g->pogs[0]->M, g->pogs[0]->X, g->pogs[0]->I, g->pogs[0]->D, -1, mem_cache);
	xs[0]->qb += qb;
	xs[0]->qe += qb;
	xs[0]->tb += tb;
	xs[0]->te += tb;
	mcns->ret = *xs[0];
}
thread_end_loop(mcns);
free(xs[0]);
free_u8list(mem_cache);
free_u1v(mem_buffer);
thread_end_func(mcns);

int revise_joint_point(u32list *cigars, int *qe, int *te){
	u4i i, op, ln, max;
	int qq, q, tt, t;
	q = t = 0;
	qq = tt = 0;
	max = 0;
	for(i=1;i<=cigars->size;i++){
		op = cigars->buffer[cigars->size - i] & 0xF;
		ln = cigars->buffer[cigars->size - i] >> 4;
		if(op == 0){
			if(ln > max){
				qq = q; tt = t;
				max = ln;
			}
			q += ln;
			t += ln;
		} else if(op == 1){
			q += ln;
		} else {
			t += ln;
		}
	}
	if(cns_debug){
		fprintf(stderr, "qe = %d -> %d\n", *qe, (*qe) - qq);
		fprintf(stderr, "te = %d -> %d\n", *te, (*te) - tt);
		fflush(stderr);
	}
	*qe -= qq;
	*te -= tt;
	return 1;
}

int revise_joint_point2(u32list *cigars, int *qe, int *te, int overhang){
	u4i i, op, ln;
	int qq, q, tt, t;
	q = t = 0;
	qq = tt = 0;
	for(i=1;i<=cigars->size;i++){
		op = cigars->buffer[cigars->size - i] & 0xF;
		ln = cigars->buffer[cigars->size - i] >> 4;
		switch(op){
			case 1: q += ln; break;
			case 2: t += ln; break;
			default: qq = q; tt = t; q += ln; t += ln;
		}
		if(q >= overhang && t >= overhang){
			if(cns_debug){
				fprintf(stderr, "qe = %d -> %d\n", *qe, (*qe) - qq);
				fprintf(stderr, "te = %d -> %d\n", *te, (*te) - tt);
				fflush(stderr);
			}
			*qe -= qq;
			*te -= tt;
			return 1;
		}
	}
	return 0;
}

int run_cns(FileReader *fr, u4i ncpu, int use_sse, u4i seqmax, int winlen, int winmin, int fail_skip, int W, int M, int X, int I, int D, int rW, int mincnt, float minfreq, int reglen, FILE *out){
	String *tag, *seq;
	u1v *cseqs;
	u4v *cxs, *cys, *tes, *qes;
	u4i i, m, eid, beg, end, meths[2];
	int c, j, sl, b, e;
	char *ss;
	thread_preprocess(mcns);
	tag = init_string(32);
	seq = init_string(32);
	cseqs = init_u1v(32);
	cxs = init_u4v(32);
	cys = init_u4v(32);
	tes = init_u4v(32);
	qes = init_u4v(32);
	thread_beg_init(mcns, ncpu);
	mcns->eid = 0;
	mcns->g = init_tripog(winlen, winmin, fail_skip, M, X, I, D, W, use_sse, rW, mincnt, minfreq);
	mcns->seq1 = init_u1v(32);
	mcns->seq2 = init_u1v(32);
	mcns->node1 = 0;
	mcns->node2 = 0;
	mcns->reglen = reglen;
	mcns->ret = KSWX_NULL;
	mcns->cigars = init_u32list(16);
	mcns->task = 0;
	thread_end_init(mcns);
	eid = 0;
	meths[0] = meths[1] = 0;
	thread_wait_one(mcns);
	while(1){
		c = fread_table(fr);
		if(c == -1 || fr->line->string[0] == 'E' || fr->line->string[0] == '>'){
			thread_wake(mcns);
			thread_wait_one(mcns);
			if(mcns->task == MCNS_TASK_POACNS && mcns->g->seqs->nseq){
				meths[mcns->g->is_tripog] ++;
				if(cns_debug){
					fprintf(stderr, "%s_%d_N%u_N%u\t%d\t%d\t", tag->string, mcns->eid, mcns->node1, mcns->node2, mcns->g->seqs->nseq, (u4i)mcns->g->cns->size);
					println_seq_basebank(mcns->g->cns, 0, mcns->g->cns->size, stderr);
				}
				cxs->buffer[mcns->eid] = cseqs->size;
				inc_u1v(cseqs, mcns->g->cns->size);
				bitseq_basebank(mcns->g->cns, 0, mcns->g->cns->size, cseqs->buffer + cxs->buffer[mcns->eid]);
				cys->buffer[mcns->eid] = cseqs->size;
			}
			beg_tripog(mcns->g);
			mcns->task = MCNS_TASK_POACNS;
			mcns->eid = eid ++;
			push_u4v(cxs, 0);
			push_u4v(cys, 0);
			if(c != -1 && fr->line->string[0] == 'E'){
				mcns->node1 = atoll(get_col_str(fr, 2) + 1);
				mcns->node2 = atoll(get_col_str(fr, 4) + 1);
				continue;
			}
			if(tag->size){
				thread_beg_iter(mcns);
				thread_wait(mcns);
				if(mcns->task == MCNS_TASK_POACNS && mcns->g->seqs->nseq){
					meths[mcns->g->is_tripog] ++;
					if(cns_debug){
						fprintf(stderr, "%s_%d_N%u_N%u\t%d\t%d\t", tag->string, mcns->eid, mcns->node1, mcns->node2, mcns->g->seqs->nseq, (u4i)mcns->g->cns->size);
						println_seq_basebank(mcns->g->cns, 0, mcns->g->cns->size, stderr);
					}
					cxs->buffer[mcns->eid] = cseqs->size;
					inc_u1v(cseqs, mcns->g->cns->size);
					bitseq_basebank(mcns->g->cns, 0, mcns->g->cns->size, cseqs->buffer + cxs->buffer[mcns->eid]);
					cys->buffer[mcns->eid] = cseqs->size;
				}
				beg_tripog(mcns->g);
				mcns->eid = 0;
				mcns->task = 0;
				mcns->node1 = 0;
				mcns->node2 = 0;
				thread_end_iter(mcns);
				push_u4v(qes, 0);
				push_u4v(tes, 0);
				for(i=1;i<eid;i++){
					thread_wait_one(mcns);
					if(mcns->task == MCNS_TASK_OVERLAP && mcns->ret.aln > 0){
						if(cns_debug){
							fprintf(stderr, "#%s_%d\t%d\t%d\t%d", tag->string, mcns->eid - 1, (int)mcns->seq1->size, mcns->ret.qb, mcns->ret.qe);
							fprintf(stderr, "\t%s_%d\t%d\t%d\t%d", tag->string, mcns->eid, (int)mcns->seq2->size, mcns->ret.tb, mcns->ret.te);
							fprintf(stderr, "\t%d\t%d\t%d\t%d\t%d\n", mcns->ret.aln, mcns->ret.mat, mcns->ret.mis, mcns->ret.ins, mcns->ret.del);
						}
						{
							b = mcns->ret.qe;
							e = mcns->ret.te;
							if(1){
								revise_joint_point(mcns->cigars, &b, &e);
							}
							qes->buffer[mcns->eid] = b;
							tes->buffer[mcns->eid] = e;
						}
					}
					mcns->task = MCNS_TASK_OVERLAP;
					mcns->eid = i;
					clear_u1v(mcns->seq1); append_array_u1v(mcns->seq1, cseqs->buffer + cxs->buffer[i-1], cys->buffer[i-1] - cxs->buffer[i-1]);
					clear_u1v(mcns->seq2); append_array_u1v(mcns->seq2, cseqs->buffer + cxs->buffer[i], cys->buffer[i] - cxs->buffer[i]);
					mcns->ret = KSWX_NULL;
					push_u4v(qes, mcns->seq1->size);
					push_u4v(tes, 0);
					thread_wake(mcns);
				}
				push_u4v(qes, cys->buffer[eid-1] - cxs->buffer[eid-1]);
				push_u4v(tes, 0);
				thread_beg_iter(mcns);
				thread_wait(mcns);
				if(mcns->task == MCNS_TASK_OVERLAP && mcns->ret.aln > 0){
					if(cns_debug){
						fprintf(stderr, "#%s_%d\t%d\t%d\t%d", tag->string, mcns->eid - 1, (int)mcns->seq1->size, mcns->ret.qb, mcns->ret.qe);
						fprintf(stderr, "\t%s_%d\t%d\t%d\t%d", tag->string, mcns->eid, (int)mcns->seq2->size, mcns->ret.tb, mcns->ret.te);
						fprintf(stderr, "\t%d\t%d\t%d\t%d\t%d\n", mcns->ret.aln, mcns->ret.mat, mcns->ret.mis, mcns->ret.ins, mcns->ret.del);
					}
					{
						b = mcns->ret.qe;
						e = mcns->ret.te;
						if(1){
							revise_joint_point(mcns->cigars, &b, &e);
						}
						qes->buffer[mcns->eid] = b;
						tes->buffer[mcns->eid] = e;
					}
				}
				mcns->ret = KSWX_NULL;
				mcns->eid = 0;
				mcns->task = 0;
				thread_end_iter(mcns);
				// generate contig seq
				clear_string(seq);
				for(i=0;i<eid;i++){
					beg = cxs->buffer[i] + tes->buffer[i];
					end = cxs->buffer[i] + qes->buffer[i + 1];
					for(m=beg;m<end;m++){
						push_string(seq, bit_base_table[cseqs->buffer[m]]);
					}
					if(cns_debug){
						fprintf(stderr, "=%s_%d\t%d\t%d\n", tag->string, i, seq->size - (end - beg), seq->size);
					}
				}
				fprintf(out, ">%s len=%d\n", tag->string, seq->size);
				for(j=0;j<seq->size;j+=100){
					sl = num_min(j + 100, seq->size);
					char ch = seq->string[sl];
					seq->string[sl] = 0;
					fprintf(out, "%s\n", seq->string + j);
					seq->string[sl] = ch;
				}
			}
			if(c == -1) break;
			clear_string(tag);
			ss = get_col_str(fr, 0) + 1;
			sl = get_col_len(fr, 0) - 1;
			for(j=0;j<sl;j++){
				if(ss[j] == ' ') break;
				push_string(tag, ss[j]);
			}
			clear_string(seq);
			eid = 0;
			clear_u1v(cseqs);
			clear_u4v(cxs);
			clear_u4v(cys);
			clear_u4v(tes);
			clear_u4v(qes);
			thread_wait_one(mcns);
		} else if(fr->line->string[0] == 'S' || fr->line->string[0] == 's'){
			if(mcns->g->seqs->nseq >= seqmax) continue;
			ss = get_col_str(fr, 5);
			sl = get_col_len(fr, 5);
			if(UInt(sl) > POG_RDLEN_MAX){
				sl = POG_RDLEN_MAX;
			}
			push_tripog(mcns->g, ss, sl);
		}
	}
	fprintf(stderr, " -- %u %u in %s -- %s:%d --\n", meths[0], meths[1], __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
	thread_beg_close(mcns);
	free_tripog(mcns->g);
	free_u1v(mcns->seq1);
	free_u1v(mcns->seq2);
	free_u32list(mcns->cigars);
	thread_end_close(mcns);
	free_u1v(cseqs);
	free_u4v(cxs);
	free_u4v(cys);
	free_u4v(tes);
	free_u4v(qes);
	free_string(tag);
	free_string(seq);
	return 0;
}

int usage(){
	printf(
	"WTPOA-CNS: Consensuser for wtdbg using PO-MSA\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: wtpoa-cns [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -i <string> Input file(s) *.ctg.lay from wtdbg, +, [STDIN]\n"
	" -o <string> Output files, [STDOUT]\n"
	" -f          Force overwrite\n"
	" -j <int>    Expected max length of node, or say the overlap length of two adjacent units in layout file, [1500] bp\n"
	" -M <int>    Match score, [2]\n"
	" -X <int>    Mismatch score, [-5]\n"
	" -I <int>    Insertion score, [-2]\n"
	" -D <int>    Deletion score, [-4]\n"
	" -B <int>    Bandwidth, [96]\n"
	" -W <int>    Window size in the middle of the first read for fast align remaining reads, [200]\n"
	"             If $W is negative, will disable fast align, but use the abs($W) as Band align score cutoff\n"
	" -w <int>    Min size of aligned size in window, [$W * 0.5]\n"
	" -A          Abort TriPOA when any read cannot be fast aligned, then try POA\n"
	//" -S <int>    SSE speedup in graph alignment, [1]\n"
	" -R <int>    Realignment bandwidth, 0: disable, [16]\n"
	" -C <int>    Min count of bases to call a consensus base, [3]\n"
	" -F <float>  Min frequency of non-gap bases to call a consensus base, [0.5]\n"
	" -N <int>    Max number of reads in PO-MSA, [20]\n"
	"             Keep in mind that I am not going to generate high accurate consensus sequences here\n"
	" -v          Verbose\n"
	"\n");
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	cplist *infs;
	FILE *out;
	char *outf;
	int reglen, use_sse, bandwidth, rW, winlen, winmin, fail_skip, M, X, I, D, mincnt, seqmax;
	float minfreq;
	int c, ncpu, overwrite;
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	ncpu = 1;
	use_sse = 2;
	seqmax = 20;
	bandwidth = 96;
	winlen = 200;
	winmin = 0;
	fail_skip = 1;
	reglen = 1500;
	M = 2;
	X = -5;
	I = -2;
	D = -4;
	rW = 16;
	mincnt = 3;
	minfreq = 0.5;
	infs = init_cplist(4);
	outf = NULL;
	overwrite = 0;
	while((c = getopt(argc, argv, "hvt:i:o:fS:B:W:w:AM:X:I:D:R:C:F:N:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'i': push_cplist(infs, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'S': use_sse = atoi(optarg); break;
			case 'B': bandwidth = atoi(optarg); break;
			case 'W': winlen = atoi(optarg); break;
			case 'w': winmin = atoi(optarg); break;
			case 'A': fail_skip = 0; break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'I': I = atoi(optarg); break;
			case 'D': D = atoi(optarg); break;
			case 'R': rW = atoi(optarg); break;
			case 'C': mincnt = atoi(optarg); break;
			case 'F': minfreq = atof(optarg); break;
			case 'N': seqmax = atoi(optarg); break;
			case 'v': cns_debug ++; break;
			default: return usage();
		}
	}
	if(winmin <= 0){
		winmin = winlen * 0.5;
	}
	if(ncpu <= 0 && _sig_proc_deamon) ncpu = _sig_proc_deamon->ncpu;
	if(ncpu <= 0){
		fprintf(stderr, " -- Invalid cpu number '%d' in %s -- %s:%d --\n", ncpu, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	if(outf && !overwrite && file_exists(outf)){
		fprintf(stderr, "File exists! '%s'\n\n", outf);
		return usage();
	}
	if(infs->size) fr = fopen_m_filereader(infs->size, infs->buffer);
	else fr = stdin_filereader();
	if(outf){
		out = open_file_for_write(outf, NULL, 1);
	} else out = stdout;
	run_cns(fr, ncpu, use_sse, seqmax, winlen, winmin, fail_skip, bandwidth, M, X, I, D, rW, mincnt, minfreq, reglen, out);
	fclose_filereader(fr);
	if(outf) fclose(out);
	free_cplist(infs);
	END_STAT_PROC_INFO(stderr);
	return 0;
}
