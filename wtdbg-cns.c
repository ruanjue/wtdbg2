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

#include "kswx.h"
#include "dbgcns.h"
#include "dagcns.h"
#include "filereader.h"

static int cns_debug = 0;

#define MCNS_TASK_DBGCNS	1
#define MCNS_TASK_DAGCNS	2
#define MCNS_TASK_OVERLAP	3

thread_beg_def(mcns);
uint32_t eid;
CNS *cns;
u1v *seq1, *seq2;
u4i node1, node2;
int reglen;
int W, M, X, I, D, E;
int candidate_mode;
int corr_struct;
f4i pM, pX, pI, pD;
kswx_t ret;
u32list *cigars;
int task;
thread_end_def(mcns);

thread_beg_func(mcns);
CNS *cns;
DAGCNS *dag;
GEGraph *g;
bdpnodev *bnodes;
bdpedgev *bedges;
bdplinkv *linkstack;
u8list *mem_cache;
u1v    *mem_buffer;
u32list *cigars[2];
kswx_t *xs[2];
u4i i, nbeg;
blk_t *blk;
int qb, qe, tb, te;
mem_cache = init_u8list(1024);
cigars[0] = mcns->cigars;
cigars[1] = NULL;
xs[0] = malloc(sizeof(kswx_t));
xs[1] = NULL;
cns = mcns->cns;
dag = init_dagcns(mcns->W, mcns->M, mcns->X, mcns->I, mcns->D, mcns->E, mcns->pM, mcns->pX, mcns->pI, mcns->pD);
g = init_gegraph();
bnodes = init_bdpnodev(32);
bedges = init_bdpedgev(32);
linkstack = init_bdplinkv(32);
mem_buffer = init_u1v(1024);
thread_beg_loop(mcns);
if(mcns->task == MCNS_TASK_DBGCNS){
	ready_cns(cns);
	run_cns(cns, mcns->candidate_mode, mcns->corr_struct);
} else if(mcns->task == MCNS_TASK_DAGCNS){
	// force mcns->candidate_mode to 2
	// use data from DBGCNS struct
	ready_cns(cns);
	clear_u8list(dag->cns);
	blk = ref_blkv(cns->qblks, 0);
	append_array_u8list(dag->cns, cns->qseqs->buffer + blk->off, blk->len);
	gen_pregraph_dagcns(dag);
	for(i=1;i<cns->qblks->size;i++){
		blk = ref_blkv(mcns->cns->qblks, i);
		nbeg = branched_dynamic_programming_alignment(dag, mcns->cns->qseqs->buffer + blk->off, blk->len, g, bnodes, bedges, mem_buffer);
		if(nbeg == 0){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			continue;
		}
		bdpgraph2dagcns(dag, g, bnodes, bedges, nbeg, linkstack);
	}
	merge_nodes_dagcns(dag);
	gen_consensus_dagcns(dag, NULL);
	if(cns_debug){
		fprintf(stderr, "DAG%d\t%d\t", mcns->eid, (int)dag->cns->size);
		print_seq_dagcns(dag, stderr);
		fprintf(stderr, "\n");
	}
	clear_u1v(mcns->cns->cns);
	append_u1v(mcns->cns->cns, (u1v*)dag->cns);
	clear_string(mcns->cns->seq);
	for(i=0;i<dag->cns->size;i++) add_char_string(mcns->cns->seq, bit_base_table[dag->cns->buffer[i]]);
} else if(mcns->task == MCNS_TASK_OVERLAP){
	if(mcns->seq2->size == 0){
		mcns->ret = KSWX_NULL;
		continue;
	}
	qb = 0; qe = mcns->seq1->size;
	tb = 0; te = mcns->seq2->size;
	if(qe > mcns->reglen) qb = qe - mcns->reglen;
	if(te > mcns->reglen) te = mcns->reglen;
	kswx_overlap_align_core(xs, cigars, qe - qb, mcns->seq1->buffer + qb, te - tb, mcns->seq2->buffer + tb, 1, mcns->M, mcns->X, mcns->I, mcns->D, mcns->E, mem_cache);
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
free_dagcns(dag);
free_gegraph(g);
free_bdpnodev(bnodes);
free_bdpedgev(bedges);
free_bdplinkv(linkstack);
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

uint32_t run(int reglen, int ksize, int Z, int W, int M, int X, int I, int D, int E, int H, int L, int XX, int OO, int EE, int cns_model, f4i pM, f4i pX, f4i pI, f4i pD, int candidate_mode, int corr_struct, uint32_t ncpu, FileReader *fr, FILE *out){
	String *tag, *seq;
	u1v *cseqs;
	u4v *cxs, *cys, *tes, *qes;
	uint32_t i, m, eid, beg, end;
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
	mcns->cns = init_cns(ksize, Z, W, 0, X, I, D, E, H, L);
	mcns->seq1 = init_u1v(32);
	mcns->seq2 = init_u1v(32);
	mcns->node1 = 0;
	mcns->node2 = 0;
	mcns->reglen = reglen;
	mcns->W = 128;
	mcns->M = M;
	mcns->X = XX;
	mcns->I = OO;
	mcns->D = OO;
	mcns->E = EE;
	mcns->pM = pM;
	mcns->pX = pX;
	mcns->pI = pI;
	mcns->pD = pD;
	mcns->candidate_mode = candidate_mode;
	mcns->corr_struct = corr_struct;
	mcns->ret = KSWX_NULL;
	mcns->cigars = init_u32list(16);
	mcns->task = 0;
	thread_end_init(mcns);
	eid = 0;
	thread_wait_one(mcns);
	while(1){
		c = readtable_filereader(fr);
		if(c == -1 || fr->line->string[0] == 'E' || fr->line->string[0] == '>'){
			thread_wake(mcns);
			thread_wait_one(mcns);
			if(mcns->task == cns_model && mcns->cns->seq->size){
				if(cns_debug){
					fprintf(stderr, "%s_%d_N%u_N%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", tag->string, mcns->eid, mcns->node1, mcns->node2, mcns->cns->qlen, mcns->cns->seq->size, mcns->cns->max_score, mcns->cns->alns[0], mcns->cns->alns[1], mcns->cns->alns[2], mcns->cns->alns[3], mcns->cns->seq->string);
				}
				cxs->buffer[mcns->eid] = cseqs->size;
				for(j=0;j<mcns->cns->seq->size;j++) push_u1v(cseqs, base_bit_table[(int)mcns->cns->seq->string[j]]);
				cys->buffer[mcns->eid] = cseqs->size;
			}
			reset_cns(mcns->cns);
			clear_string(mcns->cns->seq);
			mcns->task = cns_model;
			mcns->eid = eid ++;
			push_u4v(cxs, 0);
			push_u4v(cys, 0);
			if(fr->line->string[0] == 'E'){
				mcns->node1 = atoll(get_col_str(fr, 2) + 1);
				mcns->node2 = atoll(get_col_str(fr, 4) + 1);
				continue;
			}
			if(tag->size){
				thread_beg_iter(mcns);
				thread_wait(mcns);
				if(mcns->task == cns_model && mcns->cns->seq->size){
					if(cns_debug){
						//fprintf(stderr, "%s_%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", tag->string, mcns->eid, mcns->cns->qlen, mcns->cns->seq->size, mcns->cns->max_score, mcns->cns->alns[0], mcns->cns->alns[1], mcns->cns->alns[2], mcns->cns->alns[3], mcns->cns->seq->string);
						fprintf(stderr, "%s_%d_N%u_N%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", tag->string, mcns->eid, mcns->node1, mcns->node2, mcns->cns->qlen, mcns->cns->seq->size, mcns->cns->max_score, mcns->cns->alns[0], mcns->cns->alns[1], mcns->cns->alns[2], mcns->cns->alns[3], mcns->cns->seq->string);
						if(0){
							u4i j, x;
							for(j=x=0;j<mcns->cns->cigars->size;j++){
								if(mcns->cns->cigars->buffer[j] == DBGCNS_PATH_I){
									fputc('-', stderr);
								} else {
									fputc("ACGT"[mcns->cns->cns->buffer[x]], stderr); x ++;
								}
							}
							fputc('\n', stderr);
							for(j=0;j<mcns->cns->cigars->size;j++){
								fputc("MXID"[mcns->cns->cigars->buffer[j]], stderr);
							}
							fputc('\n', stderr);
						}
					}
					cxs->buffer[mcns->eid] = cseqs->size;
					for(j=0;j<mcns->cns->seq->size;j++) push_u1v(cseqs, base_bit_table[(int)mcns->cns->seq->string[j]]);
					cys->buffer[mcns->eid] = cseqs->size;
				}
				reset_cns(mcns->cns);
				clear_string(mcns->cns->seq);
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
			ss = get_col_str(fr, 5);
			sl = get_col_len(fr, 5);
			add_seq_cns(mcns->cns, ss, sl, (fr->line->string[0] == 'S'));
		}
	}
	thread_beg_close(mcns);
	free_cns(mcns->cns);
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
	"WTDBG-CNS: Consensuser for wtdbg\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.1\n"
	"Usage: wtdbg-cns [options]\n"
	"Options:\n"
	" -t <int>    Number of threads, [1]\n"
	" -i <string> Input file(s) *.utg.cns from wtdbg, +, [STDIN]\n"
	" -o <string> Output files, [STDOUT]\n"
	" -f          Force overwrite\n"
	" -j <int>    Expected length of node, or say the overlap length of two adject units in layout file, [1000] bp\n"
	"-----------------BEG DBG options---------------------------------\n"
	" -k <int>    Kmer size for long reads, [15]\n"
	" -Z <int>    Z-cutoff, drop the lower  (score / <-X>), [4]\n"
	" -W <int>    W-cutoff, drop the lagger (position), [48]\n"
	"             In DAG correction, -W set the bandwidth of alignment\n"
	" -H <int>    High coverage bonus, [1]\n"
	" -L <int>    High coverage cutoff = avg_cov / <-L> [10]\n"
	" -c <int>    Candidate strategy, 0: best-kmers, 1: median length, 2: first (include), 3: first (exclude), 4: longest, 5, shortest, [0]\n"
	"             In DAG correction, force to use strategy 2\n"
	"-----------------END DBG options---------------------------------\n"
	" -M <int>    Match score, [2]\n"
	" -X <int>    Mismatch score, [-7]\n"
	" -I <int>    Insertion score, [-3]\n"
	" -D <int>    Deletion score, [-4]\n"
	" -E <int>    Gap extension score, [-2]\n"
	" -m <int>    1: DBG correction; 2: DAG correction, [1]\n"
	" -S <int>    whether to correct structure before error correction, [1]\n"
	" -v          Verbose\n"
	" -V          Print version information and then exit\n"
	"\n");
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	cplist *infs;
	FILE *out;
	char *outf;
	int c, ncpu, overwrite, reglen, ksize, Z, W, C, M, X, I, D, E, H, L, XX, OO, EE;
	int candidate_mode, cns_model, corr_struct;
	f4i pM, pX, pI, pD;
	ncpu = 1;
	reglen = 1000;
	ksize = 15;
	Z = 4;
	W = 48;
	C = 1;
	M = 2;
	X = -7;
	I = -3;
	D = -4;
	E = -2;
	H = 1;
	L = 10;
	XX = -4;
	OO = -2;
	EE = -1;
	candidate_mode = 0;
	cns_model = 1;
	corr_struct = 1;
	pM = log(0.85);
	pX = log(0.10);
	pI = log(0.03);
	pD = log(0.02);
	infs = init_cplist(4);
	outf = NULL;
	overwrite = 0;
	while((c = getopt(argc, argv, "hvVt:k:i:o:fj:Z:W:C:M:X:I:D:E:H:L:m:c:S:")) != -1){
		switch(c){
			case 'h': return usage();
			case 't': ncpu = atoi(optarg); break;
			case 'k': ksize = atoi(optarg); break;
			case 'i': push_cplist(infs, optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'j': reglen = atoi(optarg); break;
			case 'Z': Z = atoi(optarg); break;
			case 'W': W = atoi(optarg); break;
			case 'C': C = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 'X': X = atoi(optarg); break;
			case 'I': I = atoi(optarg); break;
			case 'D': D = atoi(optarg); break;
			case 'E': E = atoi(optarg); break;
			case 'H': H = atoi(optarg); break;
			case 'L': L = atoi(optarg); break;
			case 'm': cns_model = atoi(optarg); break;
			case 'c': candidate_mode = atoi(optarg); break;
			case 'S': corr_struct = atoi(optarg); break;
			case 'v': cns_debug ++; break;
			case 'V': fprintf(stdout, "wtdbg-cns 1.1\n"); return 0;
			default: return usage();
		}
	}
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	if(cns_model != MCNS_TASK_DBGCNS && cns_model != MCNS_TASK_DAGCNS){
		return usage();
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
	if(cns_debug > 1) DBGCNS_DEBUG = 1;
	if(infs->size) fr = open_all_filereader(infs->size, infs->buffer, 1);
	else fr = open_filereader(NULL, 1);
	if(outf){
		out = open_file_for_write(outf, NULL, 1);
	} else out = stdout;
	run(reglen, ksize, Z, W, M, X, I, D, E, H, L, XX, OO, EE, cns_model, pM, pX, pI, pD, candidate_mode, corr_struct, ncpu, fr, out);
	close_filereader(fr);
	if(outf) fclose(out);
	free_cplist(infs);
	END_STAT_PROC_INFO(stderr);
	return 0;
}

