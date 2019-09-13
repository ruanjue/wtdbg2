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

#include "kbm.h"
#include "kbmpoa.h"
#include <regex.h>

#ifndef VERSION
#define VERSION 0.0
#endif
#ifndef RELEASE
#define RELEASE 19830203
#endif

int kbm_usage(){
	fprintf(stdout, "Program: kbm is a simple instance which implemented kmer-binmap\n");
	fprintf(stdout, "         it maps query sequence against reference by kmer matching\n");
	fprintf(stdout, "         matched kmer-pairs are bined (256bp) and counted in a matrix\n");
	fprintf(stdout, "         dynamic programming is used to search the best path\n");
	fprintf(stdout, "Version: %s (%s)\n", TOSTR(VERSION), TOSTR(RELEASE));
	fprintf(stdout, "Author: Jue Ruan <ruanjue@gmail.com>\n");
	fprintf(stdout, "Usage: kbm <options> [start|list|stop]\n");
	fprintf(stdout, "Options:\n");
	fprintf(stdout, " -i <string> File(s) of query sequences, +, [STDIN]\n");
	fprintf(stdout, " -d <string> File(s) of reference sequences, +, [<-i>]\n");
	fprintf(stdout, " -L <int>    Choose the longest subread and drop reads shorter than <int> (5000 recommended for PacBio) [0]\n");
	fprintf(stdout, "             Negative integer indicate keeping read names, e.g. -5000.\n");
	fprintf(stdout, " -o <string> Output file, [STDOUT]\n");
	fprintf(stdout, " -I          Interactive mode\n");
	fprintf(stdout, "             e.g. `mkfifo pipe` then `while true; do cat pipe && sleep 1; done | kbm -t 8 -I -d ref.fa -i - -Hk 21 -S 4`\n");
	fprintf(stdout, "             then `cat 1.fq >pipe; cat 2.fq >pipe`, fastq format is better in interaction\n");
	fprintf(stdout, " -f          Force overwrite\n");
	fprintf(stdout, " -t <int>    Number of threads, 0: all cores, [1]\n");
	fprintf(stdout, " -k <int>    Kmer-f size, <= %d, [0]\n", KBM_MAX_KSIZE);
	fprintf(stdout, " -p <int>    Kmer-p size, <= %d, [21]\n", KBM_MAX_KSIZE);
	fprintf(stdout, " -K <float>  Filter high frequency kmers, maybe repetitive, [1000]\n");
	fprintf(stdout, "             if K >= 1, take the integer value as cutoff, MUST <= 65535\n");
	fprintf(stdout, "             else, mask the top fraction part high frequency kmers\n");
	fprintf(stdout, " -E <int>    Min kmer frequency, [1]\n");
	fprintf(stdout, " -O <int>    Filter low complexity bins (#indexed_kmer less than <-O>), [2]\n");
	fprintf(stdout, " -S <float>  Subsampling kmers, 1/(<-S>) kmers are indexed, [4.00]\n");
	fprintf(stdout, "             -S is very useful in saving memeory and speeding up\n");
	fprintf(stdout, "             please note that subsampling kmers will have less matched length\n");
	fprintf(stdout, " -B <int>    Select no more than n seeds in a query bin, [256]\n");
	// Obsolete
	//fprintf(stdout, " -G <int>    Recognize error kmers in a bin when be aligned >= <-G> times, [0]\n");
	fprintf(stdout, "             If you are using shared kbmidx by other process using -D too, it will bring wrong behavior\n");
	fprintf(stdout, " -D <int>    Strand of alignment, 1: forward, 2: reverse, 3: both, [3]\n");
	fprintf(stdout, " -X <int>    Max number of bin(256bp) in one gap, [4]\n");
	fprintf(stdout, " -Y <int>    Max number of bin(256bp) in one deviation, [4]\n");
	fprintf(stdout, " -Z <float>  Max fraction of gapped BINs / aligned BINs, [0.6]\n");
	fprintf(stdout, " -x <int>    penalty for BIN gap, [-7]\n");
	fprintf(stdout, " -y <int>    penalty for BIN deviation, [-21]\n");
	fprintf(stdout, " -z <int>    Enable refine alignment with -p <-z> [0]\n");
	fprintf(stdout, " -l <int>    Min alignment length, [2048]\n");
	fprintf(stdout, " -m <int>    Min matched length, [200]\n");
	fprintf(stdout, " -s <float>  Min similarity, calculated by kmer matched length / aligned length, [0.05]\n");
	fprintf(stdout, " -r <float>  Max length variation of two aligned fragments, [0.25]\n");
	fprintf(stdout, " -c          Insist to query contained reads against all\n");
	fprintf(stdout, " -C          Chainning alignments\n");
	fprintf(stdout, " -n <int>    Max hits per query, [1000]\n");
#ifdef TEST_MODE
	fprintf(stdout, " -T <int>    For debug, [0]\n");
#endif
	fprintf(stdout, " -W <string> Dump kbm index to file, [NULL]\n");
	fprintf(stdout, " -R <string> Load kbm index from file, [NULL]\n");
	fprintf(stdout, " -q          Quiet\n");
	fprintf(stdout, " -V          Print version information and then exit\n");
#if __DEBUG__
	fprintf(stdout, " -v          Verbose, +\n");
#endif
	fprintf(stdout, "Server start: {kbm -R <wt.fa.kbmidx> start}, will mmap wt.fa.kbmidx into mmeory\n");
	fprintf(stdout, "Server  list: {kbm -R <wt.fa.kbmidx> list [10]}, will list the object tree in file\n");
	fprintf(stdout, "Server  stop: {kbm -R <wt.fa.kbmidx> stop},  will remove the mmap object\n");
	return 1;
}

thread_beg_def(maln);
CTGCNS *cc;
KBMAux *aux;
String *rdtag;
BaseBank *rdseqs;
u4i qidx;
u8i rdoff;
u4i rdlen;
int corr_mode;
float corr_cov;
u4i corr_min, corr_max;
FILE *out, *lay;
int chainning;
int interactive;
int refine;
thread_end_def(maln);

thread_beg_func(maln);
KBMPar *rpar;
KBM *rkbm;
KBMAux *raux;
kbm_map_t HIT;
u4v *tidxs;
{
	rpar = init_kbmpar();
	rpar->ksize = 0;
	rpar->psize = maln->refine;
	rpar->min_bin_degree = 0;
	rpar->kmin = 1;
	rpar->kmax = 1000;
	rpar->kmer_mod = KBM_N_HASH;
	rkbm = init_kbm(rpar);
	raux = init_kbmaux(rkbm);
}
tidxs = init_u4v(16);
thread_beg_loop(maln);
if(maln->rdlen == 0) break;
if(maln->corr_mode){
	if(map_kbmpoa(maln->cc, maln->aux, maln->rdtag->size? maln->rdtag->string : NULL, maln->qidx, maln->rdseqs, maln->rdoff, maln->rdlen, maln->corr_min, maln->corr_max, maln->corr_cov, maln->lay) == 0){
		clear_kbmmapv(maln->aux->hits);
		break;
	}
} else {
	query_index_kbm(maln->aux, maln->rdtag->size? maln->rdtag->string : NULL, maln->qidx, maln->rdseqs, maln->rdoff, maln->rdlen);
	map_kbm(maln->aux);
	if(maln->refine && maln->aux->hits->size){
		kbm_read_t *rd;
		kbm_map_t *hit;
		u4i i, j, tidx;
		clear_kbm(rkbm);
		bitpush_kbm(rkbm, maln->rdtag->size? maln->rdtag->string : NULL, maln->rdtag->size, maln->rdseqs->bits, 0, maln->rdlen);
		ready_kbm(rkbm);
		simple_index_kbm(rkbm, 0, rkbm->bins->size);
		clear_u4v(tidxs);
		for(i=0;i<maln->aux->hits->size;i++){
			hit = ref_kbmmapv(maln->aux->hits, i);
			if(tidxs->size == 0 || hit->tidx != tidxs->buffer[tidxs->size - 1]){
				push_u4v(tidxs, hit->tidx);
			}
			if(KBM_LOG){
				fprintf(maln->out, "#");
				fprint_hit_kbm(maln->aux, i, maln->out);
			}
		}
		clear_kbmmapv(maln->aux->hits);
		clear_bitsvec(maln->aux->cigars);
		for(i=0;i<tidxs->size;i++){
			tidx = get_u4v(tidxs, i);
			rd = ref_kbmreadv(maln->aux->kbm->reads, tidx);
			query_index_kbm(raux, rd->tag, tidx, maln->aux->kbm->rdseqs, rd->seqoff * KBM_BIN_SIZE, rd->bincnt * KBM_BIN_SIZE);
			map_kbm(raux);
			for(j=0;j<raux->hits->size;j++){
				flip_hit_kbmaux(maln->aux, raux, j);
			}
		}
	}
}
if(maln->chainning){
	u4i idx, lst;
	for(idx=lst=0;idx<=maln->aux->hits->size;idx++){
		if(idx == maln->aux->hits->size || maln->aux->hits->buffer[lst].tidx != maln->aux->hits->buffer[idx].tidx || maln->aux->hits->buffer[idx].qdir != maln->aux->hits->buffer[lst].qdir){
			if(idx > lst + 1){
				if(simple_chain_all_maps_kbm(maln->aux->hits->buffer + lst, idx - lst, maln->aux->cigars, &HIT, maln->aux->cigars, maln->aux->par->aln_var)){
					maln->aux->hits->buffer[lst++] = HIT;
					while(lst < idx){
						maln->aux->hits->buffer[lst++].mat = 0;
					}
				}
			}
			lst = idx;
		}
	}
}
if(maln->aux->par->max_hit){
	sort_array(maln->aux->hits->buffer, maln->aux->hits->size, kbm_map_t, num_cmpgt(b.mat, a.mat));
	if(maln->aux->hits->size > maln->aux->par->max_hit) maln->aux->hits->size = maln->aux->par->max_hit;
}
if(maln->interactive){
	u4i i;
	thread_beg_syn(maln);
	for(i=0;i<maln->aux->hits->size;i++){
		fprint_hit_kbm(maln->aux, i, maln->out);
	}
	fflush(maln->out);
	thread_end_syn(maln);
}
thread_end_loop(maln);
{
	free_kbmaux(raux);
	free_kbm(rkbm);
	free_kbmpar(rpar);
}
free_u4v(tidxs);
thread_end_func(maln);

int kbm_main(int argc, char **argv){
	cplist *qrys, *refs;
	char *outf, *loadf, *dumpf;
	FILE *out, *dump;
	KBM *kbm;
	KBMPar *par;
	KBMAux *aux;
	BitVec *solids, *rdflags;
	FileReader *fr;
	BioSequence *seqs[2], *seq;
	char regtag[14];
	u8i tot_bp, max_bp, opt_flags, nhit;
	u4i qidx, i;
	int c, ncpu, buffered_read, overwrite, quiet, tidy_reads, tidy_rdtag, skip_ctn, chainning, interactive, server, tree_maxcnt;
	int solid_kmer, refine;
	float fval;
	thread_preprocess(maln);
	par = init_kbmpar();
	par->rd_len_order = 1;
	chainning = 0;
	KBM_LOG = 0;
	buffered_read = 1;
	skip_ctn = 1;
	interactive = 0;
	solid_kmer = 0;
	refine = 0;
	solids = NULL;
	rdflags = NULL;
	ncpu = 1;
	quiet = 0;
	tidy_reads = 0;
	tidy_rdtag = -1;
	max_bp = MAX_U8;
	qrys = init_cplist(4);
	refs = init_cplist(4);
	outf = NULL;
	loadf = NULL;
	dumpf = NULL;
	overwrite = 0;
	server = 0;
	tree_maxcnt = 10;
	opt_flags = 0;
	while((c = getopt(argc, argv, "hi:d:o:fIt:k:p:K:E:O:S:B:G:D:X:Y:Z:x:y:z:l:m:n:s:cr:CT:W:R:qvV")) != -1){
		switch(c){
			case 'h': return kbm_usage();
			case 'i': push_cplist(qrys, optarg); break;
			case 'd': push_cplist(refs, optarg); break;
			case 'L': tidy_reads = atoi(optarg); break;
			case 'o': outf = optarg; break;
			case 'f': overwrite = 1; break;
			case 'I': interactive = 1; break;
			case 't': ncpu = atoi(optarg); break;
			case 'k': par->ksize = atoi(optarg); opt_flags |= (1 << 1); break;
			case 'p': par->psize = atoi(optarg); opt_flags |= (1 << 0); break;
			case 'K': fval = atof(optarg); par->kmax = fval; par->ktop = fval - par->kmax; break;
			case 'E': par->kmin = atoi(optarg); break;
			case 'O': par->min_bin_degree = atoi(optarg); break;
			case 'S': par->kmer_mod = UInt(atof(optarg) * KBM_N_HASH); opt_flags |= (1 << 2); break;
			case 'B': par->ksampling = atoi(optarg); break;
			case 'G': solid_kmer = atoi(optarg); break;
			case 'D': par->strand_mask = atoi(optarg); break;
			case 'X': par->max_bgap = atoi(optarg); break;
			case 'Y': par->max_bvar = atoi(optarg); break;
			case 'Z': par->max_gap  = atof(optarg); break;
			case 'x': par->pgap = atoi(optarg); break;
			case 'y': par->pvar = atoi(optarg); break;
			case 'z': refine = atoi(optarg); break;
			case 'l': par->min_aln = atoi(optarg) / KBM_BIN_SIZE; break;
			case 'm': par->min_mat = atoi(optarg); break;
			case 'n': par->max_hit = atoi(optarg); break;
			case 's': par->min_sim = atof(optarg); break;
			case 'r': par->aln_var = atof(optarg); break;
			case 'c': skip_ctn = 0; break;
			case 'C': chainning = 1; break;
#ifdef TEST_MODE
			case 'T': par->test_mode = atoi(optarg); break;
#endif
			case 'W': dumpf = optarg; break;
			case 'R': loadf = optarg; break;
			case 'q': quiet = 1; break;
			case 'v': KBM_LOG ++; break;
			case 'V': fprintf(stdout, "kbm2 %s\n", TOSTR(VERSION)); return 0;
			default: return kbm_usage();
		}
	}
	if(quiet){
		int devnull;
		devnull = open("/dev/null", O_WRONLY);
		dup2(devnull, KBM_LOGFNO);
	}
	if(tidy_rdtag == -1){
		tidy_rdtag = (tidy_reads >= 0);
	}
	if(tidy_reads < 0) tidy_reads = - tidy_reads;
	BEG_STAT_PROC_INFO(KBM_LOGF, argc, argv);
	if(par->ksize + par->psize > KBM_MAX_KSIZE){
		fprintf(stderr, " -- Invalid kmer size %d+%d=%d > %d in %s -- %s:%d --\n", par->ksize, par->psize, par->ksize + par->psize,  KBM_MAX_KSIZE, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return 1;
	}
	out = open_file_for_write(outf, NULL, overwrite);
	if(refs->size == 0){
		if(loadf){
		} else if(qrys->size){
			append_cplist(refs, qrys);
			clear_cplist(qrys);
		} else {
			fprintf(KBM_LOGF, " ** -h to show document\n"); fflush(KBM_LOGF);
			push_cplist(refs, "-");
		}
	}
	if(qrys->size == 0){
		par->self_aln = 1;
		interactive = 0;
		if(par->kmin < 2) par->kmin = 2;
	}
	if(interactive) buffered_read = 0;
	if(ncpu <= 0){
		ncpu = _proc_deamon->ncpu;
		if(ncpu == 0) ncpu = 1; // failed to get number of cores
	}
	if(KBM_LOG > 0){
		fprintf(KBM_LOGF, "KBM_LOG_LEVEL = %d\n", KBM_LOG);
	}
	if(optind < argc){
		server = 0;
		if(strcasecmp("start", argv[optind]) == 0) server = 1;
		else if(strcasecmp("stop", argv[optind]) == 0) server = 2;
		else if(strcasecmp("list", argv[optind]) == 0){
			server = 3;
			if(optind + 1 < argc) tree_maxcnt = atoi(argv[optind + 1]);
		}
		if(loadf == NULL) server = 0;
	}
	if(loadf){
		if(server == 1){
			fprintf(KBM_LOGF, "[%s] loading kbm index from %s\n", date(), loadf);
			kbm = mem_load_obj_file(&kbm_obj_desc, loadf, NULL, NULL, NULL, NULL);
			fprintf(KBM_LOGF, "[%s] Done. %u sequences, %llu bp, parameter('-k %d -p %d -S %d')\n", date(), (u4i)kbm->reads->size, (u8i)kbm->rdseqs->size, kbm->par->ksize, kbm->par->psize, kbm->par->kmer_mod / KBM_N_HASH);
			fprintf(KBM_LOGF, "[%s] kbm-index server start\n", date());
			return 0;
		} else if(server == 2){
			if(mem_stop_obj_file(loadf)){
				fprintf(KBM_LOGF, "[%s] kbm-index server for '%s' stop\n", date(), loadf);
			} else {
				fprintf(KBM_LOGF, "[%s] unable to find kbm-index server for '%s'\n", date(), loadf);
			}
			return 0;
		} else if(server == 3){
			print_tree_obj_file(stdout, &kbm_obj_desc, loadf, tree_maxcnt, 0);
			return 0;
		} else {
			fprintf(KBM_LOGF, "[%s] loading kbm index from %s\n", date(), loadf);
			if((kbm = mem_find_obj_file(&kbm_obj_desc, loadf, NULL, NULL, NULL, NULL, 1)) == NULL){
				fprintf(KBM_LOGF, " -- cannot find mmap object %s --\n", loadf);
				fprintf(KBM_LOGF, " -- try read from file --\n");
				kbm = mem_read_obj_file(&kbm_obj_desc, loadf, NULL, NULL, NULL, NULL);
			}
		}
		fprintf(KBM_LOGF, "[%s] Done. %u sequences, %llu bp, parameter('-k %d -p %d -S %d')\n", date(), (u4i)kbm->reads->size, (u8i)kbm->rdseqs->size, kbm->par->ksize, kbm->par->psize, kbm->par->kmer_mod / KBM_N_HASH);
		// Please note that, kbm->tag2idx is not functional after mem_load
		// check KBMPar
		if((opt_flags >> 0) & 0x01){
			if(kbm->par->psize != par->psize){
				fprintf(KBM_LOGF, " ** -p is different, %d != %d\n", kbm->par->psize, par->psize); exit(1);
			}
		} else {
			par->psize = kbm->par->psize;
		}
		if((opt_flags >> 1) & 0x01){
			if(kbm->par->ksize != par->ksize){
				fprintf(KBM_LOGF, " ** -k is different, %d != %d\n", kbm->par->ksize, par->ksize); exit(1);
			}
		} else {
			par->ksize = kbm->par->ksize;
		}
		if((opt_flags >> 2) & 0x01){
			if(kbm->par->kmer_mod != par->kmer_mod){
				fprintf(KBM_LOGF, " ** -S is different, %d != %d\n", kbm->par->kmer_mod / KBM_N_HASH, par->kmer_mod / KBM_N_HASH); exit(1);
			}
		} else {
			par->kmer_mod = kbm->par->kmer_mod;
		}
		if((opt_flags >> 3) & 0x01){
			if(kbm->par->rd_len_order != par->rd_len_order){
				fprintf(KBM_LOGF, " ** par->rd_len_order is different, %d != %d\n", kbm->par->rd_len_order, par->rd_len_order); exit(1);
			}
		} else {
			par->rd_len_order = kbm->par->rd_len_order;
		}
	} else {
		kbm = init_kbm(par);
		fprintf(KBM_LOGF, "[%s] loading sequences\n", date()); fflush(KBM_LOGF);
		fr = open_all_filereader(refs->size, refs->buffer, buffered_read);
		tot_bp = 0;
		seqs[0] = init_biosequence();
		seqs[1] = init_biosequence();
		regex_t reg;
		regmatch_t mats[3];
		int z, tag_size, len;
		z = regcomp(&reg, "^(.+?)/[0-9]+_[0-9]+$", REG_EXTENDED);
		if(z){
			regerror(z, &reg, regtag, 13);
			fprintf(stderr, " -- REGCOMP: %s --\n", regtag); fflush(stderr);
			return 1;
		}
		{
			int k = 0;
			reset_biosequence(seqs[0]);
			reset_biosequence(seqs[1]);
			while(1){
				int has = readseq_filereader(fr, seqs[k]);
				if(tidy_reads){
					if(has){
						if((z = regexec(&reg, seqs[k]->tag->string, 3, mats, 0)) == 0){
							trunc_string(seqs[k]->tag, mats[1].rm_eo);
						} else if(z != REG_NOMATCH){
							regerror(z, &reg, regtag, 13);
							fprintf(stderr, " -- REGEXEC: %s --\n", regtag); fflush(stderr);
						}
						//fprintf(stderr, "1: %s len=%d\n", seqs[k]->tag->string, seqs[k]->seq->size); fflush(stderr);
						//fprintf(stderr, "2: %s len=%d\n", seqs[!k]->tag->string, seqs[!k]->seq->size); fflush(stderr);
						if(seqs[k]->tag->size == seqs[!k]->tag->size && strcmp(seqs[k]->tag->string, seqs[!k]->tag->string) == 0){
							if(seqs[k]->seq->size > seqs[!k]->seq->size){
								k = !k;
							}
							continue;
						} else {
							seq = seqs[!k];
							k = !k;
						}
					} else {
						seq = seqs[!k];
					}
					if(seq->seq->size < tidy_reads){
						if(has) continue;
						else break;
					}
					if(tidy_rdtag){
						sprintf(regtag, "S%010llu", (u8i)kbm->reads->size);
						clear_string(seq->tag);
						append_string(seq->tag, regtag, 11);
					}
				} else {
					if(has == 0) break;
					seq = seqs[k];
				}
				tag_size = seq->tag->size;
				for(i=0;i<UInt(seq->seq->size);i+=KBM_MAX_RDLEN){
					len = num_min(seq->seq->size - i, KBM_MAX_RDLEN);
					if(i){
						append_string(seq->tag, "_V", 2);
						add_int_string(seq->tag, i / KBM_MAX_RDLEN);
					}
					if(!KBM_LOG && (kbm->reads->size % 10000) == 0){ fprintf(KBM_LOGF, "\r%u", (u4i)kbm->reads->size); fflush(KBM_LOGF); }
					//fprintf(stderr, " -- %s len=%d in %s -- %s:%d --\n", seq->tag->string, seq->seq->size, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					if(kbm->reads->size >= KBM_MAX_RDCNT){
						fprintf(stderr, " -- Read Number Out of Range: %u --\n", (u4i)kbm->reads->size); fflush(stderr);
						break;
					}
					push_kbm(kbm, seq->tag->string, seq->tag->size, seq->seq->string + i, len);
					if(i){ seq->tag->size = tag_size; seq->tag->string[tag_size] = '\0'; }
				}
				tot_bp += seq->seq->size;
				if(max_bp && tot_bp >= max_bp){ break; }
				if(has == 0) break;
				if(kbm->reads->size >= KBM_MAX_RDCNT){
					fprintf(stderr, " -- Read Number Out of Range: %u --\n", (u4i)kbm->reads->size); fflush(stderr);
					break;
				}
			}
		}
		close_filereader(fr);
		regfree(&reg);
		free_biosequence(seqs[0]);
		free_biosequence(seqs[1]);
		if(!KBM_LOG){ fprintf(KBM_LOGF, "\r%u reads", (unsigned)kbm->reads->size); fflush(KBM_LOGF); }
		ready_kbm(kbm);
		fprintf(KBM_LOGF, "\n[%s] Done, %u reads, %llu bp, %u bins\n", date(), (u4i)kbm->reads->size, tot_bp, (u4i)kbm->bins->size); fflush(KBM_LOGF);
		fprintf(KBM_LOGF, "[%s] indexing, %d threads\n", date(), ncpu);
		index_kbm(kbm, 0, kbm->bins->size, ncpu, NULL);
		if(dumpf){
			u8i size;
			dump = open_file_for_write(dumpf, NULL, 1);
			size = mem_dump_obj_file(kbm, 1, &kbm_obj_desc, 1, 0, dump);
			fclose(dump);
			fprintf(KBM_LOGF, "[%s] kbm index dumped to %s\n", date(), dumpf);
		}
	}
	fprintf(KBM_LOGF, "[%s] mapping\n", date());
	if(solid_kmer && par->self_aln){
		solids = init_bitvec((kbm->rdseqs->size >> 1) + 1);
	}
	thread_beg_init(maln, ncpu);
	maln->aux = init_kbmaux(kbm);
	maln->aux->par = par; // par might be different from kbm->par
	maln->aux->solids = solids;
	//maln->aux->bmin = 611077674;
	//maln->aux->bmax = 611077674 + 172;
	maln->rdtag = init_string(64);
	maln->rdseqs = qrys->size? init_basebank() : kbm->rdseqs;
	maln->rdoff = 0;
	maln->rdlen = 0;
	maln->cc = NULL;
	maln->corr_mode = 0;
	maln->corr_cov = 0.75;
	maln->corr_min = 5;
	maln->corr_max = 10;
	maln->out = out;
	maln->lay = NULL;
	maln->chainning = chainning;
	maln->interactive = interactive;
	maln->refine = refine;
	thread_end_init(maln);
	if(qrys->size){
		int run_mode;
		fr = open_all_filereader(qrys->size, qrys->buffer, buffered_read);
		run_mode = 0;
		if(readline_filereader(fr) > 0){
			if(strcasecmp(fr->line->string, "#pairwise_test") == 0){
				run_mode = 1;
			} else if(strcasecmp(fr->line->string, "#print_exists") == 0){
				run_mode = 2;
			} else if(strcasecmp(fr->line->string, "#correct_align") == 0){
				run_mode = 3;
				cns_debug = KBM_LOG;
				thread_beg_iter(maln);
				maln->corr_mode = 1;
				KBMBlock *kb;
				POGPar par;
				kb = init_kbmblock(2048, 2048 - 512);
				par = DEFAULT_POG_PAR;
				//maln->cc = init_ctgcns(kb, iter_kbmblock, info_kbmblock, 1, 1, 1, maln->corr_max, 200, 100, 1, 96, 2, -5, -2, -4, -1, 16, 3, 0.5, 512);
				par.refmode = 1;
				maln->cc = init_ctgcns(kb, iter_kbmblock, info_kbmblock, 1, 1, maln->corr_max, 200, 100, 1, 512, &par);
				maln->lay = maln->out;
				thread_end_iter(maln);
			} else if(strcasecmp(fr->line->string, "#print_read") == 0){
				run_mode = 4;
			} else {
				rollback_filereader(fr);
			}
		}
		seq = init_biosequence();
		qidx = 0;
		nhit = 0;
		if(run_mode == 0 || run_mode == 3){
			while(readseq_filereader(fr, seq)){
				if((qidx % 100) == 0){
					fprintf(KBM_LOGF, "\r%u\t%llu", qidx, nhit); fflush(KBM_LOGF);
				}
				thread_wait_one(maln);
				if(maln->rdlen && !maln->interactive){
					{
						aux = maln->aux;
						for(i=0;i<aux->hits->size;i++){
							fprint_hit_kbm(aux, i, out);
						}
						if(run_mode == 3) fflush(out);
						nhit += aux->hits->size;
					}
				}
				trunc_string(seq->seq, kbm_cvt_length(seq->seq->size));
				clear_basebank(maln->rdseqs);
				seq2basebank(maln->rdseqs, seq->seq->string, seq->seq->size);
				clear_string(maln->rdtag);
				append_string(maln->rdtag, seq->tag->string, seq->tag->size);
				maln->qidx = qidx ++;
				maln->rdoff = 0;
				maln->rdlen = seq->seq->size;
				thread_wake(maln);
			}
		} else if(run_mode == 1){
			int nc;
			u4i tidx;
			thread_beg_operate(maln, 0);
			aux = maln->aux;
			free_basebank(maln->rdseqs);
			maln->rdseqs = kbm->rdseqs;
			while((nc = readtable_filereader(fr)) >= 0){
				if(nc < 2) continue;
				qidx = getval_cuhash(kbm->tag2idx, get_col_str(fr, 0));
				tidx = getval_cuhash(kbm->tag2idx, get_col_str(fr, 1));
				if(qidx == MAX_U4 || tidx == MAX_U4) continue;
				if(qidx > tidx){ swap_var(qidx, tidx); }
				maln->qidx = qidx;
				maln->rdoff = kbm->reads->buffer[qidx].seqoff * KBM_BIN_SIZE;
				maln->rdlen = kbm->reads->buffer[qidx].bincnt * KBM_BIN_SIZE;
				aux->bmin = kbm->reads->buffer[tidx].binoff;
				aux->bmax = kbm->reads->buffer[tidx].binoff + kbm->reads->buffer[tidx].bincnt;
				fprintf(out, "%s <-> %s\n", kbm->reads->buffer[qidx].tag, kbm->reads->buffer[tidx].tag);
				thread_wake(maln);
				thread_wait(maln);
				if(maln->rdlen && !maln->interactive){
					for(i=0;i<aux->hits->size;i++){
						fprint_hit_kbm(aux, i, out);
					}
					nhit += aux->hits->size;
				}
				fprintf(out, "//\n");
				maln->rdlen = 0;
			}
		} else if(run_mode == 2){
			kmeroffv *kmers[2];
			int nc;
			kmers[0] = adv_init_kmeroffv(32, 0, 1);
			kmers[1] = adv_init_kmeroffv(32, 0, 1);
			while((nc = readtable_filereader(fr)) >= 0){
				//fprintf(stdout, "%s\n", fr->line->string);
				if(nc < 1) continue;
				qidx = getval_cuhash(kbm->tag2idx, get_col_str(fr, 0));
				if(qidx == MAX_U4) continue;
				print_exists_index_kbm(kbm, kbm->reads->buffer[qidx].tag, kbm->rdseqs, kbm->reads->buffer[qidx].seqoff * KBM_BIN_SIZE,
					kbm->reads->buffer[qidx].bincnt * KBM_BIN_SIZE, kmers, stdout);
			}
			free_kmeroffv(kmers[0]);
			free_kmeroffv(kmers[1]);
		} else if(run_mode == 4){
			int nc;
			while((nc = readtable_filereader(fr)) >= 0){
				if(nc < 1) continue;
				if(fr->line->string[0] == '#') continue;
				qidx = getval_cuhash(kbm->tag2idx, get_col_str(fr, 0));
				if(qidx == MAX_U4){
					fprintf(out, "#%s NOT FOUND\n", get_col_str(fr, 0));
					continue;
				}
				fprintf(out, ">%s\n", kbm->reads->buffer[qidx].tag);
				println_fwdseq_basebank(kbm->rdseqs, kbm->reads->buffer[qidx].seqoff * KBM_BIN_SIZE, kbm->reads->buffer[qidx].bincnt * KBM_BIN_SIZE, out);
				fflush(out);
			}
		}
		free_filereader(fr);
		free_biosequence(seq);
		thread_beg_iter(maln);
		thread_wait(maln);
		if(maln->rdlen && !maln->interactive){
			aux = maln->aux;
			for(i=0;i<aux->hits->size;i++){
				fprint_hit_kbm(aux, i, out);
			}
			if(run_mode == 3) fflush(out);
			nhit += aux->hits->size;
			maln->rdlen = 0;
		}
		thread_end_iter(maln);
		fprintf(KBM_LOGF, "\r%u\t%llu\n", qidx, (u8i)nhit); fflush(KBM_LOGF);
	} else {
		nhit = 0;
		rdflags = init_bitvec(kbm->reads->size);
		for(qidx=0;qidx<kbm->reads->size+ncpu;qidx++){
			if(qidx < kbm->reads->size){
				thread_wait_one(maln);
			} else {
				thread_wait_next(maln);
			}
			if(maln->rdlen){
				aux = maln->aux;
				for(i=0;i<aux->hits->size;i++){
					if(par->test_mode == 0 && skip_ctn){
						// whether reads[tidx] is contained by reads[qidx]
						kbm_map_t *hit;
						int margin = KBM_BIN_SIZE;
						hit = ref_kbmmapv(aux->hits, i);
						if((hit->tb <= margin && hit->te + margin >= (int)kbm->reads->buffer[hit->tidx].bincnt * KBM_BIN_SIZE)
							&& (hit->qb > margin || hit->qe + margin < (int)kbm->reads->buffer[hit->qidx].bincnt * KBM_BIN_SIZE)){
							one_bitvec(rdflags, hit->tidx);
						}
					}
					fprint_hit_kbm(aux, i, out);
				}
				nhit += aux->hits->size;
			}
			if(qidx >= kbm->reads->size || get_bitvec(rdflags, qidx)){
				maln->rdlen = 0;
				continue;
			} else if((qidx % 100) == 0){
				fprintf(KBM_LOGF, "\r%u\t%llu", qidx, nhit); fflush(KBM_LOGF);
			}
			maln->qidx = qidx;
			maln->rdoff = kbm->reads->buffer[qidx].seqoff * KBM_BIN_SIZE;
			maln->rdlen = kbm->reads->buffer[qidx].bincnt * KBM_BIN_SIZE;
			thread_wake(maln);
		}
		fprintf(KBM_LOGF, "\r%u\t%llu\n", qidx, nhit); fflush(KBM_LOGF);
		free_bitvec(rdflags);
	}
	thread_beg_close(maln);
	free_kbmaux(maln->aux);
	free_string(maln->rdtag);
	if(maln->rdseqs && maln->rdseqs != kbm->rdseqs) free_basebank(maln->rdseqs);
	if(maln->cc){
		free_kbmblock((KBMBlock*)maln->cc->obj);
		free_ctgcns(maln->cc);
	}
	thread_end_close(maln);
	if(solids) free_bitvec(solids);
	fprintf(KBM_LOGF, "[%s] Done\n", date());
	if(outf) fclose(out);
	free_cplist(qrys);
	free_cplist(refs);
	free_kbmpar(par);
	if(loadf){
		// DON't free kbm, for it is mem_load object
	} else {
		free_kbm(kbm);
	}
	END_STAT_PROC_INFO(KBM_LOGF);
	return 0;
}

int main(int argc, char **argv){
	return kbm_main(argc, argv);
}

