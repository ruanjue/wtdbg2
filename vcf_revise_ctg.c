//#define HUGE_STRING
#include "string.h"
#include "list.h"
#include "hashset.h"
#include "bitvec.h"
#include "file_reader.h"

int usage(){
	fprintf(stdout, "Usage: vcf_revise_ctg <ctg_fa_file> <vcf_file>\n");
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	Sequence *seq;
	cuhash *hash;
	String *chrs;
	String *ctgs;
	u32list *chr_offs, *ctg_offs;
	String *cor;
	BitVec *flags;
	char *ref, *alt, *dp, *str, *chr;
	uint32_t n_seq, cid;
	int n, pos, lst, indel, dp4[4], ctg_len, ref_len, alt_len;
	if(argc < 3) return usage();
	hash = init_cuhash(1023);
	chrs = init_string(1024);
	ctgs = init_string(1024);
	chr_offs = init_u32list(1024);
	ctg_offs = init_u32list(1024);
	cor = init_string(1024);
	fr = fopen_filereader(argv[1]);
	seq = NULL;
	n_seq = 0;
	while(fread_seq(&seq, fr)){
		push_u32list(chr_offs, chrs->size);
		push_u32list(ctg_offs, ctgs->size);
		append_string(chrs, seq->name.string, seq->name.size + 1);
		append_string(ctgs, seq->seq.string, seq->seq.size + 1);
		n_seq ++;
	}
	fclose_filereader(fr);
	for(cid=0;cid<n_seq;cid++){
		chr = chrs->string + chr_offs->buffer[cid];
		kv_put_cuhash(hash, chr, cid);
	}
	flags = init_bitvec(n_seq);
	fr = fopen_filereader(argv[2]);
	cid = 0;
	lst = 0;
	chr = chrs->string + chr_offs->buffer[cid];
	str = ctgs->string + ctg_offs->buffer[cid];
	ctg_len = strlen(str);
	clear_string(cor);
	while((n = fread_table(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(n < 8) continue;
		indel = dp4[0] = dp4[1] = dp4[2] = dp4[3] = 0;
		pos = atoi(get_col_str(fr, 1));
		ref = get_col_str(fr, 3);
		alt = get_col_str(fr, 4);
		ref_len = get_col_len(fr, 3);
		for(alt_len=0;alt_len<get_col_len(fr, 4);alt_len++) if(alt[alt_len] == ',') break;
		alt[alt_len] = 0;
		if((dp = strstr(get_col_str(fr, 7), "DP4=")) == NULL) continue;
		if(strcmp(chr, get_col_str(fr, 0))){
			if(lst){
				append_string(cor, str + lst, ctg_len - lst);
				fprintf(stdout, ">%s\n", chr);
				print_pretty_seq(stdout, cor, 100); fflush(stdout);
				one_bitvec(flags, cid);
			}
			cid = kv_get_cuhash(hash, get_col_str(fr, 0));
			if(cid == 0xFFFFFFFFU){
				fprintf(stderr, " -- Cannot find %s in %s -- %s:%d --\n", get_col_str(fr, 0), __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
			}
			lst = 0;
			chr = chrs->string + chr_offs->buffer[cid];
			str = ctgs->string + ctg_offs->buffer[cid];
			ctg_len = strlen(str);
			clear_string(cor);
		}
		dp += 3;
		dp4[0] = strtol(dp + 1, &dp, 10);
		dp4[1] = strtol(dp + 1, &dp, 10);
		dp4[2] = strtol(dp + 1, &dp, 10);
		dp4[3] = strtol(dp + 1, &dp, 10);
		if(dp4[2] + dp4[3] < dp4[0] + dp4[1] + 1) continue;
		if(pos <= lst) continue;
		//fprintf(stderr, "%s\n", fr->line->string);
		//fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", get_col_str(fr, 0), pos, dp4[0], dp4[1], dp4[2], dp4[3], ref, alt);
		//fflush(stderr);
		append_string(cor, str + lst, pos - lst - 1);
		append_string(cor, alt, alt_len);
		lst = pos - 1 + ref_len;
	}
	if(lst){
		append_string(cor, str + lst, ctg_len - lst);
		fprintf(stdout, ">%s\n", chr);
		print_pretty_seq(stdout, cor, 100); fflush(stdout);
		one_bitvec(flags, cid);
	}
	fclose_filereader(fr);
	for(cid=0;cid<n_seq;cid++){
		if(get_bitvec(flags, cid)) continue;
		chr = chrs->string + chr_offs->buffer[cid];
		str = ctgs->string + ctg_offs->buffer[cid];
		ctg_len = strlen(str);
		fprintf(stdout, ">%s\n", chr);
		print_pretty_str(stdout, str, ctg_len, 100); fflush(stdout);
	}
	free_cuhash(hash);
	free_bitvec(flags);
	free_string(chrs);
	free_string(ctgs);
	free_string(cor);
	free_u32list(chr_offs);
	free_u32list(ctg_offs);
	return 0;
}
