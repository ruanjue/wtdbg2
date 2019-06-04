
#include "mem_share.h"
#include "filereader.h"
#include "list.h"
#include "hashset.h"

typedef struct {
	u4i stroff, strlen;
	u4i taglen, flag;
	u4i qlen, qb, qe;
	u4i refidx, reflen;
} lr_hit_t;
define_list(lrhitv, lr_hit_t);

int select_best_hit(lrhitv *hits, String *lines, u4i minlen, float mincov, FILE *out){
	lr_hit_t *h1, *h2;
	u4i i, j, pass;
	int x, y, ret;
	sort_array(hits->buffer, hits->size, lr_hit_t, num_cmpgt(b.qe - b.qb, a.qe - a.qb));
	ret = 0;
	for(i=0;i<hits->size;i++){
		h1 = ref_lrhitv(hits, i);
		if(h1->qe - h1->qb < minlen) break;
		if(h1->qe - h1->qb < UInt(mincov * h1->qlen)) break;
		pass = 1;
		for(j=0;j<i;j++){
			h2 = ref_lrhitv(hits, j);
			x = num_max(h1->qb, h2->qb);
			y = num_min(h1->qe, h2->qe);
			if(y - x >= (1 - mincov) * (h1->qe - h1->qb)){
				pass = 0;
				break;
			}
		}
		if(pass){
			fprintf(out, "%s\n", lines->string + h1->stroff);
			ret ++;
		}
	}
	clear_lrhitv(hits);
	clear_string(lines);
	return ret;
}

int usage(char *prog){
	printf(
	"Usage: %s [-h] [-v] [-B:retain secondary aligment] [-l min_map_len:100] [-f min_map_cov:0.70]\n"
	, prog
	);
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	FILE *out;
	cuhash *refs;
	cplist *reftags;
	u4v    *reflens;
	String *lines;
	lrhitv *hits;
	lr_hit_t *hit, HIT;
	char *str, *reftag;
	u4i minlen, i, reflen;
	float mincov;
	int c, primary_hit, verbose;
	u1i movs[256];
	minlen = 100;
	mincov = 0.70;
	primary_hit = 1;
	verbose = 0;
	out = stdout;
	while((c = getopt(argc, argv, "hvBl:f:")) != -1){
		switch(c){
			case 'l': minlen = atoi(optarg); break;
			case 'f': mincov = atof(optarg); break;
			case 'B': primary_hit = 0; break;
			case 'v': verbose = 1; break;
			default: return usage(argv[0]);
		}
	}
	fr = open_filereader(NULL, 1);
	refs = init_cuhash(13);
	reftags = init_cplist(8);
	reflens = init_u4v(8);
	hits = init_lrhitv(4);
	lines = init_string(1024);
	memset(movs, 0, 256);
	movs[(int)'M'] = 0b11;
	movs[(int)'I'] = 0b10;
	movs[(int)'D'] = 0b01;
	movs[(int)'N'] = 0b01;
	movs[(int)'S'] = 0b10;
	movs[(int)'H'] = 0b10;
	movs[(int)'P'] = 0b00;
	movs[(int)'='] = 0b11;
	movs[(int)'X'] = 0b11;
	while((c = readline_filereader(fr))){
		if(fr->line->string[0] == '@'){
			fprintf(out, "%s\n", fr->line->string);
			if(fr->line->string[1] == 'S' && fr->line->string[2] == 'Q'){
				if((c = split_line_filereader(fr, '\t')) > 2){
					reftag = NULL;
					reflen = 0;
					for(i=1;i<3;i++){
						if(get_col_len(fr, i) <= 3){
							continue;
						}
						str = get_col_str(fr, i);
						if(str[0] == 'S' && str[1] == 'N' && str[2] == ':'){
							reftag = strdup(str + 3);
						} else if(str[0] == 'L' && str[1] == 'N' && str[2] == ':'){
							reflen = atol(str + 3);
						}
					}
					if(strlen(reftag) && reflen){
						push_cplist(reftags, reftag);
						push_u4v(reflens, reflen);
						put_cuhash(refs, (cuhash_t){reftag, reftags->size - 1});
					}
				}
			}
		} else {
			hit = &HIT;
			str = index(fr->line->string, '\t');
			if(str == NULL){
				fprintf(stderr, "[WARNNING:too_few_column] %s\n", fr->line->string);
				continue;
			}
			hit->taglen = str - fr->line->string;
			if(hits->size && (hits->buffer[0].taglen != hit->taglen || strncmp(lines->string + hits->buffer[0].stroff, fr->line->string, hit->taglen))){
				select_best_hit(hits, lines, minlen, mincov, out);
			}
			hit->stroff = lines->size;
			hit->strlen = fr->line->size;
			append_string(lines, fr->line->string, fr->line->size);
			add_char_string(lines, '\0');
			if((c = split_line_filereader(fr, '\t')) < 11){
				fprintf(stderr, "[WARNNING:too_few_columns] %s\n", lines->string + hit->stroff);
				continue;
			}
			hit->taglen = get_col_len(fr, 0);
			hit->flag = atol(get_col_str(fr, 1));
			if(primary_hit && (hit->flag & 0x900)){
				continue;
			}
			if(get_col_str(fr, 2)[0] == '*'){
				continue;
			}
			hit->refidx = getval_cuhash(refs, get_col_str(fr, 2));
			if(hit->refidx == MAX_U4){
				fprintf(stderr, "[WARNNING:unknown_refname] %s\n", lines->string + hit->stroff);
				continue;
			}
			hit->reflen = reflens->buffer[hit->refidx];
			hit->qlen = 0;
			u4i tb, te, qb, qe, len, cnt, tmp;
			u4i ln;
			char op;
			tb = atoi(get_col_str(fr, 3));
			te = tb;
			qb = qe = 0;
			len = 0;
			tmp = cnt = 0;
			str = get_col_str(fr, 5); // CIGAR
			ln = 0;
			op = 0;
			while(str[0]){
				if(str[0] >= '0' && str[0] <= '9'){
					ln = ln * 10 + str[0] - '0';
				} else {
					op = movs[(int)str[0]];
					if(op & 0b01){
						te += ln;
						qe += tmp;
						tmp = 0;
						if(cnt == 0){
							qb = qe;
							cnt = 1;
						}
						if(op & 0b10){
							qe += ln;
							len += ln;
						}
					} else if(op & 0b10){
						tmp += ln;
						len += ln;
					}
					ln = 0;
				}
				str ++;
			}
			if(hit->flag & 0x10){
				tmp = len - qb;
				qb = len - qe;
				qe = tmp;
			}
			hit->qlen = len;
			hit->qb = qb;
			hit->qe = qe;
			if(verbose){
				fprintf(out, "#%s\t%u\t+\t%u\t%u\t%s\t%c\t%u\t%u\n", get_col_str(fr, 0), len, qb, qe, get_col_str(fr, 2), "+-"[(hit->flag & 0x10) >> 4], tb, te);
			}
			push_lrhitv(hits, HIT);
		}
	}
	select_best_hit(hits, lines, minlen, mincov, out);
	free_string(lines);
	free_lrhitv(hits);
	free_cuhash(refs);
	for(i=0;i<reftags->size;i++){
		free(reftags->buffer[i]);
	}
	free_cplist(reftags);
	free_u4v(reflens);
	close_filereader(fr);
	return 0;
}
