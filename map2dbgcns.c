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
#include "file_reader.h"
#include "hashset.h"
#include "list.h"

typedef struct {
	u4i rfid, rfbeg, rfend;
	u4i rdid, dir:1, rdbeg:31, rdend;
} hit_t;

define_list(hitv, hit_t);

typedef struct {
	u8i hidx:63, hdir:1;
} href_t;
define_list(hrefv, href_t);

typedef struct {
	SeqBank *rfs;
	SeqBank *rds;
	hitv    *hits;
	hrefv   *refs;
	u4i     win, step, ext;
} LAY;

LAY* init_lay(u4i win, u4i step, u4i ext){
	LAY *lay;
	lay = malloc(sizeof(LAY));
	lay->rfs = init_seqbank();
	lay->rds = init_seqbank();
	lay->hits = init_hitv(1024);
	lay->refs = init_hrefv(1024);
	lay->win = win;
	lay->step = step;
	lay->ext = ext;
	return lay;
}

void free_lay(LAY *lay){
	free_seqbank(lay->rfs);
	free_seqbank(lay->rds);
	free_hitv(lay->hits);
	free_hrefv(lay->refs);
	free(lay);
}

void load_refseqs_lay(LAY *lay, FileReader *fr){
	Sequence *seq;
	seq = NULL;
	while(fread_seq(&seq, fr)){
		push_seqbank(lay->rfs, seq->tag.string, seq->tag.size, seq->seq.string, seq->seq.size);
	}
}

void load_rdseqs_lay(LAY *lay, FileReader *fr){
	Sequence *seq;
	seq = NULL;
	while(fread_seq(&seq, fr)){
		push_seqbank(lay->rds, seq->tag.string, seq->tag.size, seq->seq.string, seq->seq.size);
	}
}

// ref_tag, ref_beg, ref_end, rd_tag, rd_dir, rd_beg, rd_end
void load_hits_lay(LAY *lay, FileReader *fr){
	hit_t HIT;
	int c;
	while((c = fread_table(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(c < 7) continue;
		if((HIT.rfid = find_seqbank(lay->rfs, get_col_str(fr, 0))) == 0xFFFFFFFFU) continue;
		if((HIT.rdid = find_seqbank(lay->rds, get_col_str(fr, 3))) == 0xFFFFFFFFU) continue;
		HIT.dir = (get_col_str(fr, 4)[0] == '-');
		HIT.rfbeg = atoi(get_col_str(fr, 1));
		HIT.rfend = atoi(get_col_str(fr, 2));
		HIT.rdbeg = atoi(get_col_str(fr, 5));
		HIT.rdend = atoi(get_col_str(fr, 6));
		push_hitv(lay->hits, HIT);
	}
}

void sort_hits_lay(LAY *lay, int ncpu){
	u8i i;
	psort_array(lay->hits->buffer, lay->hits->size, hit_t, ncpu, num_cmpgtx(a.rfid, b.rfid, a.rfbeg, b.rfbeg));
	for(i=0;i<lay->hits->size;i++){
		push_hrefv(lay->refs, (href_t){i, 0});
		push_hrefv(lay->refs, (href_t){i, 1});
	}
	psort_array(lay->refs->buffer, lay->refs->size, href_t, ncpu, num_cmpgtx(lay->hits->buffer[a.hidx].rfid, lay->hits->buffer[b.hidx].rfid, a.hdir? lay->hits->buffer[a.hidx].rfend : lay->hits->buffer[a.hidx].rfbeg, b.hdir? lay->hits->buffer[b.hidx].rfend : lay->hits->buffer[b.hidx].rfbeg));
}

void gen_lay(LAY *lay, FILE *out){
	u8i rfoff, rbeg, rend, ridx;
	hitv *hits;
	href_t *r;
	hit_t  *h;
	u4i rfidx, nidx, i;
	int rflen, beg, end, x, y, z;
	float scale;
	nidx = 0;
	rbeg = 0;
	rend = 0;
	hits = init_hitv(32);
	for(rfidx=0;rfidx<lay->rfs->nseq;rfidx++){
		rfoff = lay->rfs->rdoffs->buffer[rfidx];
		rflen = lay->rfs->rdlens->buffer[rfidx];
		fprintf(out, ">%s len=%d\n", lay->rfs->rdtags->buffer[rfidx], rflen);
		for(beg=0;beg+(int)lay->step<rflen;beg+=lay->step){
			end = beg + lay->win;
			if((int)lay->win > rflen) end = rflen;
			fprintf(out, "E\t%d\tN%d\t+\tN%d\t+\n", beg, nidx, nidx + 1);
			nidx ++;
			fprintf(out, "S\t%s_F_%d_%d\t+\t%d\t%d\t", lay->rfs->rdtags->buffer[rfidx], beg, end - beg, beg, end - beg);
			print_seq_basebank(lay->rfs->rdseqs, rfoff + beg, end - beg, out);
			fprintf(out, "\n");
			while(rbeg < lay->refs->size){
				r = ref_hrefv(lay->refs, rbeg);
				if(lay->hits->buffer[r->hidx].rfid < rfidx){
					rbeg ++;
					continue;
				} else if(lay->hits->buffer[r->hidx].rfid > rfidx){
					break;
				}
				//off = r->hdir? lay->hits->buffer[r->hidx].rfend : lay->hits->buffer[r->hidx].rfbeg;
				//if(off <= beg){
				if((int)lay->hits->buffer[r->hidx].rfend <= beg){
					rbeg ++;
					continue;
				} else {
					break;
				}
			}
			while(rend < lay->refs->size){
				r = ref_hrefv(lay->refs, rend);
				if(lay->hits->buffer[r->hidx].rfid < rfidx){
					rend ++;
					continue;
				} else if(lay->hits->buffer[r->hidx].rfid > rfidx){
					break;
				}
				//off = r->hdir? lay->hits->buffer[r->hidx].rfend : lay->hits->buffer[r->hidx].rfbeg;
				//if(off < end){
				if((int)lay->hits->buffer[r->hidx].rfbeg < end){
					rend ++;
					continue;
				} else {
					break;
				}
			}
			// remove duplicated, caused by two elements(hit->beg:0, hit->end:1) in refs for each hit
			clear_hitv(hits);
			for(ridx=rbeg;ridx<rend;ridx++){
				r = ref_hrefv(lay->refs, ridx);
				h = ref_hitv(lay->hits, r->hidx);
				push_hitv(hits, *h);
			}
			sort_array(hits->buffer, hits->size, hit_t, num_cmpgt(a.rdid, b. rdid));
			for(i=0;i<hits->size;i++){
				h = ref_hitv(hits, i);
				if(i && h->rdid == hits->buffer[i-1].rdid) continue;
				scale = 1.0 * (h->rdend - h->rdbeg) / (h->rfend - h->rfbeg);
				x = (beg - ((int)h->rfbeg)) * scale + ((int)h->rdbeg) - lay->ext;
				y = (end - ((int)h->rfbeg)) * scale + ((int)h->rdbeg) + lay->ext;
				z = lay->rds->rdlens->buffer[h->rdid];
				if(h->dir){
					x = z - x;
					y = z - y;
					swap_tmp(x, y, z);
				}
				if(x < 0) x = 0;
				if(y > z) y = z;
				if(x + (int)(lay->win / 2) >= y) continue;
				fprintf(out, "S\t%s\t%c\t%d\t%d\t", lay->rds->rdtags->buffer[h->rdid], "+-"[h->dir], x, y - x);
				if(h->dir) print_revseq_basebank(lay->rds->rdseqs, lay->rds->rdoffs->buffer[h->rdid] + x, y - x, out);
				else print_seq_basebank(lay->rds->rdseqs, lay->rds->rdoffs->buffer[h->rdid] + x, y - x, out);
				fprintf(out, "\n");
			}
		}
	}
}

int usage(){
	printf(
	"MAP2DBGCNS: Prepare layout file for wtdbg-cns from other aligner\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0\n"
	"Usage: map2dbgcns [options] <ctg_file> <rds_file> <map_file> >dbgcns.lay\n"
	"Options:\n"
	" -w <int>     window size, [2000]\n"
	" -s <int>     sliding size, [1000]\n"
	" -e <int>     extending size on reads, [200]\n"
	" -t <int>     number of threads for sorting, [8]\n"
	"Format map_file:\n"
	"<ref_name>\\t<ref_beg>\\t<ref_end>\\t<rd_name>\\t<rd_strand:+/->\\t<rd_beg>\\t<rd_end>\n"
	" all coordinates are 0-based, and referred to forward strand\n"
	);
	return 1;
}

int main(int argc, char **argv){
	LAY *lay;
	FileReader *ctgf, *rdsf, *mapf;
	int c, win, step, ext, ncpu;
	BEG_STAT_PROC_INFO(stderr, argc, argv);
	win = 2000;
	step = 1000;
	ext = 200;
	ncpu = 8;
	while((c = getopt(argc, argv, "hw:s:e:t:")) != -1){
		switch(c){
			case 'w': win = atoi(optarg); break;
			case 's': step = atoi(optarg); break;
			case 'e': ext = atoi(optarg); break;
			case 't': ncpu = atoi(optarg); break;
			default: return usage();
		}
	}
	if(optind + 3 > argc) return usage();
	ctgf = fopen_filereader(argv[optind + 0]);
	rdsf = fopen_filereader(argv[optind + 1]);
	mapf = fopen_filereader(argv[optind + 2]);
	lay = init_lay(win, step, ext);
	fprintf(stderr, "[%s] loading contigs, ", date()); fflush(stderr);
	load_refseqs_lay(lay, ctgf);
	fprintf(stderr, "%u\n", lay->rfs->nseq); fflush(stderr);
	fclose_filereader(ctgf);
	fprintf(stderr, "[%s] loading reads, ", date()); fflush(stderr);
	load_rdseqs_lay(lay, rdsf);
	fprintf(stderr, "%u\n", lay->rds->nseq); fflush(stderr);
	fclose_filereader(rdsf);
	fprintf(stderr, "[%s] loading hits, ", date()); fflush(stderr);
	load_hits_lay(lay, mapf);
	fprintf(stderr, "%llu\n", (u8i)lay->hits->size); fflush(stderr);
	fclose_filereader(mapf);
	fprintf(stderr, "[%s] sorting hits\n", date());
	sort_hits_lay(lay, ncpu);
	fprintf(stderr, "[%s] begin to output\n", date());
	gen_lay(lay, stdout);
	fprintf(stderr, "[%s] Done\n", date());
	free_lay(lay);
	END_STAT_PROC_INFO(stderr);
	return 0;
}

