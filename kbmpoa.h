/*
 * 
 * Copyright (c) 2018, Jue Ruan <ruanjue@gmail.com>
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

#ifndef KMER_BIN_MAP_PO_MSA_CNS_RJ_H
#define KMER_BIN_MAP_PO_MSA_CNS_RJ_H

#include "kbm.h"
#include "wtpoa.h"

typedef struct {
	u4i chridx, cidx, bidx;
	String   *rdtag;
	BaseBank *rdseq;
	u8i       rdoff;
	u4i       rdlen;
	KBMAux   *aux;
	layseqr  *seqs;
	u4v      *heap;
	u2i bsize, bstep; // block size, and slide steps
	u4i sidx, lidx, hidx;
} KBMBlock;

static inline KBMBlock* init_kbmblock(u2i bsize, u2i bstep){
	KBMBlock *kb;
	kb = malloc(sizeof(KBMBlock));
	kb->rdtag = init_string(1024);
	kb->rdseq = NULL;
	kb->rdoff = 0;
	kb->rdlen = 0;
	kb->aux   = NULL;
	kb->seqs  = init_layseqr(32);
	kb->heap  = init_u4v(32);
	bsize = ((bsize + KBM_BIN_SIZE - 1) / KBM_BIN_SIZE) * KBM_BIN_SIZE;
	if(bsize == 0) bsize = KBM_BIN_SIZE;
	bstep = ((bstep + KBM_BIN_SIZE - 1) / KBM_BIN_SIZE) * KBM_BIN_SIZE;
	if(bstep == 0) bstep = KBM_BIN_SIZE;
	kb->bsize = bsize;
	kb->bstep = bstep;
	kb->chridx = 0;
	kb->cidx   = 0;
	kb->bidx   = 0;
	kb->sidx  = MAX_U4;
	kb->lidx  = MAX_U4;
	kb->hidx  = 0;
	return kb;
}

static inline void free_kbmblock(KBMBlock *kb){
	free_string(kb->rdtag);
	free_layseqr(kb->seqs);
	free(kb);
}

static inline void reset_kbmblock(KBMBlock *kb, char *rdtag, u4i chridx, BaseBank *rdseq, u8i rdoff, u4i rdlen, KBMAux *aux){
	lay_seq_t *stop;
	kb->chridx = chridx;
	kb->cidx   = chridx;
	kb->bidx = 0;
	clear_string(kb->rdtag);
	if(rdtag){
		append_string(kb->rdtag, rdtag, strlen(rdtag));
	} else {
		append_string(kb->rdtag, "anonymous", strlen("anonymous"));
	}
	kb->rdseq = rdseq;
	kb->rdoff = rdoff;
	kb->rdlen = rdlen;
	kb->aux   = aux;
	recyc_all_layseqr(kb->seqs);
	clear_u4v(kb->heap);
	kb->sidx = MAX_U4;
	kb->lidx  = MAX_U4;
	sort_array(aux->hits->buffer, aux->hits->size, kbm_map_t, num_cmpgt(a.qb, b.qb));
	kb->hidx  = 0;
	stop = pop_layseqr(kb->seqs);
	stop->chridx = chridx + 1;
	stop->bidx   = 0;
}

static inline lay_seq_t* _push_padding_ref_kbmblock(KBMBlock *kb, lay_seq_t *sq){
	lay_seq_t *st;
	u8i off;
	u4i sqidx, len;
	sqidx = offset_layseqr(kb->seqs, sq);
	while(kb->chridx == sq->chridx && kb->bidx <= sq->bidx){
		st = pop_layseqr(kb->seqs);
		st->chridx = kb->cidx;
		st->bidx  = kb->bidx;
		st->rdidx = 0;
		st->rddir = 0;
		st->rdoff = st->bidx * kb->bstep;
		st->rbeg  = st->rdoff;
		off = kb->rdoff;
		len = kb->rdlen;
		st->rend = num_min(st->rbeg + kb->bsize, len);
		st->rdtag = kb->rdtag->string;
		clear_basebank(st->seq);
		fast_fwdbits2basebank(st->seq, kb->rdseq->bits, off + st->rdoff, st->rend - st->rbeg);
		array_heap_push(kb->heap->buffer, kb->heap->size, kb->heap->cap, u4i, offset_layseqr(kb->seqs, st),
			num_cmpx(ref_layseqr(kb->seqs, a)->bidx, ref_layseqr(kb->seqs, b)->bidx, ref_layseqr(kb->seqs, a)->rdidx, ref_layseqr(kb->seqs, b)->rdidx));
		if(st->rend >= len){
			kb->chridx = MAX_U4;
			kb->bidx = 0;
			break;
		} else {
			kb->bidx ++;
		}
		sq = ref_layseqr(kb->seqs, sqidx);
	}
	sq = ref_layseqr(kb->seqs, sqidx);
	return sq;
}

static inline lay_seq_t* iter_kbmblock(void *obj){
	KBMBlock *kb;
	kbm_map_t *hit;
	lay_seq_t *sc, *sl;
	u8i coff, tsoff;
	u4i off, rdoff, rdlen, nxt, val, len, clen, bt;
	kb = (KBMBlock*)obj;
	if(kb->sidx != MAX_U4){
		recyc_layseqr(kb->seqs, kb->sidx);
		kb->sidx = MAX_U4;
	}
	sc = NULL;
	sl = NULL;
	do {
		if(kb->heap->size == 0){
			kb->lidx = MAX_U4;
		}
		while(kb->lidx == MAX_U4){
			while(kb->hidx < kb->aux->hits->size){
				if(kb->aux->hits->buffer[kb->hidx].mat == 0) kb->hidx ++;
				else break;
			}
			if(kb->hidx >= kb->aux->hits->size){
				sc = ref_layseqr(kb->seqs, 0);
				sc = _push_padding_ref_kbmblock(kb, sc);
				break;
			}
			hit = ref_kbmmapv(kb->aux->hits, kb->hidx ++);
			tsoff = kb->aux->kbm->reads->buffer[hit->tidx].seqoff * KBM_BIN_SIZE;
			rdlen = kb->aux->kbm->reads->buffer[hit->tidx].bincnt * KBM_BIN_SIZE;
			off   = hit->qb;
			rdoff = hit->qdir? Int(rdlen - hit->te) : hit->tb;
			{
				encap_layseqr(kb->seqs, 2);
				sc = pop_layseqr(kb->seqs);
				sc->chridx = kb->cidx;
				sc->bidx = (off / kb->bstep);
				sc->rdidx = hit->tidx;
				sc->rddir = hit->qdir;
				sc->rdoff = rdoff;
				sc->rbeg = off;
				sc->rend = 0;
				sc->rdtag  = kb->aux->kbm->reads->buffer[hit->tidx].tag;
				clear_basebank(sc->seq);
			}
			sl = NULL;
			{
				nxt = (off / kb->bstep) * kb->bstep;
				if(nxt && off - nxt < UInt(kb->bsize - kb->bstep)){
					if(sc->bidx){
						sl = pop_layseqr(kb->seqs);
						sl->chridx = kb->cidx;
						sl->bidx   = sc->bidx - 1;
						sl->rdidx  = hit->tidx;
						sl->rddir  = hit->qdir;
						sl->rdoff  = rdoff;
						sl->rbeg   = off;
						sl->rend   = 0;
						sl->rdtag  = kb->aux->kbm->reads->buffer[hit->tidx].tag;
						clear_basebank(sl->seq);
					}
					nxt += kb->bsize - kb->bstep;
				} else {
					nxt += kb->bstep;
				}
			}
			coff = hit->cgoff;
			clen = hit->cglen;
			while(clen){
				bt  = get_bitsvec(kb->aux->cigars, coff + clen - 1);
				clen --;
				val = KBM_BIN_SIZE;
				bt  = (bt & 0x03)? (bt & 0x03) : 3;
				while(val){
					if(bt & 0b001){
						len = num_min(val, nxt - off);
						off += len;
					} else {
						len = val;
					}
					val -= len;
					if(bt & 0b010){
						{
							if(sc){
								if(sc->seq->size == 0){
									sc->rbeg = (bt & 0b001)? off - len : off;
								}
								if(hit->qdir){
									revbits2basebank(sc->seq, kb->aux->kbm->rdseqs->bits, tsoff + rdlen - (rdoff + len), len);
								} else {
									fwdbits2basebank(sc->seq, kb->aux->kbm->rdseqs->bits, tsoff + rdoff, len);
								}
								sc->rend = off;
							}
							if(sl){
								if(sl->seq->size == 0){
									sl->rbeg = (bt & 0b001)? off - len : off;
								}
								if(hit->qdir){
									revbits2basebank(sl->seq, kb->aux->kbm->rdseqs->bits, tsoff + rdlen - (rdoff + len), len);
								} else {
									fwdbits2basebank(sl->seq, kb->aux->kbm->rdseqs->bits, tsoff + rdoff, len);
								}
								sl->rend = off;
							}
						}
						rdoff += len;
					}
					if(off == nxt){
						if(sl){
							if(sl->rend > sl->rbeg){
								u4i scidx;
								scidx = offset_layseqr(kb->seqs, sc);
								sl = _push_padding_ref_kbmblock(kb, sl);
								sc = ref_layseqr(kb->seqs, scidx);
								if(kb->lidx == MAX_U4){
									kb->lidx = offset_layseqr(kb->seqs, sl);
								}
								array_heap_push(kb->heap->buffer, kb->heap->size, kb->heap->cap, u4i, offset_layseqr(kb->seqs, sl),
									num_cmpx(ref_layseqr(kb->seqs, a)->bidx, ref_layseqr(kb->seqs, b)->bidx, ref_layseqr(kb->seqs, a)->rdidx, ref_layseqr(kb->seqs, b)->rdidx));
							}
							sl = NULL;
							nxt += 2 * kb->bstep - kb->bsize;
						} else {
							u4i scidx;
							scidx = offset_layseqr(kb->seqs, sc);
							encap_layseqr(kb->seqs, 1);
							sl = ref_layseqr(kb->seqs, scidx);
							sc = pop_layseqr(kb->seqs);
							sc->chridx = kb->cidx;
							sc->bidx  = sl->bidx + 1;
							sc->rdidx = hit->tidx;
							sc->rddir = hit->qdir;
							sc->rdoff = rdoff;
							sc->rbeg = off;
							sc->rend = 0;
							sc->rdtag  = kb->aux->kbm->reads->buffer[hit->tidx].tag;
							clear_basebank(sc->seq);
							nxt += kb->bsize - kb->bstep;
						}
					}
				}
			}
			if(sl && sl->rend > sl->rbeg){
				u4i scidx;
				scidx = offset_layseqr(kb->seqs, sc);
				sl = _push_padding_ref_kbmblock(kb, sl);
				sc = ref_layseqr(kb->seqs, scidx);
				if(kb->lidx == MAX_U4){
					kb->lidx = offset_layseqr(kb->seqs, sl);
				}
				array_heap_push(kb->heap->buffer, kb->heap->size, kb->heap->cap, u4i, offset_layseqr(kb->seqs, sl),
					num_cmpx(ref_layseqr(kb->seqs, a)->bidx, ref_layseqr(kb->seqs, b)->bidx, ref_layseqr(kb->seqs, a)->rdidx, ref_layseqr(kb->seqs, b)->rdidx));
			}
			if(sc && sc->rend > sc->rbeg){
				u4i scidx;
				scidx = offset_layseqr(kb->seqs, sl);
				sc = _push_padding_ref_kbmblock(kb, sc);
				sl = ref_layseqr(kb->seqs, scidx);
				if(kb->lidx == MAX_U4){
					kb->lidx = offset_layseqr(kb->seqs, sc);
				}
				array_heap_push(kb->heap->buffer, kb->heap->size, kb->heap->cap, u4i, offset_layseqr(kb->seqs, sc),
					num_cmpx(ref_layseqr(kb->seqs, a)->bidx, ref_layseqr(kb->seqs, b)->bidx, ref_layseqr(kb->seqs, a)->rdidx, ref_layseqr(kb->seqs, b)->rdidx));
			}
		}
		if(kb->heap->size == 0) break;
		if(kb->lidx != MAX_U4){
			sl = ref_layseqr(kb->seqs, kb->lidx);
			kb->sidx = kb->heap->buffer[0];
			sc = ref_layseqr(kb->seqs, kb->sidx);
			if(sc->chridx == sl->chridx && sc->bidx + 1 >= sl->bidx){
				kb->lidx = MAX_U4;
				continue;
			}
			array_heap_remove(kb->heap->buffer, kb->heap->size, kb->heap->cap, u4i, 0,
				num_cmpx(ref_layseqr(kb->seqs, a)->bidx, ref_layseqr(kb->seqs, b)->bidx, ref_layseqr(kb->seqs, a)->rdidx, ref_layseqr(kb->seqs, b)->rdidx));
			sc->rbeg -= sc->bidx * kb->bstep;
			sc->rend -= sc->bidx * kb->bstep;
		} else {
			kb->sidx = kb->heap->buffer[0];
			array_heap_remove(kb->heap->buffer, kb->heap->size, kb->heap->cap, u4i, 0,
				num_cmpx(ref_layseqr(kb->seqs, a)->bidx, ref_layseqr(kb->seqs, b)->bidx, ref_layseqr(kb->seqs, a)->rdidx, ref_layseqr(kb->seqs, b)->rdidx));
			sc = ref_layseqr(kb->seqs, kb->sidx);
			sc->rbeg -= sc->bidx * kb->bstep;
			sc->rend -= sc->bidx * kb->bstep;
		}
		return sc;
	} while(1);
	kb->sidx = MAX_U4;
	return NULL;
}

static inline void info_kbmblock(void *obj, lay_seq_t *sq, lay_blk_t *bk){
	KBMBlock *kb;
	kb = (KBMBlock*)obj;
	if(sq == NULL) return;
	bk->node1 = sq->bidx;
	bk->node2 = sq->bidx + 1;
	bk->reftag = kb->rdtag->string;
	bk->reflen = kb->rdlen;
	bk->refoff = sq->bidx * kb->bstep;
}

static inline int map_kbmpoa(CTGCNS *cc, KBMAux *aux, char *rdtag, u4i qidx, BaseBank *rdseq, u8i seqoff, u4i seqlen, u4i corr_min, u4i corr_max, float corr_cov, FILE *layf){
	ctg_cns_t *ctg;
	KBMBlock *kb;
	u4i i;
	int self_aln, max_hit, min_aln, min_mat;
	kb = (KBMBlock*)cc->obj;
	reset_ctgcns(cc, kb, iter_kbmblock, info_kbmblock);
	seqlen = kbm_cvt_length(seqlen);
	if(seqlen < 4 * KBM_BIN_SIZE + UInt(aux->par->min_aln)) return 0;
	if(rdseq && rdseq != aux->kbm->rdseqs){
		self_aln = 0;
	} else {
		self_aln = 1;
		rdseq = aux->kbm->rdseqs;
	}
	max_hit  = aux->par->max_hit;
	min_aln  = aux->par->min_aln;
	min_mat  = aux->par->min_mat;
	aux->par->self_aln = 0;
	aux->par->max_hit = corr_max;
	aux->par->min_aln = num_max(seqlen * corr_cov, min_aln);
	aux->par->min_mat = num_max(aux->par->min_aln * aux->par->min_sim, min_mat);
	query_index_kbm(aux, rdtag, qidx, rdseq, seqoff, seqlen);
	map_kbm(aux);
	aux->par->self_aln = 0;
	aux->par->max_hit  = max_hit;
	aux->par->min_aln  = min_aln;
	aux->par->min_mat  = min_mat;
	if(KBM_LOG){
		fprintf(KBM_LOGF, ">>> Contained alignments\n"); fflush(KBM_LOGF);
		for(i=0;i<aux->hits->size;i++){
			fprint_hit_kbm(aux, i, KBM_LOGF);
		}
	}
	if(aux->hits->size < corr_min){
		return 0;
	}
	if(self_aln){
		for(i=0;i<aux->hits->size;i++){
			if(aux->hits->buffer[i].tidx == qidx){
				aux->hits->buffer[i].mat = 0;
			}
		}
	}
	reset_kbmblock(kb, rdtag, qidx, rdseq, seqoff, seqlen, aux);
	if(layf){
		print_lays_ctgcns(cc, layf);
		fflush(layf);
		reset_kbmblock(kb, rdtag, qidx, rdseq, seqoff, seqlen, aux);
		reset_ctgcns(cc, kb, iter_kbmblock, info_kbmblock);
	}
	if((ctg = iter_ctgcns(cc))){
		if(KBM_LOG){
			fprintf(KBM_LOGF, ">%s corrected\n", ctg->tag->string);
			print_lines_basebank(ctg->cns, 0, ctg->cns->size, KBM_LOGF, 100);
			fflush(KBM_LOGF);
		}
	} else {
		return 0;
	}
	clear_kbmmapv(aux->hits);
	if(ctg->cns->size > seqlen){
		ctg->cns->size = seqlen;
		normalize_basebank(ctg->cns);
	} else if(ctg->cns->size < seqlen){
		ctg->cns->size = kbm_cvt_length(ctg->cns->size);
		normalize_basebank(ctg->cns);
	}
	if(ctg->cns->size == 0){
		repay_ctgcns(cc, ctg);
		return 0;
	}
	query_index_kbm(aux, rdtag, qidx, ctg->cns, 0, ctg->cns->size);
	map_kbm(aux);
	repay_ctgcns(cc, ctg); // Please make sure ctg is not used unless this function return
	return 1;
}

#endif
