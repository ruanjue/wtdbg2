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

#ifndef __PGZF_RJ_H
#define __PGZF_RJ_H

#include <zlib.h>
#include "mem_share.h"
#include "list.h"
#include "thread.h"

#define PGZF_DEFAULT_BUFF_SIZE	(1U << 24) // 16MB
#define PGZF_MAX_BUFF_SIZE	(1U << 28)
#define PGZF_HEAD_SIZE	30
#define PGZF_HEAD_ZS_OFFSET	16
#define PGZF_HEAD_ZX_OFFSET	24
#define PGZF_TAIL_SIZE	8
#define PGZF_INDEX_BIN	64

#define PGZF_MODE_W	1 // pgzf writer
#define PGZF_MODE_R	2 // pgzf reader
#define PGZF_MODE_R_GZ	3 // gz reader
#define PGZF_MODE_R_UNKNOWN	4 // unknown file type

#define PGZF_FILETYPE_UNKNOWN	0
#define PGZF_FILETYPE_GZ	1
#define PGZF_FILETYPE_PGZF	2

struct PGZF;

#define PGZF_TASK_NULL	0
#define PGZF_TASK_DEFLATE	1
#define PGZF_TASK_INFLATE	2

thread_beg_def(pgz);
struct PGZF *pz;
u4i zsval, soff;
u8i zxval, doff;
u1v *dst, *src;
u4i token;
int level;
int task, zerosize;
thread_end_def(pgz);

typedef struct PGZF {
	u4i ncpu, ridx, widx;
	int rw_mode, seekable;
	u4i bufsize; // MUST be multiple of 1MB
	u8i xsize; // total uncompressed size
	u8v *boffs;
	u8v *xoffs;
	u8i tot_in, tot_out;
	u1v **dsts, **srcs, *tmp;
	z_stream *z;
	u8i offset;
	FILE *file;
	thread_def_shared_vars(pgz);
	int step; // used in decompress gzip file
	int eof, error;
} PGZF;

static inline void _num2bytes_pgzf(u1i *bs, u1i bl, u8i val){
	u1i i;
	for(i=0;i<bl;i++){
		bs[i] = (u1i)val;
		val >>= 8;
	}
}

static inline u8i _bytes2num_pgzf(u1i *bs, u1i bl){
	u8i val;
	u1i i;
	val = 0;
	for(i=0;i<bl;i++){
		val = val | (((u8i)bs[i]) << (8 * i));
	}
	return val;
}

/**
 * Please see https://tools.ietf.org/html/rfc1952
 */

static inline void _gen_pgzf_header(u1i bs[30], u4i z_size){
	bs[0] = 0x1f; // GZIP ID1
	bs[1] = 0x8b; // GZIP ID2
	bs[2] = 8; // CM = deflate
	bs[3] = 0x4; // FLG = 0b00000100, FEXTRA is ture
	bs[4] = 0; // MTIME
	bs[5] = 0; // MTIME
	bs[6] = 0; // MTIME
	bs[7] = 0; // MTIME
	bs[8] = 0; // XFL
	bs[9] = 3; // OS = unix
	bs[10] = 18; // XLEN
	bs[11] = 0; // = 18
	bs[12] = 'Z'; // TAG ZS
	bs[13] = 'S'; // compressed size
	bs[14] = 4; // TAG LEN
	bs[15] = 0; // 4
	bs[16] = z_size >>  0; // compressed block size
	bs[17] = z_size >>  8; //
	bs[18] = z_size >> 16; //
	bs[19] = z_size >> 24; // = z_size
	bs[20] = 'Z'; //TAG ZX
	bs[21] = 'X'; // every 64 block size
	bs[22] = 6; // TAG LEN
	bs[23] = 0; // 6
	bs[24] = 0; // reserve to store the random access index
	bs[25] = 0; //
	bs[26] = 0; //
	bs[27] = 0; //
	bs[28] = 0; //
	bs[29] = 0; //
}

static inline void _gen_pgzf_tailer(u1i bs[8], u4i crc, u4i u_size){
	_num2bytes_pgzf(bs + 0, 4, crc);
	_num2bytes_pgzf(bs + 4, 4, u_size);
}

static inline u4i _zlib_raw_deflate_all(u1i *dst, u4i dlen, u1i *src, u4i slen, int level){
	z_stream Z, *z;
	u4i ret;
	z = &Z;
	ZEROS(z);
	deflateInit2(z, level, Z_DEFLATED, -15, 9, Z_DEFAULT_STRATEGY);
	z->avail_in = slen;
	z->next_in  = src;
	z->avail_out = dlen;
	z->next_out  = dst;
	deflate(z, Z_FINISH);
	ret = dlen - z->avail_out;
	deflateEnd(z);
	return ret;
}

static inline u4i _pgzf_deflate(u1v *dst, u1v *src, int level){
	u4i z_size;
	uLong crc;
	clear_u1v(dst);
	//if(src->size == 0) return 0;
	if(src->size >= MAX_U4) return 0;
	z_size = compressBound(src->size);
	encap_u1v(dst, z_size + PGZF_HEAD_SIZE + PGZF_TAIL_SIZE);
	z_size = _zlib_raw_deflate_all(dst->buffer + PGZF_HEAD_SIZE, z_size, src->buffer, src->size, level);
	_gen_pgzf_header(dst->buffer, z_size + PGZF_HEAD_SIZE + PGZF_TAIL_SIZE);
	crc = crc32(0L, Z_NULL, 0);
	crc = crc32(crc, src->buffer, src->size);
	_gen_pgzf_tailer(dst->buffer + PGZF_HEAD_SIZE + z_size, crc, src->size);
	dst->size = PGZF_HEAD_SIZE + z_size + PGZF_TAIL_SIZE;
	return dst->size;
}

static inline int _read_pgzf_header(FILE *in, u1v *src, u4i *hoff, u4i *zsval, u8i *zxval){
	u4i off, val, sl, end;
	int ch, is_pgzf, xflag;
	char si1, si2;
	is_pgzf = 0;
	off = *hoff;
	*zsval = 0;
	*zxval  = 0;
	if(src->size < off + 10){
		encap_u1v(src, 10);
		src->size += fread(src->buffer + src->size, 1, off + 10 - src->size, in);
	}
	// At least give 10 bytes
	if(src->size < off + 10) return PGZF_FILETYPE_UNKNOWN;
	if(src->buffer[off + 0] != 0x1f) return PGZF_FILETYPE_UNKNOWN;
	if(src->buffer[off + 1] != 0x8b) return PGZF_FILETYPE_UNKNOWN;
	if((src->buffer[off + 2] != 0x08) || (src->buffer[off + 2] & 0xE0)) return PGZF_FILETYPE_UNKNOWN;
	xflag = src->buffer[off + 3];
	off += 10;
	if(xflag & 0x04){
		if(src->size < off + 2){
			encap_u1v(src, 2);
			sl = fread(src->buffer + src->size, 1, off + 2 - src->size, in);
			src->size += sl;
		}
		if(src->size < off + 2) return PGZF_FILETYPE_UNKNOWN;
		val = _bytes2num_pgzf(src->buffer + off, 2);
		off += 2;
		end = off + val;
		if(val > 0 && val < 4) return PGZF_FILETYPE_UNKNOWN;
		if(src->size < off + val){
			encap_u1v(src, val);
			sl = fread(src->buffer + src->size, 1, off + val - src->size, in);
			src->size += sl;
			if(src->size < off + val) return PGZF_FILETYPE_UNKNOWN;
		}
		//parse TAGs
		while(off < end){
			si1 = src->buffer[off + 0];
			si2 = src->buffer[off + 1];
			sl = _bytes2num_pgzf(src->buffer + off + 2, 2);
			off += 4;
			if(si1 == 'Z' && si2 == 'S' && sl == 4){
				is_pgzf = 1;
				*zsval = _bytes2num_pgzf(src->buffer + off, 4);
			} else if(is_pgzf && si1 == 'Z' && si2 == 'X' && sl == 6){
				*zxval = _bytes2num_pgzf(src->buffer + off, 6);
			}
			off += sl;
		}
	}
	if(xflag & 0x08){
		do {
			if(off < src->size){
				ch = src->buffer[off];
			} else {
				ch = getc(in);
				if(ch == -1){
					return PGZF_FILETYPE_UNKNOWN;
				}
				push_u1v(src, ch);
			}
			off ++;
		} while(ch);
	}
	if(xflag & 0x10){
		do {
			if(off < src->size){
				ch = src->buffer[off];
			} else {
				ch = getc(in);
				if(ch == -1){
					return PGZF_FILETYPE_UNKNOWN;
				}
				push_u1v(src, ch);
			}
			off ++;
		} while(ch);
	}
	if(xflag & 0x02){
		if(src->size < off + 2){
			encap_u1v(src, 2);
			sl = fread(src->buffer + src->size, 1, off + 2 - src->size, in);
			src->size += sl;
		}
		off += 2;
		if(src->size < off) return PGZF_FILETYPE_UNKNOWN;
	}
	*hoff = off;
	return is_pgzf? PGZF_FILETYPE_PGZF : PGZF_FILETYPE_GZ;
}

int pgzf_inflate_raw_core(z_stream *z, u1i *dst, u4i *dlen, u1i *src, u4i *slen, int flush){
	u4i dl, sl;
	int ret;
	ret = Z_OK;
	dl = *dlen;
	sl = *slen;
	z->avail_in = sl;
	z->next_in  = src;
	z->avail_out = dl;
	z->next_out  = dst;
	ret = inflate(z, flush);
	*dlen = dl - z->avail_out;
	*slen = sl - z->avail_in;
	return ret;
}

// src start just after gz_header, and include fz_tailer
int pgzf_inflate_core(u1i *dst, u4i *dlen, u1i *src, u4i slen, int check){
	z_stream Z, *z;
	u4i soff, dsz;
	uLong crc, rcr;
	int ret;
	z = &Z;
	ZEROS(z);
	inflateInit2(z, -15);
	z->avail_in  = slen - PGZF_TAIL_SIZE;
	z->next_in   = src;
	z->avail_out = *dlen;
	z->next_out  = dst;
	ret = inflate(z, Z_FINISH);
	*dlen -= z->avail_out;
	soff = slen - PGZF_TAIL_SIZE - z->avail_in;
	inflateEnd(z);
	if(check){
		if(soff + 8 > slen){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			return 0;
		}
		rcr = _bytes2num_pgzf(src + soff, 4);
		dsz = _bytes2num_pgzf(src + soff + 4, 4);
		if(dsz != *dlen){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			return 0;
		}
		crc = crc32(0L, Z_NULL, 0);
		crc = crc32(crc, dst, *dlen);
		if(crc != rcr){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			return 0;
		}
	}
	return 1;
}

thread_beg_func(pgz);
PGZF *pz;
u1v *dst, *src;
u4i bufsize, hsize, rsize, dsz, ssz, next;
int ret;
pz  = pgz->pz;
dst = pgz->dst;
src = pgz->src;
thread_beg_loop(pgz);
if(pgz->task == PGZF_TASK_DEFLATE){
	if(src->size == 0 && !pgz->zerosize) continue;
	clear_u1v(dst);
	_pgzf_deflate(dst, src, pgz->level);
	while(pz->ridx != pgz->token){
		nano_sleep(1);
	}
	{
		pz->tot_out += pgz->dst->size;
		push_u8v(pz->boffs, pz->tot_out);
		fwrite(pgz->dst->buffer, 1, pgz->dst->size, pz->file);
		clear_u1v(pgz->dst);
		clear_u1v(pgz->src);
	}
	pz->ridx ++;
} else if(pgz->task == PGZF_TASK_INFLATE){
	pgz->doff = 0;
	clear_u1v(pgz->dst);
	while((pz->ridx % pz->ncpu) != UInt(pgz->t_idx)){
		nano_sleep(10);
		if(pz->error) break;
		//if(pz->eof){
			//if(pz->rw_mode != PGZF_MODE_R_GZ){
				//break;
			//}
		//}
		if(pgz->running == 0){
			break;
		}
	}
	if(pz->error) break;
	if(pz->rw_mode == PGZF_MODE_R){
		if(pgz->src->size){ // loaded header, had set zsval and zxval
		} else {
			pgz->zsval = pgz->zxval = 0;
			pgz->soff = pgz->src->size = 0;
			ret = _read_pgzf_header(pz->file, pgz->src, &pgz->soff, &pgz->zsval, &pgz->zxval);
			if(pgz->src->size == 0){
				pz->eof = 1;
				pz->ridx ++;
				break;
			}
			if(ret != PGZF_FILETYPE_PGZF){
				fprintf(stderr, " -- Error: not a PGZF format at %u block, in %s -- %s:%d --\n", pz->ridx, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
				pz->error = 1;
				pz->ridx ++;
				break;
			}
		}
		hsize = pgz->soff;
		encap_u1v(pgz->src, pgz->zsval - pgz->src->size);
		rsize = fread(pgz->src->buffer + hsize, 1, pgz->zsval - pgz->src->size, pz->file);
		if(rsize < pgz->zsval - pgz->src->size){
			fprintf(stderr, " -- Error: read %u < %u at %u block, in %s -- %s:%d --\n", UInt(pgz->src->size + rsize), pgz->zsval, pz->ridx, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			pz->error = 1;
			pz->ridx ++;
			break;
		}
		pgz->src->size += rsize;
		pz->tot_in += pgz->zsval;
		pz->ridx ++;
		dsz = _bytes2num_pgzf(pgz->src->buffer + pgz->zsval - 4, 4);
		encap_u1v(pgz->dst, dsz);
		pgz->soff = 0;
		if(pgzf_inflate_core(pgz->dst->buffer, &dsz, pgz->src->buffer + hsize, pgz->zsval - hsize, 1) == 0){
			clear_u1v(pgz->src);
			pz->error = 1;
			break;
		}
		pgz->dst->size = dsz;
		clear_u1v(pgz->src);
	} else if(pz->rw_mode == PGZF_MODE_R_GZ){
		u4i bsz;
		bsz = 1024 * 1024;
		bufsize = pz->bufsize? pz->bufsize : PGZF_DEFAULT_BUFF_SIZE;
		encap_u1v(pgz->dst, bufsize);
		while(!pz->error){
			if(pgz->src->size == pgz->soff){
				pgz->soff = pgz->src->size = 0;
			}
			if(pgz->src->size < bsz){
				encap_u1v(pgz->src, bsz - pgz->src->size);
				rsize = fread(pgz->src->buffer + pgz->src->size, 1, bsz - pgz->src->size, pz->file);
				if(rsize < bsz - pgz->src->size){
					pz->eof = 1;
				}
				pz->tot_in += rsize;
				pgz->src->size += rsize;
			}
			if(pgz->src->size == pgz->soff){
				break;
			}
			if(pz->step == 0){
				u4i tsz;
				tsz = pgz->src->size;
				ret = _read_pgzf_header(pz->file, pgz->src, &pgz->soff, &pgz->zsval, &pgz->zxval);
				if(ret != PGZF_FILETYPE_GZ && ret != PGZF_FILETYPE_PGZF){
					if(pgz->src->size == pgz->soff){
						pz->eof = 1;
					} else {
						fprintf(stderr, " -- failed in read gzip header, ret = %d in %s -- %s:%d --\n", ret, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
						pz->error = 1;
					}
					break;
				} else {
					pz->tot_in += pgz->src->size - tsz;
				}
				pz->step = 1;
				continue;
			} else if(pz->step == 2){
				if(pgz->src->size >= pgz->soff + PGZF_TAIL_SIZE){
					pgz->soff += PGZF_TAIL_SIZE;
					pz->step = 0;
					inflateReset(pz->z);
					continue;
				} else if(pz->eof){
					pz->error = 2;
					break;
				} else {
					memmove(pgz->src->buffer, pgz->src->buffer + pgz->soff, pgz->src->size - pgz->soff);
					pgz->src->size -= pgz->soff;
					pgz->soff = 0;
				}
			}
			while(pgz->dst->size < bufsize && pgz->soff < pgz->src->size){
				dsz = bufsize - pgz->dst->size;
				ssz = pgz->src->size - pgz->soff;
				ret = pgzf_inflate_raw_core(pz->z, pgz->dst->buffer + pgz->dst->size, &dsz, pgz->src->buffer + pgz->soff, &ssz, Z_NO_FLUSH);
				pgz->dst->size += dsz;
				pgz->soff += ssz;
				if(ret == Z_STREAM_END){
					pz->step = 2;
					break;
				} else if(ret != Z_OK){
					fprintf(stderr, " -- ZERROR: %d in %s -- %s:%d --\n", ret, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
					pz->error = 1;
					break;
				}
			}
			if(pgz->dst->size == bufsize){
				if(pgz->soff < pgz->src->size){
					if(pz->ncpu > 1){
						next = (pz->ridx + 1) % pz->ncpu;
						if(pz->srcs[next]->size != 0){
							fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
							abort();
						}
						append_array_u1v(pz->srcs[next], pgz->src->buffer + pgz->soff, pgz->src->size - pgz->soff);
					}
				}
				pgz->soff = pgz->src->size = 0;
				break;
			}
		}
		pz->ridx ++;
	} else if(pz->rw_mode == PGZF_MODE_R_UNKNOWN){
		bufsize = pz->bufsize? pz->bufsize : PGZF_DEFAULT_BUFF_SIZE;
		encap_u1v(pgz->dst, bufsize);
		if(pgz->src->size > bufsize){
			fprintf(stderr, " -- something wrong in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			pz->error = 1;
			pz->ridx ++;
			break;
		} else if(pgz->src->size){
			append_u1v(pgz->dst, pgz->src);
			clear_u1v(pgz->src);
		}
		rsize = fread(pgz->dst->buffer + pgz->dst->size, 1, bufsize - pgz->dst->size, pz->file);
		if(rsize < bufsize - pgz->dst->size){
			pz->eof = 1;
		}
		pgz->dst->size += rsize;
		pz->tot_in += pgz->dst->size;
		pz->ridx ++;
	}
}
thread_end_loop(pgz);
thread_end_func(pgz);

static inline PGZF* open_pgzf_writer(FILE *out, u4i buffer_size, int ncpu, int level){
	PGZF *pz;
	u4i i;
	b8i offset;
	thread_prepare(pgz);
	pz = malloc(sizeof(PGZF));
	if(ncpu < 1){
		get_linux_sys_info(NULL, NULL, &ncpu);
		if(ncpu < 1) ncpu = 8;
	}
	pz->ncpu = ncpu;
	pz->ridx = 0;
	pz->widx = 0;
	offset = ftell(out);
	if(offset == -1){
		pz->offset = 0;
		pz->seekable = 0;
	} else {
		pz->offset = offset;
		pz->seekable = 1;
	}
	pz->file = out;
	pz->error = 0;
	pz->eof = 0;
	pz->step = 0;
	pz->rw_mode = 1; // write
	if(buffer_size == 0) buffer_size = PGZF_DEFAULT_BUFF_SIZE;
	pz->bufsize = (buffer_size + 0xFFFFFU) & 0xFFF00000U;
	pz->xsize = 0;
	pz->boffs = init_u8v(32);
	pz->xoffs = init_u8v(32);
	pz->z = NULL;
	pz->dsts = calloc(pz->ncpu, sizeof(u1v*));
	for(i=0;i<pz->ncpu;i++){
		pz->dsts[i] = init_u1v(pz->bufsize);
	}
	pz->srcs = calloc(pz->ncpu, sizeof(u1v*));
	for(i=0;i<pz->ncpu;i++){
		pz->srcs[i] = init_u1v(pz->bufsize);
	}
	pz->tmp = init_u1v(32);
	pz->tot_in  = 0;
	pz->tot_out = 0;
	if(level == 0) level = Z_DEFAULT_COMPRESSION; // disable level 0, set to default level 6
	thread_beg_init(pgz, pz->ncpu);
	pgz->pz = pz;
	pgz->zsval = 0;
	pgz->zxval = 0;
	pgz->soff = 0;
	pgz->doff = 0;
	pgz->dst = pz->dsts[pgz->t_idx];
	pgz->src = pz->srcs[pgz->t_idx];
	pgz->token = 0;
	pgz->level = level;
	pgz->task = PGZF_TASK_NULL;
	pgz->zerosize = 0;
	thread_end_init(pgz);
	thread_export(pgz, pz);
	return pz;
}

static inline size_t write_pgzf(PGZF *pz, void *dat, size_t len){
	size_t off, cnt;
	thread_prepare(pgz);
	thread_import(pgz, pz);
	off = 0;
	while(off < len){
		thread_beg_operate(pgz, pz->widx % pz->ncpu);
		thread_wait(pgz);
		/*
		if(pgz->dst->size){
			pz->tot_out += pgz->dst->size;
			push_u8v(pz->boffs, pz->tot_out);
			fwrite(pgz->dst->buffer, 1, pgz->dst->size, pz->file);
			clear_u1v(pgz->dst);
			clear_u1v(pgz->src);
		}
		*/
		cnt = num_min(len - off, pz->bufsize - pgz->src->size);
		append_array_u1v(pgz->src, dat + off, cnt);
		off += cnt;
		if(pgz->src->size == pz->bufsize){
			pz->tot_in += pgz->src->size;
			pgz->task = PGZF_TASK_DEFLATE;
			pgz->token = pz->widx;
			thread_wake(pgz);
			pz->widx ++;
		}
	}
	thread_export(pgz, pz);
	return len;
}

static inline void _end_pgzf_writer(PGZF *pz){
	u4i i, widx;
	thread_prepare(pgz);
	thread_import(pgz, pz);
	widx = pz->widx;
	for(i=0;i<=pz->ncpu;i++){ // (pz->tidx + ncpu + 1) % ncpu
		thread_beg_operate(pgz, widx % pz->ncpu);
		thread_wait(pgz);
		/*
		if(pgz->dst->size){
			pz->tot_out += pgz->dst->size;
			push_u8v(pz->boffs, pz->tot_out);
			fwrite(pgz->dst->buffer, 1, pgz->dst->size, pz->file);
			clear_u1v(pgz->dst);
			clear_u1v(pgz->src);
		}
		*/
		if(pgz->src->size){ // will force to write un-full block
			pz->tot_in += pgz->src->size;
			pgz->task = PGZF_TASK_DEFLATE;
			pgz->token = pz->widx;
			thread_wake(pgz);
			pz->widx ++;
			thread_wait(pgz);
		}
		widx ++;
	}
	if(pz->widx == 0){
		thread_beg_operate(pgz, 0);
		pgz->zerosize = 1;
		pgz->task = PGZF_TASK_DEFLATE;
		thread_wake(pgz);
		thread_wait(pgz);
	}
	thread_export(pgz, pz);
}

static inline int write_index_pgzf(PGZF *pz){
	u8i i, x;
	u1i bs[6];
	pz->xsize = pz->tot_in;
	if(!pz->seekable) return 0;
	if(fseek(pz->file, pz->offset + PGZF_HEAD_ZX_OFFSET, SEEK_SET) == -1){
		perror("fseek error in write_index_pgzf");
		return 0;
	}
	_num2bytes_pgzf(bs, 6, pz->xsize);
	fwrite(bs, 1, 6, pz->file);
	for(i=64,x=1;i+PGZF_INDEX_BIN<=pz->boffs->size;i+=PGZF_INDEX_BIN,x++){
		push_u8v(pz->xoffs, pz->boffs->buffer[i+PGZF_INDEX_BIN]);
		_num2bytes_pgzf(bs, 6, pz->boffs->buffer[i+PGZF_INDEX_BIN]);
		if(fseek(pz->file, pz->offset + pz->boffs->buffer[x] + PGZF_HEAD_ZX_OFFSET, SEEK_SET) == -1){
			perror("fseek error in write_index_pgzf");
			return 0;
		}
		fwrite(bs, 1, 6, pz->file);
	}
	fseek(pz->file, 0, SEEK_END);
	return 1;
}

static inline PGZF* open_pgzf_reader(FILE *in, u4i bufsize, int ncpu){
	PGZF *pz;
	u8i zxval;
	b8i offset;
	u4i i, zsval, hoff;
	int ftype;
	thread_prepare(pgz);
	pz = malloc(sizeof(PGZF));
	pz->ncpu = ncpu;
	pz->ridx = 0;
	pz->widx = 0;
	offset = ftell(in);
	if(offset == -1){
		pz->offset = 0;
		pz->seekable = 0;
	} else {
		pz->offset = offset;
		pz->seekable = 1;
	}
	pz->file = in;
	pz->eof = 0;
	pz->error = 0;
	pz->step = 0;
	pz->dsts = calloc(pz->ncpu, sizeof(u1v*));
	pz->srcs = calloc(pz->ncpu, sizeof(u1v*));
	pz->tmp = init_u1v(32);
	pz->tot_in  = 0;
	pz->tot_out = 0;
	pz->boffs = init_u8v(32);
	pz->xoffs = init_u8v(32);
	// recognize PGZF
	zsval = zxval = 0;
	hoff = 0;
	pz->srcs[0] = init_u1v(1024);
	ftype = _read_pgzf_header(pz->file, pz->srcs[0], &hoff, &zsval, &zxval);
	pz->tot_in = pz->srcs[0]->size;
	switch(ftype){
		case PGZF_FILETYPE_GZ: pz->step = 1; pz->rw_mode = PGZF_MODE_R_GZ; break;
		case PGZF_FILETYPE_PGZF: pz->rw_mode = PGZF_MODE_R; break;
		default:
			fprintf(stderr, " ** WARNNING: input file is not in gzip format **\n");
			pz->rw_mode = PGZF_MODE_R_UNKNOWN; break;
	}
	if(pz->rw_mode == PGZF_MODE_R){
		pz->z = NULL;
		pz->xsize = zxval;
		push_u8v(pz->boffs, zsval);
		if(pz->seekable){
			u8i foff;
			foff = ftell(pz->file);
			if(fseek(pz->file, pz->offset + zsval - 4, SEEK_SET) == -1){
				fprintf(stderr, " ** ERROR: failed to read uncompress block size in the first block ERR(1) **\n");
				return NULL;
			}
			if(fread(&pz->bufsize, 4, 1, pz->file) == 0){
				fprintf(stderr, " ** ERROR: failed to read uncompress block size in the first block ERR(2) **\n");
				return NULL;
			}
			if(fseek(pz->file, foff, SEEK_SET) == -1){
				fprintf(stderr, " ** ERROR: failed to read uncompress block size in the first block ERR(3) **\n");
				return NULL;
			}
		} else {
			pz->bufsize = bufsize;
		}
	} else if(pz->rw_mode == PGZF_MODE_R_GZ){
		pz->z = calloc(1, sizeof(z_stream));
		inflateInit2(pz->z, -15);
		pz->bufsize = bufsize;
	} else {
		pz->z = NULL;
		pz->bufsize = bufsize;
	}
	if(pz->bufsize == 0) pz->bufsize = PGZF_DEFAULT_BUFF_SIZE;
	pz->bufsize = (pz->bufsize + 0xFFFFFU) & 0xFFF00000U;
	for(i=0;i<pz->ncpu;i++){
		pz->dsts[i] = init_u1v(pz->bufsize);
	}
	if(pz->bufsize > pz->srcs[0]->size){
		encap_u1v(pz->srcs[0], pz->bufsize - pz->srcs[0]->size);
	}
	for(i=1;i<pz->ncpu;i++){
		pz->srcs[i] = init_u1v(pz->bufsize);
	}
	thread_beg_init(pgz, pz->ncpu);
	pgz->pz = pz;
	pgz->zsval = pgz->t_idx? 0 : zsval;
	pgz->zxval = pgz->t_idx? 0 : zxval;
	pgz->soff = pgz->t_idx? 0 : hoff;
	pgz->doff = 0;
	pgz->src = pz->srcs[pgz->t_idx];
	pgz->dst = pz->dsts[pgz->t_idx];
	pgz->level = Z_DEFAULT_COMPRESSION; // useless in inflating
	pgz->task = PGZF_TASK_INFLATE;
	thread_end_init(pgz);
	thread_wake_all(pgz);
	thread_export(pgz, pz);
	return pz;
}

/*
static inline _clear_pgzf_reader(PGZF *pz){
	UNUSED(pz);
}

static inline off_t seek_pgzf(PGZF *pz, u8i offset){
	u4i bidx, boff, xidx, xoff;
	if(offset > pz->xsize) return -1;
	else if(offset == pz->xsize){
		pz->eof = 1;
		return offset;
	}
	if(!pz->seekable) return -1;
	bidx = offset / pz->bufsize;
	boff = offset % pz->>bufsize;
	xidx = bidx / PGZF_INDEX_BIN;
	xoff = bidx % PGZF_INDEX_BIN;
	if(xidx > pz->xoffs->size){
		
	}
	return 0;
}
*/

static inline size_t read_pgzf(PGZF *pz, void *dat, size_t len){
	size_t off;
	u4i nrun;
	thread_prepare(pgz);
	thread_import(pgz, pz);
	nrun = 0;
	for(off=0;off<len;){
		thread_beg_operate(pgz, pz->widx % pz->ncpu);
		thread_wait(pgz);
		if(pz->error) break;
		if(len - off < pgz->dst->size - pgz->doff){
			memcpy(dat + off, pgz->dst->buffer + pgz->doff, len - off);
			pz->tot_out += len - off;
			pgz->doff += len - off;
			off = len;
			break;
		} else if(pgz->dst->size){
			memcpy(dat + off, pgz->dst->buffer + pgz->doff, pgz->dst->size - pgz->doff);
			pz->tot_out += pgz->dst->size - pgz->doff;
			off += pgz->dst->size - pgz->doff;
			pgz->doff = pgz->dst->size;
			pgz->task = PGZF_TASK_INFLATE;
			thread_wake(pgz);
			pz->widx ++;
		} else if(pz->eof){
			nrun ++;
			if(nrun >= pz->ncpu){
				break;
			}
		}
	}
	return off;
}

static inline void close_pgzf(PGZF *pz){
	thread_prepare(pgz);
	if(pz->rw_mode == PGZF_MODE_W){
		_end_pgzf_writer(pz);
	}
	thread_import(pgz, pz);
	thread_beg_close(pgz);
	free_u1v(pgz->dst);
	free_u1v(pgz->src);
	thread_end_close(pgz);
	free(pz->dsts);
	free(pz->srcs);
	free_u1v(pz->tmp);
	switch(pz->rw_mode){
		case PGZF_MODE_W: write_index_pgzf(pz); fflush(pz->file); break;
		case PGZF_MODE_R: break;
		case PGZF_MODE_R_GZ:
			if(pz->z){
				inflateEnd(pz->z);
				free(pz->z);
			}
			break;
	}
	free_u8v(pz->boffs);
	free_u8v(pz->xoffs);
	free(pz);
}

static inline size_t read_pgzf4filereader(void *obj, void *dat, size_t len){ return read_pgzf((PGZF*)obj, dat, len); }
static inline void close_pgzf4filereader(void *obj){
	PGZF *pz;
	pz = (PGZF*)obj;
	if(pz->file != stdin){
		fclose(pz->file);
	}
	return close_pgzf(pz);
}

static inline size_t write_pgzf4filewriter(void *obj, void *dat, size_t len){ return write_pgzf((PGZF*)obj, dat, len); }
static inline void close_pgzf4filewriter(void *obj){
	PGZF *pz;
	pz = (PGZF*)obj;
	return close_pgzf(pz);
}

#endif
