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

#ifndef __FILEWRITER_RJ_H
#define __FILEWRITER_RJ_H

#include "mem_share.h"
#include "thread.h"
#include "pgzf.h"

typedef size_t (*write_data_func)(void *obj, void *dat, size_t len);
typedef void (*close_output_func)(void *obj);

static inline size_t _write_data_file(void *obj, void *dat, size_t len){ return fwrite(dat, 1, len, (FILE*)obj); }
static inline void _close_output_file(void *obj){ if(obj) fclose((FILE*)obj); }

/**
 * BufferedWriter
 */
typedef struct {
	FILE *bios[2];
	FILE *out;
	void *_file;
	write_data_func _write;
	close_output_func _close;
	int bidx;
	size_t buf_size;
	char *buffs[2];
	size_t blens[2];
	size_t nbytes;
	pthread_mutex_t lock;
	pthread_t pid;
	int running, flush;
} BufferedWriter;

static inline void* _buffered_writer_thread_func(void *obj){
	BufferedWriter *bw;
	size_t bsize[2];
	int bidx, lock;
	bw = (BufferedWriter*)obj;
	bw->running = 1;
	bw->flush = 0;
	bw->nbytes = 0;
	while(bw->running){
		bidx = bw->bidx;
		bsize[0] = ftell(bw->bios[0]);
		bsize[1] = ftell(bw->bios[1]);
		if(bsize[bidx] >= bw->buf_size || (bsize[bidx] && bw->flush == 1)){
			lock = 1;
			pthread_mutex_lock(&bw->lock);
		} else {
			lock = 0;
		}
		if(bsize[!bidx]){
			fflush(bw->bios[!bidx]);
			bw->_write(bw->_file, bw->buffs[!bidx], bsize[!bidx]);
			bw->nbytes += bsize[!bidx];
			fseek(bw->bios[!bidx], 0, SEEK_SET);
		}
		if(lock){
			bw->bidx = !bidx;
			pthread_mutex_unlock(&bw->lock);
		} else if(bsize[bidx]){
			pthread_mutex_lock(&bw->lock);
			bw->bidx = !bidx;
			pthread_mutex_unlock(&bw->lock);
		}
		if(bw->flush && bsize[0] == 0 && bsize[1] == 0){
			bw->flush = 2;
			while(bw->flush == 2){
				nano_sleep(1);
			}
			bw->flush = 0;
		}
		nano_sleep(10);
	}
	{
		bsize[0] = ftell(bw->bios[0]);
		bsize[1] = ftell(bw->bios[1]);
		fflush(bw->bios[0]);
		fflush(bw->bios[1]);
		bidx = bw->bidx;
		if(bsize[!bidx]){
			bw->_write(bw->_file, bw->buffs[!bidx], bsize[!bidx]);
			bw->nbytes += bsize[!bidx];
		}
		if(bsize[bidx]){
			bw->_write(bw->_file, bw->buffs[bidx], bsize[bidx]);
			bw->nbytes += bsize[bidx];
		}
	}
	return NULL;
}

static inline BufferedWriter* open2_bufferedwriter(void *obj, write_data_func _write, close_output_func _close, size_t buf_size){
	BufferedWriter *bw;
	bw = malloc(sizeof(BufferedWriter));
	bw->_file = obj;
	bw->_write = _write;
	bw->_close = _close;
	bw->buffs[0] = NULL;
	bw->buffs[1] = NULL;
	bw->blens[0] = 0;
	bw->blens[1] = 0;
	bw->bios[0] = open_memstream(bw->buffs + 0, bw->blens + 0);
	bw->bios[1] = open_memstream(bw->buffs + 1, bw->blens + 1);
	bw->out = NULL;
	bw->bidx = 0;
	bw->buf_size = buf_size? buf_size : 4 * 1024;
	bw->nbytes = 0;
	bw->lock = (pthread_mutex_t)PTHREAD_MUTEX_INITIALIZER;
	bw->running = 0;
	bw->flush = 0;
	if(pthread_create(&bw->pid, NULL, _buffered_writer_thread_func, bw) != 0){
		fprintf(stderr, " -- Failed to create thread [%s] in %s -- %s:%d --\n", "_buffered_writer_thread_func", __FUNCTION__, __FILE__, __LINE__);
		bw->pid = 0;
	}
	while(bw->running != 1){ nano_sleep(1); }
	return bw;
}

static inline BufferedWriter* open_bufferedwriter(FILE *out, size_t buf_size){
	return open2_bufferedwriter(out, _write_data_file, NULL, buf_size);
}

static inline BufferedWriter* zopen_bufferedwriter(FILE *out, size_t buf_size, int ncpu, int level){
	PGZF *pz;
	pz = open_pgzf_writer(out, buf_size, ncpu, level);
	return open2_bufferedwriter(pz, write_pgzf4filewriter, close_pgzf4filewriter, pz->bufsize);
}

static inline int beg_bufferedwriter(BufferedWriter *bw){
	if(bw->pid){
		while(bw->flush){ nano_sleep(1); }
		pthread_mutex_lock(&bw->lock);
		bw->out = bw->bios[bw->bidx];
		return 0;
	} else {
		bw->out = NULL;
		return 1; // error
	}
}

static inline int end_bufferedwriter(BufferedWriter *bw){
	if(bw->pid){
		pthread_mutex_unlock(&bw->lock);
	}
	bw->out = NULL;
	return 0;
}

static inline size_t flush_bufferedwriter(BufferedWriter *bw){
	size_t ret;
	if(bw->pid){
		pthread_mutex_unlock(&bw->lock);
		while(bw->flush == 1){ nano_sleep(1); }
		bw->flush = 1;
		while(bw->flush == 1){
			nano_sleep(1);
		}
		pthread_mutex_lock(&bw->lock);
		bw->flush = 0;
		bw->out = bw->bios[bw->bidx];
		ret = bw->nbytes;
	} else {
		ret = 0;
	}
	return ret;
}

static inline size_t close_bufferedwriter(BufferedWriter *bw){
	size_t ret;
	if(bw->pid){
		bw->running = 0;
		pthread_join(bw->pid, NULL);
	}
	fclose(bw->bios[0]);
	fclose(bw->bios[1]);
	if(bw->buffs[0]) free(bw->buffs[0]);
	if(bw->buffs[1]) free(bw->buffs[1]);
	if(bw->_close){
		bw->_close(bw->_file);
	}
	ret = bw->nbytes;
	free(bw);
	return ret;
}

#endif
