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

typedef struct {
	FILE *bios[2];
	FILE *out, *_out_;
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
			fwrite(bw->buffs[!bidx], 1, bsize[!bidx], bw->_out_);
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
		nano_sleep(1);
	}
	{
		bsize[0] = ftell(bw->bios[0]);
		bsize[1] = ftell(bw->bios[1]);
		fflush(bw->bios[0]);
		fflush(bw->bios[1]);
		bidx = bw->bidx;
		if(bsize[!bidx]){
			fwrite(bw->buffs[!bidx], 1, bsize[!bidx], bw->_out_);
			bw->nbytes += bsize[!bidx];
		}
		if(bsize[bidx]){
			fwrite(bw->buffs[bidx], 1, bsize[bidx], bw->_out_);
			bw->nbytes += bsize[bidx];
		}
	}
	return NULL;
}

static inline BufferedWriter* open_bufferedwriter(FILE *out, size_t buf_size){
	BufferedWriter *bw;
	bw = malloc(sizeof(BufferedWriter));
	bw->_out_ = out;
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
	return bw;
}

static inline void beg_bufferedwriter(BufferedWriter *bw){
	if(bw->pid){
		while(bw->flush){ nano_sleep(1); }
		pthread_mutex_lock(&bw->lock);
		bw->out = bw->bios[bw->bidx];
	} else {
		bw->out = bw->_out_;
	}
}

static inline void end_bufferedwriter(BufferedWriter *bw){
	if(bw->pid){
		pthread_mutex_unlock(&bw->lock);
	}
	bw->out = NULL;
}

static inline size_t flush_bufferedwriter(BufferedWriter *bw){
	size_t ret;
	bw->flush = 1;
	while(bw->flush == 1){
		nano_sleep(1);
	}
	pthread_mutex_lock(&bw->lock);
	bw->flush = 0;
	ret = bw->nbytes;
	pthread_mutex_unlock(&bw->lock);
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
	if(bw->blens[0]) free(bw->buffs[0]);
	if(bw->blens[1]) free(bw->buffs[1]);
	ret = bw->nbytes;
	free(bw);
	return ret;
}


#endif
