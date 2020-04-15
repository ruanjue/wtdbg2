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

#ifndef __FILEREADER_RJ_H
#define __FILEREADER_RJ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "chararray.h"
#include "mem_share.h"
#include "list.h"
#include "thread.h"
#include "pgzf.h"

#define BIOSEQ_ATTR_NULL	0
#define BIOSEQ_ATTR_TAG		1
#define BIOSEQ_ATTR_SEQ		2
#define BIOSEQ_ATTR_QLT		4
#define BIOSEQ_ATTR_FULL	7

typedef struct {
	String *tag, *seq, *dsc, *qlt;
	u4i attr;
} BioSequence;

#define FILEREADER_TYPE_NULL	0
#define FILEREADER_TYPE_FASTA	1
#define FILEREADER_TYPE_FASTQ	2
#define FILEREADER_TYPE_TEXT	3

#define FILEREADER_ATTR_NULL	0
#define FILEREADER_ATTR_NORMAL	1
#define FILEREADER_ATTR_STDIN	2
#define FILEREADER_ATTR_PROC	3
#define FILEREADER_ATTR_TEXT	4
#define FILEREADER_ATTR_USER	5 // defined by user

typedef size_t (*read_data_func)(void *obj, void *dat, size_t len);
typedef void (*close_input_func)(void *obj);

static inline size_t _read_data_file(void *obj, void *dat, size_t len){ return fread(dat, 1, len, (FILE*)obj); }
static inline void _close_input_file(void *obj){ if(obj) fclose((FILE*)obj); }
static inline void _close_input_proc(void *obj){ if(obj) pclose((FILE*)obj); }

typedef struct {
	int file_attr;
	char *filename;
	void *_file;
	read_data_func _read;
	close_input_func _close;
} file_src_t;
define_list_core(filesrcv, file_src_t, int, 0xFF);

typedef struct {
	filesrcv *files;
	int fidx;
	char *buffer[2];
	int ridx, widx, flag;
	u8i bufmax, bufoff, bufcnt[2];
#ifdef FR_USE_SPINLOCK
	pthread_spinlock_t lock;
#else
	pthread_mutex_t lock;
#endif
	char line_breaker;
	char delimiter;
	u8i  n_char, n_line;
	String *line, *line2;
	VStrv  *tabs;
	int rollback; // line will be re-used in next readline
	// thread
	pthread_t pid;
	int running;
	int eof;
} FileReader;

static inline BioSequence* init_biosequence(){
	BioSequence *seq;
	seq = malloc(sizeof(BioSequence));
	seq->tag = init_string(32);
	seq->seq = init_string(32);
	seq->dsc = init_string(32);
	seq->qlt = init_string(32);
	seq->attr = BIOSEQ_ATTR_FULL;
	return seq;
}

static inline void reset_biosequence(BioSequence *seq){
	clear_string(seq->tag);
	clear_string(seq->seq);
	clear_string(seq->dsc);
	clear_string(seq->qlt);
}

static inline void free_biosequence(BioSequence *seq){
	free_string(seq->tag);
	free_string(seq->seq);
	free_string(seq->dsc);
	free_string(seq->qlt);
	free(seq);
}

static inline void* file_src_thread_func(void *obj){
	FileReader *fr;
	file_src_t *fc;
	void *_file;
	read_data_func _read;
	close_input_func _close;
	size_t off, cnt, len;
	fr = (FileReader*)obj;
	while(fr->running){
		if(fr->fidx >= fr->files->size){
			fr->eof = 1;
			microsleep(1);
		} else {
			fr->eof = 0;
			fc = ref_filesrcv(fr->files, fr->fidx);
			_file = NULL;
			_read = NULL;
			_close = NULL;
			switch(fc->file_attr){
				case FILEREADER_ATTR_TEXT:
					len = strlen(fc->filename);
					off = 0;
					while(fr->running && len){
						while(fr->flag == 1 && fr->running){ nano_sleep(1); }
						cnt = num_min(len, fr->bufmax);
						memcpy(fr->buffer[fr->widx], fc->filename + off, cnt);
						fr->flag = 1;
						off += cnt;
						len -= cnt;
						fr->widx = !fr->widx;
					}
					break;
				case FILEREADER_ATTR_STDIN:
					if(_file == NULL){
						_file = fc->_file = stdin;
						_read = fc->_read = _read_data_file;
						_close = fc->_close = NULL;
					}
					// fall through
				case FILEREADER_ATTR_PROC:
					if(_file == NULL){
						_file = fc->_file = popen(fc->filename, "r");
						_read = fc->_read = _read_data_file;
						_close = fc->_close = _close_input_proc;
					}
					// fall through
				case FILEREADER_ATTR_USER:
					if(_file == NULL){
						_file = fc->_file;
						_read = fc->_read;
						_close = fc->_close;
					}
					// fall through
				default:
					if(_file == NULL){
						_file = fc->_file = open_file_for_read(fc->filename, NULL);
						_read = fc->_read = _read_data_file;
						_close = fc->_close = _close_input_file;
					}
					while(fr->running){
						while(fr->flag == 1){
							nano_sleep(1);
							if(fr->running == 0){
								break;
							}
						}
						if(fr->flag == 1) break;
						fr->bufcnt[fr->widx] = _read(_file, fr->buffer[fr->widx], fr->bufmax);
						fr->widx = !fr->widx;
						fr->flag = 1;
						if(fr->bufcnt[!fr->widx] == 0) break;
					}
			}
			if(_file && _close){
				_close(_file);
			}
			fr->fidx ++;
		}
	}
	return NULL;
}

static inline FileReader* init_filereader(){
	FileReader *fr;
	fr = malloc(sizeof(FileReader));
	fr->files = init_filesrcv(4);
	fr->fidx = 0;
	fr->bufmax = 128 * 1024;
	fr->bufoff = 0;
	fr->bufcnt[0] = 0;
	fr->bufcnt[1] = 0;
	fr->ridx = 0;
	fr->widx = 1;
	fr->flag = 0;
#ifdef FR_USE_SPINLOCK
	pthread_spin_init(&fr->lock, 0);
#else
	pthread_mutex_init(&fr->lock, NULL);
#endif
	fr->buffer[0] = malloc(fr->bufmax);
	fr->buffer[1] = malloc(fr->bufmax);
	fr->line_breaker = '\n';
	fr->delimiter = '\t';
	fr->n_char = 0;
	fr->n_line = 0;
	fr->line = init_string(32);
	fr->line2 = init_string(32);
	fr->tabs = init_VStrv(16);
	fr->rollback = 0;
	fr->pid  = 0;
	fr->running = 1;
	fr->eof = 0;
	return fr;
}

static inline void beg_asyn_filereader(FileReader *fr){
	if(pthread_create(&fr->pid, NULL, file_src_thread_func, fr) != 0){
		fprintf(stderr, " -- Failed to create thread [%s] in %s -- %s:%d --\n", "file_src_thread_func", __FUNCTION__, __FILE__, __LINE__);
		fr->pid = 0; // switch to directed read
	}
}

static inline void reset_filereader(FileReader *fr){
	if(fr->pid){
		fr->running = 0;
		pthread_join(fr->pid, NULL);
	}
	fr->fidx = 0;
	fr->bufoff = 0;
	fr->bufcnt[0] = 0;
	fr->bufcnt[1] = 0;
	fr->ridx = 0;
	fr->widx = 0;
	fr->flag = 0;
#ifdef FR_USE_SPINLOCK
	pthread_spin_destroy(&fr->lock);
	pthread_spin_init(&fr->lock, 0);
#else
	pthread_mutex_destroy(&fr->lock);
	pthread_mutex_init(&fr->lock, NULL);
#endif
	clear_string(fr->line);
	clear_VStrv(fr->tabs);
	fr->rollback = 0;
	fr->n_line = 0;
	fr->n_char = 0;
	fr->running = 1;
	fr->eof = 0;
	if(fr->pid){
		fr->pid = 0;
		beg_asyn_filereader(fr);
	}
}

static inline void free_filereader(FileReader *fr){
	file_src_t *f;
	int i;
	if(fr->pid){
		fr->running = 0;
		pthread_join(fr->pid, NULL);
	}
	for(i=0;i<fr->files->size;i++){
		f = ref_filesrcv(fr->files, i);
		if(f->filename) free(f->filename);
	}
#ifdef FR_USE_SPINLOCK
	pthread_spin_destroy(&fr->lock);
#else
	pthread_mutex_destroy(&fr->lock);
#endif
	free(fr->buffer[0]);
	free(fr->buffer[1]);
	free_filesrcv(fr->files);
	free_string(fr->line);
	free_string(fr->line2);
	free_VStrv(fr->tabs);
	free(fr);
}

static inline int push_filereader(FileReader *fr, char *filename){
	file_src_t *f;
	int len;
	f = next_ref_filesrcv(fr->files);
	f->_file = NULL;
	f->_read = NULL;
	f->_close = NULL;
	len = filename? strlen(filename) : 0;
	while(len && filename[len-1] == ' ') len --;
	if(len == 0 || strcmp(filename, "-") == 0){
		f->filename = NULL;
		f->file_attr = FILEREADER_ATTR_STDIN;
	} else if(filename[len-1] == '|'){
		f->filename = malloc(len);
		strncpy(f->filename, filename, len - 1);
		f->file_attr = FILEREADER_ATTR_PROC;
	} else if(len > 3 && strcmp(filename + len - 3, ".gz") == 0){
		//f->filename = malloc(len + 20);
		//sprintf(f->filename, "gzip -dc %s", filename);
		//f->file_attr = FILEREADER_ATTR_PROC;
		f->filename = strdup(filename);
		f->file_attr = FILEREADER_ATTR_USER;
		f->_file = open_pgzf_reader(open_file_for_read(f->filename, NULL), 0, 4);
		f->_read = read_pgzf4filereader;
		f->_close = close_pgzf4filereader;
	} else if(len > 5 && strcmp(filename + len - 5, ".pgzf") == 0){
		f->filename = strdup(filename);
		f->file_attr = FILEREADER_ATTR_USER;
		f->_file = open_pgzf_reader(open_file_for_read(f->filename, NULL), 0, 4);
		f->_read = read_pgzf4filereader;
		f->_close = close_pgzf4filereader;
	} else {
		f->filename = strdup(filename);
		f->file_attr = FILEREADER_ATTR_NORMAL;
	}
	return f->file_attr;
}

static inline int push_text_filereader(FileReader *fr, char *str, size_t len){
	file_src_t *f;
	UNUSED(len);
	f = next_ref_filesrcv(fr->files);
	f->_file = NULL;
	f->_read = NULL;
	f->_close = NULL;
	f->filename = str;
	f->file_attr = FILEREADER_ATTR_TEXT;
	return f->file_attr;
}

static inline int push_user_filereader(FileReader *fr, void *_file, read_data_func _read, close_input_func _close){
	file_src_t *f;
	f = next_ref_filesrcv(fr->files);
	f->_file = _file;
	f->_read = _read;
	f->_close = _close;
	f->filename = NULL;
	f->file_attr = FILEREADER_ATTR_USER;
	return f->file_attr;
}

static inline void push_all_filereader(FileReader *fr, int nfile, char **filenames){
	int i;
	for(i=0;i<nfile;i++) push_filereader(fr, filenames[i]);
}

// asyn: asynchronous reading
static inline FileReader* open_filereader(char *filename, int asyn){
	FileReader *fr;
	fr = init_filereader();
	push_filereader(fr, filename);
	if(asyn) beg_asyn_filereader(fr);
	return fr;
}

static inline FileReader* string_filereader(char *str, int asyn){
	int len;
	FileReader *fr;
	len = str? strlen(str) : 0;
	fr = init_filereader();
	push_text_filereader(fr, str, len);
	if(asyn) beg_asyn_filereader(fr);
	return fr;
}

static inline FileReader* open_all_filereader(int nfile, char **filenames, int asyn){
	FileReader *fr;
	fr = init_filereader();
	push_all_filereader(fr, nfile, filenames);
	if(asyn) beg_asyn_filereader(fr);
	return fr;
}

#define close_filereader(fr) free_filereader(fr)

static inline int asyn_readline_filereader(FileReader *fr, String *line){
	char *buffer;
	u8i i, nc;
	int ret;
	if(fr->rollback){
		fr->rollback = 0;
		return line->size + 1; // in case of end of file and not terminated by line_breaker, the return value is bigger by 1
	} else if(fr->eof && fr->bufoff == fr->bufcnt[fr->ridx]){
		return 0;
	} else {
		clear_string(line);
		nc = fr->n_char;
		while(1){
			buffer = fr->buffer[fr->ridx];
			ret = 0;
			for(i=fr->bufoff;i<fr->bufcnt[fr->ridx];){
				if(buffer[i++] == fr->line_breaker){
					ret = 1;
					break;
				}
			}
			fr->n_char += i - fr->bufoff;
			encap_string(line, i - fr->bufoff);
			append_string(line, buffer + fr->bufoff, i - fr->bufoff - ret);
			fr->bufoff = i;
			if(ret){
				return fr->n_char - nc;
			} else if(fr->eof){
				return fr->n_char - nc;
			}
			fr->bufoff = 0;
			fr->bufcnt[fr->ridx] = 0;
			while(fr->flag == 0){
				nano_sleep(1);
				if(fr->eof){
					if(fr->flag) break;
					else {
						return fr->n_char - nc;
					}
				}
			}
			fr->flag = 0;
			fr->ridx = !fr->ridx;
		}
		return 0;
	}
}

static inline int directed_readline_filereader(FileReader *fr, String *line){
	file_src_t *fc;
	void *_file;
	read_data_func _read;
	close_input_func _close;
	u8i i, nc;
	int ch;
	int ret;
	if(fr->eof) return 0;
	else if(fr->rollback){
		fr->rollback = 0;
		return line->size + 1; // in case of end of file and not terminated by line_breaker, the return value is bigger by 1
	}
	clear_string(line);
	nc = fr->n_char;
	while(fr->fidx < fr->files->size){
		fc = ref_filesrcv(fr->files, fr->fidx);
		_file = NULL;
		_read = NULL;
		_close = NULL;
		if(fr->flag == 0){
			switch(fc->file_attr){
				case FILEREADER_ATTR_TEXT:
					break;
				case FILEREADER_ATTR_STDIN:
					_file = fc->_file = stdin;
					_read = fc->_read = _read_data_file;
					_close = fc->_close = NULL;
					break;
				case FILEREADER_ATTR_PROC:
					_file = fc->_file = popen(fc->filename, "r");
					_read = fc->_read = _read_data_file;
					_close = fc->_close = _close_input_proc;
					break;
				case FILEREADER_ATTR_USER:
					_file = fc->_file;
					_read = fc->_read;
					_close = fc->_close;
					break;
				default:
					_file = fc->_file = open_file_for_read(fc->filename, NULL);
					_read = fc->_read = _read_data_file;
					_close = fc->_close = _close_input_file;
					break;
			}
			fr->flag = 1;
			fr->bufoff = 0;
			fr->bufcnt[0] = fr->bufcnt[1] = 0;
		} else {
			_file = fc->_file;
			_read = fc->_read;
			_close = fc->_close;
		}
		switch(fc->file_attr){
			case FILEREADER_ATTR_TEXT:
				ret = 0;
				for(i=fr->bufoff;fc->filename[i];){
					if(fc->filename[i++] == fr->line_breaker){
						ret = 1;
						break;
					}
				}
				fr->n_char += i - fr->bufoff;
				encap_string(line, i - fr->bufoff);
				append_string(line, fc->filename + fr->bufoff, i - fr->bufoff - ret);
				fr->bufoff = i;
				if(ret){
					break;
				}
				break;
			case FILEREADER_ATTR_STDIN:
				while((ch = fgetc(stdin)) != EOF){
					fr->n_char ++;
					if(ch == fr->line_breaker){
						break;
					}
					add_char_string(line, ch);
				}
				break;
			default:
				while(1){
					if(fr->bufoff >= fr->bufcnt[0]){
						fr->bufoff = 0;
						fr->bufcnt[0] = _read(_file, fr->buffer[0], fr->bufmax);
						if(fr->bufcnt[0] == 0) break;
					}
					ret = 0;
					for(i=fr->bufoff;i<fr->bufcnt[0];){
						if(fr->buffer[0][i++] == fr->line_breaker){
							ret = 1;
							break;
						}
					}
					fr->n_char += i - fr->bufoff;
					encap_string(line, i - fr->bufoff);
					append_string(line, fr->buffer[0] + fr->bufoff, i - fr->bufoff - ret);
					fr->bufoff = i;
					if(ret){
						break;
					}
				}
				break;
		}
		if(fr->n_char > nc){
			return fr->n_char - nc;
		} else {
			if(_file && _close){
				_close(_file);
			}
			fr->flag = 0;
			fr->fidx ++;
		}
	}
	fr->eof = 1;
	return 0;
}

int readline_filereader(FileReader *fr){
	int ret;
	ret = ((fr)->pid? asyn_readline_filereader(fr, (fr)->line) : directed_readline_filereader(fr, (fr)->line));
	if(ret > 0){
		fr->n_line ++;
	}
	return ret;
}

static inline void rollback_filereader(FileReader *fr){
	fr->rollback = 1;
	fr->n_line --;
}

static inline int split_line_filereader(FileReader *fr, char delimiter){
	VString *vs;
	int i;
	clear_VStrv(fr->tabs);
	vs = next_ref_VStrv(fr->tabs);
	vs->string = fr->line->string;
	vs->size = 0;
	for(i=0;i<fr->line->size;i++){
		if(fr->line->string[i] == delimiter){
			vs->size = fr->line->string + i - vs->string;
			vs = next_ref_VStrv(fr->tabs);
			vs->string = fr->line->string + i + 1;
			vs->size = 0;
		}
	}
	vs->size = fr->line->string + fr->line->size - vs->string;
	return (int)fr->tabs->size;
}

static inline int readtable_filereader(FileReader *fr){
	if(readline_filereader(fr) == 0) return -1;
	return split_line_filereader(fr, fr->delimiter);
}

static inline int get_col_len(FileReader *fr, int col){
	return fr->tabs->buffer[col].size;
}

static inline char* get_col_str(FileReader *fr, int col){
	VString *vs;
	vs = ref_VStrv(fr->tabs, col);
	vs->string[vs->size] = '\0';
	return vs->string;
}

static inline char* get_line_str(FileReader *fr){
	int i;
	for(i=0;i<fr->line->size;i++){
		if(fr->line->string[i] == 0){
			fr->line->string[i] = fr->delimiter;
		}
	}
	return fr->line->string;
}

// @return FILEREADER_TYPE_NULL (end of files), _FASTA, _FASTQ, or _TEXT (cannot parse sequence type)
static inline int readseq_filereader(FileReader *fr, BioSequence *seq){
	int n, i;
	do {
		if((n = readline_filereader(fr)) == 0) return FILEREADER_TYPE_NULL;
	} while(n == 0);
	reset_biosequence(seq);
	if(fr->line->string[0] == '>'){
		if(seq->attr & BIOSEQ_ATTR_TAG){
			for(i=1;i<fr->line->size;i++){
				if(fr->line->string[i] == ' ' || fr->line->string[i] == '\t') break;
			}
			append_string(seq->tag, fr->line->string + 1, i - 1);
			append_string(seq->dsc, fr->line->string + i, fr->line->size - i);
		}
		while((n = readline_filereader(fr))){
			if(fr->line->string[0] == '>'){
				rollback_filereader(fr);
				break;
			} else if(seq->attr & BIOSEQ_ATTR_SEQ){
				append_string(seq->seq, fr->line->string, fr->line->size);
			}
		}
		return FILEREADER_TYPE_FASTA;
	} else if(fr->line->string[0] == '@'){
		if(seq->attr & BIOSEQ_ATTR_TAG){
			for(i=1;i<fr->line->size;i++){
				if(fr->line->string[i] == ' ' || fr->line->string[i] == '\t') break;
			}
			append_string(seq->tag, fr->line->string + 1, i - 1);
			append_string(seq->dsc, fr->line->string + i, fr->line->size - i);
		}
		if((n = readline_filereader(fr))){
			if(seq->attr & BIOSEQ_ATTR_SEQ) append_string(seq->seq, fr->line->string, fr->line->size);
		} else {
			return FILEREADER_TYPE_FASTQ;
		}
		if((n = readline_filereader(fr))){
			// expected '+'
		} else {
			return FILEREADER_TYPE_FASTQ;
		}
		if((n = readline_filereader(fr))){
			if(seq->attr & BIOSEQ_ATTR_QLT) append_string(seq->qlt, fr->line->string, fr->line->size);
		} else {
			return FILEREADER_TYPE_FASTQ;
		}
		return FILEREADER_TYPE_FASTQ;
	} else {
		append_string(seq->dsc, fr->line->string, fr->line->size);
		return FILEREADER_TYPE_TEXT;
	}
}

#endif
