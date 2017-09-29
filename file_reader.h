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

#ifndef __FILE_READER_RJ_H
#define __FILE_READER_RJ_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "string.h"
#include "list.h"

/**
 * Sequence IO
 */

typedef struct {
	union { String name; String tag; };
	String header;
	String seq;
	String qual;
} Sequence;

typedef struct {
	FILE *file;
	char *filename;
	int is_proc;
} fr_file_t;
define_list(fr_filev, fr_file_t);

typedef struct {
	fr_filev *files;
	uint32_t fidx;
	char *buffer;
	int size;
	int capacity;
	int ptr;
	int last_brk;
	char line_breaker;
	char delimiter;
	unsigned long long n_line;
	String *line;
	String *vline;
	VStrv  *tabs;
	int seq_type;
} FileReader;

#define free_sequence(sequence) { if(sequence->name.string) free(sequence->name.string);\
	if(sequence->header.string) free(sequence->header.string);\
	if(sequence->seq.string) free(sequence->seq.string);\
	if(sequence->qual.string) free(sequence->qual.string);\
	free(sequence); }

FileReader* fopen_filereader(char *filename);

FileReader* fopen_filereader2(char *prefix, char *postfix);

FileReader* fopen_m_filereader(int n_file, char **file_names);

FileReader* stdin_filereader();

/**
 * Read characters from a copy of string
 */

FileReader* string_filereader(char *string);

void fclose_filereader(FileReader *fr);

size_t fseek_filereader(FileReader *fr, size_t pos);

int reset_filereader(FileReader *fr);

int fread_line(String *line, FileReader *fr);
int froll_back(FileReader *fr);

int fread_table(FileReader *fr);
#define get_col_str(fr, col) ref_VStrv((fr)->tabs, col)->string
#define get_col_len(fr, col) ref_VStrv((fr)->tabs, col)->size
#define col_length(fr, col) get_col_len(fr, (col) - 1)
#define col_string(fr, col) get_col_str(fr, (col) - 1)

typedef struct {
	int is_fq;
	int avg_seq_len;
	int min_seq_len;
	int max_seq_len;
} SeqFileAttr;

int guess_seq_file_type(FileReader *fr);

void guess_seq_file(FileReader *fr, SeqFileAttr *attr);

#define FASTA_FLAG_NORMAL		0
#define FASTA_FLAG_NO_NAME		1
#define FASTA_FLAG_NO_SEQ		2

int fread_fasta_adv(Sequence **seq, FileReader *fr, int flag);

#define fread_fasta(seq, fr) fread_fasta_adv(seq, fr, FASTA_FLAG_NORMAL)

#define FASTQ_FLAG_NORMAL		0
#define FASTQ_FLAG_NO_NAME		1
#define FASTQ_FLAG_NO_SEQ		2
#define FASTQ_FLAG_NO_QUAL		4

int fread_fastq_adv(Sequence **seq, FileReader *fr, int flag);

#define fread_fastq(seq, fr) fread_fastq_adv(seq, fr, FASTQ_FLAG_NORMAL)

#define SEQ_FLAG_NORMAL	0
#define SEQ_FLAG_NO_NAME	1
#define SEQ_FLAG_NO_SEQ	2
#define SEQ_FLAG_NO_QUAL	4

int fread_seq_adv(Sequence **seq, FileReader *fr, int flag);
#define fread_seq(seq, fr) fread_seq_adv(seq, fr, SEQ_FLAG_NORMAL)

char * fread_all(FileReader *fr);

static inline void print_pretty_seq(FILE *out, String *seq, int line_width){
	char c;
	int i, j;
	i = 0;
	while(i < seq->size){
		j = i + line_width;
		if(j > seq->size) j = seq->size;
		c  = seq->string[j];
		seq->string[j] = '\0';
		fprintf(out, "%s\n", seq->string + i);
		seq->string[j] = c;
		i = j;
	}
}

static inline void print_pretty_str(FILE *out, char *seq, int size, int line_width){
	char c;
	int i, j;
	i = 0;
	while(i < size){
		j = i + line_width;
		if(j > size) j = size;
		c  = seq[j];
		seq[j] = '\0';
		fprintf(out, "%s\n", seq + i);
		seq[j] = c;
		i = j;
	}
}

#endif
