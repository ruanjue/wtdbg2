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

#ifndef __MEM_SHARE_RJ_H
#define __MEM_SHARE_RJ_H

//#ifndef _GNU_SOURCE
//#define _GNU_SOURCE
//#endif
#if defined(__APPLE__) && defined(__MACH__)
#include <machine/endian.h>
#else
#include <endian.h>
#endif
#include <sys/stat.h>
#include <sys/mman.h>
//#include <sys/times.h>
#include <sys/time.h>
#include <sys/signal.h>
//#include <sys/resource.h>
//#include <sys/sysinfo.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <alloca.h>
#ifndef __HEADER_NO_EXECINFO
#include <execinfo.h>
#endif
#include <time.h>
#include "thread.h"

/**
 * Data types
 */
#define MAX_VALUE_U1	0xFFU
#define MAX_U1	MAX_VALUE_U1
typedef uint8_t	u1i;
#define MAX_VALUE_U2	0xFFFFU
#define MAX_U2	MAX_VALUE_U2
typedef uint16_t	u2i;
#define MAX_VALUE_U4	0xFFFFFFFFU
#define MAX_U4	MAX_VALUE_U4
typedef uint32_t	u4i;
#define MAX_VALUE_U8	0xFFFFFFFFFFFFFFFFLLU
#define MAX_U8	MAX_VALUE_U8
typedef unsigned long long u8i;
#define MAX_VALUE_B1	0x7F
#define MAX_B1	MAX_VALUE_B1
typedef int8_t	b1i;
#define MAX_VALUE_B2	0x7FFF
#define MAX_B2	MAX_VALUE_B2
typedef int16_t	b2i;
#define MAX_VALUE_B4	0x7FFFFFFF
#define MAX_B4	MAX_VALUE_B4
typedef int32_t	b4i;
#define MAX_VALUE_B8	0x7FFFFFFFFFFFFFFFLL
#define MAX_B8	MAX_VALUE_B8
typedef long long b8i;

typedef float f4i;
typedef long double f8i;

#define Int(x) ((int)(x))
#define UInt(x) ((unsigned)(x))
#define Int32(x) ((b4i)(x))
#define Int64(x) ((b8i)(x))
#define UInt64(x) ((u8i)(x))

/**
 * Useful functions
 */

#define num_min(n1, n2) (((n1) < (n2))? (n1) : (n2))
#define num_max(n1, n2) (((n1) > (n2))? (n1) : (n2))
#define num_diff(n1, n2) (((n1) < (n2))? ((n2) - (n1)) : ((n1) - (n2)))
#define num_cmp(a, b) (((a) > (b))? 1 : (((a) < (b))? -1 : -0))
#define num_cmpgt(a, b) ((a) > (b))
#define num_cmpx(a, b, c, d) (((a) > (b))? 1 : (((a) < (b))? -1 : (((c) > (d))? 1 : (((c) < (d))? -1 : 0))))
#define num_cmpgtx(a, b, c, d) (((a) > (b))? 1 : (((a) < (b))? 0 : (((c) > (d)))))
#define num_cmpxx(a, b, c, d, e, f) (((a) > (b))? 1 : (((a) < (b))? -1 : (((c) > (d))? 1 : (((c) < (d))? -1 : (((e) > (f))? 1 : (((e) < (f))? -1 : 0))))))
#define num_cmpgtxx(a, b, c, d, e, f) (((a) > (b))? 1 : (((a) < (b))? 0 : (((c) > (d))? 1 : (((c) < (d))? 0 : ((e) > (f))))))
#define num_abs(n) ((n) < 0? -(n) : (n))

#ifndef SWAP_TMP
#define SWAP_TMP
#define swap_tmp(a, b, t) { t = a; a = b; b = t; }
#endif

#define swap_var(a, b) { typeof(a) __swap_tmp__ = (a); (a) = (b); (b) = __swap_tmp__; }

#define UNUSED(x) (void)(x)

#define EXPR(...) __VA_ARGS__

#define _QUOTE_STR(x) #x
#define TOSTR(x) _QUOTE_STR(x)

#define uc(ch) (((ch) >= 'a' && (ch) <= 'z')? (ch) + 'A' - 'a' : (ch))
#define lc(ch) (((ch) >= 'A' && (ch) <= 'Z')? (ch) + 'a' - 'A' : (ch))

#define get_bit8(bits, idx) ((((bits)[(idx) >> 3]) >> ((idx) & 0x07)) & 0x01)
#define get_bit16(bits, idx) ((((bits)[(idx) >> 4]) >> ((idx) & 0x0F)) & 0x01)
#define get_bit32(bits, idx) ((((bits)[(idx) >> 5]) >> ((idx) & 0x1F)) & 0x01)
#define get_bit64(bits, idx) ((((bits)[(idx) >> 6]) >> ((idx) & 0x3F)) & 0x01)

#define get_2bit8(bits, idx) ((((bits)[(idx) >> 2]) >> (((idx) & 0x03) << 1)) & 0x03)
#define get_2bit16(bits, idx) ((((bits)[(idx) >> 3]) >> (((idx) & 0x07) << 1)) & 0x03)
#define get_2bit32(bits, idx) ((((bits)[(idx) >> 4]) >> (((idx) & 0x0F) << 1)) & 0x03)
#define get_2bit64(bits, idx) ((((bits)[(idx) >> 5]) >> (((idx) & 0x1F) << 1)) & 0x03)
#define inc_2bit64(bits, idx) (((bits)[(idx) >> 5]) += (1LLU << (((idx) & 0x1F) << 1)))

#define get_4bit8(bits, idx) ((((bits)[(idx) >> 1]) >> (((idx) & 0x01) << 2)) & 0x0F)
#define get_4bit16(bits, idx) ((((bits)[(idx) >> 2]) >> (((idx) & 0x03) << 2)) & 0x0F)
#define get_4bit32(bits, idx) ((((bits)[(idx) >> 3]) >> (((idx) & 0x07) << 2)) & 0x0F)
#define get_4bit64(bits, idx) ((((bits)[(idx) >> 4]) >> (((idx) & 0x0F) << 2)) & 0x0F)

#define get_8bit16(bits, idx) ((((bits)[(idx) >> 1]) >> (((idx) & 0x01) << 3)) & 0xFF)
#define get_8bit32(bits, idx) ((((bits)[(idx) >> 2]) >> (((idx) & 0x03) << 3)) & 0xFF)
#define get_8bit64(bits, idx) ((((bits)[(idx) >> 3]) >> (((idx) & 0x07) << 3)) & 0xFF)

#ifndef __HEADER_NO_EXECINFO
static inline void print_backtrace(FILE *out, int max_frame){
	void **buffer;
	int frames;
	if(max_frame < 1) max_frame = 1;
	buffer = malloc(sizeof(void*) * max_frame);
	frames = backtrace(buffer, max_frame);
	backtrace_symbols_fd(buffer, frames, fileno(out));
	free(buffer);
}
#else
static inline void print_backtrace(FILE *out, int max_frame){
	UNUSED(max_frame);
	fprintf(out, "-- ** Cannot back trace frames ** --\n");
}
#endif

static inline void num2bits(u8i num, char bits[64]){
	int i;
	for(i=0;i<64;i++){
		bits[i] = '0' + ((num >> (63 - i)) & 0x01);
	}
}

static inline size_t roundup_power2(size_t v){
	if(v == 0) return 0;
	v --;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v + 1;
}

static inline size_t roundup_times(size_t v, size_t base){
	return (v + base - 1) / base * base;
}

static inline void* malloc16(size_t size){
	u1i *p, *q;
	p = malloc(size + 16);
	if(p == NULL) return NULL;
	q = (u1i*)(((u8i)(p + 16)) & (~0xFLLU));
	*(q - 1) = q - p;
	return q;
}

static inline void free16(void *ptr){
	u1i *p, *q;
	q = (u1i*)ptr;
	p = q - *(q - 1);
	free(p);
}

static inline size_t encap_list(void **buffer, size_t e_size, size_t size, size_t cur_cap, size_t inc, int mem_zeros, size_t n_head){
	void *ptr;
	size_t cap;
	if(cur_cap - size >= inc) return cur_cap;
	if(MAX_U8 - inc <= size){
		fprintf(stderr, " -- Overflow(64bits) %llu + %llu, in %s -- %s:%d --\n", (u8i)size, (u8i)inc, __FUNCTION__, __FILE__, __LINE__);
		print_backtrace(stderr, 20);
		abort();
	}
	if(MAX_U8 - inc < 0x3FFFFFFFLLU){
		fprintf(stderr, " -- Overflow(64bits) %llu + %llu, in %s -- %s:%d --\n", (u8i)size, (u8i)inc, __FUNCTION__, __FILE__, __LINE__);
		print_backtrace(stderr, 20);
		abort();
	}
	if(size + inc < 0x3FFFFFFFLLU){
		cap = roundup_power2(size + inc);
	} else {
		cap = (size + inc + 0x3FFFFFFFLLU) & (MAX_U8 << 30);
	}
	ptr = realloc((*buffer) - n_head * e_size, e_size * (cap + n_head));
	if(ptr == NULL){
		fprintf(stderr, " -- Out of memory, try to allocate %llu bytes, old size %llu, old addr %p in %s -- %s:%d --\n", (unsigned long long)(e_size * (cap + n_head)), (unsigned long long)(e_size * (cur_cap + n_head)), *buffer, __FUNCTION__, __FILE__, __LINE__);
		print_backtrace(stderr, 20);
		abort();
	}
	*buffer = ptr + n_head * e_size;
	if(mem_zeros) memset((*buffer) + (cur_cap * e_size), 0, e_size * (cap - cur_cap));
	return cap;
}


#define ZEROS(e) memset((void*)(e), 0, sizeof(*(e)))

#ifndef __USE_GNU

# define TEMP_FAILURE_RETRY(expression) \
	({	\
		long int __result;	\
		do __result = (long int) (expression);	\
		while (__result == -1L && errno == EINTR);	\
		__result;	\
	})

#endif


static inline void nano_sleep(u8i nsec){
	struct timespec timeout;
	timeout.tv_sec  = nsec/1000000000;
	timeout.tv_nsec = nsec%1000000000;
	nanosleep(&timeout, NULL);
}

#define micro_sleep(usec) nano_sleep(((u8i)(usec)) * 1000LLU)
#define microsleep(usec) micro_sleep(usec)

static inline long long microtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((long long)tv.tv_sec) * 1000000 + tv.tv_usec;
}

static inline char* date(){
	time_t tm;
	char *dstr, *p;
	tm = time(NULL);
	p = dstr = asctime(localtime(&tm));
	while(*p){
		if(*p == '\n'){ *p = 0; break; }
		p ++;
	}
	return dstr;
}

static inline int replace_char(char *str, char src, char dst, int max){
	int ret;
	ret = 0;
	while(*str){
		if((*str) == src){
			(*str) = dst;
			ret ++;
			if(ret == max) break;
		}
		str ++;
	}
	return ret;
}

static inline int file_exists(const char *filename){
	char *rpath;
	struct stat s;
	//realpath = canonicalize_file_name(filename);
	rpath = realpath(filename, NULL);
	if(rpath == NULL) return 0;
	if(stat(rpath, &s) == -1){ free(rpath); return 0; }
	free(rpath);
	switch(s.st_mode & S_IFMT){
		//case S_IFBLK:
		//case S_IFCHR:
		//case S_IFDIR:
		//case S_IFIFO:
		//case S_IFSOCK:
		case S_IFLNK:
		case S_IFREG: return 1;
		default: return 0;
	}
}

static inline int dir_exists(const char *filename){
	char *rpath;
	struct stat s;
	//realpath = canonicalize_file_name(filename);
	rpath = realpath(filename, NULL);
	if(rpath == NULL) return 0;
	if(stat(rpath, &s) == -1){ free(rpath); return 0; }
	free(rpath);
	switch(s.st_mode & S_IFMT){
		//case S_IFBLK:
		//case S_IFCHR:
		//case S_IFREG:
		//case S_IFIFO:
		//case S_IFSOCK:
		case S_IFLNK:
		case S_IFDIR: return 1;
		default: return 0;
	}
}

static inline char* relative_filename(char *filename){
	char *ptr;
	if(filename == NULL) return NULL;
	ptr = filename + strlen(filename);
	while(ptr >= filename){
		if(*ptr == '/') break;
		ptr --;
	}
	return strdup(ptr + 1);
}

static inline char* absolute_filename(char *filename){
	char *path, *cwd, *ptr;
	int x, y, z, i;
	if(filename == NULL) return NULL;
	if(filename[0] == '/') return strdup(filename);
	cwd = getcwd(NULL, 0);
	path = malloc(strlen(cwd) + strlen(filename) + 2);
	x = 0;
	y = 0;
	z = 0;
	ptr = filename;
	while(*ptr){
		if((*ptr) == '.'){
			z ++;
			if(z > 2) break;
		} else if((*ptr) == '/'){
			y = ptr - filename + 1;
			if(z == 2) x ++;
			z = 0;
		} else {
			break;
		}
		ptr ++;
	}
	z = 0;
	i = strlen(cwd);
	while(x && i){
		while(i){
			i --;
			if(cwd[i] == '/') break;
		}
		x --;
	}
	if(x || i == 0){
		fprintf(stderr, " -- BAD File name: '%s'. PWD = '%s' --\n", filename, cwd); fflush(stderr);
		free(cwd);
		return NULL;
	}
	strncpy(path, cwd, i);
	free(cwd);
	path[i] = '/';
	strncpy(path + i + 1, filename + y, strlen(filename) - y);
	path[i + 1 + strlen(filename) - y] = 0;
	return path;
}

static inline int exists_file(char *dir, char *filename){
	char *realpath, *fullname;
	int ret;
	realpath = absolute_filename(dir? dir : ".");
	if(!dir_exists(realpath)){ free(realpath); return 0; }
	fullname = malloc(strlen(realpath) + strlen(filename) + 3);
	sprintf(fullname, "%s/%s", realpath, filename);
	free(realpath);
	ret = file_exists(fullname);
	free(fullname);
	return ret;
}

static inline FILE* open_file_for_read(char *name, char *suffix){
	char *full_name;
	FILE *file;
	if(name == NULL && suffix == NULL){
		full_name = "-";
	} else if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	if(strcmp(full_name, "-") == 0){
		file = stdin;
	} else {
		file = fopen(full_name, "r");
	}
	if(file == NULL){
		fprintf(stderr, "Cannot open file for read: %s\n", full_name);
		perror(NULL);
		exit(1);
	}
	return file;
}

static inline FILE* open_file_for_write(char *name, char *suffix, int overwrite){
	char *full_name;
	FILE *file;
	if(name == NULL && suffix == NULL){
		full_name = "-";
	} else if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	if(strcmp(full_name, "-") == 0){
		file = stdout;
	} else if(!overwrite && file_exists(full_name)){
		fprintf(stderr, "File exists: %s\n", full_name); exit(1);
	} else {
		file = fopen(full_name, "w+");
	}
	if(file == NULL){
		fprintf(stderr, "Cannot open file for write: %s\n", full_name);
		perror(NULL);
		exit(1);
	}
	return file;
}

static inline FILE* open_file_for_append(char *name, char *suffix){
	char *full_name;
	FILE *file;
	if(name == NULL && suffix == NULL){
		full_name = "-";
	} else if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	if(strcmp(full_name, "-") == 0){
		file = stdout;
	} else {
		file = fopen(full_name, "a+");
	}
	if(file == NULL){
		fprintf(stderr, "Cannot open file for append: %s\n", full_name);
		perror(NULL);
		exit(1);
	}
	return file;
}

static inline void close_file(FILE *file){
	if(file == NULL) return;
	if(file == stdin || file == stdout || file == stderr) return;
	if(fclose(file)) perror("Error on close file");
}

#define STEP_WISE_BUFF_SIZE	4096
static inline size_t fread_stepwise(void *_ptr, size_t size, size_t ne, FILE *stream){
	void *ptr;
	size_t ret, len, n, sz;
	ret = 0;
	ptr = _ptr;
	sz = STEP_WISE_BUFF_SIZE / size;
	if(sz < 1) sz = 1;
	while(ne){
		len = num_min(ne, sz);
		if((n = fread(ptr, size, len, stream)) == 0){
			break;
		}
		ret += n;
		ptr += n;
		ne  -= n;
	}
	return ret;
}

//static inline u8i total_sys_memory(){ return getpagesize() * get_phys_pages(); }
//static inline u8i avail_sys_memory(){ return getpagesize() * get_avphys_pages(); }
//static inline u8i proc_maxrss_memory(){ struct rusage ru; if(getrusage(RUSAGE_SELF, &ru) != -1){ return ru.ru_maxrss * 1024; } else { return 0; } }

/*
* runtime memory usage, adopt from Heng's runlib
*/
static inline int get_linux_sys_info(u8i *memtotal, u8i *memavail, int *ncpu){
	FILE *fp;
	char buffer[64];
	u8i freed, buffered, cached;
	int core;
	if(memtotal || memavail){
		freed = buffered = cached = 0;
		fp = open_file_for_read("/proc/meminfo", NULL);
		while((fscanf(fp, "%s", buffer)) > 0){
			if (strstr(buffer, "MemTotal") == buffer) {
				fscanf(fp, "%llu", memtotal);
				(*memtotal) *= 1024;
			} else if (strstr(buffer, "MemFree") == buffer) {
				fscanf(fp, "%llu", &freed);
			} else if (strstr(buffer, "Buffers") == buffer) {
				fscanf(fp, "%llu", &buffered);
			} else if (strstr(buffer, "Cached") == buffer) {
				fscanf(fp, "%llu", &cached);
			}
		}
		fclose(fp);
		if(memavail){
			*memavail = (freed + buffered + cached) * 1024;
		}
	}
	if(ncpu){
		fp = open_file_for_read("/proc/cpuinfo", NULL);
		core = 0;
		while((fscanf(fp, "%s", buffer)) > 0){
			if (strstr(buffer, "processor") == buffer) {
				core ++;
			}
		}
		fclose(fp);
		*ncpu = core;
	}
	return 0;
}

static inline int get_linux_proc_info(u8i *rss, u8i *vsize, double *utime, double *stime){
	int c, n_spc;
	char str[64];
	FILE *fp;
	unsigned long tmp, tmp2;
	size_t page_size;
	*rss = *vsize = 0;
	*utime = *stime = 0;
	//page_size = getpagesize();
	page_size = sysconf(_SC_PAGESIZE);
	sprintf(str, "/proc/%u/stat", getpid());
	fp = open_file_for_read(str, NULL);
	n_spc = 0;
	while ((c = fgetc(fp)) != EOF) {
		if (c == ' ') ++n_spc;
		if (n_spc == 13) break;
	}
	fscanf(fp, "%lu%lu", &tmp, &tmp2);
	*utime = tmp / 100.0;
	*stime = tmp2 / 100.0;
	++n_spc;
	while ((c = fgetc(fp)) != EOF) {
		if (c == ' ') ++n_spc;
		if (n_spc == 22) break;
	}
	fscanf(fp, "%lu%lu", &tmp, &tmp2);
	fclose(fp);
	*vsize = tmp;
	*rss = tmp2 * page_size;
	return 0;
}

thread_beg_def(_proc_deamon);
u8i memtotal, memavail;
int ncpu;
u8i max_rss, max_vsz;
u8i rss_limit;
double utime, stime, rtime;
double rtime_limit;
int interval;
thread_end_def(_proc_deamon);

static struct _proc_deamon_struct* _sig_proc_deamon = NULL;
static inline void print_proc_stat_info(int signum){
	FILE *log;
	if(_sig_proc_deamon == NULL) return;
	log = stderr;
	thread_beg_syn_read(_sig_proc_deamon);
	fprintf(log, "** PROC_STAT(%d) **: real %.3f sec, user %.3f sec, sys %.3f sec, maxrss %.1f kB, maxvsize %.1f kB\n", signum, _sig_proc_deamon->rtime, _sig_proc_deamon->utime, _sig_proc_deamon->stime,  _sig_proc_deamon->max_rss / 1024.0,  _sig_proc_deamon->max_vsz / 1024.0);
	thread_end_syn_read(_sig_proc_deamon);
	//if(signum > 0) fprintf(log, "-- noticed by signal %d\n", signum);
	//else fprintf(log, "-- noticed by program\n");
	//fprintf(log, "-- real         %16.3f sec\n", _sig_proc_deamon->rtime);
	//fprintf(log, "-- user         %16.3f sec\n", _sig_proc_deamon->utime);
	//fprintf(log, "-- sys          %16.3f sec\n", _sig_proc_deamon->stime);
	//fprintf(log, "-- maxrss       %16.3f kB\n",  _sig_proc_deamon->max_rss / 1024.0);
	//fprintf(log, "-- maxvsize     %16.3f kB\n",  _sig_proc_deamon->max_vsz / 1024.0);
	//fprintf(log, "--\n");
	fflush(log);
}

static inline void _deamon_config_proc_limit(int signum){
	char *val;
	u8i rss_limit, rtime_limit;
	UNUSED(signum);
	if(_sig_proc_deamon == NULL) return;
	thread_beg_syn_write(_sig_proc_deamon);
	val = getenv("LIMIT_RSS");
	if(val){
		rss_limit = atol(val) * 1024 * 1024;
	} else {
		rss_limit = 0;
	}
	val = getenv("LIMIT_RTIME");
	if(val){
		rtime_limit = atol(val);
	} else {
		rtime_limit = 0;
	}
	_sig_proc_deamon->rss_limit = rss_limit;
	_sig_proc_deamon->rtime_limit = rtime_limit;
	fprintf(stderr, "** PROC_LIMIT: max rss = %llu MB rtime = %llu secs. **\n", rss_limit / 1024 * 1024, rtime_limit);
	thread_end_syn_write(_sig_proc_deamon);
}

thread_beg_func(_proc_deamon);
u8i rss, vsz;
u8i cb, ce;
cb = microtime();
get_linux_sys_info((u8i*)&_proc_deamon->memtotal, (u8i*)&_proc_deamon->memavail, (int*)&_proc_deamon->ncpu);
_proc_deamon->max_rss = 0;
_proc_deamon->max_vsz = 0;
_proc_deamon->rtime = 0;
_proc_deamon->rss_limit = 0;
_proc_deamon->rtime_limit = 0;
_proc_deamon->interval = 100000; // 0.1 sec
get_linux_proc_info(&rss, &vsz, (double*)&_proc_deamon->utime, (double*)&_proc_deamon->stime);
_sig_proc_deamon = (struct _proc_deamon_struct*)_proc_deamon;
signal(SIGUSR1, _deamon_config_proc_limit);
signal(SIGUSR2, print_proc_stat_info);
_proc_deamon->once = 0; // Don't set ->state to 0 at thread_end_loop
thread_beg_loop(_proc_deamon);
thread_beg_syn_write(_proc_deamon);
get_linux_proc_info(&rss, &vsz, (double*)&_proc_deamon->utime, (double*)&_proc_deamon->stime);
if(rss > _proc_deamon->max_rss) _proc_deamon->max_rss = rss;
if(vsz > _proc_deamon->max_vsz) _proc_deamon->max_vsz = vsz;
ce = microtime();
_proc_deamon->rtime = (ce - cb) / 1000000.0;
thread_end_syn_write(_proc_deamon);
if(_proc_deamon->rss_limit && rss > _proc_deamon->rss_limit){
	fprintf(stderr, "-- Exceed memory limit: %llu > %llu\n", rss, _proc_deamon->rss_limit);
	abort();
}
if(_proc_deamon->rtime_limit){
	if(_proc_deamon->rtime > _proc_deamon->rtime_limit){
		fprintf(stderr, "-- Timeout: %16.3f > %16.3f\n", _proc_deamon->rtime, _proc_deamon->rtime_limit);
		abort();
	}
}
microsleep(_proc_deamon->interval); // 0.1 sec
thread_end_loop(_proc_deamon);
get_linux_proc_info(&rss, &vsz, (double*)&_proc_deamon->utime, (double*)&_proc_deamon->stime);
if(rss > _proc_deamon->max_rss) _proc_deamon->max_rss = rss;
if(vsz > _proc_deamon->max_vsz) _proc_deamon->max_vsz = vsz;
ce = microtime();
_proc_deamon->rtime = (ce - cb) / 1000000.0;
thread_end_func(_proc_deamon);

#if defined(__linux__) || defined(__unix__) || defined(__CYGWIN__)
#define BEG_STAT_PROC_INFO(log, argc, argv)	\
thread_preprocess(_proc_deamon);	\
thread_beg_init(_proc_deamon, 1);	\
thread_end_init(_proc_deamon);	\
thread_wake_all(_proc_deamon);	\
if(log){	\
	fprintf(log, "--\n");	\
	fprintf(log, "-- total memory %16.1f kB\n", _proc_deamon->memtotal / 1024.0);	\
	fprintf(log, "-- available    %16.1f kB\n", _proc_deamon->memavail / 1024.0);	\
	fprintf(log, "-- %d cores\n", _proc_deamon->ncpu);	\
	int i;	\
	char **opts;	\
	opts = (char**)argv;	\
	fprintf(log, "-- Starting program:");	\
	for(i=0;i<argc;i++) fprintf(log, " %s", opts[i]);	\
	fprintf(log, "\n");	\
	fprintf(log, "-- pid          %16d\n", getpid());	\
	fprintf(log, "-- date         %s\n", date());	\
	fprintf(log, "--\n");	\
	fflush(log);	\
}
#else
#define BEG_STAT_PROC_INFO(log, argc, argv)	\
if(log){	\
	int i;	\
	char **opts;	\
	opts = (char**)argv;	\
	fprintf(log, "-- Starting program:");	\
	for(i=0;i<argc;i++) fprintf(log, " %s", opts[i]);	\
	fprintf(log, "\n");	\
	fflush(log);	\
}
#endif

#if defined(__linux__) || defined(__unix__) || defined(__CYGWIN__)
#define SET_PROC_LIMIT(memory, time)	\
thread_beg_iter(_proc_deamon);	\
_proc_deamon->rss_limit = memory;	\
_proc_deamon->rtime_limit = time;	\
thread_end_iter(_proc_deamon)
#else
#define SET_PROC_LIMIT(memory, time)
#endif

#if defined(__linux__) || defined(__unix__) || defined(__CYGWIN__)
#define END_STAT_PROC_INFO(log)	\
_proc_deamon->once = 1;	\
thread_beg_close(_proc_deamon);	\
if(log){	\
	fprintf(log, "** PROC_STAT(TOTAL) **: real %.3f sec, user %.3f sec, sys %.3f sec, maxrss %.1f kB, maxvsize %.1f kB\n---\n", _sig_proc_deamon->rtime, _sig_proc_deamon->utime, _sig_proc_deamon->stime,  _sig_proc_deamon->max_rss / 1024.0,  _sig_proc_deamon->max_vsz / 1024.0);	\
}	\
thread_end_close(_proc_deamon)
#else
#define END_STAT_PROC_INFO(log)	\
if(log){	\
	fprintf(log, "-- proc stat failed: only support 'linux, unix, cygwin' OS --\n");	\
}
#endif



// Object I/O

/**
 * Careful: when two pointers refer to the same addr of a object, we cannot make sure of it after mem_load
 */

#ifndef OBJ_DESC_MAX_CHILD
#define OBJ_DESC_MAX_CHILD 64
#endif

static inline size_t mem_size_round(size_t size){ return (size + 7) & 0xFFFFFFFFFFFFFFF8LLU; }

static inline uint8_t mem_size_gap(size_t size){ return (size & 0x07U)? 8 - (size & 0x07U) : 0; }

static inline size_t mem_dump(void *mem, size_t len, FILE *out){
	size_t size;
	uint8_t i, v;
	if(mem == NULL) return 0;
	size = mem_size_round(len);
	if(out){
		fwrite(mem, 1, len, out);
		v = 0;
		for(i=0;i<mem_size_gap(len);i++) fwrite(&v, 1, 1, out);
	}
	return size;
}

typedef size_t (*mem_array_count)(void *obj, int idx);
typedef void (*mem_load_post)(void *obj, size_t aux_data);

struct obj_desc_t;

#define MEM_PTR_TYPE_DUMP	1 // whether it take bytes
#define MEM_PTR_TYPE_POINTER	2 // take desc->size or sizeof(void*) bytes

typedef struct obj_desc_t {
	const char *tag;
	size_t size; // If size = 0, it is virtual, mem_type will be OR with MEM_PTR_TYPE_DUMP. See OBJ_DESC_CHAR_ARRAY
	int n_child; // <= OBJ_DESC_MAX_CHILD.
	uint8_t mem_type[OBJ_DESC_MAX_CHILD];
	off_t  addr[OBJ_DESC_MAX_CHILD]; // offsetof(type, field)
	const struct obj_desc_t *desc[OBJ_DESC_MAX_CHILD];
	mem_array_count cnt;
	mem_load_post post;
} obj_desc_t;

// Basic obj_desc_t, size = 1 byte
static const struct obj_desc_t OBJ_DESC_DATA = {"OBJ_DESC_DATA", 1, 0, {}, {}, {}, NULL, NULL};
// Special obj_desc_t for string, set mem_type=0 and addr=0 to call the _char_array_obj_desc_cnt on itself
// so that we know the length of string, then set size=0 to indicate that it is an virtual reference, program should add MEM_PTR_TYPE_DUMP to its mem_type
static inline size_t _char_array_obj_desc_cnt(void *obj, int idx){ if(idx == 0 && obj) return strlen((char*)obj) + 1; else return 0; }
static const struct obj_desc_t OBJ_DESC_CHAR_ARRAY = {"OBJ_DESC_CHAR_ARRAY", 0, 1, {0}, {0}, {&OBJ_DESC_DATA}, _char_array_obj_desc_cnt, NULL};

static inline size_t mem_size_obj(void *obj, uint8_t mem_type, const obj_desc_t *desc, size_t size, size_t cnt){
	size_t m;
	void *ref;
	int i;
	if(desc == NULL) return size;
	if(obj == NULL) return size;
	switch(mem_type){
		case 3: size += mem_size_round(sizeof(void*) * cnt);
		// fall through
		case 2:
			for(m=0;m<cnt;m++) if(((void**)obj)[m]) size += mem_size_round(desc->size);
			break;
		case 1: size += mem_size_round(cnt * desc->size); // TODO: if desc == &OBJ_DESC_DATA, mem_size_round may waste many memory
		// fall through
		case 0: break;
	}
	if(desc->n_child == 0) return size;
	for(m=0;m<cnt;m++){
		if(mem_type & 0x02){
			ref = ((void**)obj)[m];
			if(ref == NULL) continue;
		} else {
			ref = obj + m * desc->size;
		}
		for(i=0;i<desc->n_child;i++){
			if(desc->mem_type[i] & 0x01){
				size += mem_size_obj(*((void**)(ref + desc->addr[i])), desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], 0, desc->cnt? desc->cnt(ref, i) : 1);
			} else {
				size += mem_size_obj(ref + desc->addr[i], desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], 0, desc->cnt? desc->cnt(ref, i) : 1);
			}
		}
	}
	return size;
}

static inline void mem_dup_obj(void **ret, void *obj, size_t aux_data, uint8_t mem_type, const obj_desc_t *desc, size_t cnt){
	size_t m;
	void *ref, *chd;
	int i;
	if(desc == NULL || obj == NULL){ *ret = NULL; return; }
	if(mem_type & 0x01){
		if(mem_type & 0x02){
			*ret = calloc(cnt, sizeof(void*) * cnt);
		} else if(desc->size){
			*ret = calloc(cnt, desc->size);
			memcpy(*ret, obj, desc->size * cnt);
		} else {
			// OBJ_DESC_CHAR_ARRAY
		}
	}
	if(desc->n_child == 0 && desc->post == NULL){
	//if((mem_type & 0x02) == 0 && desc->n_child == 0){
	} else {
		for(m=0;m<cnt;m++){
			if(mem_type & 0x02){
				ref = ((void**)obj)[m];
				if(ref == NULL) continue;
				chd = malloc(desc->size);
				memcpy(chd, ref, desc->size);
				ret[m] = chd;
			} else if(desc->size){
				ref = obj + m * desc->size;
				chd = (*ret) + m * desc->size;
			} else {
				ref = obj;
				chd = ret;
			}
			for(i=0;i<desc->n_child;i++){
				if(desc->mem_type[i] & 0x01){
					mem_dup_obj((void**)(chd + desc->addr[i]), *((void**)(ref + desc->addr[i])), aux_data, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], desc->cnt? desc->cnt(ref, i) : 1);
				} else {
					mem_dup_obj((void**)(chd + desc->addr[i]), ref + desc->addr[i], aux_data, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], desc->cnt? desc->cnt(ref, i) : 1);
				}
			}
			if(desc->post) desc->post(chd, aux_data);
		}
	}
}

static inline size_t mem_dump_obj(void *obj, uint8_t mem_type, const obj_desc_t *desc, size_t offset, size_t cnt, FILE *out, int mem_free){
	void *ref;
	size_t size, m;
	int i;
	if(obj == NULL) return offset;
	size = offset;
	if(mem_type & 0x01){
		if(mem_type & 0x02){
			size += mem_dump(obj, cnt * sizeof(void*), out);
		} else {
			size += mem_dump(obj, cnt * desc->size, out);
		}
	}
	if((mem_type & 0x02) == 0 && desc->n_child == 0){
	} else {
		for(m=0;m<cnt;m++){
			if(mem_type & 0x02){
				ref = ((void**)obj)[m];
				if(ref == NULL) continue;
				size += mem_dump(ref, desc->size, out);
			} else {
				ref = obj + m * desc->size;
			}
			for(i=0;i<desc->n_child;i++){
				if(desc->mem_type[i] & 0x01){
					size = mem_dump_obj(*((void**)(ref + desc->addr[i])), desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], size, desc->cnt? desc->cnt(ref, i) : 1, out, mem_free);
				} else {
					size = mem_dump_obj(ref + desc->addr[i], desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], size, desc->cnt? desc->cnt(ref, i) : 1, out, mem_free);
				}
			}
		}
	}
	if(mem_free){
		if((mem_type & 0x02) && desc->size) for(m=0;m<cnt;m++) free(((void**)obj)[m]);
		if((mem_type & 0x01) && ((mem_type & 0x02) | desc->size)) free(obj);
	}
	return size;
}

static inline size_t mem_load_obj(void *obj, size_t aux_data, uint8_t mem_type, const obj_desc_t *desc, size_t addr_beg, size_t cnt){
	size_t addr, m;
	int i;
	void *ref, **ptr;
	if(obj == NULL) return 0;
	addr = addr_beg? : (size_t)obj;
	if(mem_type & 0x01){
		if(mem_type & 0x02){
			addr += mem_size_round(cnt * sizeof(void*));
		} else {
			addr += mem_size_round(cnt * desc->size);
		}
	}
	if(desc->n_child == 0 && desc->post == NULL){
		if(mem_type == 2 || mem_type == 3){
			for(m=0;m<cnt;m++){
				ptr = ((void**)obj) + m;
				if(*ptr == NULL) continue;
				*ptr = (void*)addr;
				addr += mem_size_round(desc->size);
			}
		}
		return addr;
	}
	for(m=0;m<cnt;m++){
		if(mem_type & 0x02){
			ptr = ((void**)obj) + m;
			if(*ptr == NULL) continue;
			ref = *ptr = (void*)addr;
			addr += mem_size_round(desc->size);
		} else {
			ref = obj + m * desc->size;
		}
		for(i=0;i<desc->n_child;i++){
			ptr = ref + desc->addr[i];
			if(desc->mem_type[i] & 0x01){
				if(*ptr == NULL) continue;
				*ptr = (void*)addr;
				addr = mem_load_obj(*ptr, aux_data, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], addr, desc->cnt? desc->cnt(ref, i) : 1);
			} else {
				addr = mem_load_obj((void*)ptr, aux_data, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], addr, desc->cnt? desc->cnt(ref, i) : 1);
			}
		}
		if(desc->post) desc->post(ref, aux_data);
	}
	return addr;
}

static inline size_t mem_tree_obj(FILE *out, void *obj, size_t mem_type, const obj_desc_t *desc, size_t addr, size_t cnt, size_t max_cnt, int level, int max_level){
	size_t m;
	int i, j;
	void *ref, **ptr;
	if(obj == NULL) return 0;
	if(addr == 0) addr = (size_t)obj;
	if(max_level <= 0 || level <= max_level){
		for(j=0;j<level;j++) fputc('-', out);
		fprintf(out, "%s:%p:%llu:%llu:", desc->tag, obj, (u8i)mem_type, (u8i)cnt);
	}
	if(mem_type & 0x01){
		if(mem_type & 0x02){
			addr += mem_size_round(cnt * sizeof(void*));
		} else {
			addr += mem_size_round(cnt * desc->size);
		}
	}
	if(desc->n_child == 0){
		if(mem_type == 3){
			for(m=0;m<cnt;m++){
				ptr = ((void**)obj) + m;
				if(*ptr == NULL) continue;
				addr += mem_size_round(desc->size);
			}
		}
		if(max_level <= 0 || level <= max_level) fprintf(out, "%llu\n", (u8i)(addr - (size_t)obj));
		return addr;
	}
	if(max_level <= 0 || level <= max_level) fprintf(out, "%llu\n", (u8i)(addr - (size_t)obj));
	for(m=0;m<cnt;m++){
		if(cnt > 1 && (max_level <= 0 || level <= max_level) && (max_cnt == 0 || m < max_cnt)){
			for(j=0;j<level;j++) fputc('-', out);
			fprintf(out, "[%llu/%llu]\n", (u8i)m, (u8i)cnt);
		}
		if(mem_type & 0x02){
			ptr = ((void**)obj) + m;
			if(*ptr == NULL) continue;
			ref = (void*)addr;
			addr += mem_size_round(desc->size);
		} else {
			ref = obj + m * desc->size;
		}
		for(i=0;i<desc->n_child;i++){
			ptr = ref + desc->addr[i];
			if(desc->mem_type[i] & 0x01){
				if(*ptr == NULL) continue;
				if(max_cnt && m >= max_cnt){
					addr = mem_tree_obj(out, (void*)addr, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], addr, desc->cnt? desc->cnt(ref, i) : 1, max_cnt, 2, 1);
				} else {
					addr = mem_tree_obj(out, (void*)addr, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], addr, desc->cnt? desc->cnt(ref, i) : 1, max_cnt, level + 1, max_level);
				}
			} else {
				if(max_cnt && m >= max_cnt){
					addr = mem_tree_obj(out, (void*)ptr, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], addr, desc->cnt? desc->cnt(ref, i) : 1, max_cnt, 2, 1);
				} else {
					addr = mem_tree_obj(out, (void*)ptr, desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP), desc->desc[i], addr, desc->cnt? desc->cnt(ref, i) : 1, max_cnt, level + 1, max_level);
				}
			}
		}
	}
	return addr;
}

static inline const obj_desc_t* mem_locate_obj(void *obj, size_t *_mem_type, const obj_desc_t *desc, u4i *trace_childs, size_t *trace_cnts, u4i ntrace, size_t *addr, size_t *cnt, int selected){
	const obj_desc_t *ret;
	size_t m, mem_type;
	u4i ns;
	int i, nc;
	void *ref, **ptr;
	if(obj == NULL) return desc;
	if(*addr == 0) *addr = (size_t)obj;
	if(selected && ntrace == 0) return desc;
	mem_type = *_mem_type;
	if(mem_type & 0x01){
		if(mem_type & 0x02){
			*addr += mem_size_round((*cnt) * sizeof(void*));
		} else {
			*addr += mem_size_round((*cnt) * desc->size);
		}
	}
	if(selected){
		if(trace_childs[0] >= (u4i)desc->n_child){
			fprintf(stderr, " -- Illegal trace child (%d >= %d) in OBJ_DESC[%s] in %s -- %s:%d --\n", trace_childs[0], desc->n_child, desc->tag, __FUNCTION__, __FILE__, __LINE__);
			return desc;
		}
		if(trace_cnts[0] >= *cnt){
			fprintf(stderr, " -- Illegal trace cont (%llu >= %llu) in OBJ_DESC[%s] in %s -- %s:%d --\n", (u8i)trace_cnts[0], (u8i)*cnt, desc->tag, __FUNCTION__, __FILE__, __LINE__);
			return desc;
		}
	}
	if(desc->n_child == 0){
		switch(mem_type){
			case 2:
			case 3:
				*addr += (*cnt) * mem_size_round(desc->size);
			default: return desc;
		}
	}
	ret = desc;
	ns = selected? trace_cnts[0] + 1 : *cnt;
	for(m=0;m<ns;m++){
		if(mem_type & 0x02){
			ptr = ((void**)obj) + m;
			if(*ptr == NULL){
				if(selected && m + 1 == ns){
					*_mem_type = desc->mem_type[trace_childs[0]] | (desc->size? 0 : MEM_PTR_TYPE_DUMP);
					*cnt = 0;
				}
				continue;
			}
			ref = (void*)*addr;
			*addr += mem_size_round(desc->size);
		} else {
			ref = obj + m * desc->size;
		}
		nc = (selected && m + 1 == ns)? trace_childs[0] + 1 : (u4i)desc->n_child;
		for(i=0;i<nc;i++){
			ptr = ref + desc->addr[i];
			if(desc->mem_type[i] & 0x01){
				*_mem_type = desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP);
				if(*ptr == NULL){
					*cnt = 0;
					continue;
				}
				*cnt = desc->cnt? desc->cnt(ref, i) : 1;
				ret = mem_locate_obj((void*)*addr, _mem_type, desc->desc[i], trace_childs + 1, trace_cnts + 1, ntrace - 1, addr, cnt, (selected && m + 1 == ns && i == nc));
			} else {
				*_mem_type = desc->mem_type[i] | (desc->size? 0 : MEM_PTR_TYPE_DUMP);
				*cnt = desc->cnt? desc->cnt(ref, i) : 1;
				ret = mem_locate_obj((void*)ptr, _mem_type, desc->desc[i], trace_childs + 1, trace_cnts + 1, ntrace - 1, addr, cnt, (selected && m + 1 == ns && i == nc));
			}
		}
	}
	return ret;
}

static inline size_t mem_dump_obj_file(void *obj, size_t mem_type, const obj_desc_t *desc, size_t cnt, size_t aux_data, FILE *out){
	size_t size;
	if(desc == NULL) return 0;
	if((mem_type & 0x01) == 0){
		fprintf(stderr, " -- Illegal mem_type (%u) to call mem_dump, object should have standalone memory in %s -- %s:%d --\n", (int)mem_type, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	size = mem_size_obj(obj, mem_type, desc, 0, cnt);
	if(out){
		fwrite(&size, sizeof(size_t), 1, out);
		fwrite(&mem_type, sizeof(size_t), 1, out);
		fwrite(&cnt, sizeof(size_t), 1, out);
		fwrite(&aux_data, sizeof(size_t), 1, out);
	}
	size = 4 * sizeof(size_t);
	size = mem_dump_obj(obj, mem_type, desc, size, cnt, out, 0);
	if(out) fflush(out);
	return size;
}

static inline size_t mem_dump_free_obj_file(void *obj, size_t mem_type, const obj_desc_t *desc, size_t cnt, size_t aux_data, FILE *out){
	size_t size;
	if(desc == NULL) return 0;
	if((mem_type & 0x01) == 0){
		fprintf(stderr, " -- Illegal mem_type (%u) to call mem_dump, object should have standalone memory in %s -- %s:%d --\n", (int)mem_type, __FUNCTION__, __FILE__, __LINE__);
		exit(1);
	}
	size = mem_size_obj(obj, mem_type, desc, 0, cnt);
	if(out){
		fwrite(&size, sizeof(size_t), 1, out);
		fwrite(&mem_type, sizeof(size_t), 1, out);
		fwrite(&cnt, sizeof(size_t), 1, out);
		fwrite(&aux_data, sizeof(size_t), 1, out);
	}
	size = 4 * sizeof(size_t);
	size = mem_dump_obj(obj, mem_type, desc, size, cnt, out, 1);
	if(out) fflush(out);
	return size;
}

static char *mem_share_locks = NULL;
static int mem_share_lock_size = 0;

static inline void cleanup_mem_share_file_locks(){
	int off;
	off = 0;
	while(off < mem_share_lock_size){
		shm_unlink(mem_share_locks + off);
		off += strlen(mem_share_locks + off) + 1;
	}
	if(mem_share_locks) free(mem_share_locks);
	mem_share_lock_size = 0;
	mem_share_locks = NULL;
}

#ifndef sighandler_t
typedef void (*sighandler_t)(int sig);
#endif
static sighandler_t sig_term = SIG_IGN;
static sighandler_t sig_int  = SIG_IGN;
static sighandler_t sig_hup  = SIG_IGN;
static sighandler_t sig_kill = SIG_IGN;
static volatile sig_atomic_t cleanup_mem_share_in_progress = 0;

static inline void sig_cleanup_mem_share_file_locks(int sig){
	if(cleanup_mem_share_in_progress) raise(sig);
	cleanup_mem_share_in_progress = 1;
	cleanup_mem_share_file_locks();
	signal(SIGTERM, sig_term);
	signal(SIGINT , sig_int);
	signal(SIGHUP, sig_hup);
	signal(SIGKILL, sig_kill);
	raise(sig);
}

static inline void register_mem_share_file_lock(char *file){
	int len;
	if(mem_share_lock_size == 0){
		if((sig_term = signal(SIGTERM, sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGTERM, SIG_IGN);
		if((sig_int  = signal(SIGINT , sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGINT , SIG_IGN);
		if((sig_hup  = signal(SIGHUP , sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGHUP , SIG_IGN);
		if((sig_kill = signal(SIGKILL, sig_cleanup_mem_share_file_locks)) == SIG_IGN) signal(SIGKILL, SIG_IGN);
		atexit(cleanup_mem_share_file_locks);
	}
	len = strlen(file);
	mem_share_locks = realloc(mem_share_locks, mem_share_lock_size + len + 1);
	strcpy(mem_share_locks + mem_share_lock_size, file);
	mem_share_lock_size += len + 1;
}

static inline void print_tree_obj_file(FILE *out, const obj_desc_t *desc, char *path, size_t max_cnt, int max_level){
	FILE *file;
	void *mem;
	size_t psize, *size, *mem_type, *cnt, *aux_data;
	if(desc == NULL) return;
	if((file = fopen(path, "r")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return;
	}
	size = alloca(sizeof(size_t));
	mem_type = alloca(sizeof(size_t));
	cnt = alloca(sizeof(size_t));
	aux_data = alloca(sizeof(size_t));
	fread(size, sizeof(size_t), 1, file);
	fread(mem_type, sizeof(size_t), 1, file);
	fread(cnt, sizeof(size_t), 1, file);
	fread(aux_data, sizeof(size_t), 1, file);
	psize = getpagesize();
	mem = mmap(0, (((*size) + 4 * sizeof(size_t)) + psize - 1) / psize * psize, PROT_READ, MAP_PRIVATE, fileno(file), 0);
	if(mem == MAP_FAILED){
		perror("Cannot mmap");
		return;
	}
	fprintf(out, "OBJ_TREE[%s]{\n", path);
	mem_tree_obj(out, mem + 4 * sizeof(size_t), *mem_type, desc, 0, *cnt, max_cnt, 1, max_level);
	fprintf(out, "}\n");
	munmap(mem, (((*size) + 4 * sizeof(size_t)) + psize - 1) / psize * psize);
	fclose(file);
}

// Directly read from file, don't share this object
static inline void* mem_read_obj_file(const obj_desc_t *desc, char *path, size_t *size, size_t *mem_type, size_t *cnt, size_t *aux_data){
	void *mem;
	size_t nin;
	FILE *file;
	if(desc == NULL) return NULL;
	if((file = fopen(path, "r")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return NULL;
	}
	if(size == NULL) size = alloca(sizeof(size_t));
	if(mem_type == NULL) mem_type = alloca(sizeof(size_t));
	if(cnt == NULL) cnt = alloca(sizeof(size_t));
	if(aux_data == NULL) aux_data = alloca(sizeof(size_t));
	fread(size, sizeof(size_t), 1, file);
	fread(mem_type, sizeof(size_t), 1, file);
	fread(cnt, sizeof(size_t), 1, file);
	fread(aux_data, sizeof(size_t), 1, file);
	mem = malloc(*size);
	if(mem == NULL){
		fprintf(stderr, " -- Cannot alloc %llu bytes memory for %s in %s -- %s:%d --\n", (unsigned long long)*size, path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return NULL;
	}
	if((nin = fread_stepwise(mem, 1, *size, file)) != *size){
		fprintf(stderr, " -- Read %llu bytes, not %llu bytes in %s -- %s:%d --\n", (unsigned long long)nin, (unsigned long long)*size, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return NULL;
	}
	fclose(file);
	mem_load_obj(mem, *aux_data, *mem_type, desc, 0, *cnt);
	//if(desc->post) desc->post(mem, *aux_data);
	return mem;
}

// only read a child object specified by indexs and nidx
static inline void* mem_read_sub_obj_file(const obj_desc_t *desc, u4i *trace_childs, size_t *trace_cnts, u4i ntrace, char *path, size_t *size){
	const obj_desc_t *sub;
	void *mem;
	size_t *mem_type, *cnt, *aux_data, psize, addr, sz, nin;
	FILE *file;
	if(desc == NULL) return NULL;
	if((file = fopen(path, "r")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	if(size == NULL) size = alloca(sizeof(size_t));
	mem_type = alloca(sizeof(size_t));
	cnt = alloca(sizeof(size_t));
	aux_data = alloca(sizeof(size_t));
	fread(size, sizeof(size_t), 1, file);
	fread(mem_type, sizeof(size_t), 1, file);
	fread(cnt, sizeof(size_t), 1, file);
	fread(aux_data, sizeof(size_t), 1, file);
	psize = getpagesize();
	mem = mmap(0, (((*size) + 4 * sizeof(size_t)) + psize - 1) / psize * psize, PROT_READ, MAP_PRIVATE, fileno(file), 0);
	if(mem == MAP_FAILED){
		perror("Cannot mmap");
		return NULL;
	}
	addr = 0;
	sub = mem_locate_obj(mem + 4 * sizeof(size_t), mem_type, desc, trace_childs, trace_cnts, ntrace, &addr, cnt, 1);
	sz  = mem_size_obj((void*)addr, *mem_type, sub, 0, *cnt);
	munmap(mem, (((*size) + 4 * sizeof(size_t)) + psize - 1) / psize * psize);
	mem = malloc(sz);
	fseek(file, addr - ((size_t)mem) + 4 * sizeof(size_t), SEEK_SET);
	if((nin = fread_stepwise(mem, 1, sz, file)) != sz){
		fprintf(stderr, " -- Read %llu bytes, not %llu bytes in %s -- %s:%d --\n", (unsigned long long)nin, (unsigned long long)sz, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		return NULL;
	}
	fclose(file);
	mem_load_obj(mem, *aux_data, *mem_type, sub, 0, *cnt);
	//if(sub->post) sub->post(mem, *aux_data);
	*size = sz;
	return mem;
}

// WARNNING ** the content of shared memory object can be modified by any mmaped process, thus is not multiple-process-safety
static inline void* mem_load_obj_file(const obj_desc_t *desc, char *_path, size_t *size, size_t *mem_type, size_t *cnt, size_t *aux_data){
	void *mem;
	size_t sz, nin, psize, *msg;
	char *lock, *shadow, *path, *shmp;
	char hostname[65];
	FILE *file;
	int fd, ret;
	if(desc == NULL) return NULL;
	path = absolute_filename(_path);
	shmp = strdup(path);
	replace_char(shmp, '/', '_', 0);
	if((file = fopen(path, "r+")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", path, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	if(size == NULL) size = alloca(sizeof(size_t));
	if(mem_type == NULL) mem_type = alloca(sizeof(size_t));
	if(cnt == NULL) cnt = alloca(sizeof(size_t));
	if(aux_data == NULL) aux_data = alloca(sizeof(size_t));
	fread(size, sizeof(size_t), 1, file);
	fread(mem_type, sizeof(size_t), 1, file);
	fread(cnt, sizeof(size_t), 1, file);
	fread(aux_data, sizeof(size_t), 1, file);
	psize = getpagesize();
	sz = (*size) + 4 * sizeof(size_t);
	//fd = fileno(file);
	//mem = mmap(0, (size + 4 * sizeof(size_t) + psize - 1) / psize * psize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	gethostname(hostname, 64);
	shadow = alloca(strlen(shmp) + strlen(hostname) + 40);
	sprintf(shadow, "%s.mem_share.%s.%ld.shm", shmp, hostname, gethostid());
	fd = shm_open(shadow, O_CREAT | O_RDWR, 0777);
	if(fd == -1){
		fprintf(stderr, " -- shm_open failed: %s in %s -- %s:%d --\n", shadow, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	ret = ftruncate(fd, (sz + psize - 1) / psize * psize);
	if(ret == -1){
		fprintf(stderr, " -- ftruncate failed: %s in %s -- %s:%d --\n", shadow, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	//register_mem_share_file_lock(shadow);
	mem = mmap(0, (sz + psize - 1) / psize * psize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if(mem == MAP_FAILED){
		perror("Cannot mmap");
		return NULL;
	}
	fseek(file, 0, SEEK_SET);
	if((nin = fread_stepwise(mem, 1, sz, file)) != sz){
		fprintf(stderr, " -- Read %llu bytes, not %llu bytes in %s -- %s:%d --\n", (unsigned long long)nin, (unsigned long long)sz, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	fclose(file);
	fprintf(stderr, " -- Read %llu bytes from %s --\n", (unsigned long long)nin, path); fflush(stderr);
	mem_load_obj(mem + 4 * sizeof(size_t), *aux_data, *mem_type, desc, 0, *cnt);
	//if(desc->post) desc->post(mem + 4 * sizeof(size_t), *aux_data);
	lock = alloca(strlen(shmp) + strlen(hostname) + 20);
	sprintf(lock, "%s.mem_share.%s.%ld", shmp, hostname, gethostid());
	if((fd = shm_open(lock, O_CREAT | O_RDWR, 0777)) == -1){
		fprintf(stderr, " -- shm_open failed: %s in %s -- %s:%d --\n", lock, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	if((ret = ftruncate(fd, psize)) == -1){
		fprintf(stderr, " -- ftruncate failed: %s in %s -- %s:%d --\n", lock, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); exit(1);
	}
	msg = mmap(0, psize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if(msg == MAP_FAILED){
		perror("Cannot mmap");
		exit(1);
	}
	msg[0] = (size_t)mem;
	msg[1] = *size;
	fprintf(stderr, "-- mem_share '%s' at %p:%llu --\n", path, mem, (u8i)*size); fflush(stderr);
	free(path);
	free(shmp);
	//register_mem_share_file_lock(lock);
	return mem + 4 * sizeof(size_t);
}

static inline int mem_stop_obj_file(char *_path){
	char *lock, *shadow, *path, *shmp;
	char hostname[65];
	path = absolute_filename(_path);
	shmp = strdup(path);
	replace_char(shmp, '/', '_', 0);
	gethostname(hostname, 64);
	shadow = alloca(strlen(shmp) + strlen(hostname) + 40);
	sprintf(shadow, "%s.mem_share.%s.%ld.shm", shmp, hostname, gethostid());
	if(shm_unlink(shadow) == -1){
		fprintf(stderr, " -- Failed to remove mmap object %s --\n", shadow); fflush(stderr);
		return 0;
	}
	lock = alloca(strlen(shmp) + strlen(hostname) + 20);
	sprintf(lock, "%s.mem_share.%s.%ld", shmp, hostname, gethostid());
	if(shm_unlink(lock) == -1){
		fprintf(stderr, " -- Failed to remove mmap object %s --\n", lock); fflush(stderr);
		return 0;
	}
	free(path);
	free(shmp);
	return 1;
}

// wr: whether to obtain write permission the mmap object
static inline void* mem_find_obj_file(const obj_desc_t *desc, char *_path, size_t *size, size_t *mem_type, size_t *cnt, size_t *aux_data, int wr){
	char *lock, *shadow, *path, *shmp;
	char hostname[65];
	void *addr, *mem;
	size_t psize, *msg;
	int fd, prot;
	UNUSED(desc);
	gethostname(hostname, 64);
	psize = getpagesize();
	path = absolute_filename(_path);
	shmp = strdup(path);
	replace_char(shmp, '/', '_', 0);
	lock = alloca(strlen(shmp) + strlen(hostname) + 32);
	sprintf(lock, "%s.mem_share.%s.%ld", shmp, hostname, gethostid());
	if(size == NULL) size = alloca(sizeof(size_t));
	if(mem_type == NULL) mem_type = alloca(sizeof(size_t));
	if(cnt == NULL) cnt = alloca(sizeof(size_t));
	if(aux_data == NULL) aux_data = alloca(sizeof(size_t));
	if((fd = shm_open(lock, O_RDWR, 0777)) == -1){
		return NULL;
	}
	msg = mmap(0, psize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if(msg == MAP_FAILED){
		perror("Cannot mmap");
		return NULL;
	}
	addr = (void*)msg[0];
	*size = msg[1];
	munmap(msg, psize);
	shadow = alloca(strlen(shmp) + strlen(hostname) + 40);
	sprintf(shadow, "%s.mem_share.%s.%ld.shm", shmp, hostname, gethostid());
	fd = shm_open(shadow, O_RDWR, 0777);
	prot = PROT_READ;
	if(wr) prot |= PROT_WRITE;
	mem = mmap(addr, ((*size) + 4 * sizeof(size_t) + psize - 1) / psize * psize, prot, MAP_SHARED | MAP_FIXED, fd, 0);
	if(mem == MAP_FAILED){
		perror("Cannot map shared object");
		return NULL;
	}
	*mem_type = ((size_t*)addr)[1];
	*cnt      = ((size_t*)addr)[2];
	*aux_data = ((size_t*)addr)[3];
	//fclose(file);
	//mem_load_obj(mem + 4 * sizeof(size_t), *mem_type, desc, 0, *cnt);
	fprintf(stderr, "-- mem_map '%s' at %p:%llu --\n", path, addr, (u8i)*size); fflush(stderr);
	free(path);
	free(shmp);
	return mem + 4 * sizeof(size_t);
}

#endif

/*
 * An example to use mem_dump
*/

/*

#include "mem_share.h"

typedef struct {
	char *str;
	int val;
} Type1;

size_t type1_count(void *obj, int idx){ if(idx == 0){ return strlen(((Type1*)obj)->str) + 1; } else return 0; }
//const obj_desc_t type1_obj_desc = {"TYPE1", sizeof(Type1), 1, {1}, {offsetof(Type1, str)}, {&OBJ_DESC_DATA}, type1_count, NULL};
const obj_desc_t type1_obj_desc = {"TYPE1", sizeof(Type1), 1, {1}, {offsetof(Type1, str)}, {&OBJ_DESC_CHAR_ARRAY}, NULL, NULL};

typedef struct {
	int a, b, c;
	Type1 d1, d2[10], *d3, *d4[10], **d5;
	char **strs;
	int d3len, d5len;
} Type2;

size_t type2_count(void *obj, int idx){
	switch(idx){
		case 0: return 1;
		case 1: return 10;
		case 2: return ((Type2*)obj)->d3len;
		case 3: return 10;
		case 4: return ((Type2*)obj)->d5len;
		default: return 10;
	}
}

const obj_desc_t type2_obj_desc = {"TYPE2", sizeof(Type2), 6, {0, 0, 1, 2, 3, 3}, {offsetof(Type2, d1), offsetof(Type2, d2), offsetof(Type2, d3), offsetof(Type2, d4), offsetof(Type2, d5), offsetof(Type2, strs)}, {&type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &OBJ_DESC_CHAR_ARRAY}, type2_count, NULL};

int main(){
	Type2 *t2, *t3;
	t2 = calloc(1, sizeof(Type2));
	int idx = 0;
	t2->d1.val = idx ++;
	t2->d1.str = strdup("d1");
	int i;
	for(i=0;i<10;i++){
		t2->d2[i].val = idx ++;
		t2->d2[i].str = strdup("d2");
	}
	t2->d3len = 10;
	t2->d3 = malloc(sizeof(Type1) * t2->d3len);
	for(i=0;i<t2->d3len;i++){
		t2->d3[i].val = idx ++;
		t2->d3[i].str = strdup("d3");
	}
	for(i=0;i<10;i++){
		t2->d4[i] = malloc(sizeof(Type1));
		t2->d4[i]->val = idx ++;
		t2->d4[i]->str = strdup("d4");
	}
	t2->d5len = 10;
	t2->d5 = malloc(sizeof(Type1*) * t2->d5len);
	for(i=0;i<t2->d5len;i++){
		t2->d5[i] = malloc(sizeof(Type1));
		t2->d5[i]->val = idx ++;
		t2->d5[i]->str = strdup("d5");
	}
	t2->strs = malloc(sizeof(char*) * 10);
	for(i=0;i<10;i++){
		t2->strs[i] = malloc(32);
		sprintf(t2->strs[i], "strs[%d,%d]", i, idx ++);
	}
	size_t aux_data, size, cnt, mem_type;
	FILE *file;
	size = mem_size_obj(t2, 1, &type2_obj_desc, 0, 1);
	fprintf(stdout, " -- size = %d in %s -- %s:%d --\n", (int)size, __FUNCTION__, __FILE__, __LINE__);
	aux_data = 1000999900;
	file = fopen("test.mem_share", "w");
	size = mem_dump_free_obj_file(t2, 1, &type2_obj_desc, 1, aux_data, file);
	fclose(file);
	fprintf(stdout, " -- size = %d in %s -- %s:%d --\n", (int)(size - 4 * sizeof(size_t)), __FUNCTION__, __FILE__, __LINE__);
	t3 = mem_read_obj_file(&type2_obj_desc, "test.mem_share", &mem_type, &cnt, &aux_data);
	fprintf(stdout, " -- aux_data = %d in %s -- %s:%d --\n", (int)aux_data, __FUNCTION__, __FILE__, __LINE__);
	printf("%d %s\n", t3->d1.val, t3->d1.str);
	for(i=0;i<10;i++) printf("%d %s\n", t3->d2[i].val, t3->d2[i].str);
	for(i=0;i<t3->d3len;i++) printf("%d %s\n", t3->d3[i].val, t3->d3[i].str);
	for(i=0;i<10;i++) printf("%d %s\n", t3->d4[i]->val, t3->d4[i]->str);
	for(i=0;i<t3->d5len;i++) printf("%d %s\n", t3->d5[i]->val, t3->d5[i]->str);
	for(i=0;i<10;i++) printf("%s\n", t3->strs[i]);
	free(t3);
	return 0;
}

*/
