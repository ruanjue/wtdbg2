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
 
#ifndef __TIMER_RJ_H
#define __TIMER_RJ_H
#include <time.h>
//#include <sys/times.h>
#include <sys/time.h>
#include <sys/signal.h>
#include <errno.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#define def_counter(name) size_t name##_cnt; int name##_newline
#define beg_counter(name) name##_cnt = 0; name##_newline = 1
#define run_counter(name, tick, out, newline) {	\
	if(((name##_cnt) % tick) == 0){	\
		if(newline) fprintf(out, "[%s] %llu\n", date(), (unsigned long long)name##_cnt);	\
		else        fprintf(out, "\r[%s] %llu", date(), (unsigned long long)name##_cnt);	\
		fflush(out);	\
		name##_newline = newline;	\
	}	\
	name##_cnt ++;	\
}
#define end_counter(name, out) {	\
	if(name##_newline) fprintf(out, "[%s] %llu\n", date(), (unsigned long long)name##_cnt);	\
	else             fprintf(out, "\r[%s] %llu\n", date(), (unsigned long long)name##_cnt);	\
}

#define def_clock(name)	clock_t name##_t0, name##_t1, name##_t2
#define beg_clock(name)	name##_t0 = name##_t1 = name##_t2 = clock()
#define read_clock(name)	(name##_t1 = name##_t2, name##_t2 = clock(), name##_t2 - name##_t1)
#define count_clock(name)	(name##_t1 = name##_t2, name##_t2 = clock(), name##_t2 - name##_t0)
#define end_clock(name)		(name##_t1 = name##_t2, name##_t2 = clock(), name##_t2 - name##_t0)

typedef void (*timer_handler)(int signal, void *userdata);

#define TIMER_SIG_NULL	-1
#define TIMER_SIG_QUIT	0x7FFFFFFF

typedef struct {
	uint64_t elapsed_usec;
	struct timeval  tv;
	int pipes[2];
	void *userdata;
	timer_handler **handlers;
	int n_handler;
} Timer;

static inline Timer* init_timer(){
	Timer *timer;
	timer = malloc(sizeof(Timer));
	timer->elapsed_usec = 0;
	if(pipe(timer->pipes)){ free(timer); return NULL; }
	timer->userdata     = NULL;
	timer->handlers     = NULL;
	timer->n_handler    = 0;
	return timer;
}

static inline void free_timer(Timer *timer){
	close(timer->pipes[0]);
	close(timer->pipes[1]);
	if(timer->n_handler) free(timer->handlers);
	free(timer);
}

static inline void clear_listener_timer(Timer *timer){
	if(timer->handlers){
		free(timer->handlers);
		timer->handlers = NULL;
	}
	timer->n_handler = 0;
}

static inline void add_listener_timer(Timer *timer, timer_handler *handler){
	if(handler == NULL) return;
	timer->handlers = realloc(timer->handlers, sizeof(timer_handler*) * (timer->n_handler + 1));
	timer->handlers[timer->n_handler] = handler;
	timer->n_handler ++;
}

static inline void notify_listeners_timer(Timer *timer, int signal){
	int i;
	timer_handler *handler;
	for(i=0;i<timer->n_handler;i++){
		handler = timer->handlers[i];
		(*handler)(signal, timer->userdata);
	}
}

// RETURN: signal
static inline int once_timer(Timer *timer, uint32_t timeout){
	fd_set set;
	int signal;
	FD_ZERO(&set);
	FD_SET(timer->pipes[0], &set);
	timer->tv.tv_sec  = timeout / 1000000;
	timer->tv.tv_usec = timeout % 1000000;
	TEMP_FAILURE_RETRY(select(FD_SETSIZE, &set, NULL, NULL, timeout ? &timer->tv : NULL));
	if(read(timer->pipes[0], &signal, sizeof(int)) != sizeof(int)) signal = -1;
	if(signal != TIMER_SIG_QUIT) notify_listeners_timer(timer, signal);
	return signal;
}

static inline void notify_timer(Timer *timer, int signal){
	if(signal < 0) return;
	else write(timer->pipes[1], &signal, sizeof(int));
}

static inline int run_timer(Timer *timer, uint32_t interval, int times){
	fd_set set;
	int i;
	int signal;
	FD_ZERO(&set);
	FD_SET(timer->pipes[0], &set);
	timer->tv.tv_sec  = interval / 1000000;
	timer->tv.tv_usec = interval % 1000000;
	i = 0;
	while(i++ != times){
		TEMP_FAILURE_RETRY(select(FD_SETSIZE, &set, NULL, NULL, interval? &timer->tv : NULL));
		if(read(timer->pipes[0], &signal, sizeof(int)) != sizeof(int)) signal = -1;
		if(signal == TIMER_SIG_QUIT) break;
		notify_listeners_timer(timer, signal);
	}
	return i;
}

static inline void reset_timer(Timer *timer){
	int signal;
	while(read(timer->pipes[0], &signal, sizeof(int)) == sizeof(int));
}

static inline void stop_timer(Timer *timer){ notify_timer(timer, TIMER_SIG_QUIT); }

#endif
