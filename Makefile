#VERSION=2.8
#MINOR_VER=20160621
#MINOR_VER=20170330
#MINOR_VER=20170707

CC := gcc

ifeq (0, ${MAKELEVEL})
TIMESTAMP=$(shell date)
endif

ifeq (1, ${DEBUG})
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O0 -DTIMESTAMP="$(TIMESTAMP)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -msse4.2
else
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O4 -DTIMESTAMP="$(TIMESTAMP)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -msse4.2
endif

INSTALLDIR=/usr/local/bin
GLIBS=-lm -lrt -lpthread
GENERIC_SRC=mem_share.h string.h filereader.h bitvec.h bit2vec.h bitsvec.h hashset.h sort.h list.h dna.h thread.h

PROGS=kbm wtdbg2 wtdbg-cns wtpoa-cns

all: $(PROGS)

kbm: $(GENERIC_SRC) kbm.c kbm.h
	$(CC) $(CFLAGS) -o $@ kbm.c $(GLIBS)

wtdbg2: $(GENERIC_SRC) wtdbg.c kbm.h
	$(CC) $(CFLAGS) -o $@ wtdbg.c $(GLIBS)

wtdbg-cns: $(GENERIC_SRC) wtdbg-cns.c kswx.h ksw.h ksw.c dbgcns.h dagcns.h queue.h general_graph.h
	$(CC) $(CFLAGS) -o wtdbg-cns wtdbg-cns.c ksw.c $(GLIBS)

wtpoa-cns: $(GENERIC_SRC) wtpoa-cns.c poacns.h tripoa.h ksw.h ksw.c
	$(CC) $(CFLAGS) -o $@ wtpoa-cns.c ksw.c $(GLIBS)

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out $(PROGS)

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out

install:
	cp -fvu $(PROGS) $(INSTALLDIR)
