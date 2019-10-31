VERSION=2.5
RELEASE=20190621

CC  := gcc
BIN := /usr/local/bin

ifeq (0, ${MAKELEVEL})
TIMESTAMP=$(shell date)
endif

ifeq (1, ${DEBUG})
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O0 -DDEBUG=1 -DVERSION="$(VERSION)" -DRELEASE="$(RELEASE)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -msse4.2
else
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O4 -DVERSION="$(VERSION)" -DRELEASE="$(RELEASE)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -msse4.2
endif

GLIBS=-lm -lrt -lpthread -lz
GENERIC_SRC=mem_share.h chararray.h sort.h list.h pgzf.h sort.h list.h dna.h thread.h filereader.h filewriter.h bitvec.h bit2vec.h bitsvec.h hashset.h

PROGS=kbm2 wtdbg2 wtdbg-cns wtpoa-cns pgzf

all: $(PROGS)

kbm2: $(GENERIC_SRC) kbm.c kbm.h kbmpoa.h wtpoa.h tripoa.h poacns.h kswx.h ksw.h ksw.c
	$(CC) $(CFLAGS) -o $@ kbm.c ksw.c $(GLIBS)

wtdbg2: $(GENERIC_SRC) wtdbg.c wtdbg-graph.h wtdbg.h kbm.h kswx.h ksw.h ksw.c kbmpoa.h wtpoa.h tripoa.h poacns.h
	$(CC) $(CFLAGS) -o $@ wtdbg.c ksw.c $(GLIBS)

wtdbg-cns: $(GENERIC_SRC) wtdbg-cns.c kswx.h ksw.h ksw.c dbgcns.h dagcns.h queue.h general_graph.h
	$(CC) $(CFLAGS) -o wtdbg-cns wtdbg-cns.c ksw.c $(GLIBS)

wtpoa-cns: $(GENERIC_SRC) wtpoa.h wtpoa-cns.c poacns.h tripoa.h ksw.h ksw.c
	$(CC) $(CFLAGS) -o $@ wtpoa-cns.c ksw.c $(GLIBS)

pgzf: mem_share.h sort.h list.h thread.h pgzf.h pgzf.c
	$(CC) $(CFLAGS) -o $@ pgzf.c $(GLIBS)

best_sam_hits4longreads: $(GENERIC_SRC) best_sam_hits4longreads.c
	$(CC) $(CFLAGS) -o $@ best_sam_hits4longreads.c $(GLIBS)

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out $(PROGS)

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out

install: $(PROGS)
	mkdir -p $(BIN) && cp -fvu $(PROGS) wtdbg2.pl $(BIN)
