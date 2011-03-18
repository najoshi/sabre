PROGRAM_NAME = sabre
VERSION = 1.00
CC = gcc
CFLAGS = -Wall -pedantic -DVERSION=$(VERSION)
DEBUG = -g
OPT = -O3
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)
LDFLAGS = -lz
SDIR = src

.PHONY: clean default build distclean dist debug

default: build

barcode.o: $(SDIR)/barcode.c
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

demulti_single.o: $(SDIR)/demulti_single.c $(SDIR)/sabre.h $(SDIR)/kseq.h
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

demulti_paired.o: $(SDIR)/demulti_paired.c $(SDIR)/sabre.h $(SDIR)/kseq.h
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

sabre.o: $(SDIR)/sabre.c $(SDIR)/sabre.h
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

clean:
	rm -rf *.o $(SDIR)/*.gch ./sabre

distclean: clean
	rm -rf *.tar.gz

dist:
	tar -zcf $(ARCHIVE).tar.gz src Makefile

build: barcode.o demulti_single.o demulti_paired.o sabre.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(OPT) $? -o $(PROGRAM_NAME)

debug:
	$(MAKE) build "CFLAGS=-Wall -pedantic -g -DDEBUG"

