# Source, Executable, Includes, Library Defines
VERSION = 1.00
CC = gcc
INCL = kseq.h sabre.h
#SRC = demulti_paired.c demulti_single.c sabre.c utils.c
SRC = sabre.c usage.c demultiplex.c utils.c
OBJ = $(SRC:.c=.o)
DSRC=src

CFLAGS = -Wall -O2 -std=c99 -pedantic -DVERSION=$(VERSION)
CFLAGSDEV = -Wall -O0 -ggdb -std=c99 -pedantic -DVERSION=$(VERSION)

LDFLAGS = -lz -lpthread
GPROF = -pg
EXE = sabre

.PHONY: default

default: build
# a smarter way to have an if statement here instead of explicit grpof target
# have a look at gcc -M

%.o: %.c
	$(CC) -c $(CFLAGS) $(SRC)

usage.o: usage.h sabre.h kseq.h
utils.o: utils.h sabre.h kseq.h
demultiplex.o: demultiplex.h utils.h sabre.h kseq.h
sabre.o: sabre.h

build: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LDFLAGS)
	#ln -sf $(DSRC)/$(EXE) ..

dev: $(OBJ)
	$(CC) $(CFLAGSDEV) $(OBJ) -o $(EXE)-dev $(LDFLAGS)

metrics:
	$(CC) $(CFLAGSDEV) -o metrics metrics.c $(LDFLAGS)

gprof:
	$(CC) $(CFLAGS) $(GPROF) $(SRC) -o $(EXE).gprof $(LDFLAGS)

clean:
	$(RM) $(OBJ) $(EXE) core

clean-all:
	$(RM) $(OBJ) $(EXE) $(EXE)-dev $(EXE).gprof core gmon.out metrics
